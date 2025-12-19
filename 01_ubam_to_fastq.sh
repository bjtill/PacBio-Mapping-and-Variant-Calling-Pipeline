#!/bin/bash
# Script: 01_ubam_to_fastq.sh
# Purpose: Convert unmapped BAM (uBAM) to FASTQ format
# Usage: ./01_ubam_to_fastq.sh -i input.bam -o output.fastq.gz [options]

set -euo pipefail

# Source utility functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/utils.sh"

# Source configuration
if [[ -f "${SCRIPT_DIR}/../configs/pipeline_config.sh" ]]; then
    source "${SCRIPT_DIR}/../configs/pipeline_config.sh"
fi

# ==============================================================================
# USAGE
# ==============================================================================

usage() {
    cat << EOF
Usage: $0 -i INPUT_BAM -o OUTPUT_FASTQ [OPTIONS]

Required arguments:
    -i, --input         Input unmapped BAM file
    -o, --output        Output FASTQ file (will be gzipped)

Optional arguments:
    -t, --threads       Number of threads (default: 4)
    -s, --samtools      Path to samtools (default: samtools)
    -l, --log           Log file path (default: output_dir/ubam_to_fastq.log)
    -h, --help          Show this help message

Description:
    Converts unmapped BAM files to FASTQ format. The output will be
    automatically compressed with gzip.

Example:
    $0 -i sample.ubam -o sample.fastq.gz -t 8

EOF
    exit 1
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

INPUT_BAM=""
OUTPUT_FASTQ=""
THREADS=4
SAMTOOLS_BIN="${SAMTOOLS:-samtools}"
LOG_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_BAM="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FASTQ="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -s|--samtools)
            SAMTOOLS_BIN="$2"
            shift 2
            ;;
        -l|--log)
            LOG_FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# ==============================================================================
# VALIDATION
# ==============================================================================

# Check required arguments
[[ -z "$INPUT_BAM" ]] && die "Input BAM file not specified (-i)"
[[ -z "$OUTPUT_FASTQ" ]] && die "Output FASTQ file not specified (-o)"

# Check if input file exists
check_file "$INPUT_BAM"

# Check if samtools is available
check_command "$SAMTOOLS_BIN"

# Set up logging
OUTPUT_DIR=$(dirname "$OUTPUT_FASTQ")
ensure_dir "$OUTPUT_DIR"

if [[ -z "$LOG_FILE" ]]; then
    LOG_FILE="${OUTPUT_DIR}/ubam_to_fastq.log"
fi

# Set up trap for cleanup
trap cleanup EXIT

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

log_info "========================================="
log_info "uBAM to FASTQ Conversion"
log_info "========================================="
log_info "Input BAM: $INPUT_BAM"
log_info "Output FASTQ: $OUTPUT_FASTQ"
log_info "Threads: $THREADS"
log_info "Log file: $LOG_FILE"
log_info "========================================="

start_timer

# Validate BAM file
log_info "Validating input BAM file..."
if ! validate_bam "$INPUT_BAM"; then
    die "BAM validation failed"
fi

# Get BAM statistics
log_info "Getting BAM statistics..."
READ_COUNT=$("$SAMTOOLS_BIN" view -c "$INPUT_BAM")
log_info "Total reads in BAM: $READ_COUNT"

# Check available disk space (estimate: FASTQ is ~2x BAM size when compressed)
INPUT_SIZE_GB=$(du -BG "$INPUT_BAM" | cut -f1 | sed 's/G//')
REQUIRED_SPACE=$((INPUT_SIZE_GB * 3))
check_disk_space "$OUTPUT_DIR" "$REQUIRED_SPACE"

# Perform conversion
log_info "Converting uBAM to FASTQ..."

# Check if BAM is paired-end or single-end
PAIRED_COUNT=$("$SAMTOOLS_BIN" view -c -f 1 "$INPUT_BAM" || echo "0")

if [[ "$PAIRED_COUNT" -gt 0 ]]; then
    log_info "Detected paired-end reads"
    
    # For paired-end, create two output files
    OUTPUT_R1="${OUTPUT_FASTQ%.fastq.gz}_R1.fastq.gz"
    OUTPUT_R2="${OUTPUT_FASTQ%.fastq.gz}_R2.fastq.gz"
    
    log_info "Output R1: $OUTPUT_R1"
    log_info "Output R2: $OUTPUT_R2"
    
    "$SAMTOOLS_BIN" fastq \
        -@ "$THREADS" \
        -1 "$OUTPUT_R1" \
        -2 "$OUTPUT_R2" \
        -0 /dev/null \
        -s /dev/null \
        -n \
        "$INPUT_BAM" 2>&1 | tee -a "$LOG_FILE"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        die "FASTQ conversion failed"
    fi
    
    # Validate output files
    log_info "Validating output files..."
    validate_fastq "$OUTPUT_R1" || die "Output R1 validation failed"
    validate_fastq "$OUTPUT_R2" || die "Output R2 validation failed"
    
    # Count reads in output
    R1_READS=$(($(zcat "$OUTPUT_R1" | wc -l) / 4))
    R2_READS=$(($(zcat "$OUTPUT_R2" | wc -l) / 4))
    log_info "Reads in R1: $R1_READS"
    log_info "Reads in R2: $R2_READS"
    
    # Report file sizes
    log_info "Output R1 size: $(get_file_size "$OUTPUT_R1")"
    log_info "Output R2 size: $(get_file_size "$OUTPUT_R2")"
    
else
    log_info "Detected single-end reads"
    
    "$SAMTOOLS_BIN" fastq \
        -@ "$THREADS" \
        -0 "$OUTPUT_FASTQ" \
        -n \
        "$INPUT_BAM" 2>&1 | tee -a "$LOG_FILE"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        die "FASTQ conversion failed"
    fi
    
    # Validate output file
    log_info "Validating output file..."
    validate_fastq "$OUTPUT_FASTQ" || die "Output validation failed"
    
    # Count reads in output
    OUTPUT_READS=$(($(zcat "$OUTPUT_FASTQ" | wc -l) / 4))
    log_info "Reads in output: $OUTPUT_READS"
    
    # Verify read count matches
    if [[ "$OUTPUT_READS" -ne "$READ_COUNT" ]]; then
        log_warning "Read count mismatch: Input=$READ_COUNT, Output=$OUTPUT_READS"
    fi
    
    # Report file size
    log_info "Output size: $(get_file_size "$OUTPUT_FASTQ")"
fi

stop_timer

log_info "========================================="
log_info "Conversion completed successfully!"
log_info "========================================="

# YAD-friendly output
yad_success "Conversion completed: $OUTPUT_FASTQ"

exit 0
