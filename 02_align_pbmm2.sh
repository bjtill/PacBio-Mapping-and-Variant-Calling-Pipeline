#!/bin/bash
# Script: 02_align_pbmm2.sh
# Purpose: Align PacBio reads to reference genome using pbmm2
# Usage: ./02_align_pbmm2.sh -i input.fastq.gz -r reference.fasta -o output.bam [options]

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
Usage: $0 -i INPUT -r REFERENCE -o OUTPUT_BAM [OPTIONS]

Required arguments:
    -i, --input         Input file (FASTQ, FASTQ.gz, or uBAM)02_align_pbmm2.sh
    -r, --reference     Reference genome FASTA file
    -o, --output        Output BAM file

Optional arguments:
    -t, --threads       Number of threads (default: $ALIGN_THREADS)
    -p, --preset        Alignment preset (default: $ALIGN_PRESET)
                        Options: HiFi, SUBREAD, CCS, ISOSEQ, UNROLLED
    -s, --sample        Sample name for read group (default: from filename)
    -m, --memory        Sort memory per thread (default: $ALIGN_SORT_MEMORY)
    --pbmm2             Path to pbmm2 (default: pbmm2)
    --no-sort           Don't sort output (not recommended)
    --no-index          Don't create BAM index
    -l, --log           Log file path
    -h, --help          Show this help message

Description:
    Aligns PacBio reads to a reference genome using pbmm2 (minimap2 wrapper).
    Output is sorted and indexed by default.

Presets:
    HiFi     - PacBio HiFi reads (CCS, >Q20)
    SUBREAD  - PacBio subreads (CLR)
    CCS      - PacBio CCS reads (same as HiFi)
    ISOSEQ   - PacBio Iso-Seq full-length transcripts
    UNROLLED - Unrolled reads from SMRTbell libraries

Example:
    $0 -i sample.fastq.gz -r hg38.fasta -o sample.bam -t 16 -p HiFi

EOF
    exit 1
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

INPUT=""
REFERENCE="${REFERENCE_GENOME:-}"
OUTPUT_BAM=""
THREADS="${ALIGN_THREADS:-16}"
PRESET="${ALIGN_PRESET:-HiFi}"
SAMPLE_NAME=""
SORT_MEMORY="${ALIGN_SORT_MEMORY:-4G}"
PBMM2_BIN="${PBMM2:-pbmm2}"
DO_SORT=true
DO_INDEX=true
LOG_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_BAM="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -p|--preset)
            PRESET="$2"
            shift 2
            ;;
        -s|--sample)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        -m|--memory)
            SORT_MEMORY="$2"
            shift 2
            ;;
        --pbmm2)
            PBMM2_BIN="$2"
            shift 2
            ;;
        --no-sort)
            DO_SORT=false
            shift
            ;;
        --no-index)
            DO_INDEX=false
            shift
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
[[ -z "$INPUT" ]] && die "Input file not specified (-i)"
[[ -z "$REFERENCE" ]] && die "Reference genome not specified (-r)"
[[ -z "$OUTPUT_BAM" ]] && die "Output BAM file not specified (-o)"

# Check if files exist
check_file "$INPUT"
check_file "$REFERENCE"

# Check if pbmm2 is available
check_command "$PBMM2_BIN"

# Validate preset
case "$PRESET" in
    HiFi|SUBREAD|CCS|ISOSEQ|UNROLLED) ;;
    *)
        die "Invalid preset: $PRESET. Must be one of: HiFi, SUBREAD, CCS, ISOSEQ, UNROLLED"
        ;;
esac

# Set up logging
OUTPUT_DIR=$(dirname "$OUTPUT_BAM")
ensure_dir "$OUTPUT_DIR"

if [[ -z "$LOG_FILE" ]]; then
    LOG_FILE="${OUTPUT_DIR}/alignment.log"
fi

# Set sample name from filename if not provided
if [[ -z "$SAMPLE_NAME" ]]; then
    SAMPLE_NAME=$(basename "$INPUT" | sed 's/\.\(fastq\|fq\|bam\).*//')
fi

# Set up trap for cleanup
trap cleanup EXIT

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

log_info "========================================="
log_info "PacBio Alignment with pbmm2"
log_info "========================================="
log_info "Input: $INPUT"
log_info "Reference: $REFERENCE"
log_info "Output: $OUTPUT_BAM"
log_info "Sample: $SAMPLE_NAME"
log_info "Preset: $PRESET"
log_info "Threads: $THREADS"
log_info "Log file: $LOG_FILE"
log_info "========================================="

start_timer

# Validate reference
log_info "Validating reference genome..."
validate_reference "$REFERENCE"

# Check input file type and validate
log_info "Validating input file..."
if [[ "$INPUT" == *.bam ]]; then
    validate_bam "$INPUT" || die "Input BAM validation failed"
    INPUT_TYPE="BAM"
elif [[ "$INPUT" == *.fastq.gz || "$INPUT" == *.fq.gz ]]; then
    validate_fastq "$INPUT" || die "Input FASTQ validation failed"
    INPUT_TYPE="FASTQ"
elif [[ "$INPUT" == *.fastq || "$INPUT" == *.fq ]]; then
    validate_fastq "$INPUT" || die "Input FASTQ validation failed"
    INPUT_TYPE="FASTQ"
else
    die "Unknown input file type. Expected .bam, .fastq, or .fastq.gz"
fi
log_info "Input type: $INPUT_TYPE"

# Estimate required disk space (BAM is ~1.5x input size)
INPUT_SIZE_GB=$(du -BG "$INPUT" | cut -f1 | sed 's/G//')
REQUIRED_SPACE=$((INPUT_SIZE_GB * 2))
check_disk_space "$OUTPUT_DIR" "$REQUIRED_SPACE"

# Build pbmm2 command
log_info "Starting alignment..."

PBMM2_CMD=("$PBMM2_BIN" align)
PBMM2_CMD+=(--preset "$PRESET")
PBMM2_CMD+=(--sample "$SAMPLE_NAME")
PBMM2_CMD+=(-j "$THREADS")
PBMM2_CMD+=(-J "$THREADS")  # Sort threads

if [[ "$DO_SORT" == true ]]; then
    PBMM2_CMD+=(--sort)
    PBMM2_CMD+=(--sort-memory "$SORT_MEMORY")
fi

PBMM2_CMD+=("$REFERENCE")
PBMM2_CMD+=("$INPUT")
PBMM2_CMD+=("$OUTPUT_BAM")

# Log the command
log_debug "Command: ${PBMM2_CMD[*]}"

# Run pbmm2
"${PBMM2_CMD[@]}" 2>&1 | tee -a "$LOG_FILE"

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    die "pbmm2 alignment failed"
fi

# Create index if requested
if [[ "$DO_INDEX" == true ]]; then
    log_info "Creating BAM index..."
    samtools index -@ "$THREADS" "$OUTPUT_BAM" 2>&1 | tee -a "$LOG_FILE"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        die "BAM indexing failed"
    fi
fi

# Generate alignment statistics
log_info "Generating alignment statistics..."
STATS_FILE="${OUTPUT_BAM%.bam}.stats.txt"
samtools stats -@ "$THREADS" "$OUTPUT_BAM" > "$STATS_FILE" 2>&1

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    log_warning "Failed to generate stats file"
else
    log_info "Stats file: $STATS_FILE"
    
    # Extract key statistics
    TOTAL_READS=$(grep "^SN" "$STATS_FILE" | grep "sequences:" | cut -f3)
    MAPPED_READS=$(grep "^SN" "$STATS_FILE" | grep "reads mapped:" | cut -f3)
    
    if [[ -n "$TOTAL_READS" && -n "$MAPPED_READS" ]]; then
        MAPPING_RATE=$(echo "scale=2; $MAPPED_READS * 100 / $TOTAL_READS" | bc)
        log_info "Total reads: $TOTAL_READS"
        log_info "Mapped reads: $MAPPED_READS"
        log_info "Mapping rate: ${MAPPING_RATE}%"
    fi
fi

# Generate flagstat
log_info "Generating flagstat..."
FLAGSTAT_FILE="${OUTPUT_BAM%.bam}.flagstat.txt"
samtools flagstat -@ "$THREADS" "$OUTPUT_BAM" > "$FLAGSTAT_FILE" 2>&1

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    log_warning "Failed to generate flagstat file"
else
    log_info "Flagstat file: $FLAGSTAT_FILE"
fi

# Report output size
log_info "Output BAM size: $(get_file_size "$OUTPUT_BAM")"
if [[ -f "${OUTPUT_BAM}.bai" ]]; then
    log_info "Index size: $(get_file_size "${OUTPUT_BAM}.bai")"
fi

stop_timer

log_info "========================================="
log_info "Alignment completed successfully!"
log_info "========================================="

# YAD-friendly output
yad_success "Alignment completed: $OUTPUT_BAM"

exit 0
