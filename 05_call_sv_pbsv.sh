#!/bin/bash
# Script: 05_call_sv_pbsv.sh
# Purpose: Call structural variants using pbsv (PacBio Structural Variant caller)
# Usage: ./05_call_sv_pbsv.sh -i input.bam -r reference.fasta -o output.vcf [options]

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
Usage: $0 -i INPUT_BAM -r REFERENCE -o OUTPUT_VCF [OPTIONS]

Required arguments:
    -i, --input         Input BAM file (must be sorted and indexed)
    -r, --reference     Reference genome FASTA file
    -o, --output        Output VCF file

Optional arguments:
    -t, --threads       Number of threads (default: $PBSV_THREADS)
    --min-length        Minimum SV length (default: $SV_MIN_LENGTH)
    --max-length        Maximum SV length (default: $SV_MAX_LENGTH)
    --tandem-repeats    Tandem repeat annotation BED file (optional)
    --sample            Sample name (default: from BAM)
    --types             SV types to call (comma-separated)
                        Options: DEL,INS,INV,DUP,BND (default: all)
    --ccs               Input is CCS/HiFi data (recommended)
    --pbsv              Path to pbsv (default: pbsv)
    -l, --log           Log file path
    -h, --help          Show this help message

Description:
    Calls structural variants (SVs) using pbsv, PacBio's official SV caller.
    pbsv works in two stages:
        1. discover - Find SV signatures in aligned reads
        2. call - Genotype and output SVs
    
    SV Types detected:
        DEL - Deletions
        INS - Insertions  
        INV - Inversions
        DUP - Duplications
        BND - Breakends (translocations)

Example:
    $0 -i sample.bam -r hg38.fasta -o sample.sv.vcf -t 16 --ccs

EOF
    exit 1
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

INPUT_BAM=""
REFERENCE="${REFERENCE_GENOME:-}"
OUTPUT_VCF=""
THREADS="${PBSV_THREADS:-16}"
MIN_SV_LENGTH="${SV_MIN_LENGTH:-50}"
MAX_SV_LENGTH="${SV_MAX_LENGTH:-100000}"
TANDEM_REPEATS=""
SAMPLE_NAME=""
SV_TYPES=""
IS_CCS=false
PBSV_BIN="${PBSV:-pbsv}"
LOG_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_BAM="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_VCF="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --min-length)
            MIN_SV_LENGTH="$2"
            shift 2
            ;;
        --max-length)
            MAX_SV_LENGTH="$2"
            shift 2
            ;;
        --tandem-repeats)
            TANDEM_REPEATS="$2"
            shift 2
            ;;
        --sample)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        --types)
            SV_TYPES="$2"
            shift 2
            ;;
        --ccs)
            IS_CCS=true
            shift
            ;;
        --pbsv)
            PBSV_BIN="$2"
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
[[ -z "$REFERENCE" ]] && die "Reference genome not specified (-r)"
[[ -z "$OUTPUT_VCF" ]] && die "Output VCF file not specified (-o)"

# Check if files exist
check_file "$INPUT_BAM"
check_file "$REFERENCE"

# Check for BAM index
if [[ ! -f "${INPUT_BAM}.bai" ]]; then
    die "BAM index not found: ${INPUT_BAM}.bai"
fi

# Check if tandem repeats file exists (if specified)
if [[ -n "$TANDEM_REPEATS" ]]; then
    check_file "$TANDEM_REPEATS"
fi

# Check if pbsv is available
check_command "$PBSV_BIN"

# Set up logging
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
ensure_dir "$OUTPUT_DIR"

if [[ -z "$LOG_FILE" ]]; then
    LOG_FILE="${OUTPUT_DIR}/pbsv.log"
fi

# Get sample name from BAM if not provided
if [[ -z "$SAMPLE_NAME" ]]; then
    SAMPLE_NAME=$(samtools view -H "$INPUT_BAM" | grep "^@RG" | head -1 | sed 's/.*SM:\([^\t]*\).*/\1/' || basename "$INPUT_BAM" .bam)
fi

# Set up trap for cleanup
TEMP_FILES=""
trap cleanup EXIT

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

log_info "========================================="
log_info "pbsv Structural Variant Calling"
log_info "========================================="
log_info "Input BAM: $INPUT_BAM"
log_info "Reference: $REFERENCE"
log_info "Output VCF: $OUTPUT_VCF"
log_info "Sample: $SAMPLE_NAME"
log_info "CCS/HiFi data: $IS_CCS"
log_info "SV length range: $MIN_SV_LENGTH - $MAX_SV_LENGTH bp"
log_info "Threads: $THREADS"
if [[ -n "$TANDEM_REPEATS" ]]; then
    log_info "Tandem repeats: $TANDEM_REPEATS"
fi
if [[ -n "$SV_TYPES" ]]; then
    log_info "SV types: $SV_TYPES"
fi
log_info "Log file: $LOG_FILE"
log_info "========================================="

start_timer

# Validate inputs
log_info "Validating inputs..."
validate_bam "$INPUT_BAM" || die "BAM validation failed"
validate_reference "$REFERENCE" || die "Reference validation failed"

# Estimate disk space
INPUT_SIZE_GB=$(du -BG "$INPUT_BAM" | cut -f1 | sed 's/G//')
REQUIRED_SPACE=$((INPUT_SIZE_GB * 2))
check_disk_space "$OUTPUT_DIR" "$REQUIRED_SPACE"

# Step 1: Discover SV signatures
log_info "Step 1/2: Discovering SV signatures..."

SVSIG_FILE="${OUTPUT_VCF%.vcf}.svsig.gz"

DISCOVER_CMD=("$PBSV_BIN" discover)
DISCOVER_CMD+=(--log-level INFO)

if [[ -n "$TANDEM_REPEATS" ]]; then
    DISCOVER_CMD+=(--tandem-repeats "$TANDEM_REPEATS")
fi

DISCOVER_CMD+=("$INPUT_BAM")
DISCOVER_CMD+=("$SVSIG_FILE")

# Log command
log_debug "Command: ${DISCOVER_CMD[*]}"

# Run pbsv discover
"${DISCOVER_CMD[@]}" 2>&1 | tee -a "$LOG_FILE"

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    die "pbsv discover failed"
fi

log_info "SV signature file created: $SVSIG_FILE"
log_info "Signature file size: $(get_file_size "$SVSIG_FILE")"

# Step 2: Call and genotype SVs
log_info "Step 2/2: Calling and genotyping SVs..."

CALL_CMD=("$PBSV_BIN" call)
CALL_CMD+=(--log-level INFO)
CALL_CMD+=(-j "$THREADS")
CALL_CMD+=(--min-sv-length "$MIN_SV_LENGTH")

if [[ "$IS_CCS" == true ]]; then
    CALL_CMD+=(--ccs)
fi

if [[ -n "$SV_TYPES" ]]; then
    CALL_CMD+=(--types "$SV_TYPES")
fi

CALL_CMD+=("$REFERENCE")
CALL_CMD+=("$SVSIG_FILE")
CALL_CMD+=("$OUTPUT_VCF")

# Log command
log_debug "Command: ${CALL_CMD[*]}"

# Run pbsv call
"${CALL_CMD[@]}" 2>&1 | tee -a "$LOG_FILE"

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    die "pbsv call failed"
fi

# Compress and index VCF if not already compressed
if [[ "$OUTPUT_VCF" != *.gz ]]; then
    log_info "Compressing VCF..."
    bgzip -@ "$THREADS" "$OUTPUT_VCF"
    OUTPUT_VCF="${OUTPUT_VCF}.gz"
fi

log_info "Indexing VCF..."
if command -v tabix &> /dev/null; then
    tabix -p vcf "$OUTPUT_VCF" 2>&1 | tee -a "$LOG_FILE"
else
    log_warning "tabix not found, VCF index not created"
fi

# Generate SV statistics
log_info "Generating SV statistics..."
STATS_FILE="${OUTPUT_VCF%.vcf.gz}.stats.txt"

if command -v bcftools &> /dev/null; then
    bcftools stats "$OUTPUT_VCF" > "$STATS_FILE" 2>&1
    
    # Extract key statistics
    TOTAL_SVS=$(grep "^SN" "$STATS_FILE" | grep "number of records:" | cut -f4 || echo "N/A")
    
    log_info "SV Statistics:"
    log_info "  Total SVs: $TOTAL_SVS"
    
    # Count by type
    if [[ -f "$OUTPUT_VCF" ]]; then
        for svtype in DEL INS INV DUP BND; do
            COUNT=$(zgrep -v "^#" "$OUTPUT_VCF" | grep "SVTYPE=$svtype" | wc -l || echo "0")
            if [[ "$COUNT" -gt 0 ]]; then
                log_info "  ${svtype}: $COUNT"
            fi
        done
    fi
    
    # Length distribution
    log_info "  Size distribution:"
    zgrep -v "^#" "$OUTPUT_VCF" | grep -o "SVLEN=[^;]*" | cut -d= -f2 | \
        awk '{
            if ($1 < 0) $1 = -$1;
            if ($1 < 100) count_0_100++;
            else if ($1 < 1000) count_100_1k++;
            else if ($1 < 10000) count_1k_10k++;
            else count_10k_plus++;
        }
        END {
            print "    50-100bp: " count_0_100;
            print "    100bp-1kb: " count_100_1k;
            print "    1kb-10kb: " count_1k_10k;
            print "    >10kb: " count_10k_plus;
        }' | tee -a "$LOG_FILE"
else
    log_warning "bcftools not found, skipping detailed statistics"
    
    # Basic count
    TOTAL_SVS=$(zgrep -v "^#" "$OUTPUT_VCF" | wc -l)
    log_info "Total SVs: $TOTAL_SVS"
fi

# Report output sizes
log_info "Output VCF size: $(get_file_size "$OUTPUT_VCF")"

stop_timer

log_info "========================================="
log_info "SV calling completed successfully!"
log_info "========================================="
log_info "Output files:"
log_info "  VCF: $OUTPUT_VCF"
log_info "  Signature: $SVSIG_FILE"
if [[ -f "$STATS_FILE" ]]; then
    log_info "  Stats: $STATS_FILE"
fi
log_info "========================================="

# YAD-friendly output
yad_success "SV calling completed: $OUTPUT_VCF"

exit 0
