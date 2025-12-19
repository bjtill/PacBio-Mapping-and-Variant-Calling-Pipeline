#!/bin/bash
# Script: 03_mark_duplicates.sh
# Purpose: Mark PCR and optical duplicates in aligned PacBio BAM files
# Usage: ./03_mark_duplicates.sh -i input.bam -o output.bam [options]

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
Usage: $0 -i INPUT_BAM -o OUTPUT_BAM [OPTIONS]

Required arguments:
    -i, --input         Input sorted BAM file
    -o, --output        Output BAM file with duplicates marked

Optional arguments:
    -t, --threads       Number of threads (default: 8)
    -m, --metrics       Duplicate metrics output file
    --remove            Remove duplicates instead of marking
    --tool              Tool to use: pbmarkdup or samtools (default: pbmarkdup)
    -l, --log           Log file path
    -h, --help          Show this help message

Description:
    Marks or removes PCR and optical duplicates from aligned BAM files.
    
    pbmarkdup (recommended for PacBio):
        - PacBio-specific duplicate marking
        - Handles HiFi kinetics tags
        - More accurate for long reads
    
    samtools markdup (alternative):
        - Standard duplicate marking
        - Faster for some workflows
        - Works well for most cases

Example:
    $0 -i aligned.bam -o marked.bam -t 16 -m duplicates.txt

EOF
    exit 1
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

INPUT_BAM=""
OUTPUT_BAM=""
THREADS=8
METRICS_FILE=""
REMOVE_DUPS=false
TOOL="${PBMARKDUP:-pbmarkdup}"
TOOL_CHOICE="pbmarkdup"
LOG_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_BAM="$2"
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
        -m|--metrics)
            METRICS_FILE="$2"
            shift 2
            ;;
        --remove)
            REMOVE_DUPS=true
            shift
            ;;
        --tool)
            TOOL_CHOICE="$2"
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
[[ -z "$OUTPUT_BAM" ]] && die "Output BAM file not specified (-o)"

# Check if input file exists
check_file "$INPUT_BAM"

# Validate tool choice
case "$TOOL_CHOICE" in
    pbmarkdup)
        check_command "pbmarkdup"
        TOOL="pbmarkdup"
        ;;
    samtools)
        check_command "samtools"
        TOOL="samtools"
        ;;
    *)
        die "Invalid tool: $TOOL_CHOICE. Must be 'pbmarkdup' or 'samtools'"
        ;;
esac

# Set up logging
OUTPUT_DIR=$(dirname "$OUTPUT_BAM")
ensure_dir "$OUTPUT_DIR"

if [[ -z "$LOG_FILE" ]]; then
    LOG_FILE="${OUTPUT_DIR}/mark_duplicates.log"
fi

if [[ -z "$METRICS_FILE" ]]; then
    METRICS_FILE="${OUTPUT_BAM%.bam}.duplicates.txt"
fi

# Set up trap for cleanup
TEMP_FILES=""
trap cleanup EXIT

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

log_info "========================================="
log_info "Mark Duplicates"
log_info "========================================="
log_info "Input BAM: $INPUT_BAM"
log_info "Output BAM: $OUTPUT_BAM"
log_info "Tool: $TOOL_CHOICE"
log_info "Remove duplicates: $REMOVE_DUPS"
log_info "Threads: $THREADS"
log_info "Metrics file: $METRICS_FILE"
log_info "Log file: $LOG_FILE"
log_info "========================================="

start_timer

# Validate BAM file
log_info "Validating input BAM file..."
validate_bam "$INPUT_BAM" || die "BAM validation failed"

# Check if BAM is sorted
log_info "Checking if BAM is sorted..."
if ! samtools view -H "$INPUT_BAM" | grep -q "@HD.*SO:coordinate"; then
    log_warning "BAM file may not be coordinate sorted"
    log_warning "Sorting is recommended before marking duplicates"
fi

# Estimate disk space (output is similar size to input)
INPUT_SIZE_GB=$(du -BG "$INPUT_BAM" | cut -f1 | sed 's/G//')
REQUIRED_SPACE=$((INPUT_SIZE_GB * 2))
check_disk_space "$OUTPUT_DIR" "$REQUIRED_SPACE"

# Mark duplicates based on tool choice
if [[ "$TOOL_CHOICE" == "pbmarkdup" ]]; then
    log_info "Marking duplicates with pbmarkdup..."
    
    # Build pbmarkdup command
    PBMARKDUP_CMD=("pbmarkdup")
    PBMARKDUP_CMD+=(-j "$THREADS")
    
    if [[ "$REMOVE_DUPS" == true ]]; then
        log_info "Duplicates will be REMOVED (not just marked)"
    fi
    
    PBMARKDUP_CMD+=("$INPUT_BAM")
    PBMARKDUP_CMD+=("$OUTPUT_BAM")
    
    # Log command
    log_debug "Command: ${PBMARKDUP_CMD[*]}"
    
    # Run pbmarkdup
    "${PBMARKDUP_CMD[@]}" 2>&1 | tee -a "$LOG_FILE"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        die "pbmarkdup failed"
    fi
    
    # Extract duplicate metrics from pbmarkdup output
    if [[ -f "$LOG_FILE" ]]; then
        log_info "Extracting duplicate statistics..."
        grep -i "duplicate" "$LOG_FILE" > "$METRICS_FILE" || true
    fi
    
else  # samtools markdup
    log_info "Marking duplicates with samtools markdup..."
    
    # Build samtools markdup command
    MARKDUP_CMD=("samtools" "markdup")
    MARKDUP_CMD+=(-@ "$THREADS")
    MARKDUP_CMD+=(-s)  # Report stats
    
    if [[ "$REMOVE_DUPS" == true ]]; then
        MARKDUP_CMD+=(-r)  # Remove duplicates
        log_info "Duplicates will be REMOVED (not just marked)"
    fi
    
    MARKDUP_CMD+=(-f "$METRICS_FILE")  # Metrics file
    MARKDUP_CMD+=("$INPUT_BAM")
    MARKDUP_CMD+=("$OUTPUT_BAM")
    
    # Log command
    log_debug "Command: ${MARKDUP_CMD[*]}"
    
    # Run samtools markdup
    "${MARKDUP_CMD[@]}" 2>&1 | tee -a "$LOG_FILE"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        die "samtools markdup failed"
    fi
fi

# Index output BAM
log_info "Indexing output BAM..."
samtools index -@ "$THREADS" "$OUTPUT_BAM" 2>&1 | tee -a "$LOG_FILE"

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    die "BAM indexing failed"
fi

# Generate statistics
log_info "Generating alignment statistics..."
STATS_FILE="${OUTPUT_BAM%.bam}.stats.txt"
samtools stats -@ "$THREADS" "$OUTPUT_BAM" > "$STATS_FILE"

FLAGSTAT_FILE="${OUTPUT_BAM%.bam}.flagstat.txt"
samtools flagstat -@ "$THREADS" "$OUTPUT_BAM" > "$FLAGSTAT_FILE"

# Extract and report duplicate statistics
log_info "Duplicate Statistics:"
if [[ -f "$METRICS_FILE" ]]; then
    log_info "Metrics file: $METRICS_FILE"
    
    # Try to extract duplicate counts
    if [[ "$TOOL_CHOICE" == "samtools" ]]; then
        DUP_COUNT=$(grep "^SN" "$STATS_FILE" | grep "duplicates:" | cut -f3 || echo "N/A")
        TOTAL_READS=$(grep "^SN" "$STATS_FILE" | grep "sequences:" | cut -f3 || echo "N/A")
        
        if [[ "$DUP_COUNT" != "N/A" && "$TOTAL_READS" != "N/A" ]]; then
            DUP_RATE=$(echo "scale=2; $DUP_COUNT * 100 / $TOTAL_READS" | bc)
            log_info "  Total reads: $TOTAL_READS"
            log_info "  Duplicate reads: $DUP_COUNT"
            log_info "  Duplication rate: ${DUP_RATE}%"
        fi
    fi
fi

# Compare input and output sizes
INPUT_SIZE=$(du -h "$INPUT_BAM" | cut -f1)
OUTPUT_SIZE=$(du -h "$OUTPUT_BAM" | cut -f1)
log_info "Input size: $INPUT_SIZE"
log_info "Output size: $OUTPUT_SIZE"

stop_timer

log_info "========================================="
log_info "Duplicate marking completed successfully!"
log_info "========================================="

# YAD-friendly output
yad_success "Duplicate marking completed: $OUTPUT_BAM"

exit 0
