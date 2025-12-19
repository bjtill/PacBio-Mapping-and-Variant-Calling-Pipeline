#!/bin/bash
# Script: 04_call_variants_deepvariant.sh
# Purpose: Call small variants (SNPs/Indels) using DeepVariant
# Usage: ./04_call_variants_deepvariant.sh -i input.bam -r reference.fasta -o output.vcf.gz [options]

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
    -o, --output        Output VCF file (will be gzipped)

Optional arguments:
    -t, --threads       Number of threads (default: $DEEPVARIANT_THREADS)
    -m, --model         Model type (default: $DEEPVARIANT_MODEL)
                        Options: PACBIO, WGS, WES, HYBRID_PACBIO_ILLUMINA
    --regions           BED file with regions to call (optional)
    --sample            Sample name (default: from BAM)
    --intermediate      Keep intermediate files (default: false)
    --docker            Use Docker instead of native installation
    --singularity       Use Singularity instead of native installation
    -l, --log           Log file path
    -h, --help          Show this help message

Description:
    Calls small variants (SNPs and indels) using Google DeepVariant.
    DeepVariant uses deep learning and is highly accurate for PacBio HiFi data.
    
    The output includes:
        - VCF file with variants
        - gVCF file with all sites (optional)
        - Visual report in HTML format

Models:
    PACBIO                  - For PacBio HiFi/CCS reads (recommended)
    WGS                     - For Illumina whole genome sequencing
    WES                     - For exome sequencing
    HYBRID_PACBIO_ILLUMINA  - For hybrid assemblies

Example:
    $0 -i sample.bam -r hg38.fasta -o sample.vcf.gz -t 16 -m PACBIO

EOF
    exit 1
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

INPUT_BAM=""
REFERENCE="${REFERENCE_GENOME:-}"
OUTPUT_VCF=""
THREADS="${DEEPVARIANT_THREADS:-16}"
MODEL="${DEEPVARIANT_MODEL:-PACBIO}"
REGIONS_BED=""
SAMPLE_NAME=""
KEEP_INTERMEDIATE=false
USE_DOCKER=false
USE_SINGULARITY=false
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
        -m|--model)
            MODEL="$2"
            shift 2
            ;;
        --regions)
            REGIONS_BED="$2"
            shift 2
            ;;
        --sample)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        --intermediate)
            KEEP_INTERMEDIATE=true
            shift
            ;;
        --docker)
            USE_DOCKER=true
            shift
            ;;
        --singularity)
            USE_SINGULARITY=true
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

# Check if regions file exists (if specified)
if [[ -n "$REGIONS_BED" ]]; then
    check_file "$REGIONS_BED"
fi

# Validate model
case "$MODEL" in
    PACBIO|WGS|WES|HYBRID_PACBIO_ILLUMINA) ;;
    *)
        die "Invalid model: $MODEL"
        ;;
esac

# Check for DeepVariant installation
if [[ "$USE_DOCKER" == false && "$USE_SINGULARITY" == false ]]; then
    if ! command -v run_deepvariant &> /dev/null; then
        die "DeepVariant not found. Install it or use --docker or --singularity"
    fi
fi

# Set up logging
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
ensure_dir "$OUTPUT_DIR"

if [[ -z "$LOG_FILE" ]]; then
    LOG_FILE="${OUTPUT_DIR}/deepvariant.log"
fi

# Get sample name from BAM if not provided
if [[ -z "$SAMPLE_NAME" ]]; then
    SAMPLE_NAME=$(samtools view -H "$INPUT_BAM" | grep "^@RG" | head -1 | sed 's/.*SM:\([^\t]*\).*/\1/' || basename "$INPUT_BAM" .bam)
fi

# Set up trap for cleanup
TEMP_DIR="${OUTPUT_DIR}/tmp_deepvariant"
TEMP_FILES="$TEMP_DIR"
trap cleanup EXIT

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

log_info "========================================="
log_info "DeepVariant Small Variant Calling"
log_info "========================================="
log_info "Input BAM: $INPUT_BAM"
log_info "Reference: $REFERENCE"
log_info "Output VCF: $OUTPUT_VCF"
log_info "Sample: $SAMPLE_NAME"
log_info "Model: $MODEL"
log_info "Threads: $THREADS"
log_info "Log file: $LOG_FILE"
if [[ -n "$REGIONS_BED" ]]; then
    log_info "Regions: $REGIONS_BED"
fi
log_info "========================================="

start_timer

# Validate inputs
log_info "Validating inputs..."
validate_bam "$INPUT_BAM" || die "BAM validation failed"
validate_reference "$REFERENCE" || die "Reference validation failed"

# Create temporary directory for intermediate files
ensure_dir "$TEMP_DIR"

# Estimate disk space (VCF is typically small, but intermediate files can be large)
INPUT_SIZE_GB=$(du -BG "$INPUT_BAM" | cut -f1 | sed 's/G//')
REQUIRED_SPACE=$((INPUT_SIZE_GB * 3))
check_disk_space "$OUTPUT_DIR" "$REQUIRED_SPACE"

# Prepare output filenames
OUTPUT_BASE="${OUTPUT_VCF%.vcf.gz}"
OUTPUT_BASE="${OUTPUT_BASE%.vcf}"
GVCF_OUTPUT="${OUTPUT_BASE}.g.vcf.gz"
VISUAL_REPORT="${OUTPUT_BASE}.visual_report.html"

# Build DeepVariant command
log_info "Running DeepVariant..."

if [[ "$USE_DOCKER" == true ]]; then
    log_info "Using Docker..."
    
    DOCKER_IMAGE="${DEEPVARIANT_DOCKER:-google/deepvariant:1.6.0}"
    
    # Prepare Docker volume mounts
    INPUT_BAM_DIR=$(dirname "$(readlink -f "$INPUT_BAM")")
    REF_DIR=$(dirname "$(readlink -f "$REFERENCE")")
    OUTPUT_DIR_ABS=$(readlink -f "$OUTPUT_DIR")
    
    docker run --rm \
        -v "$INPUT_BAM_DIR:$INPUT_BAM_DIR:ro" \
        -v "$REF_DIR:$REF_DIR:ro" \
        -v "$OUTPUT_DIR_ABS:$OUTPUT_DIR_ABS" \
        "$DOCKER_IMAGE" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type="$MODEL" \
        --ref="$REFERENCE" \
        --reads="$INPUT_BAM" \
        --output_vcf="${OUTPUT_BASE}.vcf.gz" \
        --output_gvcf="$GVCF_OUTPUT" \
        --num_shards="$THREADS" \
        --intermediate_results_dir="$TEMP_DIR" \
        $(if [[ -n "$REGIONS_BED" ]]; then echo "--regions=$REGIONS_BED"; fi) \
        2>&1 | tee -a "$LOG_FILE"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        die "DeepVariant (Docker) failed"
    fi
    
elif [[ "$USE_SINGULARITY" == true ]]; then
    log_info "Using Singularity..."
    
    SINGULARITY_IMAGE="${DEEPVARIANT_SINGULARITY:-docker://google/deepvariant:1.6.0}"
    
    singularity exec "$SINGULARITY_IMAGE" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type="$MODEL" \
        --ref="$REFERENCE" \
        --reads="$INPUT_BAM" \
        --output_vcf="${OUTPUT_BASE}.vcf.gz" \
        --output_gvcf="$GVCF_OUTPUT" \
        --num_shards="$THREADS" \
        --intermediate_results_dir="$TEMP_DIR" \
        $(if [[ -n "$REGIONS_BED" ]]; then echo "--regions=$REGIONS_BED"; fi) \
        2>&1 | tee -a "$LOG_FILE"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        die "DeepVariant (Singularity) failed"
    fi
    
else
    log_info "Using native installation..."
    
    run_deepvariant \
        --model_type="$MODEL" \
        --ref="$REFERENCE" \
        --reads="$INPUT_BAM" \
        --output_vcf="${OUTPUT_BASE}.vcf.gz" \
        --output_gvcf="$GVCF_OUTPUT" \
        --num_shards="$THREADS" \
        --intermediate_results_dir="$TEMP_DIR" \
        $(if [[ -n "$REGIONS_BED" ]]; then echo "--regions=$REGIONS_BED"; fi) \
        2>&1 | tee -a "$LOG_FILE"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        die "DeepVariant failed"
    fi
fi

# Generate VCF statistics
log_info "Generating VCF statistics..."
STATS_FILE="${OUTPUT_BASE}.vcf.stats.txt"

if command -v bcftools &> /dev/null; then
    bcftools stats "${OUTPUT_BASE}.vcf.gz" > "$STATS_FILE" 2>&1
    
    # Extract key statistics
    SNP_COUNT=$(grep "^SN" "$STATS_FILE" | grep "number of SNPs:" | cut -f4 || echo "N/A")
    INDEL_COUNT=$(grep "^SN" "$STATS_FILE" | grep "number of indels:" | cut -f4 || echo "N/A")
    TSTV_RATIO=$(grep "^TSTV" "$STATS_FILE" | cut -f5 || echo "N/A")
    
    log_info "Variant Statistics:"
    log_info "  SNPs: $SNP_COUNT"
    log_info "  Indels: $INDEL_COUNT"
    log_info "  Ts/Tv ratio: $TSTV_RATIO"
else
    log_warning "bcftools not found, skipping VCF statistics"
fi

# Clean up intermediate files if requested
if [[ "$KEEP_INTERMEDIATE" == false ]]; then
    log_info "Cleaning up intermediate files..."
    rm -rf "$TEMP_DIR"
fi

# Report output sizes
log_info "Output VCF size: $(get_file_size "${OUTPUT_BASE}.vcf.gz")"
if [[ -f "$GVCF_OUTPUT" ]]; then
    log_info "Output gVCF size: $(get_file_size "$GVCF_OUTPUT")"
fi

stop_timer

log_info "========================================="
log_info "Variant calling completed successfully!"
log_info "========================================="
log_info "Output files:"
log_info "  VCF: ${OUTPUT_BASE}.vcf.gz"
log_info "  gVCF: $GVCF_OUTPUT"
if [[ -f "$STATS_FILE" ]]; then
    log_info "  Stats: $STATS_FILE"
fi
log_info "========================================="

# YAD-friendly output
yad_success "Variant calling completed: ${OUTPUT_BASE}.vcf.gz"

exit 0
