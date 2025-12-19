#!/bin/bash
# Example YAD GUI wrapper for PacBio pipeline
# This demonstrates how to integrate the pipeline scripts with YAD

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# Source configuration
if [[ -f "${PIPELINE_DIR}/configs/pipeline_config.sh" ]]; then
    source "${PIPELINE_DIR}/configs/pipeline_config.sh"
fi

# ==============================================================================
# MAIN MENU
# ==============================================================================

show_main_menu() {
    CHOICE=$(yad --list \
        --title="PacBio Pipeline" \
        --text="Select a pipeline step:" \
        --width=600 --height=400 \
        --column="Step" --column="Description" \
        "1" "Convert uBAM to FASTQ" \
        "2" "Align reads (pbmm2)" \
        "3" "Mark duplicates" \
        "4" "Call small variants (DeepVariant)" \
        "5" "Call structural variants (pbsv)" \
        "6" "QC and coverage analysis" \
        "7" "Run complete pipeline" \
        --button="Exit:1" --button="Continue:0")
    
    if [[ $? -ne 0 ]]; then
        exit 0
    fi
    
    STEP=$(echo "$CHOICE" | cut -d'|' -f1)
    
    case "$STEP" in
        1) run_ubam_to_fastq ;;
        2) run_alignment ;;
        3) run_mark_duplicates ;;
        4) run_deepvariant ;;
        5) run_pbsv ;;
        6) run_qc ;;
        7) run_complete_pipeline ;;
        *) show_main_menu ;;
    esac
}

# ==============================================================================
# STEP 1: uBAM to FASTQ
# ==============================================================================

run_ubam_to_fastq() {
    # Get input file
    INPUT=$(yad --file \
        --title="Select Input uBAM" \
        --file-filter="BAM files|*.bam" \
        --file-filter="All files|*")
    
    if [[ -z "$INPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get output file
    DEFAULT_OUTPUT="${INPUT%.bam}.fastq.gz"
    OUTPUT=$(yad --file --save \
        --title="Save Output FASTQ" \
        --filename="$DEFAULT_OUTPUT" \
        --file-filter="FASTQ files|*.fastq.gz *.fq.gz")
    
    if [[ -z "$OUTPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get parameters
    PARAMS=$(yad --form \
        --title="uBAM to FASTQ Parameters" \
        --field="Threads:NUM" "16!1..64!1" \
        --button="Cancel:1" --button="Run:0")
    
    if [[ $? -ne 0 ]]; then
        show_main_menu
        return
    fi
    
    THREADS=$(echo "$PARAMS" | cut -d'|' -f1 | cut -d',' -f1)
    
    # Run conversion
    LOG_FILE="/tmp/ubam_to_fastq_$$.log"
    
    "${SCRIPT_DIR}/01_ubam_to_fastq.sh" \
        -i "$INPUT" \
        -o "$OUTPUT" \
        -t "$THREADS" \
        -l "$LOG_FILE" \
        2>&1 | tee "$LOG_FILE" | \
        yad --progress --pulsate --auto-close \
            --title="Converting uBAM to FASTQ" \
            --text="Processing: $(basename "$INPUT")"
    
    RESULT=$?
    
    if [[ $RESULT -eq 0 ]]; then
        yad --info --title="Success" \
            --text="Conversion completed!\n\nOutput: $OUTPUT" \
            --button="View Log:bash -c 'xdg-open $LOG_FILE'" \
            --button="OK:0"
    else
        yad --error --title="Error" \
            --text="Conversion failed!\n\nCheck log: $LOG_FILE"
    fi
    
    show_main_menu
}

# ==============================================================================
# STEP 2: Alignment
# ==============================================================================

run_alignment() {
    # Get input file
    INPUT=$(yad --file \
        --title="Select Input File" \
        --file-filter="FASTQ/BAM files|*.fastq *.fastq.gz *.fq *.fq.gz *.bam" \
        --file-filter="All files|*")
    
    if [[ -z "$INPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get reference
    REFERENCE=$(yad --file \
        --title="Select Reference Genome" \
        --filename="${REFERENCE_GENOME}" \
        --file-filter="FASTA files|*.fasta *.fa *.fna")
    
    if [[ -z "$REFERENCE" ]]; then
        show_main_menu
        return
    fi
    
    # Get output file
    DEFAULT_OUTPUT="${INPUT%.*}.bam"
    OUTPUT=$(yad --file --save \
        --title="Save Output BAM" \
        --filename="$DEFAULT_OUTPUT" \
        --file-filter="BAM files|*.bam")
    
    if [[ -z "$OUTPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get parameters
    PARAMS=$(yad --form \
        --title="Alignment Parameters" \
        --field="Sample Name" "$(basename "${INPUT%.*}")" \
        --field="Preset:CB" "HiFi!SUBREAD!CCS!ISOSEQ!UNROLLED" \
        --field="Threads:NUM" "16!1..64!1" \
        --field="Sort Memory (GB):NUM" "4!1..32!1" \
        --button="Cancel:1" --button="Run:0")
    
    if [[ $? -ne 0 ]]; then
        show_main_menu
        return
    fi
    
    SAMPLE=$(echo "$PARAMS" | cut -d'|' -f1)
    PRESET=$(echo "$PARAMS" | cut -d'|' -f2)
    THREADS=$(echo "$PARAMS" | cut -d'|' -f3 | cut -d',' -f1)
    MEMORY=$(echo "$PARAMS" | cut -d'|' -f4 | cut -d',' -f1)
    
    # Run alignment
    LOG_FILE="/tmp/alignment_$$.log"
    
    "${SCRIPT_DIR}/02_align_pbmm2.sh" \
        -i "$INPUT" \
        -r "$REFERENCE" \
        -o "$OUTPUT" \
        -s "$SAMPLE" \
        -p "$PRESET" \
        -t "$THREADS" \
        -m "${MEMORY}G" \
        -l "$LOG_FILE" \
        2>&1 | tee "$LOG_FILE" | \
        yad --progress --pulsate --auto-close \
            --title="Aligning Reads" \
            --text="Processing: $(basename "$INPUT")"
    
    RESULT=$?
    
    if [[ $RESULT -eq 0 ]]; then
        # Show statistics
        STATS_FILE="${OUTPUT%.bam}.stats.txt"
        if [[ -f "$STATS_FILE" ]]; then
            MAPPED=$(grep "reads mapped:" "$STATS_FILE" | cut -f3)
            TOTAL=$(grep "sequences:" "$STATS_FILE" | cut -f3)
            RATE=$(echo "scale=2; $MAPPED * 100 / $TOTAL" | bc)
            
            yad --info --title="Success" \
                --text="Alignment completed!\n\nOutput: $OUTPUT\nMapping rate: ${RATE}%" \
                --button="View Stats:bash -c 'xdg-open $STATS_FILE'" \
                --button="View Log:bash -c 'xdg-open $LOG_FILE'" \
                --button="OK:0"
        else
            yad --info --title="Success" \
                --text="Alignment completed!\n\nOutput: $OUTPUT"
        fi
    else
        yad --error --title="Error" \
            --text="Alignment failed!\n\nCheck log: $LOG_FILE"
    fi
    
    show_main_menu
}

# ==============================================================================
# STEP 3: Mark Duplicates
# ==============================================================================

run_mark_duplicates() {
    # Get input file
    INPUT=$(yad --file \
        --title="Select Input BAM" \
        --file-filter="BAM files|*.bam")
    
    if [[ -z "$INPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get output file
    DEFAULT_OUTPUT="${INPUT%.bam}.markdup.bam"
    OUTPUT=$(yad --file --save \
        --title="Save Output BAM" \
        --filename="$DEFAULT_OUTPUT" \
        --file-filter="BAM files|*.bam")
    
    if [[ -z "$OUTPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get parameters
    PARAMS=$(yad --form \
        --title="Mark Duplicates Parameters" \
        --field="Tool:CB" "pbmarkdup!samtools" \
        --field="Remove duplicates:CHK" "FALSE" \
        --field="Threads:NUM" "8!1..32!1" \
        --button="Cancel:1" --button="Run:0")
    
    if [[ $? -ne 0 ]]; then
        show_main_menu
        return
    fi
    
    TOOL=$(echo "$PARAMS" | cut -d'|' -f1)
    REMOVE=$(echo "$PARAMS" | cut -d'|' -f2)
    THREADS=$(echo "$PARAMS" | cut -d'|' -f3 | cut -d',' -f1)
    
    # Build command
    CMD=("${SCRIPT_DIR}/03_mark_duplicates.sh" \
        -i "$INPUT" \
        -o "$OUTPUT" \
        --tool "$TOOL" \
        -t "$THREADS")
    
    if [[ "$REMOVE" == "TRUE" ]]; then
        CMD+=(--remove)
    fi
    
    LOG_FILE="/tmp/markdup_$$.log"
    CMD+=(-l "$LOG_FILE")
    
    # Run mark duplicates
    "${CMD[@]}" 2>&1 | tee "$LOG_FILE" | \
        yad --progress --pulsate --auto-close \
            --title="Marking Duplicates" \
            --text="Processing: $(basename "$INPUT")"
    
    RESULT=$?
    
    if [[ $RESULT -eq 0 ]]; then
        yad --info --title="Success" \
            --text="Duplicate marking completed!\n\nOutput: $OUTPUT" \
            --button="View Log:bash -c 'xdg-open $LOG_FILE'" \
            --button="OK:0"
    else
        yad --error --title="Error" \
            --text="Duplicate marking failed!\n\nCheck log: $LOG_FILE"
    fi
    
    show_main_menu
}

# ==============================================================================
# STEP 4: DeepVariant
# ==============================================================================

run_deepvariant() {
    # Get input file
    INPUT=$(yad --file \
        --title="Select Input BAM" \
        --file-filter="BAM files|*.bam")
    
    if [[ -z "$INPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get reference
    REFERENCE=$(yad --file \
        --title="Select Reference Genome" \
        --filename="${REFERENCE_GENOME}" \
        --file-filter="FASTA files|*.fasta *.fa *.fna")
    
    if [[ -z "$REFERENCE" ]]; then
        show_main_menu
        return
    fi
    
    # Get output file
    DEFAULT_OUTPUT="${INPUT%.bam}.vcf.gz"
    OUTPUT=$(yad --file --save \
        --title="Save Output VCF" \
        --filename="$DEFAULT_OUTPUT" \
        --file-filter="VCF files|*.vcf.gz")
    
    if [[ -z "$OUTPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get parameters
    PARAMS=$(yad --form \
        --title="DeepVariant Parameters" \
        --field="Model:CB" "PACBIO!WGS!WES!HYBRID_PACBIO_ILLUMINA" \
        --field="Threads:NUM" "16!1..64!1" \
        --field="Use Docker:CHK" "FALSE" \
        --field="Regions BED (optional):FL" "" \
        --button="Cancel:1" --button="Run:0")
    
    if [[ $? -ne 0 ]]; then
        show_main_menu
        return
    fi
    
    MODEL=$(echo "$PARAMS" | cut -d'|' -f1)
    THREADS=$(echo "$PARAMS" | cut -d'|' -f2 | cut -d',' -f1)
    USE_DOCKER=$(echo "$PARAMS" | cut -d'|' -f3)
    REGIONS=$(echo "$PARAMS" | cut -d'|' -f4)
    
    # Build command
    CMD=("${SCRIPT_DIR}/04_call_variants_deepvariant.sh" \
        -i "$INPUT" \
        -r "$REFERENCE" \
        -o "$OUTPUT" \
        -m "$MODEL" \
        -t "$THREADS")
    
    if [[ "$USE_DOCKER" == "TRUE" ]]; then
        CMD+=(--docker)
    fi
    
    if [[ -n "$REGIONS" ]]; then
        CMD+=(--regions "$REGIONS")
    fi
    
    LOG_FILE="/tmp/deepvariant_$$.log"
    CMD+=(-l "$LOG_FILE")
    
    # Run DeepVariant
    "${CMD[@]}" 2>&1 | tee "$LOG_FILE" | \
        yad --progress --pulsate --auto-close \
            --title="Calling Variants (DeepVariant)" \
            --text="Processing: $(basename "$INPUT")\n\nThis may take a while..."
    
    RESULT=$?
    
    if [[ $RESULT -eq 0 ]]; then
        yad --info --title="Success" \
            --text="Variant calling completed!\n\nOutput: $OUTPUT" \
            --button="View Log:bash -c 'xdg-open $LOG_FILE'" \
            --button="OK:0"
    else
        yad --error --title="Error" \
            --text="Variant calling failed!\n\nCheck log: $LOG_FILE"
    fi
    
    show_main_menu
}

# ==============================================================================
# STEP 5: pbsv
# ==============================================================================

run_pbsv() {
    # Get input file
    INPUT=$(yad --file \
        --title="Select Input BAM" \
        --file-filter="BAM files|*.bam")
    
    if [[ -z "$INPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get reference
    REFERENCE=$(yad --file \
        --title="Select Reference Genome" \
        --filename="${REFERENCE_GENOME}" \
        --file-filter="FASTA files|*.fasta *.fa *.fna")
    
    if [[ -z "$REFERENCE" ]]; then
        show_main_menu
        return
    fi
    
    # Get output file
    DEFAULT_OUTPUT="${INPUT%.bam}.sv.vcf"
    OUTPUT=$(yad --file --save \
        --title="Save Output VCF" \
        --filename="$DEFAULT_OUTPUT" \
        --file-filter="VCF files|*.vcf")
    
    if [[ -z "$OUTPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get parameters
    PARAMS=$(yad --form \
        --title="pbsv Parameters" \
        --field="Min SV Length (bp):NUM" "50!10..1000!10" \
        --field="Max SV Length (bp):NUM" "100000!1000..1000000!1000" \
        --field="Threads:NUM" "16!1..64!1" \
        --field="CCS/HiFi data:CHK" "TRUE" \
        --button="Cancel:1" --button="Run:0")
    
    if [[ $? -ne 0 ]]; then
        show_main_menu
        return
    fi
    
    MIN_LEN=$(echo "$PARAMS" | cut -d'|' -f1 | cut -d',' -f1)
    MAX_LEN=$(echo "$PARAMS" | cut -d'|' -f2 | cut -d',' -f1)
    THREADS=$(echo "$PARAMS" | cut -d'|' -f3 | cut -d',' -f1)
    IS_CCS=$(echo "$PARAMS" | cut -d'|' -f4)
    
    # Build command
    CMD=("${SCRIPT_DIR}/05_call_sv_pbsv.sh" \
        -i "$INPUT" \
        -r "$REFERENCE" \
        -o "$OUTPUT" \
        --min-length "$MIN_LEN" \
        --max-length "$MAX_LEN" \
        -t "$THREADS")
    
    if [[ "$IS_CCS" == "TRUE" ]]; then
        CMD+=(--ccs)
    fi
    
    LOG_FILE="/tmp/pbsv_$$.log"
    CMD+=(-l "$LOG_FILE")
    
    # Run pbsv
    "${CMD[@]}" 2>&1 | tee "$LOG_FILE" | \
        yad --progress --pulsate --auto-close \
            --title="Calling Structural Variants (pbsv)" \
            --text="Processing: $(basename "$INPUT")"
    
    RESULT=$?
    
    if [[ $RESULT -eq 0 ]]; then
        yad --info --title="Success" \
            --text="SV calling completed!\n\nOutput: $OUTPUT" \
            --button="View Log:bash -c 'xdg-open $LOG_FILE'" \
            --button="OK:0"
    else
        yad --error --title="Error" \
            --text="SV calling failed!\n\nCheck log: $LOG_FILE"
    fi
    
    show_main_menu
}

# ==============================================================================
# STEP 6: QC
# ==============================================================================

run_qc() {
    # Get input file
    INPUT=$(yad --file \
        --title="Select Input BAM" \
        --file-filter="BAM files|*.bam")
    
    if [[ -z "$INPUT" ]]; then
        show_main_menu
        return
    fi
    
    # Get reference
    REFERENCE=$(yad --file \
        --title="Select Reference Genome" \
        --filename="${REFERENCE_GENOME}" \
        --file-filter="FASTA files|*.fasta *.fa *.fna")
    
    if [[ -z "$REFERENCE" ]]; then
        show_main_menu
        return
    fi
    
    # Get output directory
    DEFAULT_DIR="$(dirname "$INPUT")/qc_$(basename "${INPUT%.bam}")"
    OUTPUT_DIR=$(yad --file --directory \
        --title="Select Output Directory" \
        --filename="$DEFAULT_DIR")
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        show_main_menu
        return
    fi
    
    # Get parameters
    PARAMS=$(yad --form \
        --title="QC Parameters" \
        --field="Threads:NUM" "8!1..32!1" \
        --field="Min MAPQ:NUM" "20!0..60!1" \
        --field="Min Coverage:NUM" "10!1..100!1" \
        --field="Regions BED (optional):FL" "" \
        --button="Cancel:1" --button="Run:0")
    
    if [[ $? -ne 0 ]]; then
        show_main_menu
        return
    fi
    
    THREADS=$(echo "$PARAMS" | cut -d'|' -f1 | cut -d',' -f1)
    MIN_MAPQ=$(echo "$PARAMS" | cut -d'|' -f2 | cut -d',' -f1)
    MIN_COV=$(echo "$PARAMS" | cut -d'|' -f3 | cut -d',' -f1)
    BED=$(echo "$PARAMS" | cut -d'|' -f4)
    
    # Build command
    CMD=("${SCRIPT_DIR}/06_qc_coverage.sh" \
        -i "$INPUT" \
        -r "$REFERENCE" \
        -o "$OUTPUT_DIR" \
        -t "$THREADS" \
        --min-mapq "$MIN_MAPQ" \
        --min-coverage "$MIN_COV")
    
    if [[ -n "$BED" ]]; then
        CMD+=(-b "$BED")
    fi
    
    LOG_FILE="/tmp/qc_$$.log"
    CMD+=(-l "$LOG_FILE")
    
    # Run QC
    "${CMD[@]}" 2>&1 | tee "$LOG_FILE" | \
        yad --progress --pulsate --auto-close \
            --title="Running QC Analysis" \
            --text="Processing: $(basename "$INPUT")"
    
    RESULT=$?
    
    if [[ $RESULT -eq 0 ]]; then
        REPORT="${OUTPUT_DIR}/$(basename "${INPUT%.bam}")_qc_report.html"
        yad --info --title="Success" \
            --text="QC analysis completed!\n\nOutput directory: $OUTPUT_DIR" \
            --button="View Report:bash -c 'xdg-open $REPORT'" \
            --button="OK:0"
    else
        yad --error --title="Error" \
            --text="QC analysis failed!\n\nCheck log: $LOG_FILE"
    fi
    
    show_main_menu
}

# ==============================================================================
# STEP 7: Complete Pipeline
# ==============================================================================

run_complete_pipeline() {
    yad --info --title="Complete Pipeline" \
        --text="This feature will run the entire pipeline sequentially.\n\nNot yet implemented in this example."
    
    show_main_menu
}

# ==============================================================================
# START APPLICATION
# ==============================================================================

show_main_menu
