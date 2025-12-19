#!/bin/bash
# Script: 06_qc_coverage.sh
# Purpose: Generate QC metrics and coverage statistics for aligned BAM files
# Usage: ./06_qc_coverage.sh -i input.bam -r reference.fasta -o output_dir [options]

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
Usage: $0 -i INPUT_BAM -r REFERENCE -o OUTPUT_DIR [OPTIONS]

Required arguments:
    -i, --input         Input BAM file (must be sorted and indexed)
    -r, --reference     Reference genome FASTA file
    -o, --output        Output directory for QC reports

Optional arguments:
    -t, --threads       Number of threads (default: 8)
    -b, --bed           BED file with regions of interest (optional)
    --min-mapq          Minimum mapping quality (default: $MIN_MAPPING_QUALITY)
    --min-coverage      Minimum coverage threshold (default: $MIN_COVERAGE)
    -l, --log           Log file path
    -h, --help          Show this help message

Description:
    Generates comprehensive QC metrics including:
        - Alignment statistics (samtools stats, flagstat, idxstats)
        - Coverage depth analysis (per-base and per-region)
        - Read length distribution
        - Mapping quality distribution
        - Insert size distribution (if paired-end)
        - Coverage uniformity metrics
    
    Output includes:
        - Text-based statistics files
        - Coverage plots (if plotting tools available)
        - HTML summary report

Example:
    $0 -i sample.bam -r hg38.fasta -o qc_output -t 16

EOF
    exit 1
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

INPUT_BAM=""
REFERENCE="${REFERENCE_GENOME:-}"
OUTPUT_DIR=""
THREADS=8
BED_FILE=""
MIN_MAPQ="${MIN_MAPPING_QUALITY:-20}"
MIN_COV="${MIN_COVERAGE:-10}"
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
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -b|--bed)
            BED_FILE="$2"
            shift 2
            ;;
        --min-mapq)
            MIN_MAPQ="$2"
            shift 2
            ;;
        --min-coverage)
            MIN_COV="$2"
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
[[ -z "$OUTPUT_DIR" ]] && die "Output directory not specified (-o)"

# Check if files exist
check_file "$INPUT_BAM"
check_file "$REFERENCE"

# Check for BAM index
if [[ ! -f "${INPUT_BAM}.bai" ]]; then
    die "BAM index not found: ${INPUT_BAM}.bai"
fi

# Check if BED file exists (if specified)
if [[ -n "$BED_FILE" ]]; then
    check_file "$BED_FILE"
fi

# Check for samtools
check_command "samtools"

# Create output directory
ensure_dir "$OUTPUT_DIR"

# Set up logging
if [[ -z "$LOG_FILE" ]]; then
    LOG_FILE="${OUTPUT_DIR}/qc.log"
fi

# Set up trap for cleanup
TEMP_FILES=""
trap cleanup EXIT

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

SAMPLE_NAME=$(basename "$INPUT_BAM" .bam)

log_info "========================================="
log_info "Quality Control and Coverage Analysis"
log_info "========================================="
log_info "Input BAM: $INPUT_BAM"
log_info "Reference: $REFERENCE"
log_info "Output directory: $OUTPUT_DIR"
log_info "Sample name: $SAMPLE_NAME"
log_info "Threads: $THREADS"
log_info "Min MAPQ: $MIN_MAPQ"
log_info "Min coverage: $MIN_COV"
if [[ -n "$BED_FILE" ]]; then
    log_info "Regions BED: $BED_FILE"
fi
log_info "Log file: $LOG_FILE"
log_info "========================================="

start_timer

# Validate inputs
log_info "Validating inputs..."
validate_bam "$INPUT_BAM" || die "BAM validation failed"
validate_reference "$REFERENCE" || die "Reference validation failed"

# ==============================================================================
# ALIGNMENT STATISTICS
# ==============================================================================

log_info "Generating alignment statistics..."

# samtools stats
STATS_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.stats.txt"
log_info "Running samtools stats..."
samtools stats -@ "$THREADS" "$INPUT_BAM" > "$STATS_FILE" 2>&1

# samtools flagstat
FLAGSTAT_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.flagstat.txt"
log_info "Running samtools flagstat..."
samtools flagstat -@ "$THREADS" "$INPUT_BAM" > "$FLAGSTAT_FILE" 2>&1

# samtools idxstats
IDXSTATS_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.idxstats.txt"
log_info "Running samtools idxstats..."
samtools idxstats "$INPUT_BAM" > "$IDXSTATS_FILE" 2>&1

# Extract key metrics
log_info "Extracting key metrics..."

TOTAL_READS=$(grep "^SN" "$STATS_FILE" | grep "sequences:" | cut -f3)
MAPPED_READS=$(grep "^SN" "$STATS_FILE" | grep "reads mapped:" | cut -f3)
MAPPED_BASES=$(grep "^SN" "$STATS_FILE" | grep "bases mapped:" | cut -f3)
AVG_LENGTH=$(grep "^SN" "$STATS_FILE" | grep "average length:" | cut -f3)
AVG_QUALITY=$(grep "^SN" "$STATS_FILE" | grep "average quality:" | cut -f3)
ERROR_RATE=$(grep "^SN" "$STATS_FILE" | grep "error rate:" | cut -f3)

if [[ -n "$TOTAL_READS" && -n "$MAPPED_READS" ]]; then
    MAPPING_RATE=$(echo "scale=2; $MAPPED_READS * 100 / $TOTAL_READS" | bc)
else
    MAPPING_RATE="N/A"
fi

log_info "Basic Statistics:"
log_info "  Total reads: $TOTAL_READS"
log_info "  Mapped reads: $MAPPED_READS"
log_info "  Mapping rate: ${MAPPING_RATE}%"
log_info "  Mapped bases: $MAPPED_BASES"
log_info "  Average read length: $AVG_LENGTH"
log_info "  Average quality: $AVG_QUALITY"
log_info "  Error rate: $ERROR_RATE"

# ==============================================================================
# COVERAGE ANALYSIS
# ==============================================================================

log_info "Analyzing coverage depth..."

# Generate per-base coverage
COVERAGE_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.coverage.txt"
log_info "Calculating per-base coverage..."

if [[ -n "$BED_FILE" ]]; then
    samtools depth -@ "$THREADS" -q "$MIN_MAPQ" -b "$BED_FILE" "$INPUT_BAM" > "$COVERAGE_FILE" 2>&1
else
    samtools depth -@ "$THREADS" -q "$MIN_MAPQ" "$INPUT_BAM" > "$COVERAGE_FILE" 2>&1
fi

# Calculate coverage statistics
log_info "Calculating coverage statistics..."
COVERAGE_STATS="${OUTPUT_DIR}/${SAMPLE_NAME}.coverage_stats.txt"

awk -v min_cov="$MIN_COV" '
BEGIN {
    total_bases = 0;
    covered_bases = 0;
    sum_coverage = 0;
    min = 999999999;
    max = 0;
    bases_above_min = 0;
}
{
    total_bases++;
    cov = $3;
    
    if (cov > 0) {
        covered_bases++;
    }
    
    sum_coverage += cov;
    
    if (cov < min) min = cov;
    if (cov > max) max = cov;
    
    if (cov >= min_cov) {
        bases_above_min++;
    }
    
    # Store for percentile calculation
    coverage[total_bases] = cov;
}
END {
    if (total_bases > 0) {
        mean = sum_coverage / total_bases;
        pct_covered = (covered_bases / total_bases) * 100;
        pct_above_min = (bases_above_min / total_bases) * 100;
        
        # Calculate median
        asort(coverage);
        if (total_bases % 2 == 0) {
            median = (coverage[total_bases/2] + coverage[total_bases/2 + 1]) / 2;
        } else {
            median = coverage[int(total_bases/2) + 1];
        }
        
        print "Total bases analyzed: " total_bases;
        print "Covered bases (depth > 0): " covered_bases;
        print "Percent covered: " sprintf("%.2f%%", pct_covered);
        print "Mean coverage: " sprintf("%.2f", mean);
        print "Median coverage: " median;
        print "Min coverage: " min;
        print "Max coverage: " max;
        print "Bases with coverage >= " min_cov ": " bases_above_min;
        print "Percent >= " min_cov "x: " sprintf("%.2f%%", pct_above_min);
    }
}' "$COVERAGE_FILE" > "$COVERAGE_STATS"

log_info "Coverage Statistics:"
cat "$COVERAGE_STATS" | while read line; do
    log_info "  $line"
done

# Calculate coverage per chromosome
log_info "Calculating per-chromosome coverage..."
CHROM_COVERAGE="${OUTPUT_DIR}/${SAMPLE_NAME}.chromosome_coverage.txt"

awk '{chr[$1] += $3; count[$1]++} 
     END {
         print "Chromosome\tMean_Coverage\tBases";
         for (c in chr) {
             print c "\t" chr[c]/count[c] "\t" count[c];
         }
     }' "$COVERAGE_FILE" | sort -k1,1V > "$CHROM_COVERAGE"

log_info "Per-chromosome coverage written to: $CHROM_COVERAGE"

# ==============================================================================
# READ LENGTH DISTRIBUTION
# ==============================================================================

log_info "Analyzing read length distribution..."
LENGTH_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.read_lengths.txt"

samtools view -@ "$THREADS" "$INPUT_BAM" | awk '{print length($10)}' | \
    sort -n | uniq -c | awk '{print $2 "\t" $1}' > "$LENGTH_FILE"

# Calculate length statistics
LENGTH_STATS=$(awk '
    {len = $1; count = $2; sum += len * count; total += count; lengths[len] = count}
    END {
        mean = sum / total;
        
        # Find median
        cumsum = 0;
        median = 0;
        for (l in lengths) {
            cumsum += lengths[l];
            if (cumsum >= total/2 && median == 0) {
                median = l;
                break;
            }
        }
        
        print "Mean: " int(mean) "bp";
        print "Median: " median "bp";
        print "Total reads: " total;
    }' "$LENGTH_FILE")

log_info "Read Length Distribution:"
echo "$LENGTH_STATS" | while read line; do
    log_info "  $line"
done

# ==============================================================================
# MAPPING QUALITY DISTRIBUTION
# ==============================================================================

log_info "Analyzing mapping quality distribution..."
MAPQ_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.mapq_distribution.txt"

samtools view -@ "$THREADS" "$INPUT_BAM" | awk '{print $5}' | \
    sort -n | uniq -c | awk '{print $2 "\t" $1}' > "$MAPQ_FILE"

HIGH_MAPQ=$(awk -v min="$MIN_MAPQ" '$1 >= min {sum += $2} END {print sum}' "$MAPQ_FILE")
TOTAL_MAPQ=$(awk '{sum += $2} END {print sum}' "$MAPQ_FILE")

if [[ -n "$TOTAL_MAPQ" && "$TOTAL_MAPQ" -gt 0 ]]; then
    PCT_HIGH_MAPQ=$(echo "scale=2; $HIGH_MAPQ * 100 / $TOTAL_MAPQ" | bc)
    log_info "Reads with MAPQ >= $MIN_MAPQ: $HIGH_MAPQ (${PCT_HIGH_MAPQ}%)"
fi

# ==============================================================================
# GENERATE HTML REPORT
# ==============================================================================

log_info "Generating HTML summary report..."
HTML_REPORT="${OUTPUT_DIR}/${SAMPLE_NAME}_qc_report.html"

cat > "$HTML_REPORT" << EOF
<!DOCTYPE html>
<html>
<head>
    <title>QC Report - $SAMPLE_NAME</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; margin-top: 30px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { text-align: left; padding: 12px; border-bottom: 1px solid #ddd; }
        th { background-color: #3498db; color: white; }
        tr:hover { background-color: #f5f5f5; }
        .metric { font-weight: bold; color: #2980b9; }
        .warning { color: #e74c3c; }
        .good { color: #27ae60; }
    </style>
</head>
<body>
    <h1>Quality Control Report</h1>
    <p><strong>Sample:</strong> $SAMPLE_NAME</p>
    <p><strong>Date:</strong> $(date '+%Y-%m-%d %H:%M:%S')</p>
    
    <h2>Alignment Statistics</h2>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
        <tr><td>Total Reads</td><td>$TOTAL_READS</td></tr>
        <tr><td>Mapped Reads</td><td>$MAPPED_READS</td></tr>
        <tr><td>Mapping Rate</td><td class="$([ ${MAPPING_RATE%.*} -ge 90 ] 2>/dev/null && echo "good" || echo "warning")">${MAPPING_RATE}%</td></tr>
        <tr><td>Mapped Bases</td><td>$MAPPED_BASES</td></tr>
        <tr><td>Average Read Length</td><td>$AVG_LENGTH bp</td></tr>
        <tr><td>Average Quality</td><td>$AVG_QUALITY</td></tr>
        <tr><td>Error Rate</td><td>$ERROR_RATE</td></tr>
    </table>
    
    <h2>Coverage Statistics</h2>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
EOF

# Add coverage stats to HTML
while IFS= read -r line; do
    metric=$(echo "$line" | cut -d: -f1)
    value=$(echo "$line" | cut -d: -f2- | sed 's/^ //')
    echo "        <tr><td>$metric</td><td>$value</td></tr>" >> "$HTML_REPORT"
done < "$COVERAGE_STATS"

cat >> "$HTML_REPORT" << EOF
    </table>
    
    <h2>Output Files</h2>
    <ul>
        <li><a href="${SAMPLE_NAME}.stats.txt">Samtools Stats</a></li>
        <li><a href="${SAMPLE_NAME}.flagstat.txt">Flagstat</a></li>
        <li><a href="${SAMPLE_NAME}.coverage.txt">Per-base Coverage</a></li>
        <li><a href="${SAMPLE_NAME}.coverage_stats.txt">Coverage Statistics</a></li>
        <li><a href="${SAMPLE_NAME}.chromosome_coverage.txt">Per-chromosome Coverage</a></li>
        <li><a href="${SAMPLE_NAME}.read_lengths.txt">Read Length Distribution</a></li>
        <li><a href="${SAMPLE_NAME}.mapq_distribution.txt">MAPQ Distribution</a></li>
    </ul>
</body>
</html>
EOF

log_info "HTML report created: $HTML_REPORT"

stop_timer

log_info "========================================="
log_info "QC analysis completed successfully!"
log_info "========================================="
log_info "Output files in: $OUTPUT_DIR"
log_info "HTML report: $HTML_REPORT"
log_info "========================================="

# YAD-friendly output
yad_success "QC analysis completed: $HTML_REPORT"

exit 0
