#!/bin/bash
# Utility functions for PacBio pipeline scripts
# Source this file in other scripts: source utils.sh

# ==============================================================================
# LOGGING FUNCTIONS
# ==============================================================================

# Colors for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Log levels
LOG_LEVEL=${LOG_LEVEL:-"INFO"}

# Timestamp function
timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

# Logging function
log() {
    local level=$1
    shift
    local message="$@"
    
    case $level in
        DEBUG)
            if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
                echo -e "${BLUE}[$(timestamp)] [DEBUG]${NC} $message"
            fi
            ;;
        INFO)
            if [[ "$LOG_LEVEL" == "DEBUG" || "$LOG_LEVEL" == "INFO" ]]; then
                echo -e "${GREEN}[$(timestamp)] [INFO]${NC} $message"
            fi
            ;;
        WARNING)
            if [[ "$LOG_LEVEL" != "ERROR" ]]; then
                echo -e "${YELLOW}[$(timestamp)] [WARNING]${NC} $message"
            fi
            ;;
        ERROR)
            echo -e "${RED}[$(timestamp)] [ERROR]${NC} $message" >&2
            ;;
    esac
    
    # Also write to log file if LOG_FILE is set
    if [[ -n "$LOG_FILE" ]]; then
        echo "[$(timestamp)] [$level] $message" >> "$LOG_FILE"
    fi
}

# Convenience functions
log_debug() { log DEBUG "$@"; }
log_info() { log INFO "$@"; }
log_warning() { log WARNING "$@"; }
log_error() { log ERROR "$@"; }

# ==============================================================================
# ERROR HANDLING
# ==============================================================================

# Exit with error message
die() {
    log_error "$@"
    exit 1
}

# Check if command exists
check_command() {
    local cmd=$1
    if ! command -v "$cmd" &> /dev/null; then
        die "Required command not found: $cmd"
    fi
}

# Check if file exists
check_file() {
    local file=$1
    if [[ ! -f "$file" ]]; then
        die "Required file not found: $file"
    fi
}

# Check if directory exists
check_dir() {
    local dir=$1
    if [[ ! -d "$dir" ]]; then
        die "Required directory not found: $dir"
    fi
}

# Create directory if it doesn't exist
ensure_dir() {
    local dir=$1
    if [[ ! -d "$dir" ]]; then
        mkdir -p "$dir" || die "Failed to create directory: $dir"
        log_debug "Created directory: $dir"
    fi
}

# ==============================================================================
# RESOURCE CHECKING
# ==============================================================================

# Check available memory (in GB)
check_memory() {
    local required_gb=$1
    local available_gb=$(free -g | awk '/^Mem:/{print $7}')
    
    if (( available_gb < required_gb )); then
        log_warning "Available memory (${available_gb}GB) is less than required (${required_gb}GB)"
        return 1
    fi
    return 0
}

# Check available disk space (in GB)
check_disk_space() {
    local directory=$1
    local required_gb=$2
    local available_gb=$(df -BG "$directory" | awk 'NR==2 {print $4}' | sed 's/G//')
    
    if (( available_gb < required_gb )); then
        log_warning "Available disk space (${available_gb}GB) is less than required (${required_gb}GB)"
        return 1
    fi
    return 0
}

# ==============================================================================
# FILE OPERATIONS
# ==============================================================================

# Get file size in human-readable format
get_file_size() {
    local file=$1
    du -h "$file" | cut -f1
}

# Count lines in file (handles compressed files)
count_lines() {
    local file=$1
    if [[ "$file" == *.gz ]]; then
        zcat "$file" | wc -l
    else
        wc -l < "$file"
    fi
}

# ==============================================================================
# PROGRESS TRACKING
# ==============================================================================

# Initialize progress tracking
init_progress() {
    local total_steps=$1
    echo "0/$total_steps" > "$PROGRESS_FILE"
}

# Update progress
update_progress() {
    local current_step=$1
    local total_steps=$2
    local description=$3
    
    echo "$current_step/$total_steps" > "$PROGRESS_FILE"
    
    local percent=$(( current_step * 100 / total_steps ))
    log_info "Progress: $percent% - $description"
    
    # Output for YAD progress bar
    echo "$percent"
    echo "# $description"
}

# ==============================================================================
# VALIDATION FUNCTIONS
# ==============================================================================

# Validate BAM file
validate_bam() {
    local bam=$1
    log_info "Validating BAM file: $bam"
    
    if ! samtools quickcheck "$bam" 2>/dev/null; then
        log_error "Invalid BAM file: $bam"
        return 1
    fi
    
    log_debug "BAM file is valid"
    return 0
}

# Validate FASTQ file
validate_fastq() {
    local fastq=$1
    log_info "Validating FASTQ file: $fastq"
    
    # Check if file can be read and has proper format
    local first_char
    if [[ "$fastq" == *.gz ]]; then
        first_char=$(zcat "$fastq" | head -c 1)
    else
        first_char=$(head -c 1 "$fastq")
    fi
    
    if [[ "$first_char" != "@" ]]; then
        log_error "Invalid FASTQ file (doesn't start with @): $fastq"
        return 1
    fi
    
    log_debug "FASTQ file appears valid"
    return 0
}

# Validate reference genome
validate_reference() {
    local ref=$1
    check_file "$ref"
    
    # Check for index
    if [[ ! -f "${ref}.fai" ]]; then
        log_warning "Reference index not found. Creating index..."
        samtools faidx "$ref" || die "Failed to index reference"
    fi
    
    log_debug "Reference genome is valid"
    return 0
}

# ==============================================================================
# CLEANUP FUNCTIONS
# ==============================================================================

# Cleanup function for trap
cleanup() {
    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        log_error "Script failed with exit code: $exit_code"
    fi
    
    # Remove temp files if specified
    if [[ "$KEEP_INTERMEDIATE" == "false" && -n "$TEMP_FILES" ]]; then
        log_info "Cleaning up temporary files..."
        rm -f $TEMP_FILES
    fi
}

# ==============================================================================
# TIMING FUNCTIONS
# ==============================================================================

# Start timer
start_timer() {
    TIMER_START=$(date +%s)
}

# Stop timer and report elapsed time
stop_timer() {
    local end=$(date +%s)
    local elapsed=$((end - TIMER_START))
    local hours=$((elapsed / 3600))
    local minutes=$(((elapsed % 3600) / 60))
    local seconds=$((elapsed % 60))
    
    log_info "Elapsed time: ${hours}h ${minutes}m ${seconds}s"
}

# ==============================================================================
# YAD OUTPUT FUNCTIONS
# ==============================================================================

# Output for YAD list dialog
yad_list_output() {
    local item=$1
    local value=$2
    echo "$item|$value"
}

# Output for YAD form
yad_form_output() {
    local field=$1
    local value=$2
    echo "$field:$value"
}

# Success message for YAD
yad_success() {
    local message=$1
    echo "SUCCESS|$message"
}

# Error message for YAD
yad_error() {
    local message=$1
    echo "ERROR|$message"
}
