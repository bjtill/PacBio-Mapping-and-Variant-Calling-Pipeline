#!/bin/bash
# Script: check_dependencies.sh
# Purpose: Check if all required tools are installed and accessible
# Usage: ./check_dependencies.sh

set -euo pipefail

# Source utility functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/utils.sh"

# ==============================================================================
# TOOL CHECKS
# ==============================================================================

echo "========================================"
echo "PacBio Pipeline Dependency Checker"
echo "========================================"
echo ""

MISSING_TOOLS=()
OPTIONAL_MISSING=()

# Check required tools
echo "Checking required tools:"
echo ""

check_tool() {
    local tool=$1
    local required=$2
    local min_version=$3
    
    if command -v "$tool" &> /dev/null; then
        echo -e "${GREEN}✓${NC} $tool found: $(command -v "$tool")"
        
        # Try to get version
        case "$tool" in
            samtools|bcftools|tabix|bgzip)
                version=$("$tool" --version 2>&1 | head -1 || echo "version unknown")
                echo "    Version: $version"
                ;;
            pbmm2|pbsv|pbmarkdup)
                version=$("$tool" --version 2>&1 || echo "version unknown")
                echo "    Version: $version"
                ;;
            *)
                echo "    Version: (unable to determine)"
                ;;
        esac
    else
        if [[ "$required" == "true" ]]; then
            echo -e "${RED}✗${NC} $tool NOT FOUND"
            MISSING_TOOLS+=("$tool")
        else
            echo -e "${YELLOW}⚠${NC} $tool NOT FOUND (optional)"
            OPTIONAL_MISSING+=("$tool")
        fi
    fi
    echo ""
}

# Required tools
check_tool "samtools" "true" "1.17"
check_tool "pbmm2" "true" "1.0"
check_tool "pbsv" "true" "2.0"
check_tool "bcftools" "true" "1.17"

# Optional but recommended tools
check_tool "pbmarkdup" "false" "1.0"
check_tool "tabix" "false" "1.17"
check_tool "bgzip" "false" "1.17"

# Check for DeepVariant
echo "Checking DeepVariant:"
echo ""
if command -v run_deepvariant &> /dev/null; then
    echo -e "${GREEN}✓${NC} DeepVariant found (native installation)"
    echo "    Location: $(command -v run_deepvariant)"
elif command -v docker &> /dev/null; then
    echo -e "${GREEN}✓${NC} Docker found - can use DeepVariant Docker"
    echo "    Location: $(command -v docker)"
    
    # Check if DeepVariant image exists
    if docker images | grep -q "deepvariant"; then
        echo "    DeepVariant Docker image: found"
    else
        echo -e "    ${YELLOW}⚠${NC} DeepVariant Docker image: not found"
        echo "    Run: docker pull google/deepvariant:1.6.0"
    fi
elif command -v singularity &> /dev/null; then
    echo -e "${GREEN}✓${NC} Singularity found - can use DeepVariant Singularity"
    echo "    Location: $(command -v singularity)"
else
    echo -e "${RED}✗${NC} DeepVariant NOT FOUND"
    echo "    Install DeepVariant, Docker, or Singularity"
    MISSING_TOOLS+=("deepvariant/docker/singularity")
fi
echo ""

# ==============================================================================
# SYSTEM CHECKS
# ==============================================================================

echo "========================================"
echo "System Resources:"
echo "========================================"
echo ""

# Check CPU cores
CPU_CORES=$(nproc)
echo "CPU cores: $CPU_CORES"
if [[ $CPU_CORES -ge 16 ]]; then
    echo -e "  ${GREEN}✓${NC} Sufficient cores for pipeline"
else
    echo -e "  ${YELLOW}⚠${NC} Recommended: 16+ cores"
fi
echo ""

# Check memory
TOTAL_MEM_GB=$(free -g | awk '/^Mem:/{print $2}')
AVAIL_MEM_GB=$(free -g | awk '/^Mem:/{print $7}')
echo "Total memory: ${TOTAL_MEM_GB}GB"
echo "Available memory: ${AVAIL_MEM_GB}GB"
if [[ $TOTAL_MEM_GB -ge 128 ]]; then
    echo -e "  ${GREEN}✓${NC} Sufficient memory for pipeline"
else
    echo -e "  ${YELLOW}⚠${NC} Recommended: 128GB+ total memory"
fi
echo ""

# Check disk space
HOME_SPACE=$(df -BG "$HOME" | awk 'NR==2 {print $4}' | sed 's/G//')
echo "Available disk space in $HOME: ${HOME_SPACE}GB"
if [[ $HOME_SPACE -ge 100 ]]; then
    echo -e "  ${GREEN}✓${NC} Sufficient disk space"
else
    echo -e "  ${YELLOW}⚠${NC} Warning: Low disk space"
    echo "  Recommended: 500GB+ free space for analysis"
fi
echo ""

# ==============================================================================
# CONFIGURATION CHECK
# ==============================================================================

echo "========================================"
echo "Configuration:"
echo "========================================"
echo ""

CONFIG_FILE="${SCRIPT_DIR}/../configs/pipeline_config.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    echo -e "${GREEN}✓${NC} Configuration file found"
    echo "    Location: $CONFIG_FILE"
    
    # Source config
    source "$CONFIG_FILE"
    
    echo ""
    echo "Reference genome:"
    if [[ -n "$REFERENCE_GENOME" && -f "$REFERENCE_GENOME" ]]; then
        echo -e "  ${GREEN}✓${NC} $REFERENCE_GENOME"
        
        # Check for index
        if [[ -f "${REFERENCE_GENOME}.fai" ]]; then
            echo -e "  ${GREEN}✓${NC} Reference index found"
        else
            echo -e "  ${YELLOW}⚠${NC} Reference index not found"
            echo "    Run: samtools faidx $REFERENCE_GENOME"
        fi
    else
        echo -e "  ${YELLOW}⚠${NC} Not configured or file not found"
        echo "    Edit $CONFIG_FILE to set REFERENCE_GENOME"
    fi
else
    echo -e "${YELLOW}⚠${NC} Configuration file not found"
    echo "    Expected: $CONFIG_FILE"
fi
echo ""

# ==============================================================================
# SUMMARY
# ==============================================================================

echo "========================================"
echo "Summary:"
echo "========================================"
echo ""

if [[ ${#MISSING_TOOLS[@]} -eq 0 ]]; then
    echo -e "${GREEN}✓${NC} All required tools are installed!"
    echo ""
    echo "You can now use the pipeline:"
    echo "  ${SCRIPT_DIR}/01_ubam_to_fastq.sh --help"
    echo "  ${SCRIPT_DIR}/02_align_pbmm2.sh --help"
    echo "  ..."
    echo ""
    echo "Or use the YAD GUI wrapper:"
    echo "  ${SCRIPT_DIR}/yad_wrapper_example.sh"
    EXIT_CODE=0
else
    echo -e "${RED}✗${NC} Missing required tools:"
    for tool in "${MISSING_TOOLS[@]}"; do
        echo "  - $tool"
    done
    echo ""
    echo "Installation instructions:"
    echo ""
    echo "Using conda (recommended):"
    echo "  conda create -n pacbio python=3.10"
    echo "  conda activate pacbio"
    echo "  conda install -c bioconda pbmm2 pbsv samtools bcftools"
    echo ""
    echo "For DeepVariant:"
    echo "  # Option 1: Docker (easiest)"
    echo "  sudo apt install docker.io"
    echo "  docker pull google/deepvariant:1.6.0"
    echo ""
    echo "  # Option 2: Conda"
    echo "  conda install -c bioconda deepvariant"
    EXIT_CODE=1
fi

if [[ ${#OPTIONAL_MISSING[@]} -gt 0 ]]; then
    echo ""
    echo "Optional tools not found (pipeline will still work):"
    for tool in "${OPTIONAL_MISSING[@]}"; do
        echo "  - $tool"
    done
fi

echo ""
echo "========================================"

exit $EXIT_CODE
