#!/bin/bash
# PacBio Pipeline Configuration File
# Edit this file to set default paths and parameters

# ==============================================================================
# REFERENCE GENOME SETTINGS
# ==============================================================================
REFERENCE_GENOME="/path/to/reference/genome.fasta"
REFERENCE_INDEX="/path/to/reference/genome.fasta.fai"

# ==============================================================================
# TOOL PATHS (Leave empty to use system PATH)
# ==============================================================================
PBMM2=""              # Leave empty if in PATH
MINIMAP2=""           # Leave empty if in PATH
SAMTOOLS=""           # Leave empty if in PATH
PBMARKDUP=""          # Leave empty if in PATH
DEEPVARIANT_DIR=""    # Directory containing DeepVariant (if not using Docker)
BCFTOOLS=""           # Leave empty if in PATH
PBSV=""               # Leave empty if in PATH
SNIFFLES2=""          # Leave empty if in PATH

# ==============================================================================
# DEFAULT PARAMETERS
# ==============================================================================

# Alignment parameters
ALIGN_THREADS=16
ALIGN_PRESET="HiFi"      # Options: HiFi, SUBREAD, CCS, ISOSEQ
ALIGN_SORT_MEMORY="4G"   # Memory per thread for sorting

# Variant calling parameters
DEEPVARIANT_MODEL="PACBIO"
DEEPVARIANT_THREADS=16
MIN_MAPPING_QUALITY=20
MIN_BASE_QUALITY=20

# SV calling parameters
SV_MIN_LENGTH=50
SV_MAX_LENGTH=100000
PBSV_THREADS=16
SNIFFLES_THREADS=16

# Quality control parameters
MIN_READ_LENGTH=1000
MIN_READ_QUALITY=20
MIN_COVERAGE=10

# ==============================================================================
# OUTPUT SETTINGS
# ==============================================================================
LOG_LEVEL="INFO"         # Options: DEBUG, INFO, WARNING, ERROR
KEEP_INTERMEDIATE=true   # Keep intermediate files (true/false)
COMPRESS_OUTPUT=true     # Compress output files (true/false)

# ==============================================================================
# DOCKER/SINGULARITY SETTINGS (if using containers)
# ==============================================================================
USE_DOCKER=false
DEEPVARIANT_DOCKER="google/deepvariant:1.6.0"
SINGULARITY_CACHEDIR="/tmp/singularity_cache"

# ==============================================================================
# RESOURCE LIMITS
# ==============================================================================
MAX_MEMORY="256G"        # Maximum memory to use
MAX_THREADS=16           # Maximum threads to use
TEMP_DIR="/tmp"          # Temporary directory for intermediate files
