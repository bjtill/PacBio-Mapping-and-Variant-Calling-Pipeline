# PacBio-Mapping-and-Variant-Calling-Pipeline
A modular pipeline for PacBio HiFi data analysis.
_____________________________________________________________________________________________________________________________________________________________________________

## OVERVIEW  

A PacBio mapping and variant calling tool designed for your Ubuntu workstations with 16+ cores and 256GB+ RAM. The pipeline is modular, and designed for easy integration with YAD GUI applications. 

### Core Pipeline Scripts (6 main steps)

1. **01_ubam_to_fastq.sh**
   - Converts unmapped BAM files to FASTQ format
   - Handles both single-end and paired-end reads
   - Validates input/output and provides detailed statistics

2. **02_align_pbmm2.sh**
   - Aligns PacBio reads using pbmm2 (minimap2 wrapper)
   - Supports multiple presets (HiFi, SUBREAD, CCS, ISOSEQ)
   - Produces sorted, indexed BAM with alignment statistics

3. **03_mark_duplicates.sh**
   - Marks PCR and optical duplicates
   - Supports pbmarkdup (PacBio-specific) or samtools
   - Option to remove duplicates if needed

4. **04_call_variants_deepvariant.sh**
   - Calls small variants (SNPs/Indels) using DeepVariant
   - Supports native, Docker, or Singularity execution
   - Outputs VCF, gVCF, and statistics

5. **05_call_sv_pbsv.sh**
   - Calls structural variants using pbsv
   - Detects DEL, INS, INV, DUP, and BND
   - Two-stage process: discover signatures, then call/genotype

6. **06_qc_coverage.sh**
   - Comprehensive quality control and coverage analysis
   - Generates alignment stats, coverage metrics, read distributions
   - Creates HTML summary report

### Supporting Infrastructure

**utils.sh**
- Shared utility functions library
- Logging (with color-coded output)
- Error handling and validation
- Progress tracking
- Resource checking (memory, disk space)
- File operations and validation
- YAD-compatible output functions

**pipeline_config.sh**
- Central configuration file
- Reference genome paths
- Tool paths and parameters
- Resource limits
- Docker/Singularity settings

**check_dependencies.sh**
- Validates all required tools are installed
- Checks system resources
- Verifies configuration
- Provides installation instructions for missing tools

**yad_wrapper_example.sh**
- Complete example GUI wrapper
- Demonstrates all pipeline steps
- File selection dialogs
- Parameter configuration forms
- Progress indicators and result notifications

### Documentation

**QUICKSTART.md** (7.3KB)
- Quick installation guide
- Basic usage examples
- Common workflows
- Performance tips
- Troubleshooting basics

## File Organization

```
pacbio_pipeline/
├── QUICKSTART.md              # Quick start guide
├── configs/
│   └── pipeline_config.sh     # Central configuration
├── docs/
│   └── README.md              # Comprehensive documentation
└── scripts/
    ├── utils.sh               # Shared utilities
    ├── 01_ubam_to_fastq.sh   # Step 1: Convert uBAM
    ├── 02_align_pbmm2.sh     # Step 2: Alignment
    ├── 03_mark_duplicates.sh  # Step 3: Mark duplicates
    ├── 04_call_variants_deepvariant.sh  # Step 4: Small variants
    ├── 05_call_sv_pbsv.sh    # Step 5: Structural variants
    ├── 06_qc_coverage.sh     # Step 6: QC analysis
    ├── check_dependencies.sh  # Dependency checker
    └── yad_wrapper_example.sh # GUI wrapper example
```

## Installation on Your Systems

### 1. Copy Pipeline to Your Workstations
```bash
# Extract or copy the pacbio_pipeline directory
scp -r pacbio_pipeline/ user@workstation:/opt/
```

### 2. Install Dependencies
```bash
# Using conda (recommended)
conda create -n pacbio python=3.10
conda activate pacbio
conda install -c bioconda pbmm2 pbsv samtools bcftools deepvariant

# Or install Docker for DeepVariant
sudo apt install docker.io
docker pull google/deepvariant:1.6.0
```

### 3. Configure
```bash
cd /opt/pacbio_pipeline
nano configs/pipeline_config.sh
# Set REFERENCE_GENOME path and other parameters
```

### 4. Verify Installation
```bash
cd scripts
./check_dependencies.sh
```

## Usage Examples

### Command Line (Basic)
```bash
# Complete workflow for one sample
cd /opt/pacbio_pipeline/scripts

./02_align_pbmm2.sh -i sample.fastq.gz -r /data/hg38.fasta -o sample.bam -t 16 -p HiFi
./03_mark_duplicates.sh -i sample.bam -o sample.markdup.bam -t 8
./04_call_variants_deepvariant.sh -i sample.markdup.bam -r /data/hg38.fasta -o sample.vcf.gz -t 16 --docker
./05_call_sv_pbsv.sh -i sample.markdup.bam -r /data/hg38.fasta -o sample.sv.vcf -t 16 --ccs
./06_qc_coverage.sh -i sample.markdup.bam -r /data/hg38.fasta -o qc_results -t 8
```

### Using YAD GUI
```bash
# Install YAD if not present
sudo apt install yad

# Launch GUI
cd /opt/pacbio_pipeline/scripts
./yad_wrapper_example.sh
```

### Batch Processing
```bash
# Process multiple samples in parallel
for sample in /data/*.fastq.gz; do
    (
        name=$(basename "$sample" .fastq.gz)
        ./02_align_pbmm2.sh -i "$sample" -r /data/hg38.fasta -o "${name}.bam" -t 8
    ) &
done
wait
```
