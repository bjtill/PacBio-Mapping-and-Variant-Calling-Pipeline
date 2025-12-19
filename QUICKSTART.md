# Quick Start Guide - PacBio Pipeline

## Installation

### 1. Extract the Pipeline

```bash
# If you received a tar.gz file:
tar -xzf pacbio_pipeline.tar.gz
cd pacbio_pipeline

# Or if you cloned from git:
cd pacbio_pipeline
```

### 2. Install Dependencies

The easiest way to install all dependencies is using conda:

```bash
# Create a new conda environment
conda create -n pacbio python=3.10
conda activate pacbio

# Install PacBio tools and samtools
conda install -c bioconda pbmm2 pbsv pbmarkdup samtools bcftools tabix

# Install DeepVariant (optional - you can also use Docker)
conda install -c bioconda deepvariant
```

**Alternative: Using Docker for DeepVariant**

If you prefer Docker (recommended for DeepVariant):

```bash
# Install Docker
sudo apt update
sudo apt install docker.io
sudo usermod -aG docker $USER
# Log out and back in for group changes to take effect

# Pull DeepVariant image
docker pull google/deepvariant:1.6.0
```

### 3. Check Installation

Run the dependency checker to verify everything is installed:

```bash
cd scripts
./check_dependencies.sh
```

This will report:
- ✓ Tools that are installed correctly
- ✗ Tools that are missing
- ⚠ Optional tools that are missing but not required

### 4. Configure the Pipeline

Edit the configuration file:

```bash
nano configs/pipeline_config.sh
```

**Essential settings to configure:**

```bash
# Set your reference genome path
REFERENCE_GENOME="/path/to/your/reference/genome.fasta"

# Adjust thread counts based on your hardware
ALIGN_THREADS=16
DEEPVARIANT_THREADS=16
PBSV_THREADS=16

# If using Docker for DeepVariant
USE_DOCKER=true
DEEPVARIANT_DOCKER="google/deepvariant:1.6.0"
```

### 5. Test the Pipeline

Run a test with small data:

```bash
# Test alignment
cd scripts
./02_align_pbmm2.sh --help

# Try with your data
./02_align_pbmm2.sh \
    -i /path/to/your/sample.fastq.gz \
    -r /path/to/reference.fasta \
    -o test_output.bam \
    -t 16
```

## Using the Pipeline

### Command Line Usage

Each script can be run independently:

```bash
# Step 1: Convert uBAM to FASTQ (if needed)
./01_ubam_to_fastq.sh -i sample.ubam -o sample.fastq.gz -t 8

# Step 2: Align reads
./02_align_pbmm2.sh \
    -i sample.fastq.gz \
    -r reference.fasta \
    -o sample.bam \
    -t 16 \
    -p HiFi

# Step 3: Mark duplicates
./03_mark_duplicates.sh \
    -i sample.bam \
    -o sample.markdup.bam \
    -t 8

# Step 4: Call small variants
./04_call_variants_deepvariant.sh \
    -i sample.markdup.bam \
    -r reference.fasta \
    -o sample.vcf.gz \
    -t 16 \
    -m PACBIO \
    --docker

# Step 5: Call structural variants
./05_call_sv_pbsv.sh \
    -i sample.markdup.bam \
    -r reference.fasta \
    -o sample.sv.vcf \
    -t 16 \
    --ccs

# Step 6: Generate QC report
./06_qc_coverage.sh \
    -i sample.markdup.bam \
    -r reference.fasta \
    -o qc_results \
    -t 8
```

### Using the YAD GUI Wrapper

If you have YAD installed:

```bash
# Install YAD (if not already installed)
sudo apt install yad

# Launch the GUI
./yad_wrapper_example.sh
```

The GUI provides:
- File selection dialogs
- Parameter configuration forms
- Progress indicators
- Result notifications

### Batch Processing

To process multiple samples:

```bash
#!/bin/bash
# process_all.sh

REFERENCE="/path/to/reference.fasta"
THREADS=16

for fastq in /path/to/samples/*.fastq.gz; do
    sample=$(basename "$fastq" .fastq.gz)
    echo "Processing $sample..."
    
    # Align
    ./02_align_pbmm2.sh -i "$fastq" -r "$REFERENCE" -o "${sample}.bam" -t $THREADS
    
    # Mark duplicates
    ./03_mark_duplicates.sh -i "${sample}.bam" -o "${sample}.markdup.bam" -t 8
    
    # Call variants
    ./04_call_variants_deepvariant.sh \
        -i "${sample}.markdup.bam" \
        -r "$REFERENCE" \
        -o "${sample}.vcf.gz" \
        -t $THREADS --docker
    
    # QC
    ./06_qc_coverage.sh \
        -i "${sample}.markdup.bam" \
        -r "$REFERENCE" \
        -o "qc_${sample}" \
        -t 8
done
```

## Troubleshooting

### "Command not found" errors

Make sure conda environment is activated:
```bash
conda activate pacbio
```

Or add tools to your PATH:
```bash
export PATH="/path/to/conda/envs/pacbio/bin:$PATH"
```

### "BAM index not found" errors

The pipeline should create indexes automatically, but if needed:
```bash
samtools index -@ 8 yourfile.bam
```

### DeepVariant memory errors

Reduce the number of threads or use Docker:
```bash
./04_call_variants_deepvariant.sh ... -t 8 --docker
```

### Out of disk space

Check space before running:
```bash
df -h /path/to/output/directory
```

Estimate required space: ~5-10x the input file size.

### pbmm2 alignment fails

Make sure reference is indexed:
```bash
samtools faidx /path/to/reference.fasta
```

Check that preset matches your data type (HiFi for CCS reads).

## Getting Help

For help with any script:
```bash
./script_name.sh --help
```

Check log files (created in output directories):
- `alignment.log`
- `deepvariant.log`
- `pbsv.log`
- `qc.log`

## Directory Organization

Recommended directory structure:

```
project/
├── raw/              # Original uBAM or FASTQ files
├── aligned/          # Aligned BAM files
├── variants/         # VCF files
│   ├── snps/        # Small variants
│   └── svs/         # Structural variants
├── qc/              # QC reports
└── logs/            # Log files
```

## Performance Tips

### For 16-core, 256GB systems:

- **Alignment**: Use `-t 16` (all cores)
- **DeepVariant**: Use `-t 16` (benefits from more cores)
- **pbsv**: Use `-t 16`
- **QC**: Use `-t 8` (less intensive)

### Process multiple samples in parallel:

```bash
# Process 2 samples simultaneously
sample1.sh & sample2.sh & wait
```

But be careful not to exceed memory limits!

### For large genomes:

- Consider processing by chromosome/region
- Use BED files to specify regions
- Increase sort memory: `-m 8G`

## Common Workflows

### Workflow 1: Human WGS (HiFi)

```bash
./02_align_pbmm2.sh -i sample.fastq.gz -r hg38.fasta -o sample.bam -t 16 -p HiFi
./03_mark_duplicates.sh -i sample.bam -o sample.markdup.bam -t 8
./04_call_variants_deepvariant.sh -i sample.markdup.bam -r hg38.fasta -o sample.vcf.gz -t 16 -m PACBIO --docker
./05_call_sv_pbsv.sh -i sample.markdup.bam -r hg38.fasta -o sample.sv.vcf -t 16 --ccs
./06_qc_coverage.sh -i sample.markdup.bam -r hg38.fasta -o qc -t 8
```

### Workflow 2: Targeted/Exome (HiFi)

```bash
./02_align_pbmm2.sh -i sample.fastq.gz -r hg38.fasta -o sample.bam -t 16 -p HiFi
./03_mark_duplicates.sh -i sample.bam -o sample.markdup.bam -t 8
./04_call_variants_deepvariant.sh -i sample.markdup.bam -r hg38.fasta -o sample.vcf.gz -t 16 -m PACBIO --regions targets.bed --docker
./06_qc_coverage.sh -i sample.markdup.bam -r hg38.fasta -o qc -t 8 -b targets.bed
```

## Next Steps

1. Read the full documentation: `docs/README.md`
2. Customize the YAD wrapper for your specific needs
3. Set up automated batch processing
4. Integrate with your existing pipelines

## Updates and Support

Keep your tools updated:
```bash
conda activate pacbio
conda update --all
```

For issues with specific tools, consult their documentation:
- pbmm2: https://github.com/PacificBiosciences/pbmm2
- DeepVariant: https://github.com/google/deepvariant
- pbsv: https://github.com/PacificBiosciences/pbsv
