# BAM to FASTQ Conversion - Local Docker Setup

This directory contains Dockerfiles and scripts for converting BAM files to FASTQ format locally, with a path to cloud deployment.

## Quick Start (Recommended: samtools)

### 1. Build the Docker image

```bash
docker build -f Dockerfile.samtools -t bam2fq-samtools:latest .
```

### 2. Run the conversion

```bash
chmod +x bam_to_fastq_samtools_simple.sh
./bam_to_fastq_samtools_simple.sh -i /path/to/input.bam -o /path/to/output/sample -t 8
```

This will create:
- `sample_1.fastq.gz` (R1 reads)
- `sample_2.fastq.gz` (R2 reads)
- `sample_unpaired.fastq.gz` (singleton reads)

## Methods Available

### Method 1: samtools (RECOMMENDED)

**Pros:**
- Fast and lightweight (~20MB image)
- Single-step conversion
- Perfect for most use cases
- Lower resource requirements

**When to use:** Default choice for BAM to FASTQ conversion

```bash
# Build
docker build -f Dockerfile.samtools -t bam2fq-samtools:latest .

# Run
./bam_to_fastq_samtools_simple.sh -i input.bam -o output_prefix -t 8
```

### Method 2: GATK (For special cases)

**Pros:**
- Restores original quality scores (OQ tags)
- Data sanitization and validation
- Better for BQSR-processed BAMs

**Cons:**
- Larger image (~500MB+)
- Slower (two-step process)
- More resource intensive

**When to use:** 
- You need to restore pre-BQSR quality scores
- You need data validation/sanitization
- You're in a GATK-heavy pipeline

```bash
# Build
docker build -f Dockerfile.gatk -t bam2fq-gatk:latest .

# Run
./bam_to_fastq_gatk.sh -i input.bam -o output_prefix
```

## Script Options

### bam_to_fastq_samtools_simple.sh

```
Required:
  -i INPUT_BAM       Path to input BAM file
  -o OUTPUT_PREFIX   Output prefix for FASTQ files

Optional:
  -t THREADS         Number of threads (default: 4)
  -d DOCKER_IMAGE    Docker image to use (default: bam2fq-samtools:latest)
```

### bam_to_fastq_gatk.sh

```
Required:
  -i INPUT_BAM       Path to input BAM file
  -o OUTPUT_PREFIX   Output prefix for FASTQ files

Optional:
  -r REFERENCE       Reference FASTA (rarely needed for BAM)
  -d DOCKER_IMAGE    Docker image to use (default: bam2fq-gatk:latest)
```

## Path to Cloud Deployment

These Docker images and scripts are designed to work both locally and in the cloud:

### For Google Cloud (GCP)

1. **Push image to Container Registry:**
```bash
# Tag for GCR
docker tag bam2fq-samtools:latest gcr.io/YOUR-PROJECT/bam2fq-samtools:latest

# Push to GCR
docker push gcr.io/YOUR-PROJECT/bam2fq-samtools:latest
```

2. **Update your WDL to use the custom image:**
```wdl
runtime {
    docker: "gcr.io/YOUR-PROJECT/bam2fq-samtools:latest"
    memory: "4G"
    cpu: 4
    disks: "local-disk 50 HDD"
}
```

### For AWS

1. **Push to ECR:**
```bash
# Authenticate
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin YOUR-ACCOUNT.dkr.ecr.us-east-1.amazonaws.com

# Tag and push
docker tag bam2fq-samtools:latest YOUR-ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/bam2fq-samtools:latest
docker push YOUR-ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/bam2fq-samtools:latest
```

### For Azure

1. **Push to ACR:**
```bash
# Login
az acr login --name YOUR-REGISTRY

# Tag and push
docker tag bam2fq-samtools:latest YOUR-REGISTRY.azurecr.io/bam2fq-samtools:latest
docker push YOUR-REGISTRY.azurecr.io/bam2fq-samtools:latest
```

## Performance Tips

1. **Use appropriate thread count**: Set `-t` to the number of available CPU cores
2. **Use SSD storage**: Significantly faster for I/O-heavy operations
3. **Monitor memory**: BAM files can require substantial memory for processing
4. **Consider file size**: For very large BAM files (>100GB), consider:
   - Using more memory
   - Processing in chunks
   - Using cloud compute with high I/O

## Troubleshooting

### Docker permission issues
```bash
# Add user to docker group
sudo usermod -aG docker $USER
# Then log out and back in
```

### Out of memory errors
Increase Docker's memory allocation in Docker Desktop settings or add more swap space.

### Slow performance
- Increase thread count with `-t`
- Use faster storage (SSD)
- Ensure Docker has adequate CPU/memory resources

## Testing

Test with a small BAM file first:

```bash
# Create a small test BAM (if you have samtools installed)
samtools view -h -s 0.01 large.bam -o small_test.bam

# Run conversion
./bam_to_fastq_samtools_simple.sh -i small_test.bam -o test_output -t 4
```

## Resource Requirements

### samtools method
- **Memory**: 2-4GB for most BAM files
- **CPU**: Scales well with thread count
- **Disk**: 2-3x input BAM size (for output FASTQ.gz)
- **Time**: ~10-30 minutes per 30GB BAM (8 threads)

### GATK method
- **Memory**: 4-8GB for most BAM files
- **CPU**: 2-4 cores
- **Disk**: 3-4x input BAM size (intermediate + output files)
- **Time**: ~20-60 minutes per 30GB BAM

## Support

For issues or questions, contact Shadi Zaheri.

