#!/bin/bash

set -euo pipefail

# BAM to FASTQ conversion using samtools - streaming from S3 to S3
# This script streams the BAM file from S3, converts to FASTQ, uploads to S3, and cleans up local files

usage() {
    cat << EOF
Usage: $0 -i S3_BAM_PATH -o S3_OUTPUT_PREFIX [-t THREADS] [-d DOCKER_IMAGE]

Required:
  -i S3_BAM_PATH        S3 path to input BAM file (e.g., s3://bucket/path/file.bam)
  -o S3_OUTPUT_PREFIX   S3 output prefix for FASTQ files (e.g., s3://bucket/path/sample1)

Optional:
  -t THREADS            Number of threads (default: 4)
  -d DOCKER_IMAGE       Docker image to use (default: bam2fq-samtools:latest)
  -h                    Show this help message

Example:
  $0 -i s3://mybucket/sample.bam -o s3://mybucket/output/sample -t 8

Output files (uploaded to S3):
  {S3_OUTPUT_PREFIX}_1.fastq.gz     - Forward reads (R1)
  {S3_OUTPUT_PREFIX}_2.fastq.gz     - Reverse reads (R2)
  {S3_OUTPUT_PREFIX}_unpaired.fastq.gz - Singleton reads

EOF
    exit 1
}

# Default values
THREADS=4
DOCKER_IMAGE="bam2fq-samtools:latest"

# Parse arguments
while getopts "i:o:t:d:h" opt; do
    case $opt in
        i) S3_BAM_PATH="$OPTARG" ;;
        o) S3_OUTPUT_PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) DOCKER_IMAGE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [ -z "${S3_BAM_PATH:-}" ] || [ -z "${S3_OUTPUT_PREFIX:-}" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate S3 paths
if [[ ! "$S3_BAM_PATH" =~ ^s3:// ]]; then
    echo "Error: Input path must be an S3 path (starting with s3://)"
    exit 1
fi

if [[ ! "$S3_OUTPUT_PREFIX" =~ ^s3:// ]]; then
    echo "Error: Output prefix must be an S3 path (starting with s3://)"
    exit 1
fi

# Check if S3 input file exists
echo "Checking S3 input file..."
if ! aws s3 ls "$S3_BAM_PATH" > /dev/null 2>&1; then
    echo "Error: S3 BAM file not found or not accessible: $S3_BAM_PATH"
    exit 1
fi

# Extract bucket and path from S3 output prefix
S3_OUTPUT_DIR=$(dirname "$S3_OUTPUT_PREFIX")
OUTPUT_BASE=$(basename "$S3_OUTPUT_PREFIX")

# Create temporary directory for local output
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "Starting BAM to FASTQ conversion (S3 to S3)..."
echo "Input S3 BAM: $S3_BAM_PATH"
echo "Output S3 prefix: $S3_OUTPUT_PREFIX"
echo "Temporary local dir: $TEMP_DIR"
echo "Threads: $THREADS"
echo "Docker image: $DOCKER_IMAGE"

# Stream from S3 and convert to FASTQ
echo "Step 1/2: Converting BAM to FASTQ..."
aws s3 cp "$S3_BAM_PATH" - | \
docker run --rm -i \
    -v "$TEMP_DIR:/output" \
    --entrypoint /bin/bash \
    "$DOCKER_IMAGE" \
    -c "
        set -euo pipefail
        samtools bam2fq \
            -@ $THREADS \
            -1 >(pigz -p $THREADS > /output/${OUTPUT_BASE}_1.fastq.gz) \
            -2 >(pigz -p $THREADS > /output/${OUTPUT_BASE}_2.fastq.gz) \
            -s >(pigz -p $THREADS > /output/${OUTPUT_BASE}_unpaired.fastq.gz) \
            -0 /dev/null \
            -
    "

echo "Step 2/2: Uploading FASTQ files to S3..."
aws s3 cp "$TEMP_DIR/${OUTPUT_BASE}_1.fastq.gz" "${S3_OUTPUT_PREFIX}_1.fastq.gz"
aws s3 cp "$TEMP_DIR/${OUTPUT_BASE}_2.fastq.gz" "${S3_OUTPUT_PREFIX}_2.fastq.gz"
aws s3 cp "$TEMP_DIR/${OUTPUT_BASE}_unpaired.fastq.gz" "${S3_OUTPUT_PREFIX}_unpaired.fastq.gz"

echo "Cleaning up temporary files..."
# Cleanup happens automatically via trap

echo "Conversion complete!"
echo "Output files uploaded to S3:"
echo "  ${S3_OUTPUT_PREFIX}_1.fastq.gz"
echo "  ${S3_OUTPUT_PREFIX}_2.fastq.gz"
echo "  ${S3_OUTPUT_PREFIX}_unpaired.fastq.gz"
