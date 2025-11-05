#!/bin/bash

set -euo pipefail

# BAM to FASTQ conversion using samtools - streaming from S3
# This script streams the BAM file directly from S3 without downloading

usage() {
    cat << EOF
Usage: $0 -i S3_BAM_PATH -o OUTPUT_PREFIX [-t THREADS] [-d DOCKER_IMAGE]

Required:
  -i S3_BAM_PATH     S3 path to input BAM file (e.g., s3://bucket/path/file.bam)
  -o OUTPUT_PREFIX   Output prefix for FASTQ files (e.g., sample1)

Optional:
  -t THREADS         Number of threads (default: 4)
  -d DOCKER_IMAGE    Docker image to use (default: bam2fq-samtools:latest)
  -h                 Show this help message

Example:
  $0 -i s3://mybucket/sample.bam -o /data/sample -t 8

Output files:
  {OUTPUT_PREFIX}_1.fastq.gz     - Forward reads (R1)
  {OUTPUT_PREFIX}_2.fastq.gz     - Reverse reads (R2)
  {OUTPUT_PREFIX}_unpaired.fastq.gz - Singleton reads

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
        o) OUTPUT_PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) DOCKER_IMAGE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [ -z "${S3_BAM_PATH:-}" ] || [ -z "${OUTPUT_PREFIX:-}" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate S3 path
if [[ ! "$S3_BAM_PATH" =~ ^s3:// ]]; then
    echo "Error: Input path must be an S3 path (starting with s3://)"
    exit 1
fi

# Check if S3 file exists
echo "Checking S3 file..."
if ! aws s3 ls "$S3_BAM_PATH" > /dev/null 2>&1; then
    echo "Error: S3 BAM file not found or not accessible: $S3_BAM_PATH"
    exit 1
fi

# Get absolute path for output
OUTPUT_DIR=$(dirname "$(realpath "$OUTPUT_PREFIX")")
OUTPUT_BASE=$(basename "$OUTPUT_PREFIX")

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Starting BAM to FASTQ conversion (streaming from S3)..."
echo "Input S3 BAM: $S3_BAM_PATH"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Threads: $THREADS"
echo "Docker image: $DOCKER_IMAGE"

# Stream from S3 and convert to FASTQ
aws s3 cp "$S3_BAM_PATH" - | \
docker run --rm -i \
    -v "$OUTPUT_DIR:/output" \
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

echo "Conversion complete!"
echo "Output files:"
echo "  ${OUTPUT_PREFIX}_1.fastq.gz"
echo "  ${OUTPUT_PREFIX}_2.fastq.gz"
echo "  ${OUTPUT_PREFIX}_unpaired.fastq.gz"
