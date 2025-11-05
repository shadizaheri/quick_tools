#!/bin/bash

set -euo pipefail

# Simplified BAM to FASTQ conversion using samtools
# RECOMMENDED: Simpler approach with direct file output

usage() {
    cat << EOF
Usage: $0 -i INPUT_BAM -o OUTPUT_PREFIX [-t THREADS] [-d DOCKER_IMAGE]

Required:
  -i INPUT_BAM       Path to input BAM file
  -o OUTPUT_PREFIX   Output prefix for FASTQ files (e.g., sample1)

Optional:
  -t THREADS         Number of threads (default: 4)
  -d DOCKER_IMAGE    Docker image to use (default: bam2fq-samtools:latest)
  -h                 Show this help message

Example:
  $0 -i /data/sample.bam -o /data/sample -t 8

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
        i) INPUT_BAM="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) DOCKER_IMAGE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [ -z "${INPUT_BAM:-}" ] || [ -z "${OUTPUT_PREFIX:-}" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input file exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file not found: $INPUT_BAM"
    exit 1
fi

# Get absolute paths
INPUT_BAM=$(realpath "$INPUT_BAM")
OUTPUT_DIR=$(dirname "$(realpath "$OUTPUT_PREFIX")")
OUTPUT_BASE=$(basename "$OUTPUT_PREFIX")

echo "Starting BAM to FASTQ conversion..."
echo "Input BAM: $INPUT_BAM"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Threads: $THREADS"
echo "Docker image: $DOCKER_IMAGE"

# Run samtools bam2fq in Docker with direct file output
docker run --rm \
    -v "$OUTPUT_DIR:/output" \
    -v "$(dirname "$INPUT_BAM"):/input" \
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
            /input/$(basename "$INPUT_BAM")
    "

echo "Conversion complete!"
echo "Output files:"
echo "  ${OUTPUT_PREFIX}_1.fastq.gz"
echo "  ${OUTPUT_PREFIX}_2.fastq.gz"
echo "  ${OUTPUT_PREFIX}_unpaired.fastq.gz"


