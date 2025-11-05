#!/bin/bash

set -euo pipefail

# BAM to FASTQ conversion using GATK
# Use this if you need OQ tag restoration or data sanitization

usage() {
    cat << EOF
Usage: $0 -i INPUT_BAM -o OUTPUT_PREFIX [-r REFERENCE] [-d DOCKER_IMAGE]

Required:
  -i INPUT_BAM       Path to input BAM file
  -o OUTPUT_PREFIX   Output prefix for FASTQ files (e.g., sample1)

Optional:
  -r REFERENCE       Reference FASTA (not typically needed for BAM)
  -d DOCKER_IMAGE    Docker image to use (default: bam2fq-gatk:latest)
  -h                 Show this help message

Example:
  $0 -i /data/sample.bam -o /data/sample

Output files:
  {OUTPUT_PREFIX}_1.fastq.gz     - Forward reads (R1)
  {OUTPUT_PREFIX}_2.fastq.gz     - Reverse reads (R2)
  {OUTPUT_PREFIX}_unpaired.fastq.gz - Singleton reads

EOF
    exit 1
}

# Default values
DOCKER_IMAGE="bam2fq-gatk:latest"
REFERENCE=""

# Parse arguments
while getopts "i:o:r:d:h" opt; do
    case $opt in
        i) INPUT_BAM="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
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

echo "Starting BAM to FASTQ conversion (GATK method)..."
echo "Input BAM: $INPUT_BAM"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Docker image: $DOCKER_IMAGE"

# Create temp directory for intermediate files
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Step 1: RevertSam - Convert BAM to uBAM
echo "Step 1/2: Reverting BAM to unmapped BAM..."
docker run --rm \
    -v "$(dirname "$INPUT_BAM"):/input" \
    -v "$TEMP_DIR:/temp" \
    "$DOCKER_IMAGE" \
    RevertSam \
    -I "/input/$(basename "$INPUT_BAM")" \
    -O /temp/unmapped.bam \
    --SANITIZE true \
    --REMOVE_ALIGNMENT_INFORMATION true \
    --RESTORE_ORIGINAL_QUALITIES true

# Step 2: SamToFastq - Extract FASTQ from uBAM
echo "Step 2/2: Extracting FASTQ files..."
docker run --rm \
    -v "$TEMP_DIR:/temp" \
    -v "$OUTPUT_DIR:/output" \
    --entrypoint /bin/bash \
    "$DOCKER_IMAGE" \
    -c "
        set -euo pipefail
        
        # Extract FASTQ files
        gatk SamToFastq \
            -I /temp/unmapped.bam \
            -F /temp/${OUTPUT_BASE}_1.fastq \
            -F2 /temp/${OUTPUT_BASE}_2.fastq \
            -FU /temp/${OUTPUT_BASE}_unpaired.fastq
        
        # Compress with pigz
        pigz -c /temp/${OUTPUT_BASE}_1.fastq > /output/${OUTPUT_BASE}_1.fastq.gz
        pigz -c /temp/${OUTPUT_BASE}_2.fastq > /output/${OUTPUT_BASE}_2.fastq.gz
        pigz -c /temp/${OUTPUT_BASE}_unpaired.fastq > /output/${OUTPUT_BASE}_unpaired.fastq.gz
    "

echo "Conversion complete!"
echo "Output files:"
echo "  ${OUTPUT_PREFIX}_1.fastq.gz"
echo "  ${OUTPUT_PREFIX}_2.fastq.gz"
echo "  ${OUTPUT_PREFIX}_unpaired.fastq.gz"


