#!/bin/bash
# save as: trim.sh

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Get absolute paths
INPUT_DIR=$(realpath "$INPUT_DIR")
OUTPUT_DIR=$(realpath -m "$OUTPUT_DIR")

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Go to input directory
cd "$INPUT_DIR" || { echo "Cannot cd to $INPUT_DIR"; exit 1; }

# Find and process all samples
for file in *_sub_R1.fq; do
    # Extract sample name (remove _sub_R1.fq)
    sample="${file%_sub_R1.fq}"
    
    echo "Processing: $sample"
    
    # Run cutadapt with FULL paths
    cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
             -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
             -m 215 -M 285 --discard-untrimmed \
             -o "${OUTPUT_DIR}/${sample}_R1_trim.fq.gz" \
             -p "${OUTPUT_DIR}/${sample}_R2_trim.fq.gz" \
             "${INPUT_DIR}/${sample}_sub_R1.fq" \
             "${INPUT_DIR}/${sample}_sub_R2.fq"
done

echo "Done! Output in: $OUTPUT_DIR"
