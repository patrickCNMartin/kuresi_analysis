#!/bin/bash

# Output file
OUTPUT_FILE="data/fastqs/md5_fastq.txt"

# Define your file paths here
FILE_PATHS=(data/fastqs/*)

# Clear or create the output file
> "$OUTPUT_FILE"

# Write header
printf "%-40s %s\n" "FILENAME" "MD5_HASH" >> "$OUTPUT_FILE"
printf "%-40s %s\n" "--------" "--------" >> "$OUTPUT_FILE"

# Loop over files and compute md5sum
for FILEPATH in "${FILE_PATHS[@]}"; do
    if [ -f "$FILEPATH" ]; then
        FILENAME=$(basename "$FILEPATH")
        HASH=$(md5sum "$FILEPATH" | awk '{print $1}')
        printf "%-40s %s\n" "$FILENAME" "$HASH" >> "$OUTPUT_FILE"
    else
        echo "WARNING: File not found: $FILEPATH" >&2
    fi
done

echo "Done! Results written to $OUTPUT_FILE"