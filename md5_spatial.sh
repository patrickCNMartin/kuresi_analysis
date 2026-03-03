#!/bin/bash

SPATIAL_DIR="data/spatial"
OUTPUT_FILE="data/spatial/md5_spatial.txt"

> "$OUTPUT_FILE"
printf "%-60s %s\n" "FILENAME" "MD5_HASH" >> "$OUTPUT_FILE"
printf "%-60s %s\n" "--------" "--------" >> "$OUTPUT_FILE"

for SAMPLE_DIR in "$SPATIAL_DIR"/*/; do
    echo "Found: $SAMPLE_DIR"
    SAMPLE=$(basename "$SAMPLE_DIR")
    for FILEPATH in "$SAMPLE_DIR"*; do
        if [ -f "$FILEPATH" ]; then
            FILENAME="${SAMPLE}/$(basename "$FILEPATH")"
            HASH=$(md5sum "$FILEPATH" | awk '{print $1}')
            printf "%-60s %s\n" "$FILENAME" "$HASH" >> "$OUTPUT_FILE"
        fi
    done
done

echo "Done! Results written to $OUTPUT_FILE"