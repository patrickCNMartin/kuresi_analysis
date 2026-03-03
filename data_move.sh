#!/bin/bash
YAML_FILE="config/config.yaml"
SPATIAL_DIR="data/spatial"

IN_INPUT_FILES=false

while IFS= read -r line; do

    if [[ "$line" =~ ^input_files:[[:space:]]*$ ]]; then
        IN_INPUT_FILES=true
        continue
    fi

    if [[ "$IN_INPUT_FILES" == true && "$line" =~ ^[a-zA-Z] ]]; then
        IN_INPUT_FILES=false
    fi

    if [[ "$IN_INPUT_FILES" == false ]]; then
        continue
    fi

    if [[ "$line" =~ ^[[:space:]]{2}([a-zA-Z0-9_]+):[[:space:]]*$ ]]; then
        SAMPLE="${BASH_REMATCH[1]}"
    fi

    if [[ "$line" =~ ^[[:space:]]+(coordinates|counts|image|scale)[[:space:]]*:[[:space:]]*\"(.+)\" ]]; then
        FILE_PATH="${BASH_REMATCH[2]}"
        DEST_DIR="$SPATIAL_DIR/$SAMPLE"
        mkdir -p "$DEST_DIR"

        if [ -f "$FILE_PATH" ]; then
            cp "$FILE_PATH" "$DEST_DIR/"
            echo "Copied: $FILE_PATH -> $DEST_DIR/"
        else
            echo "WARNING: File not found: $FILE_PATH" >&2
        fi
    fi

done < "$YAML_FILE"

echo "Done!"