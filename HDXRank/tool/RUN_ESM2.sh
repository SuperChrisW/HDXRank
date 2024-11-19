#!/bin/bash

# Define the paths
MODEL="esm2_t6_8M_UR50D"
FASTA_DIR="/home/lwang/models/HDX_LSTM/data/Latest_test/fasta_files"
OUTPUT_DIR="/home/lwang/models/HDX_LSTM/data/Latest_test/esm2_t6_8M_UR50D"

mkdir -p "$OUTPUT_DIR"

# Define the filename pattern to match (use "*" to match all files)
FILE_NAME="*"

for FASTA_FILE in "$FASTA_DIR"/*.fasta; do
  BASENAME=$(basename "$FASTA_FILE" .fasta)
  
  if [[ "$FILE_NAME" == "*" || "$BASENAME" == "$FILE_NAME" ]]; then
    OUTPUT_PATH="$OUTPUT_DIR"

    # Construct the esm-extract command
    COMMAND="esm-extract $MODEL $FASTA_FILE $OUTPUT_PATH --include per_tok"

    if eval "$COMMAND"; then
      echo "Processed $FASTA_FILE successfully"
    else
      echo "Error processing $FASTA_FILE" >&2
    fi
  fi
done

echo "All matching FASTA files have been processed."
