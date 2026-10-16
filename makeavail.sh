#!/bin/bash

TARGET_DIR="${1:-.}" # Defaults to current directory
PREFIX="${2:-EA}"    # Defaults to "EA"
OUTPUT_FILE="AVAILABLE"

# Setup header
printf "XXXXXX EMPTY XXXXXXXX\n" >"$OUTPUT_FILE"
printf "XXXXXX EMPTY XXXXXXXX\n" >>"$OUTPUT_FILE"
printf "%-8s %-6s      %-15s      %s\n" "YYYYMMDD" "HHMMSS" "name" "ON DISK" >>"$OUTPUT_FILE"

# Search path
for filepath in "${TARGET_DIR%/}/"$PREFIX*; do

  if [ -f "$filepath" ]; then
    filename=$(basename "$filepath")

    if [ "$filename" = "$OUTPUT_FILE" ]; then
      continue
    fi

    # Strip prefix so we can get the date
    digits="${filename#$PREFIX}"

    # Format assumed: YYMMDDHH (8 characters total left after prefix)
    YY=${digits:0:2}
    MM=${digits:2:2}
    DD=${digits:4:2}
    HH=${digits:6:2}

    YYYYMMDD="20${YY}${MM}${DD}"
    HHMMSS="${HH}0000"
    printf "%-8s %-6s      %-15s      ON DISK\n" "$YYYYMMDD" "$HHMMSS" "$filename" >>"$OUTPUT_FILE"

  fi
done

echo "AVAILABLE file successfully generated: $OUTPUT_FILE"
