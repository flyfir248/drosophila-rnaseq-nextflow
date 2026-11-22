#!/usr/bin/env bash
set -euo pipefail

MANIFEST="logs/corrupt_fastq_manifest.tsv"
mkdir -p logs data/raw

: > "$MANIFEST"

while IFS=$'\t' read -r run sizes checksums urls; do
    [[ "$run" == "run_accession" ]] && continue

    IFS=';' read -ra size_array <<< "$sizes"
    IFS=';' read -ra md5_array <<< "$checksums"
    IFS=';' read -ra url_array <<< "$urls"

    for i in "${!url_array[@]}"; do
        filename="${url_array[$i]##*/}"
        file="data/raw/$filename"
        expected_md5="${md5_array[$i]}"

        if [[ -f "$file" ]] &&
           [[ "$(md5sum "$file" | awk '{print $1}')" == "$expected_md5" ]] &&
           gzip -t "$file" 2>/dev/null; then
            echo "VALID: $filename"
        else
            printf '%s\t%s\t%s\t%s\n' \
                "$run" "https://${url_array[$i]}" \
                "${size_array[$i]}" "$expected_md5" >> "$MANIFEST"
        fi
    done
done < data/metadata/ena_fastq_checksums.tsv

echo ""
echo "Files requiring repair: $(wc -l < "$MANIFEST")"
echo "Starting up to two fresh downloads concurrently."

xargs -P 2 -n 4 bash -c '
    set -euo pipefail

    run="$1"
    url="$2"
    expected_size="$3"
    expected_md5="$4"

    filename="${url##*/}"
    output="data/raw/$filename"
    part="${output}.part"
    logfile="logs/${filename}.repair.log"

    echo "REPAIR: $run / $filename"

    rm -f "$part"

    curl -fLsS \
        --retry 8 \
        --retry-all-errors \
        --retry-delay 15 \
        --connect-timeout 30 \
        -A "Mozilla/5.0" \
        --output "$part" \
        "$url" 2>> "$logfile"

    actual_size=$(stat -c%s "$part")
    [[ "$actual_size" -eq "$expected_size" ]] || {
        echo "SIZE ERROR: $filename expected=$expected_size actual=$actual_size"
        exit 1
    }

    echo "$expected_md5  $part" | md5sum -c -
    gzip -t "$part"

    mv "$part" "$output"
    echo "REPAIRED: $filename"
' _ < "$MANIFEST"

echo ""
echo "Repair downloads completed."
