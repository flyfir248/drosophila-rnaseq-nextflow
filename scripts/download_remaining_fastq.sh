#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw logs

MANIFEST="logs/remaining_fastq_manifest.tsv"

awk -F '\t' '
    NR > 1 && $1 != "SRR1576457" {
        n = split($3, sizes, ";")
        split($7, urls, ";")
        for (i = 1; i <= n; i++)
            print $1 "\t" "https://" urls[i] "\t" sizes[i]
    }
' data/metadata/ena_run_metadata.tsv > "$MANIFEST"

echo "Prepared $(wc -l < "$MANIFEST") download tasks."
echo "Maximum concurrent downloads: 2"

xargs -P 2 -n 3 bash -c '
    set -euo pipefail

    run="$1"
    url="$2"
    expected="$3"
    filename="${url##*/}"
    output="data/raw/$filename"
    logfile="logs/$filename.log"

    if [[ -f "$output" ]]; then
        actual=$(stat -c%s "$output")

        if [[ "$actual" -eq "$expected" ]]; then
            echo "SKIP: $filename already complete"
            exit 0
        fi

        if [[ "$actual" -gt "$expected" ]]; then
            echo "RESTART: $filename is oversized ($actual bytes)"
            rm -f "$output"
        fi
    fi

    echo "START: $run / $filename"

    curl -fL \
        --retry 8 \
        --retry-all-errors \
        --retry-delay 15 \
        --connect-timeout 30 \
        --max-time 0 \
        -A "Mozilla/5.0" \
        --continue-at - \
        --output "$output" \
        "$url" 2>> "$logfile"

    actual=$(stat -c%s "$output")

    if [[ "$actual" -ne "$expected" ]]; then
        echo "ERROR: $filename expected=$expected actual=$actual"
        exit 1
    fi

    echo "COMPLETE: $filename ($actual bytes)"
' _ < "$MANIFEST"

rm -f "$MANIFEST"

echo ""
echo "All remaining FASTQ files downloaded and size-verified."
