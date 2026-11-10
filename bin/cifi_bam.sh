#!/usr/bin/env bash
# CiFi reads -> one unaligned BAM (as in cifiasm): FASTQ is imported, several BAMs are concatenated.
# usage: cifi_bam.sh OUT.bam THREADS file1 [file2 ...]
set -euo pipefail
out=$1; t=$2; shift 2
bams=()
for f in "$@"; do
  if [[ "$f" == *.bam ]]; then
    bams+=("$f")
  else
    b=$(basename "$f"); ub="${b%%.f*}.import.bam"
    samtools import -@ "$t" -0 "$f" -o "$ub"
    bams+=("$ub")
  fi
done
if (( ${#bams[@]} > 1 )); then
  samtools cat -@ "$t" -o "$out" "${bams[@]}"
else
  mv "${bams[0]}" "$out"      # single imported FASTQ
fi
