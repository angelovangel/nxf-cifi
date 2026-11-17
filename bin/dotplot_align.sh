#!/usr/bin/env bash
# Whole-assembly alignment for a dot plot: QUERY vs TARGET -> NAME.paf.gz, q.len, t.len
# usage: dotplot_align.sh NAME QUERY.fa[.gz] TARGET.fa[.gz] PRESET THREADS CHUNK_BP
#
# Query sequences are cut into CHUNK_BP pieces before mapping and the PAF is lifted back to the
# original sequences. Chromosome-scale queries (100s of Mb) in repetitive genomes otherwise make
# minimap2 build chains with > 2^31 anchors (assertion `n_v < INT32_MAX` in mg_chain_backtrack)
# and use a lot of memory. Alignments are split at chunk boundaries, which a dot plot does not show.
set -euo pipefail
name=$1; qfa=$2; tfa=$3; preset=$4; threads=$5; chunk=$6

seqlen() { zcat -f "$1" | awk '/^>/{if(n!="")print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)} END{if(n!="")print n"\t"l}'; }
seqlen "$qfa" > q.len
seqlen "$tfa" > t.len

# chunks are named <seq>::<offset>
zcat -f "$qfa" | awk -v CH="$chunk" '
    function out(s) { print ">" n "::" off; print s; off += length(s) }
    /^>/ { if (seq != "") out(seq); seq = ""; n = substr($1, 2); off = 0; next }
    { seq = seq $0; while (length(seq) >= CH) { out(substr(seq, 1, CH)); seq = substr(seq, CH + 1) } }
    END { if (seq != "") out(seq) }' > q.chunks.fa

minimap2 -x "$preset" -t "$threads" --secondary=no "$tfa" q.chunks.fa 2> minimap2.log \
  | awk -v OFS='\t' '
      NR == FNR { len[$1] = $2; next }
      { i = index($1, "::"); while ((j = index(substr($1, i + 2), "::")) > 0) i += j + 1   # last "::"
        off = substr($1, i + 2) + 0; $1 = substr($1, 1, i - 1)
        $2 = len[$1]; $3 += off; $4 += off; print }' q.len - \
  | gzip -c > "$name.paf.gz"
rm -f q.chunks.fa
