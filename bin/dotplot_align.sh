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

zcat -f "$tfa" | awk '/^>/{if(n!="")print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)} END{if(n!="")print n"\t"l}' > t.len

# one streaming pass over the query: cut every sequence into <seq>::<offset> pieces of CH bp and
# write q.len. Lines are cut where a piece ends (no string building: linear time, one line in memory).
zcat -f "$qfa" | awk -v CH="$chunk" '
    function done() { if (n != "") print n "\t" pos > "q.len" }
    /^>/ { done(); n = substr($1, 2); pos = 0; next }
    {
        L = length($0); p = 1
        while (L > 0) {
            if (pos % CH == 0) print ">" n "::" pos          # a new piece starts here
            k = CH - pos % CH; if (k > L) k = L
            print (p == 1 && k == L) ? $0 : substr($0, p, k)
            pos += k; p += k; L -= k
        }
    }
    END { done() }' > q.chunks.fa

minimap2 -x "$preset" -t "$threads" --secondary=no "$tfa" q.chunks.fa 2> minimap2.log \
  | awk -v OFS='\t' '
      NR == FNR { len[$1] = $2; next }
      { match($1, /::[0-9]+$/); off = substr($1, RSTART + 2) + 0; $1 = substr($1, 1, RSTART - 1)
        $2 = len[$1]; $3 += off; $4 += off; print }' q.len - \
  | gzip -c > "$name.paf.gz"
rm -f q.chunks.fa
