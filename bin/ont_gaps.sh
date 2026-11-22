#!/usr/bin/env bash
# ONT support of every scaffold join (N gap): do long reads cross it?
# usage: ont_gaps.sh NAME FASTA BAM   (BAM sorted + indexed)
#
#   NAME.ont.gaps.tsv  seq start end spanning broken status   (0-based gap coordinates)
#   NAME.ont.gaps.bed  the same as a coloured IGV track: green supported, red conflict, grey no evidence
#
# spanning: reads aligned >= FLANK bp into both sides of the gap on the same strand, in one alignment or in
#           several pieces (minimap2 often breaks at the N run) = the join is right.
# broken:   other reads with an alignment that stops within FLANK bp of the gap with >= CLIP bp clipped
#           (the rest of the read aligns elsewhere or not at all) = evidence against the join.
# status:   supported (spanning >= 2), conflict (spanning = 0 and broken >= 2), no_evidence (otherwise)
set -euo pipefail
name=$1; fa=$2; bam=$3
FLANK=1000; CLIP=1000; MINQ=1

# N runs: one streaming pass, a run can continue over line ends
zcat -f "$fa" | awk -v OFS='\t' '
    function flush() { if (gs != "") print n, gs, ge; gs = "" }
    /^>/ { flush(); n = substr($1, 2); pos = 0; next }
    {
        L = length($0); rest = $0; off = 0
        while (match(rest, /[Nn]+/)) {
            s = pos + off + RSTART - 1; e = s + RLENGTH
            if (gs != "" && s == ge) ge = e; else { flush(); gs = s; ge = e }
            off += RSTART + RLENGTH - 1; rest = substr(rest, RSTART + RLENGTH)
        }
        pos += L
        if (gs != "" && ge < pos) flush()
    }
    END { flush() }' > gaps.tmp

printf "seq\tstart\tend\tspanning\tbroken\tstatus\n" > "$name.ont.gaps.tsv"
printf 'track name="ONT support of joins" itemRgb="On"\n' > "$name.ont.gaps.bed"
while IFS=$'\t' read -r chr gs ge; do
    a=$(( gs - FLANK + 1 )); [ "$a" -lt 1 ] && a=1
    samtools view -F 0x104 -q $MINQ "$bam" "$chr:$a-$(( ge + FLANK ))" \
    | awk -v gs="$gs" -v ge="$ge" -v F=$FLANK -v C=$CLIP '
        {
            q = $1; st = int($2 / 16) % 2; s0 = $4 - 1; cig = $6; rl = 0; lc = 0; tc = 0; first = 1
            while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
                k = substr(cig, 1, RLENGTH - 1) + 0; op = substr(cig, RLENGTH, 1)
                if (op ~ /[MDN=X]/) rl += k
                if (op ~ /[SH]/) { if (first) lc = k; else tc = k }
                first = 0; cig = substr(cig, RLENGTH + 1)
            }
            e0 = s0 + rl
            if (s0 <= gs - F && e0 >= gs - 50) left[q, st] = 1           # reaches the gap from the left
            if (e0 >= ge + F && s0 <= ge + 50) right[q, st] = 1          # reaches it from the right
            if (tc >= C && e0 >= gs - F && e0 <= ge) brk[q] = 1           # stops at the gap, rest clipped
            if (lc >= C && s0 >= gs && s0 <= ge + F) brk[q] = 1
            seen[q] = 1
        }
        END {
            sp = 0; br = 0
            for (q in seen) {
                if ((left[q, 0] && right[q, 0]) || (left[q, 1] && right[q, 1])) sp++
                else if (brk[q]) br++
            }
            print sp + 0, br + 0
        }' | {
        read -r sp br
        st=no_evidence; rgb=150,150,150
        if [ "$sp" -ge 2 ]; then st=supported; rgb=46,154,82
        elif [ "$br" -ge 2 ]; then st=conflict; rgb=200,50,40; fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$chr" "$gs" "$ge" "$sp" "$br" "$st" >> "$name.ont.gaps.tsv"
        printf "%s\t%s\t%s\t%s: %s spanning, %s broken\t0\t.\t%s\t%s\t%s\n" \
            "$chr" "$gs" "$ge" "$st" "$sp" "$br" "$gs" "$ge" "$rgb" >> "$name.ont.gaps.bed"
    }
done < gaps.tmp
rm -f gaps.tmp
