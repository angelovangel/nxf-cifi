#!/usr/bin/env python3
"""Contiguity, gap and telomere stats for one assembly (same definitions as compare_assemblies.sh).

usage: asm_stats.py LABEL FASTA TELO_F TELO_R TELO_WINDOW TELO_MIN

writes LABEL.stats.tsv (one row), LABEL.telo.tsv, LABEL.gaps.bed
"""
import gzip
import re
import sys

label, fasta, telo_f, telo_r = sys.argv[1:5]
window, tmin = int(sys.argv[5]), int(sys.argv[6])
telo_f, telo_r = telo_f.upper(), telo_r.upper()


def records(path):
    op = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    half, c = sum(lengths) / 2, 0
    for L in lengths:
        c += L
        if c >= half:
            return L
    return 0


scaf, ctg, ngaps = [], [], 0
gap_re = re.compile(r"[Nn]+")
with open(f"{label}.telo.tsv", "w") as telo, open(f"{label}.gaps.bed", "w") as bed:
    telo.write("seq\tlength\tcopies_start\tcopies_end\ttelo_start\ttelo_end\n")
    for name, s in records(fasta):
        L = len(s)
        scaf.append(L)
        prev = 0
        for m in gap_re.finditer(s):  # contigs = pieces between N runs
            bed.write(f"{name}\t{m.start()}\t{m.end()}\n")
            ngaps += 1
            if m.start() > prev:
                ctg.append(m.start() - prev)
            prev = m.end()
        if L > prev:
            ctg.append(L - prev)
        a, b = s[:window].upper(), s[-window:].upper()
        cs, ce = a.count(telo_f) + a.count(telo_r), b.count(telo_f) + b.count(telo_r)
        telo.write(f"{name}\t{L}\t{cs}\t{ce}\t{int(cs >= tmin)}\t{int(ce >= tmin)}\n")

size = sum(scaf)
big = [L for L in scaf if L >= 1_000_000]
t2t = one = 0
with open(f"{label}.telo.tsv") as fh:
    next(fh)
    for line in fh:
        f = line.split("\t")
        k = int(f[4]) + int(f[5])
        t2t += k == 2
        one += k == 1

cols = {
    "assembly": label, "size_bp": size, "scaffolds": len(scaf), "scaffold_N50": n50(scaf),
    "longest": max(scaf, default=0), "scaf_ge_1Mb": len(big),
    "pct_bp_in_ge_1Mb": f"{100 * sum(big) / size:.2f}" if size else "0.00",
    "contigs": len(ctg), "contig_N50": n50(ctg), "gaps": ngaps, "T2T": t2t, "one_telo": one,
}
with open(f"{label}.stats.tsv", "w") as out:
    out.write("\t".join(cols) + "\n" + "\t".join(str(v) for v in cols.values()) + "\n")
