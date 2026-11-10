#!/usr/bin/env python3
"""summary.tsv for one assembly run: per-assembly stats + compleasm + merqury.

run in a directory holding <name>.stats.tsv, <name>.compleasm.txt, haps.qv / primary.qv,
haps.completeness.stats / primary.completeness.stats.
Scaffolds (<name>.scaffolds) inherit BUSCO and QV of their contigs (same sequence).
"""
import glob
import os


def busco(path):
    out = {k: "NA" for k in "SDFM"}
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:  # e.g. "S:97.21%, 1384"
                if len(line) > 2 and line[0] in out and line[1] == ":":
                    out[line[0]] = line[2:].split(",")[0].strip()
    return out


def lookup(files, name, col):
    for path in files:
        if os.path.exists(path):
            with open(path) as fh:
                for line in fh:
                    f = line.split()
                    if f and f[0] == name and len(f) > col:
                        return f[col]
    return "NA"


qv_files = ["haps.qv", "primary.qv"]
cp_files = ["haps.completeness.stats", "primary.completeness.stats"]
names = sorted(p[: -len(".stats.tsv")] for p in glob.glob("*.stats.tsv"))

header, rows = None, []
for n in names:
    with open(f"{n}.stats.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        vals = fh.readline().rstrip("\n").split("\t")
    base = n[: -len(".scaffolds")] if n.endswith(".scaffolds") else n
    b = busco(f"{base}.compleasm.txt")
    vals += [b["S"], b["D"], b["F"], b["M"], lookup(qv_files, base, 3), lookup(cp_files, base, 4)]
    rows.append(vals)

with open("summary.tsv", "w") as out:
    out.write("\t".join(header + ["BUSCO_S", "BUSCO_D", "BUSCO_F", "BUSCO_M", "QV", "kmer_compl_pct"]) + "\n")
    for r in rows:
        out.write("\t".join(r) + "\n")

both = lookup(["haps.qv"], "Both", 3)
if both != "NA":
    print(f"Combined (hap1+hap2) merqury QV: {both}")
