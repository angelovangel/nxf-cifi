#!/usr/bin/env python3
"""IGV desktop session for one haplotype assembly (File > Open Session...).

usage: igv_session.py NAME FASTA TELO_TSV GAPS_BED AGP|none ONT(0/1) TELO_MOTIF_LEN

writes NAME.igv.xml for <outdir>/ and the tracks for <outdir>/igv/ (session paths are relative to <outdir> and never
use ../, which IGV does not resolve for every track; the output folder can be moved):
  NAME.igv.xml          genome fasta/NAME.fa + igv/<the tracks below> (+ ont_mapping/NAME.ont.bam with ONT=1)
  NAME.contigs.bed      contigs from the YaHS AGP, alternating colours, name = contig (+/-); without an AGP the
                        pieces between N gaps (seq.1, seq.2, ...; no names or orientation), if there are gaps
  NAME.gaps.bed         N runs (from asm_stats.py)
  NAME.telomeres.bed    telomere repeat at sequence ends (asm_stats.py calls), length ~ copies x motif
  NAME.fa.fai           scaffolds only (contig FASTAs already have one in fasta/): IGV needs it next to the FASTA
with ONT=1 the session also loads ont_mapping/NAME.ont.coverage.wig (1 kb bins), NAME.ont.gaps.bed (ONT support
of each join, from ont_gaps.sh) and NAME.ont.bam
"""
import os
import sys
from xml.sax.saxutils import quoteattr

name, fasta, telo, gaps, agp, ont, mlen = sys.argv[1:8]
ont, mlen = ont == "1", int(mlen)


def write_fai(path, out):
    """samtools faidx equivalent (uniform line length per sequence, as written by yahs / seqkit)"""
    rows, off, cur = [], 0, None
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if cur: rows.append(cur)
                cur = [line[1:].split()[0].decode(), 0, off + len(line), 0, 0]
            elif cur:
                b = len(line.rstrip(b"\r\n"))
                if not cur[3]: cur[3], cur[4] = b, len(line)
                cur[1] += b
            off += len(line)
    if cur: rows.append(cur)
    with open(out, "w") as fh:
        for r in rows: fh.write("\t".join(map(str, r)) + "\n")
    return {r[0]: r[1] for r in rows}


lengths = write_fai(fasta, f"{name}.fa.fai")
if not name.endswith(".scaffolds"): os.remove(f"{name}.fa.fai")    # contig FASTAs: fasta/ already has the index
gap_rows = [l.split("\t")[:3] for l in open(gaps) if l.strip() and not l.startswith(("track", "#"))]

tracks = []
if agp != "none":
    with open(agp) as fh, open(f"{name}.contigs.bed", "w") as out:
        out.write('track name="contigs (YaHS AGP)" itemRgb="On"\n')
        k = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 9 or f[4] in ("N", "U"): continue
            i = k[f[0]] = k.get(f[0], -1) + 1                     # contig index within its scaffold
            rgb = "63,120,181" if i % 2 == 0 else "156,195,230"   # as in the report's scaffold layout
            out.write(f"{f[0]}\t{int(f[1]) - 1}\t{f[2]}\t{f[5]} ({f[8]})\t0\t{f[8] if f[8] in '+-' else '.'}\t"
                      f"{int(f[1]) - 1}\t{f[2]}\t{rgb}\n")
    tracks.append(f"{name}.contigs.bed")
elif gap_rows:
    # no AGP: the contigs are the pieces between N gaps
    by = {}
    for c, a, b in gap_rows: by.setdefault(c, []).append((int(a), int(b)))
    with open(f"{name}.contigs.bed", "w") as out:
        out.write('track name="contigs (split at N gaps)" itemRgb="On"\n')
        for c, L in lengths.items():
            p, k = 0, 0
            for a, b in sorted(by.get(c, [])) + [(L, L)]:
                if a > p:
                    out.write(f"{c}\t{p}\t{a}\t{c}.{k + 1}\t0\t.\t{p}\t{a}\t{'63,120,181' if k % 2 == 0 else '156,195,230'}\n")
                    k += 1
                p = b
    tracks.append(f"{name}.contigs.bed")

with open(gaps) as fh, open(f"{name}.gaps.bed", "w") as out:
    out.write('track name="gaps (N)" color=40,40,40\n')
    for line in fh:
        if line.strip() and not line.startswith(("track", "#")): out.write(line)
tracks.append(f"{name}.gaps.bed")

with open(telo) as fh, open(f"{name}.telomeres.bed", "w") as out:
    out.write('track name="telomeres" color=46,154,82\n')
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")                       # seq length copies_start copies_end telo_start telo_end
        if len(f) < 6: continue
        L = int(f[1])
        if f[4] == "1": out.write(f"{f[0]}\t0\t{min(L, int(f[2]) * mlen)}\ttelomere ({f[2]} copies)\n")
        if f[5] == "1": out.write(f"{f[0]}\t{max(0, L - int(f[3]) * mlen)}\t{L}\ttelomere ({f[3]} copies)\n")
tracks.append(f"{name}.telomeres.bed")

res = [f"    <Resource path={quoteattr('igv/' + t)} type=\"bed\"/>" for t in tracks]
if ont:
    bam = f"ont_mapping/{name}.ont.bam"
    res.append(f"    <Resource path={quoteattr(f'ont_mapping/{name}.ont.coverage.wig')} type=\"wig\"/>")
    res.append(f"    <Resource path={quoteattr(f'ont_mapping/{name}.ont.gaps.bed')} type=\"bed\"/>")
    res.append(f"    <Resource path={quoteattr(bam)} index={quoteattr(bam + '.csi')} type=\"bam\"/>")
# the first sequence opens; the genome is the published FASTA (relative to this session file)
first = next(iter(lengths), None)
locus = f' locus="{first}"' if first else ""
with open(f"{name}.igv.xml", "w") as out:
    out.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
              f'<Session genome={quoteattr("fasta/" + name + ".fa")} hasGeneTrack="false" hasSequenceTrack="true"'
              f'{locus} version="8">\n  <Resources>\n' + "\n".join(res) + "\n  </Resources>\n</Session>\n")
