#!/usr/bin/env python3
"""Dot plot + alignment coverage from a minimap2 PAF (query vs target assembly).

usage: dotplot.py NAME PAF QUERY.len TARGET.len MIN_ALN MIN_SEQ QLABEL TLABEL

  *.len    : two columns, sequence name and length (all sequences of the assembly)
  MIN_ALN  : alignment blocks shorter than this are ignored (plot and coverage)
  MIN_SEQ  : sequences shorter than this are not drawn (coverage still uses all)

writes NAME.dotplot.png, NAME.dotplot.json (layout for the interactive html report) and NAME.dotplot.tsv (bp of query / target covered by alignments,
bp in inverted blocks relative to the main orientation of each query sequence)
Targets are sorted by length; each query is placed next to the target it aligns to best and
flipped when most of it aligns on the minus strand, so collinear sequence shows as a diagonal.
"""
import gzip
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

name, paf, qlen_f, tlen_f = sys.argv[1:5]
min_aln, min_seq = int(sys.argv[5]), int(sys.argv[6])
qlabel, tlabel = sys.argv[7], sys.argv[8]


def read_len(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                out[f[0]] = int(f[1])
    return out


qlen, tlen = read_len(qlen_f), read_len(tlen_f)
alns = []  # qname, qs, qe, strand, tname, ts, te
op = gzip.open if paf.endswith(".gz") else open
with op(paf, "rt") as fh:
    for line in fh:
        f = line.split("\t")
        if len(f) < 12 or int(f[10]) < min_aln:
            continue
        alns.append((f[0], int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8])))


def merged_bp(iv):
    tot, end = 0, -1
    for a, b in sorted(iv):
        if a > end:
            tot += b - a
            end = b
        elif b > end:
            tot += b - end
            end = b
    return tot


# per query: best target, strand balance, weighted position on the best target
best = {}
for q, qs, qe, s, t, ts, te in alns:
    d = best.setdefault(q, {"t": {}, "plus": 0, "minus": 0})
    w = qe - qs
    d["t"][t] = d["t"].get(t, 0) + w
    d["plus" if s == "+" else "minus"] += w
flip = {q: d["minus"] > d["plus"] for q, d in best.items()}

# coverage table
qiv, tiv, inv = {}, {}, 0
for q, qs, qe, s, t, ts, te in alns:
    qiv.setdefault(q, []).append((qs, qe))
    tiv.setdefault(t, []).append((ts, te))
    if (s == "-") != flip[q]:
        inv += qe - qs
qcov = sum(merged_bp(v) for v in qiv.values())
tcov = sum(merged_bp(v) for v in tiv.values())
qtot, ttot = sum(qlen.values()), sum(tlen.values())
pct = lambda a, b: f"{100 * a / b:.2f}" if b else "NA"
with open(f"{name}.dotplot.tsv", "w") as out:
    out.write("query\ttarget\tquery_bp\ttarget_bp\tquery_aligned_bp\tpct_query_aligned\t"
              "target_aligned_bp\tpct_target_aligned\tinverted_bp\talignments\n")
    out.write(f"{qlabel}\t{tlabel}\t{qtot}\t{ttot}\t{qcov}\t{pct(qcov, qtot)}\t"
              f"{tcov}\t{pct(tcov, ttot)}\t{inv}\t{len(alns)}\n")

# layout
tseq = sorted((t for t in tlen if tlen[t] >= min_seq), key=lambda t: -tlen[t])
tidx = {t: i for i, t in enumerate(tseq)}
toff, c = {}, 0
for t in tseq:
    toff[t] = c
    c += tlen[t]
tsum = c

qpos = {}
for q, qs, qe, s, t, ts, te in alns:
    if t in tidx and q in qlen and qlen[q] >= min_seq:
        bt = max(best[q]["t"], key=best[q]["t"].get)
        if t == bt:
            p = qpos.setdefault(q, [0, 0])
            p[0] += (ts + te) / 2 * (qe - qs)
            p[1] += qe - qs
qseq = []
for q, (sw, w) in qpos.items():
    bt = max(best[q]["t"], key=best[q]["t"].get)
    qseq.append((tidx.get(bt, len(tseq)), sw / w, q))
qseq = [q for *_, q in sorted(qseq)]
qoff, c = {}, 0
for q in qseq:
    qoff[q] = c
    c += qlen[q]
qsum = c

segs, cols = [], []
for q, qs, qe, s, t, ts, te in alns:
    if q not in qoff or t not in toff:
        continue
    L = qlen[q]
    pts = [(ts, qs), (te, qe)] if s == "+" else [(ts, qe), (te, qs)]
    fwd = (s == "+") != flip[q]
    segs.append([((toff[t] + x) / 1e6, (qoff[q] + (L - y if flip[q] else y)) / 1e6) for x, y in pts])
    cols.append("#1f5fa8" if fwd else "#c0392b")

# layout for the interactive report (longest alignments first, capped)
MAX_JSON = 40000
ti = {t: i for i, t in enumerate(tseq)}
qi = {q: i for i, q in enumerate(qseq)}
jsegs = sorted(([ti[t], qi[q], ts, te, qs, qe, int(s == "+")] for q, qs, qe, s, t, ts, te in alns
                if q in qi and t in ti), key=lambda r: r[5] - r[4], reverse=True)
with open(f"{name}.dotplot.json", "w") as out:
    json.dump({"q": qlabel, "t": tlabel, "min_aln": min_aln, "min_seq": min_seq,
               "qpct": pct(qcov, qtot), "tpct": pct(tcov, ttot), "shown": min(len(jsegs), MAX_JSON),
               "total": len(jsegs),
               "tseq": [[t, toff[t], tlen[t]] for t in tseq],
               "qseq": [[q, qoff[q], qlen[q], int(flip[q])] for q in qseq],
               "segs": jsegs[:MAX_JSON]}, out, separators=(",", ":"))

fig, ax = plt.subplots(figsize=(9, 9))
ax.add_collection(LineCollection(segs, colors=cols, linewidths=1.2))
for t in tseq[1:]:
    ax.axvline(toff[t] / 1e6, color="0.85", lw=0.5, zorder=0)
for q in qseq[1:]:
    ax.axhline(qoff[q] / 1e6, color="0.85", lw=0.5, zorder=0)
ax.set_xlim(0, max(tsum, 1) / 1e6)
ax.set_ylim(0, max(qsum, 1) / 1e6)
ax.set_xlabel(f"{tlabel} (Mb)")
ax.set_ylabel(f"{qlabel} (Mb)")
ax.set_title(f"{qlabel} vs {tlabel}\n{pct(qcov, qtot)}% of query / {pct(tcov, ttot)}% of target in "
             f"alignments >= {min_aln // 1000} kb (blue +, red -)", fontsize=10)
# sequence names on the top / right axes when they fit
if 0 < len(tseq) <= 40:
    top = ax.secondary_xaxis("top")
    top.set_xticks([(toff[t] + tlen[t] / 2) / 1e6 for t in tseq])
    top.set_xticklabels(tseq, rotation=90, fontsize=6)
if 0 < len(qseq) <= 40:
    right = ax.secondary_yaxis("right")
    right.set_yticks([(qoff[q] + qlen[q] / 2) / 1e6 for q in qseq])
    right.set_yticklabels([q + (" (-)" if flip[q] else "") for q in qseq], fontsize=6)
if not segs:
    ax.text(0.5, 0.5, "no alignments to show", ha="center", va="center", transform=ax.transAxes)
fig.tight_layout()
fig.savefig(f"{name}.dotplot.png", dpi=150)
print(f"{name}: {len(alns)} alignments, {pct(qcov, qtot)}% query / {pct(tcov, ttot)}% target aligned")
