#!/usr/bin/env python3
"""CiFi contact analysis of one scaffolded haplotype (after McGinty et al., Dennis lab, Nat Commun 2026).

usage: cifi_stats.py NAME SEGMENTS.tsv.gz SCAFFOLDS.agp MAPQ MAX_PAIRS MIN_SCAF BIN

  SEGMENTS.tsv.gz : primary alignments of the cifi digest segments to the contigs, one per line:
                    qname flag contig pos(1-based) mapq ref_span  (qname = <read>__CIFI_SEG__<index>)
  SCAFFOLDS.agp   : YaHS AGP, used to lift contig coordinates to scaffold coordinates
  MAPQ            : minimum segment MAPQ to take part in contacts (as cifi contacts -q)
  MAX_PAIRS       : pairwise contacts analysed (distance decay, cis/trans); reads are taken in file order
                    until this many pairs are reached, segment-level statistics use all reads
  MIN_SCAF        : scaffolds >= this length count as chromosomes for cis/trans
  BIN             : coverage bin size (bp)

Per read, every pair of usable segments is a contact (k segments -> k(k-1)/2), as in cifi contacts and
the paper. Adjacent segments (consecutive along the CiFi read) are direct ligation junctions; non-adjacent
pairs are the extra, multi-way contacts that long reads add.

writes NAME.cifi_stats.tsv (summary row), NAME.cifi_scaffolds.tsv (per chromosome-scale scaffold),
       NAME.cifi_stats.json (summary + curves for the html report), NAME.cifi_stats.png
"""
import gzip
import json
import math
import sys
from bisect import bisect_right

name, seg_tsv, agp_file = sys.argv[1:4]
min_q, max_pairs, min_scaf, bin_size = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
SEP = "__CIFI_SEG__"
MAPQ_LEVELS = [1, 10, 20, 30, 60]
RANGES = [(0, 1e3, "<1 kb"), (1e3, 1e4, "1-10 kb"), (1e4, 1e5, "10-100 kb"), (1e5, 1e6, "0.1-1 Mb"),
          (1e6, float("inf"), ">=1 Mb")]
PS_PER_DECADE, PS_MIN = 8, 100          # distance-decay bins: 8 per decade from 100 bp

# ------------------------------------------------------------------ AGP: contig -> scaffold coordinates
agp, scaf_len, order = {}, {}, []
with open(agp_file) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        obj, obeg, oend = f[0], int(f[1]), int(f[2])
        if obj not in scaf_len:
            order.append(obj)
        scaf_len[obj] = max(scaf_len.get(obj, 0), oend)
        if f[4] in ("N", "U"):
            continue
        agp.setdefault(f[5], []).append((int(f[6]), int(f[7]), obj, obeg, f[8]))
for v in agp.values():
    v.sort()
agp_starts = {c: [x[0] for x in v] for c, v in agp.items()}


def lift(ctg, p):
    v = agp.get(ctg)
    if not v:
        return None, None
    i = bisect_right(agp_starts[ctg], p) - 1
    if i < 0:
        return None, None
    cb, ce, s, sb, o = v[i]
    if p > ce:
        return None, None
    return s, (sb + (ce - p) if o == "-" else sb + (p - cb))


big = {s for s, L in scaf_len.items() if L >= min_scaf}
big_order = sorted(big, key=lambda s: -scaf_len[s])

# ------------------------------------------------------------------ counters
n_seg = n_mapped = 0
mapq_ge = dict.fromkeys(MAPQ_LEVELS, 0)
n_reads = 0
card = {}                    # usable segments per read -> reads
nscaf_hist = {}              # distinct chromosome-scale scaffolds per read (reads with >= 2 usable) -> reads
contacts_total = 0
pairs = cis = trans = other = 0
adj = {"cis": 0, "trans": 0}
nonadj = {"cis": 0, "trans": 0}
range_counts = [0] * len(RANGES)
n_ps = int(math.ceil(math.log10(max(max(scaf_len.values(), default=1), 10) / PS_MIN) * PS_PER_DECADE)) + 1
ps_all, ps_adj, ps_non = [0] * n_ps, [0] * n_ps, [0] * n_ps
per_scaf = {s: [0, 0, 0] for s in big}       # segments, cis contacts, trans contacts
cov = {s: [0] * (scaf_len[s] // bin_size + 1) for s in big}
pair_reads = 0


def ps_bin(d):
    if d < PS_MIN:
        return 0
    return min(n_ps - 1, int(math.log10(d / PS_MIN) * PS_PER_DECADE))


def flush(segs):
    """segs: [(index, scaffold, pos)] usable segments of one CiFi read"""
    global contacts_total, pairs, cis, trans, other, pair_reads
    k = len(segs)
    card[k] = card.get(k, 0) + 1
    if k < 2:
        return
    contacts_total += k * (k - 1) // 2
    nb = len({s for _, s, _ in segs if s in big})
    nscaf_hist[nb] = nscaf_hist.get(nb, 0) + 1
    if pairs >= max_pairs:
        return
    pair_reads += 1
    segs.sort()
    for i in range(k):
        _, si, pi = segs[i]
        for j in range(i + 1, k):
            _, sj, pj = segs[j]
            pairs += 1
            a = j == i + 1
            if si == sj:
                cis += 1
                d = abs(pi - pj)
                b = ps_bin(d)
                ps_all[b] += 1
                (ps_adj if a else ps_non)[b] += 1
                (adj if a else nonadj)["cis"] += 1
                for r, (lo, hi, _) in enumerate(RANGES):
                    if lo <= d < hi:
                        range_counts[r] += 1
                        break
                if si in per_scaf:
                    per_scaf[si][1] += 1
            elif si in big and sj in big:
                trans += 1
                (adj if a else nonadj)["trans"] += 1
                per_scaf[si][2] += 1
                per_scaf[sj][2] += 1
            else:
                other += 1


op = gzip.open if seg_tsv.endswith(".gz") else open
cur, segs = None, []
with op(seg_tsv, "rt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        rid, _, idx = f[0].rpartition(SEP)
        if not rid:
            rid, idx = f[0], "0"
        if rid != cur:
            if cur is not None:
                flush(segs)
                n_reads += 1
            cur, segs = rid, []
        n_seg += 1
        flag, ctg, mq = int(f[1]), f[2], int(f[4])
        if flag & 4 or ctg == "*":
            continue
        n_mapped += 1
        for q in MAPQ_LEVELS:
            if mq >= q:
                mapq_ge[q] += 1
        if mq < min_q:
            continue
        s, p = lift(ctg, int(f[3]) + int(f[5]) // 2)
        if s is None:
            continue
        segs.append((int(idx), s, p))
        if s in big:
            per_scaf[s][0] += 1
            cov[s][p // bin_size] += 1
if cur is not None:
    flush(segs)
    n_reads += 1

# ------------------------------------------------------------------ summary
pct = lambda a, b: round(100 * a / b, 2) if b else None
reads_ge2 = sum(v for k, v in card.items() if k >= 2)
reads_ge3 = sum(v for k, v in card.items() if k >= 3)


def median_from_hist(h, keep=lambda k: True):
    items = sorted((k, v) for k, v in h.items() if keep(k))
    tot = sum(v for _, v in items)
    if not tot:
        return None
    c = 0
    for k, v in items:
        c += v
        if c >= tot / 2:
            return k


def ps_curve(counts):
    tot = sum(counts)
    out = []
    for b, c in enumerate(counts):
        if not c:
            continue
        lo, hi = PS_MIN * 10 ** (b / PS_PER_DECADE), PS_MIN * 10 ** ((b + 1) / PS_PER_DECADE)
        if b == 0:
            lo = 1
        out.append([round(math.sqrt(lo * hi)), c / tot / (hi - lo)])   # contact probability per bp
    return out


def slope(curve, lo=5e4, hi=5e6):
    pts = [(math.log10(s), math.log10(p)) for s, p in curve if lo <= s <= hi and p > 0]
    if len(pts) < 3:
        return None
    mx, my = sum(x for x, _ in pts) / len(pts), sum(y for _, y in pts) / len(pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    return round(sum((x - mx) * (y - my) for x, y in pts) / sxx, 3) if sxx else None


curves = {"all": ps_curve(ps_all), "adjacent": ps_curve(ps_adj), "non_adjacent": ps_curve(ps_non)}
ct = cis + trans
summary = {
    "haplotype": name,
    "cifi_reads": n_reads,
    "segments": n_seg,
    "pct_segments_mapped": pct(n_mapped, n_seg),
    **{f"pct_segments_mapq_ge{q}": pct(mapq_ge[q], n_seg) for q in MAPQ_LEVELS},
    "pct_reads_ge2_usable": pct(reads_ge2, n_reads),
    "median_usable_segments": median_from_hist(card, lambda k: k >= 1),
    "pct_multiway_reads": pct(reads_ge3, reads_ge2),
    "contacts": contacts_total,
    "contacts_per_read": round(contacts_total / reads_ge2, 2) if reads_ge2 else None,
    "pairs_analysed": pairs,
    "pct_cis": pct(cis, ct),
    "pct_trans": pct(trans, ct),
    "pct_cis_ge10kb": pct(sum(range_counts[2:]), ct),
    "pct_cis_ge1Mb": pct(range_counts[4], ct),
    "pct_trans_adjacent": pct(adj["trans"], adj["cis"] + adj["trans"]),
    "pct_trans_non_adjacent": pct(nonadj["trans"], nonadj["cis"] + nonadj["trans"]),
    "pct_reads_multi_chrom": pct(sum(v for k, v in nscaf_hist.items() if k >= 2), reads_ge2),
    "ps_slope_50kb_5Mb": slope(curves["all"]),
    "chrom_scaffolds": len(big),
    "min_mapq": min_q,
}
with open(f"{name}.cifi_stats.tsv", "w") as out:
    out.write("\t".join(summary) + "\n" + "\t".join("NA" if v is None else str(v) for v in summary.values()) + "\n")

scaf_rows = []
for s in big_order:
    seg, c, t = per_scaf[s]
    scaf_rows.append([s, scaf_len[s], seg, round(seg / scaf_len[s] * 1e6, 1), c, t, pct(c, c + t)])
with open(f"{name}.cifi_scaffolds.tsv", "w") as out:
    out.write("scaffold\tlength\tsegments\tsegments_per_Mb\tcis_contacts\ttrans_contacts\tpct_cis\n")
    for r in scaf_rows:
        out.write("\t".join("NA" if v is None else str(v) for v in r) + "\n")

cmax = max(card, default=0)
json.dump({
    "summary": summary,
    "mapq": {str(q): pct(mapq_ge[q], n_seg) for q in MAPQ_LEVELS},
    "cardinality": [[k, card.get(k, 0)] for k in range(1, cmax + 1)],
    "chroms_per_read": sorted([k, v] for k, v in nscaf_hist.items()),
    "ranges": [[lab, pct(c, ct)] for (_, _, lab), c in zip(RANGES, range_counts)],
    "ps": curves,
    "scaffolds": scaf_rows,
    "pair_reads": pair_reads,
    "max_pairs": max_pairs,
    "min_scaf": min_scaf,
}, open(f"{name}.cifi_stats.json", "w"), separators=(",", ":"))

# ------------------------------------------------------------------ overview figure
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

fig = plt.figure(figsize=(14, 9), constrained_layout=True)
gs = fig.add_gridspec(2, 3)
ax = fig.add_subplot(gs[0, 0])
ks = list(range(1, min(cmax, 40) + 1))
ax.bar(ks, [card.get(k, 0) for k in ks], color="#2f6fb3")
ax.set_title(f"Usable segments per CiFi read (MAPQ >= {min_q})", fontsize=10)
ax.set_xlabel("segments"); ax.set_ylabel("reads")
ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
ax = fig.add_subplot(gs[0, 1])
for key, col, lab in (("all", "#1d1f21", "all pairs"), ("adjacent", "#2f6fb3", "adjacent (direct ligation)"),
                      ("non_adjacent", "#e07b24", "non-adjacent (multi-way)")):
    c = curves[key]
    if c:
        ax.loglog([x for x, _ in c], [y for _, y in c], "-o", ms=2.5, lw=1.4, color=col, label=lab)
ax.set_title(f"Contact probability P(s), slope 50 kb-5 Mb: {summary['ps_slope_50kb_5Mb']}", fontsize=10)
ax.set_xlabel("distance s (bp)"); ax.set_ylabel("P(s) per bp"); ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
ax = fig.add_subplot(gs[0, 2])
labs = [f"MAPQ>={q}" for q in MAPQ_LEVELS] + ["cis", "trans"] + [r[2] for r in RANGES]
vals = [summary[f"pct_segments_mapq_ge{q}"] for q in MAPQ_LEVELS] + [summary["pct_cis"], summary["pct_trans"]] + \
       [pct(c, ct) for c in range_counts]
ax.barh(range(len(labs)), [v or 0 for v in vals], color=["#7d5bb5"] * 5 + ["#2e9a52", "#c8413b"] + ["#8a8f96"] * 5)
ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=8); ax.invert_yaxis(); ax.set_xlim(0, 105)
for i, v in enumerate(vals):
    ax.text((v or 0) + 1, i, "NA" if v is None else f"{v:.1f}", va="center", fontsize=7)
ax.set_title("% of segments (MAPQ) and of contacts", fontsize=10)
ax = fig.add_subplot(gs[1, :])
x0, ticks, labels = 0, [], []
for i, s in enumerate(big_order):
    c = list(cov[s])
    last = scaf_len[s] - (len(c) - 1) * bin_size          # partial last bin, scaled to a full bin
    c[-1] = c[-1] * bin_size / last if last > bin_size / 10 else None
    ax.plot([x0 + (j + 0.5) * bin_size / 1e6 for j in range(len(c))], [v if v is not None else float("nan") for v in c], lw=0.6,
            color="#2f6fb3" if i % 2 == 0 else "#7fa6d6")
    ticks.append(x0 + scaf_len[s] / 2e6); labels.append(s); x0 += scaf_len[s] / 1e6
if 0 < len(big_order) <= 40:
    ax.set_xticks(ticks); ax.set_xticklabels(labels, rotation=90, fontsize=7)
ax.set_xlim(0, max(x0, 1)); ax.set_ylabel(f"segments / {bin_size // 1000} kb")
ax.set_title(f"CiFi segment coverage along scaffolds >= {min_scaf / 1e6:g} Mb (gaps and dips: low mappability, "
             "collapsed or missing sequence)", fontsize=10)
fig.suptitle(f"{name}: CiFi contacts", fontsize=12)
fig.savefig(f"{name}.cifi_stats.png", dpi=130)
print(f"{name}: {n_reads:,} reads, {n_seg:,} segments, {contacts_total:,} contacts, "
      f"cis {summary['pct_cis']}%, P(s) slope {summary['ps_slope_50kb_5Mb']}")
