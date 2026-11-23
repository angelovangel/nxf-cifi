#!/usr/bin/env bash
#
# Side-by-side comparison of nxf-cifi runs, e.g. HiFi+ONT (no CiFi) vs HiFi+ONT+CiFi
# (hifiasm --h1/--h2 phasing) vs YaHS scaffolds. Needs only bash + python3: it collects the results the
# pipeline already computed and runs no assembly or evaluation tools.
#
# Per assembly (from <run>/eval/...)
#   - contiguity : size, scaffolds, scaffold N50, longest, scaffolds >= 1 Mb (+ % bp), contigs
#                  (split at N gaps), contig N50, gaps                         [eval/stats]
#   - telomeres  : sequences with telomeres at both ends (T2T) / one end       [eval/stats]
#   - BUSCO      : compleasm S/D/F/M                                           [eval/compleasm]
#   - accuracy   : merqury QV + k-mer completeness                             [eval/merqury]
#   Scaffolds inherit BUSCO and QV of their contigs (same sequence).
# Per run (hap1 + hap2)
#   - combined merqury QV / k-mer completeness and the spectra-cn plot
# Plots (python3 + matplotlib; skipped with -P or if matplotlib is missing)
#   - nx_curves.png  : scaffold and contig Nx curves of all assemblies
#   - comparison.png : scaffold/contig N50, % bp in >= 1 Mb, T2T, BUSCO complete, QV
# HTML report (report.html, self-contained), one section per topic (runs as subsections): assemblies,
#   comparison plots, merqury spectra, dot plots, CiFi (library, contacts, per-haplotype contact maps +
#   analysis), ONT coverage; reads (faster stats, coverage, QC links) first; commands used (expandable)
#
# Usage
#   bash compare_assemblies.sh -d noCiFi=asm_noCiFi -d CiFi=asm_CiFi -d CiFiL0=asm_CiFi_l0 -o compare_3runs
#
# Options
#   -d  label=run_dir  nxf-cifi output dir (repeat per run); adds <label>_hap1/_hap2 (contigs) and,
#       if scaffolded, <label>Scaf_hap1/_hap2; with -p also <label>_primary
#   -a  label=fasta    any other assembly: contiguity + telomeres computed here, BUSCO/QV = NA
#   -p  include the primary assembly (*.p_ctg.fa) of each run
#   -o  output directory                          [compare_out]
#   -t  report title                              [Assembly comparison]
#   -m  telomere motif, for -a assemblies         [TTTAGGG]
#   -w  bp examined at each end, for -a           [10000]
#   -n  minimum motif copies = telomere, for -a   [30]
#   -P  no plots
#
# Output: <outdir>/report.html, comparison.tsv, hap_pairs.tsv, nx_curves.png, comparison.png

set -euo pipefail

SPECS=(); OUTDIR="compare_out"; PRIMARY=0; PLOTS=1; TITLE="Assembly comparison"
TELO="TTTAGGG"; WINDOW=10000; TMIN=30

usage() { awk 'NR>1 && !/^#/{exit} NR>1{sub(/^# ?/,""); print}' "$0"; exit 1; }
while getopts "d:a:po:t:m:w:n:Ph" opt; do
  case $opt in
    d|a) [[ "$OPTARG" == *=* ]] || { echo "ERROR: -$opt needs label=path"; exit 1; }
         SPECS+=("$opt=$OPTARG") ;;
    p) PRIMARY=1 ;;  o) OUTDIR=$OPTARG ;;  P) PLOTS=0 ;;  t) TITLE=$OPTARG ;;
    m) TELO=$OPTARG ;; w) WINDOW=$OPTARG ;; n) TMIN=$OPTARG ;;
    *) usage ;;
  esac
done
(( ${#SPECS[@]} >= 1 )) || { echo "ERROR: give at least one -d label=run_dir or -a label=fasta"; usage; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
mkdir -p "$OUTDIR"

export REPORT_CMD="$(basename "$0") $(printf '%q ' "$@")"
python3 - "$OUTDIR" "$PRIMARY" "$TELO" "$WINDOW" "$TMIN" "$PLOTS" "$TITLE" "${SPECS[@]}" <<'PY'
import glob, gzip, os, re, sys

outdir, primary, telo, window, tmin = sys.argv[1], sys.argv[2] == "1", sys.argv[3].upper(), int(sys.argv[4]), int(sys.argv[5])
plots, title, specs = sys.argv[6] == "1", sys.argv[7], sys.argv[8:]
rc = telo[::-1].translate(str.maketrans("ACGT", "TGCA"))
STATS = ["size_bp", "scaffolds", "scaffold_N50", "longest", "scaf_ge_1Mb", "pct_bp_in_ge_1Mb",
         "contigs", "contig_N50", "gaps", "T2T", "one_telo"]

def die(msg): sys.exit(f"ERROR: {msg}")

def one(pattern, what):
    hits = sorted(glob.glob(pattern))
    if len(hits) > 1: die(f"several {what}: {' '.join(hits)}; use -a")
    return hits[0] if hits else None

# ------------------------------------------------------------------ fallback: stats from a fasta
def n50(ls):
    ls, c = sorted(ls, reverse=True), 0
    for L in ls:
        c += L
        if c >= sum(ls) / 2: return L
    return 0

def fasta_stats(path):
    op = gzip.open if path.endswith(".gz") else open
    scaf, ctg, gaps, t2t, onet = [], [], 0, 0, 0
    def add(s):
        nonlocal gaps, t2t, onet
        scaf.append(len(s)); prev = 0
        for m in re.finditer(r"[Nn]+", s):
            gaps += 1
            if m.start() > prev: ctg.append(m.start() - prev)
            prev = m.end()
        if len(s) > prev: ctg.append(len(s) - prev)
        a, b = s[:window].upper(), s[-window:].upper()
        k = (a.count(telo) + a.count(rc) >= tmin) + (b.count(telo) + b.count(rc) >= tmin)
        t2t += k == 2; onet += k == 1
    chunks, started = [], False
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if started: add("".join(chunks))
                chunks, started = [], True
            else:
                chunks.append(line.strip())
    if started: add("".join(chunks))
    size, big = sum(scaf), [L for L in scaf if L >= 1e6]
    LENGTHS[path] = (scaf, ctg)
    return dict(zip(STATS, map(str, [size, len(scaf), n50(scaf), max(scaf, default=0), len(big),
                f"{100 * sum(big) / size:.2f}" if size else "0.00", len(ctg), n50(ctg), gaps, t2t, onet])))

LENGTHS = {}   # fasta -> (scaffold lengths, contig lengths), for the Nx plot

def lengths(fa, gaps_bed):
    """scaffold + contig lengths; uses .fai or telo.tsv + gaps.bed from the pipeline when present"""
    if fa in LENGTHS: return LENGTHS[fa]
    telo_tsv = os.path.join(os.path.dirname(os.path.dirname(gaps_bed)), "telomeres",
                            os.path.basename(gaps_bed).replace(".gaps.bed", ".telo.tsv")) if gaps_bed else ""
    size_file = fa + ".fai" if os.path.exists(fa + ".fai") else telo_tsv
    if gaps_bed and os.path.exists(size_file) and os.path.exists(gaps_bed):
        scaf = {}
        for line in open(size_file):
            f = line.split("\t")
            if len(f) > 1 and f[1].strip().isdigit(): scaf[f[0]] = int(f[1])
        gaps = {}
        for line in open(gaps_bed):
            f = line.split("\t"); gaps.setdefault(f[0], []).append((int(f[1]), int(f[2])))
        ctg = []
        for n, L in scaf.items():
            prev = 0
            for a, b in sorted(gaps.get(n, [])):
                if a > prev: ctg.append(a - prev)
                prev = b
            if L > prev: ctg.append(L - prev)
        LENGTHS[fa] = (list(scaf.values()), ctg)
    else:
        fasta_stats(fa)
    return LENGTHS[fa]

# ------------------------------------------------------------------ readers for pipeline results
def read_stats(path):
    with open(path) as fh:
        h, v = fh.readline().rstrip("\n").split("\t"), fh.readline().rstrip("\n").split("\t")
    d = dict(zip(h, v))
    return {k: d.get(k, "NA") for k in STATS}

def busco(path):
    out = {k: "NA" for k in "SDFM"}
    if path and os.path.exists(path):
        for line in open(path):
            if len(line) > 2 and line[0] in out and line[1] == ":":
                out[line[0]] = line[2:].split(",")[0].strip()
    return out

def field(paths, key, col):
    for p in paths:
        if os.path.exists(p):
            for line in open(p):
                f = line.split()
                if f and f[0] == key and len(f) > col: return f[col]
    return "NA"

# ------------------------------------------------------------------ collect assemblies
rows, pairs, labels, sources, RUNS = [], [], set(), [], []
def add_row(label, stats, b, qv, comp, fa=None, gaps_bed=""):
    sources.append((label, fa, gaps_bed))
    if label in labels: die(f"label '{label}' used twice")
    if not re.fullmatch(r"[A-Za-z0-9_]+", label): die(f"bad label '{label}' (letters, digits, _ only)")
    labels.add(label)
    rows.append([label] + [stats[k] for k in STATS] + [b[k] for k in "SDFM"] + [qv, comp])

for spec in specs:
    kind, rest = spec.split("=", 1)
    label, path = rest.split("=", 1)
    path = path.rstrip("/")
    if kind == "a":
        if not os.path.isfile(path): die(f"-a {label}: {path} not found")
        print(f"[{label}] computing contiguity + telomeres from {path}")
        add_row(label, fasta_stats(path), busco(None), "NA", "NA", path)
        continue

    fdir = os.path.join(path, "fasta")
    if not os.path.isdir(fdir): die(f"-d {label}: no {fdir} (not an nxf-cifi output dir?)")
    ev = os.path.join(path, "eval")
    qvs = [os.path.join(ev, "merqury", f) for f in ("haps.qv", "primary.qv")]
    cps = [os.path.join(ev, "merqury", f) for f in ("haps.completeness.stats", "primary.completeness.stats")]
    wanted = [("hap1", "p_ctg"), ("hap2", "p_ctg"), ("hap1", "p_ctg.scaffolds"), ("hap2", "p_ctg.scaffolds")]
    found = []
    if primary:
        fa = next((f for f in sorted(glob.glob(os.path.join(fdir, "*.p_ctg.fa")))
                   if not re.search(r"\.hap\d\.p_ctg\.fa$", f)), None)
        if fa: found.append((f"{label}_primary", fa))
    for h, kind_ in wanted:
        fa = one(os.path.join(fdir, f"*.{h}.{kind_}.fa"), f"*.{h}.{kind_}.fa in {fdir}")
        if fa: found.append(((label if kind_ == "p_ctg" else label + "Scaf") + "_" + h, fa))
    if not any(l == f"{label}_hap1" for l, _ in found): die(f"-d {label}: no *.hap1.p_ctg.fa in {fdir}")

    for lab, fa in found:
        name = os.path.basename(fa)[:-3]                      # e.g. myplant.hic.hap1.p_ctg(.scaffolds)
        base = name[:-len(".scaffolds")] if name.endswith(".scaffolds") else name
        sf = os.path.join(ev, "stats", f"{name}.stats.tsv")
        if os.path.exists(sf):
            st = read_stats(sf)
        else:
            print(f"[{lab}] {sf} missing, computing contiguity + telomeres from {fa}")
            st = fasta_stats(fa)
        cs = os.path.join(ev, "compleasm", f"{base}.compleasm.txt")
        if not os.path.exists(cs): cs = os.path.join(ev, "compleasm", base, "summary.txt")   # run_hifiasm.sh layout
        add_row(lab, st, busco(cs), field(qvs, base, 3), field(cps, base, 4),
                fa, os.path.join(ev, "stats", f"{name}.gaps.bed"))

    # images / links for the html report
    dp = []
    for png in sorted(glob.glob(os.path.join(path, "dotplots", "*.dotplot.png"))):
        tsv = png[:-4] + ".tsv"
        lines = open(tsv).read().splitlines() if os.path.exists(tsv) else []
        row = lines[1].split("\t") if len(lines) > 1 else []
        js = png[:-4] + ".json"
        dp.append((png, row, js if os.path.exists(js) else None))
    RUNS.append(dict(label=label, path=os.path.abspath(path), dotplots=dp,
        reads=sorted(glob.glob(os.path.join(path, "reads", "*.html"))),
        readstats=os.path.join(path, "reads", "reads.stats.tsv"),
        command=os.path.join(path, "pipeline_info", "command.txt"),
        spectra=sorted(f for f in glob.glob(os.path.join(ev, "merqury", "*.png"))
                       if re.search(r"^(haps|primary)\.spectra-(cn|asm)\.fl\.png$", os.path.basename(f))),
        contact=sorted(glob.glob(os.path.join(path, "scaffolds", "**", "*.png"), recursive=True)),
        agp=sorted(glob.glob(os.path.join(path, "scaffolds", "**", "*.yahs_scaffolds_final.agp"), recursive=True)),
        cifi=sorted(glob.glob(os.path.join(path, "cifi", "stats", "*.cifi_stats.json"))),
        markdup=sorted(glob.glob(os.path.join(path, "cifi", "markdup", "*.markdup.tsv"))),
        ontcov=sorted(glob.glob(os.path.join(path, "ont_mapping", "*.ont.coverage.tsv"))),
        cifi_html=sorted(glob.glob(os.path.join(path, "cifi", "**", "*.html"), recursive=True))))

    hq = os.path.join(ev, "merqury", "haps.qv")
    if os.path.exists(hq):
        png = os.path.join(ev, "merqury", "haps.spectra-cn.fl.png")
        pairs.append([label, field([hq], "Both", 3),
                      field([os.path.join(ev, "merqury", "haps.completeness.stats")], "both", 4),
                      os.path.abspath(png) if os.path.exists(png) else "NA"])

# ------------------------------------------------------------------ write + print
hdr = ["assembly"] + STATS + ["BUSCO_S", "BUSCO_D", "BUSCO_F", "BUSCO_M", "QV", "kmer_compl_pct"]
phdr = ["run", "combined_QV", "combined_kmer_compl_pct", "spectra_cn_plot"]
def write(path, h, rs):
    with open(path, "w") as fh:
        for r in [h] + rs: fh.write("\t".join(r) + "\n")
def show(h, rs):
    w = [max(len(r[i]) for r in [h] + rs) for i in range(len(h))]
    for r in [h] + rs: print("  ".join(c.ljust(w[i]) for i, c in enumerate(r)).rstrip())
write(os.path.join(outdir, "comparison.tsv"), hdr, rows)
write(os.path.join(outdir, "hap_pairs.tsv"), phdr, pairs)
print(); show(hdr, rows)
if pairs: print(); show(phdr[:3], [p[:3] for p in pairs])
print(f"\nTables: {os.path.join(outdir, 'comparison.tsv')}, {os.path.join(outdir, 'hap_pairs.tsv')}")

# ------------------------------------------------------------------ plots
def make_plots():
    if not plots: return []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not found: plots skipped (tables and html report are written)")
        return []

    group = lambda l: re.sub(r"_(hap[12]|primary)$", "", l)
    groups = list(dict.fromkeys(group(r[0]) for r in rows))
    cmap = plt.get_cmap("tab10" if len(groups) <= 10 else "tab20")
    color = {g: cmap(i % cmap.N) for i, g in enumerate(groups)}
    style = lambda l: "--" if l.endswith("_hap2") else (":" if l.endswith("_primary") else "-")

    # Nx curves
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for label, fa, gb in sources:
        if not fa: continue
        scaf, ctg = lengths(fa, gb)
        for ax, ls in zip(axes, (scaf, ctg)):
            ls = sorted(ls, reverse=True); tot = sum(ls)
            if not tot: continue
            xs, ys, c = [0], [ls[0] / 1e6], 0
            for L in ls:
                xs.append(100 * c / tot); ys.append(L / 1e6); c += L
                xs.append(100 * c / tot); ys.append(L / 1e6)
            ax.plot(xs, ys, style(label), color=color[group(label)], lw=1.4, label=label)
    for ax, t in zip(axes, ("Scaffold Nx", "Contig Nx (split at N gaps)")):
        ax.set_title(t); ax.set_xlabel("x (% of assembly)"); ax.set_xlim(0, 100); ax.set_yscale("log")
        ax.axvline(50, color="0.8", lw=0.8, zorder=0); ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel("sequence length (Mb)")
    axes[1].legend(fontsize=7, loc="upper right", ncol=1 if len(rows) < 12 else 2)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "nx_curves.png"), dpi=150); plt.close(fig)

    # metric bars
    names = [r[0] for r in rows]
    fig, axes = plt.subplots(2, 3, figsize=(16, 4 + 0.3 * len(rows)), constrained_layout=True)
    for ax, (mt, f) in zip(axes.flat, METRICS):
        vals = [f(r) for r in rows]
        known = [v for v in vals if v is not None]
        y = list(range(len(rows)))
        ax.barh(y, [v or 0 for v in vals], color=[color[group(n)] for n in names], edgecolor="white")
        ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8); ax.invert_yaxis()
        ax.set_title(mt, fontsize=10); ax.grid(axis="x", alpha=0.3)
        if known and (mt.startswith("BUSCO") or mt == "merqury QV"):   # zoom in on the differences
            hi = 101 if mt.startswith("BUSCO") else max(known) + 3
            ax.set_xlim(max(0, min(known) - 5), hi)
        else:
            ax.set_xlim(0, (max(known) if known and max(known) > 0 else 1) * 1.18)
        x0 = ax.get_xlim()[0]
        for i, v in enumerate(vals):
            ax.text(v if v is not None else x0, i, f" {v:.4g}" if v is not None else " NA",
                    va="center", fontsize=7, color="0.2" if v is not None else "0.5")
    fig.savefig(os.path.join(outdir, "comparison.png"), dpi=150); plt.close(fig)
    out = [os.path.join(outdir, f) for f in ("comparison.png", "nx_curves.png")]
    print(f"Plots:  {', '.join(out)}")
    return out

def num(v):
    try: return float(str(v).rstrip("%"))
    except ValueError: return None
col = {h: i for i, h in enumerate(hdr)}
METRICS = [("scaffold N50 (Mb)", lambda r: num(r[col["scaffold_N50"]]) and num(r[col["scaffold_N50"]]) / 1e6),
           ("contig N50 (Mb)",   lambda r: num(r[col["contig_N50"]]) and num(r[col["contig_N50"]]) / 1e6),
           ("% bp in scaffolds >= 1 Mb", lambda r: num(r[col["pct_bp_in_ge_1Mb"]])),
           ("T2T sequences",     lambda r: num(r[col["T2T"]])),
           ("BUSCO complete (S+D, %)", lambda r: None if num(r[col["BUSCO_S"]]) is None
                                         else num(r[col["BUSCO_S"]]) + (num(r[col["BUSCO_D"]]) or 0)),
           ("merqury QV",        lambda r: num(r[col["QV"]]))]
pngs = make_plots()   # static PNGs; the html report draws its own interactive charts

# ------------------------------------------------------------------ html report
# self-contained page: sortable tables, interactive SVG / canvas charts (no external libraries)
import base64, datetime, html, json
esc = html.escape

def img(path, alt):
    with open(path, "rb") as fh:
        b = base64.b64encode(fh.read()).decode()
    return f'<img src="data:image/png;base64,{b}" alt="{esc(alt)}" loading="lazy">'

def fig(path, caption):
    return f'<figure>{img(path, caption)}<figcaption>{caption}</figcaption></figure>'

def mb(v, d=2):
    x = num(v)
    return "NA" if x is None else f"{x / 1e6:,.{d}f}"

def intf(v):
    x = num(v)
    return "NA" if x is None else f"{x:,.0f}"

# one hue per run, as light / dark pairs (like the scaffold layout): contigs light, YaHS scaffolds dark
PAIRS = [("#8bbbe3", "#2f6aa8"), ("#f4b26f", "#c9661a"), ("#8fd19a", "#2f8a46"), ("#ee9a94", "#b8352f"),
         ("#bba4e0", "#6a48a8"), ("#86cfd4", "#1f858d"), ("#d2ad8e", "#7d4f32"), ("#f1a7cf", "#b84a8a"),
         ("#cfd07a", "#8a8a1c"), ("#b8bcc2", "#5d6269")]
PALETTE = [d for _, d in PAIRS]
RUN_ORDER = [sp.split("=", 2)[1] for sp in specs]       # -d / -a labels, in command-line order
def run_of(text):
    """run of a row label: the run itself, <run>_hap1, <run>Scaf_hap2, ... (longest matching run label)"""
    hits = [r for r in RUN_ORDER if text == r or text.startswith(r + "_") or text.startswith(r + "Scaf_")]
    return RUN_ORDER.index(max(hits, key=len)) if hits else None

grp = lambda l: re.sub(r"_(hap[12]|primary)$", "", l)
groups = list(dict.fromkeys(grp(r[0]) for r in rows))
runs_in_rows = {run_of(r[0]) for r in rows}
# hue = run (as in the CiFi charts), light = contigs, dark = YaHS scaffolds; hap2 / primary are dashed /
# dotted lines and primary a paler bar. Labels that match no run get the next hues.
_extra = [g for g in groups if run_of(g + "_x") is None]
gcol = {g: PAIRS[((run_of(g + "_x") if run_of(g + "_x") is not None else len(RUN_ORDER) + _extra.index(g)))
                 % len(PAIRS)][1 if g.endswith("Scaf") else 0] for g in groups}

def nx_points(ls):
    """[x_start %, length, rank] steps of the Nx curve, thinned for long tails"""
    ls = sorted(ls, reverse=True); tot = sum(ls)
    if not tot: return []
    pts, c, last = [], 0, -1.0
    for i, L in enumerate(ls):
        x0 = 100 * c / tot; c += L
        if i < 300 or x0 - last >= 0.05:
            pts.append([round(x0, 4), L, i + 1]); last = x0
    pts.append([100.0, ls[-1], len(ls)])
    return pts

DATA = {"rows": [], "nx": {}, "dp": {}, "hic": {}, "scaf": {}, "pair": {}}
for r in rows:
    m = re.search(r"_(hap[12]|primary)$", r[0])
    DATA["rows"].append({"label": r[0], "group": grp(r[0]), "color": gcol[grp(r[0])], "hap": m[1] if m else "",
                         "v": {h: num(r[col[h]]) for h in hdr[1:]}})
for label, fa, gb in sources:
    if fa:
        scaf, ctg = lengths(fa, gb)
        DATA["nx"][label] = {"scaf": nx_points(scaf), "ctg": nx_points(ctg)}

# summary table: (header, cell text, sort value, better = max / min / None)
COLS = [("Assembly", lambda r: r[0], None, None),
        ("Size (Mb)", lambda r: mb(r[col["size_bp"]], 1), lambda r: num(r[col["size_bp"]]), None),
        ("n", lambda r: intf(r[col["scaffolds"]]), lambda r: num(r[col["scaffolds"]]), None),
        ("N50 (Mb)", lambda r: mb(r[col["scaffold_N50"]]), lambda r: num(r[col["scaffold_N50"]]), max),
        ("Longest (Mb)", lambda r: mb(r[col["longest"]]), lambda r: num(r[col["longest"]]), max),
        ("≥ 1 Mb (n)", lambda r: intf(r[col["scaf_ge_1Mb"]]), lambda r: num(r[col["scaf_ge_1Mb"]]), None),
        ("% bp ≥ 1 Mb", lambda r: r[col["pct_bp_in_ge_1Mb"]], lambda r: num(r[col["pct_bp_in_ge_1Mb"]]), max),
        ("Gaps", lambda r: intf(r[col["gaps"]]), lambda r: num(r[col["gaps"]]), None),
        ("n", lambda r: intf(r[col["contigs"]]), lambda r: num(r[col["contigs"]]), None),
        ("N50 (Mb)", lambda r: mb(r[col["contig_N50"]]), lambda r: num(r[col["contig_N50"]]), max),
        ("T2T", lambda r: r[col["T2T"]], lambda r: num(r[col["T2T"]]), max),
        ("One end", lambda r: r[col["one_telo"]], lambda r: num(r[col["one_telo"]]), None),
        ("S", lambda r: r[col["BUSCO_S"]], lambda r: num(r[col["BUSCO_S"]]), max),
        ("D", lambda r: r[col["BUSCO_D"]], lambda r: num(r[col["BUSCO_D"]]), None),
        ("F", lambda r: r[col["BUSCO_F"]], lambda r: num(r[col["BUSCO_F"]]), None),
        ("M", lambda r: r[col["BUSCO_M"]], lambda r: num(r[col["BUSCO_M"]]), min),
        ("QV", lambda r: r[col["QV"]], lambda r: num(r[col["QV"]]), max),
        ("Compl. (%)", lambda r: r[col["kmer_compl_pct"]], lambda r: num(r[col["kmer_compl_pct"]]), max)]

# coloured header groups of the assemblies table (title, n columns, css class); must add up to len(COLS)
COL_GROUPS = [("", 2, ""), ("Scaffolds", 6, "g-scaf"), ("Contigs", 2, "g-ctg"), ("Telomeres", 2, "g-telo"),
              ("BUSCO (%)", 4, "g-busco"), ("Merqury", 2, "g-merq")]
assert sum(n for _, n, _ in COL_GROUPS) == len(COLS)

def row_tags(lab, run=None, hap=None):
    """filter attributes of a table row (report filter bar): run, haplotype, contigs/scaffolds, label.
    From the first cell (<run>_hap1, <run>Scaf_hap2, asm.hic.hap1.p_ctg.scaffolds, ...) or given for the table."""
    ri = run_of(lab)
    r = RUN_ORDER[ri] if ri is not None else run
    m = re.search(r"hap[12]|primary", lab)
    h = hap or (m[0] if m else "primary" if re.search(r"\.p_ctg(\.|$)", lab) else None)
    k = ("scaffolds" if re.search(r"Scaf(_|$)|\.scaffolds$", lab) else "contigs") if h and not hap else None
    a = {"run": r, "hap": h, "kind": k}
    return "".join(f' data-{n}="{esc(v)}"' for n, v in a.items() if v) + (f' data-lab="{esc(lab)}"' if r or h else "")

def table(head, body, numeric_from=1, run=None, hap=None, kind=None, groups=None):
    """body: rows of (text, css class, sort value or None); click a header to sort.
    Rows from several runs (first column = run or assembly label) get a grey background for every other run.
    run / hap / kind tag all rows for the filter bar when the first column does not say (e.g. per-scaffold tables).
    groups: [(title, n columns, css class)] covering all columns: a coloured header row over related columns."""
    cg = [""] * len(head)
    if groups:
        k = 0
        for _, n, c in groups:
            for j in range(k, k + n): cg[j] = (f"cg {c}" + (" cg0" if j == k else "")) if c else ""
            k += n
    h = "".join(f'<th class="{(("num " if i >= numeric_from else "") + cg[i]).strip()}">{esc(c)}</th>' for i, c in enumerate(head))
    if groups:
        h = "".join(f'<th colspan="{n}" class="grp{" " + c + " cg0" if c else ""}">{esc(t)}</th>' for t, n, c in groups) + "</tr><tr>" + h
    def td(i, cell):
        text, cls, v = (list(cell) + ["", None])[:3]
        cls = ("num " if i >= numeric_from else "") + cls + (" " + cg[i] if cg[i] else "")
        dv = "" if i < numeric_from else f' data-v="{"" if v is None else v}"'
        return f'<td class="{cls.strip()}"{dv}>{esc(str(text))}</td>'
    runs = [run_of(str((list(r[0]) + [""])[0])) for r in body]
    banded = len({x for x in runs if x is not None}) > 1
    def tr(r, ri):
        a = ' class="rb"' if banded and ri is not None and ri % 2 else ""
        t = row_tags(str((list(r[0]) + [""])[0]), run, hap)
        if kind: t += f' data-kind="{kind}"'
        return f"<tr{a}{t}>" + "".join(td(i, c) for i, c in enumerate(r)) + "</tr>"
    b = "".join(tr(r, ri) for r, ri in zip(body, runs))
    return f'<div class="scroll"><table class="sortable{" runs" if banded else ""}"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'

body = [[(f(r), "", key(r) if key else None) for name, f, key, _ in COLS] for r in rows]

parts = [f'<section id="assemblies"><h2>Assemblies</h2>'
         f'<p class="note">Click a column header to sort. '
         f'Scaffolds inherit BUSCO and QV from their contigs. All columns: <code>comparison.tsv</code>.</p>'
         f'{table([c[0] for c in COLS], body, groups=COL_GROUPS)}']
if pairs:
    parts.append('<h3>Haplotype pairs (hap1 + hap2)</h3>' +
                 table(["Run", "Combined QV", "Combined k-mer completeness (%)"],
                       [[(p[0], ""), (p[1], "", num(p[1])), (p[2], "", num(p[2]))] for p in pairs]))
parts.append("</section>")
_cn = (["colour = run"] if len(runs_in_rows) > 1 else []) + \
      (["light = contigs, dark = YaHS scaffolds"] if any(g.endswith("Scaf") for g in groups) else [])
COLOR_NOTE = (" " + "; ".join(_cn)[0].upper() + "; ".join(_cn)[1:] + "; grey lines separate the assemblies.") if _cn else ""
parts.append('<section id="plots"><h2>Comparison plots</h2>'
             '<p class="note">Hover a bar for details; the same assembly is highlighted in every panel. '
             'Plots show the assemblies left by the table filter (top bar).'
             + COLOR_NOTE
             + '</p>'
             '<div id="bars" class="bars"></div>'
             '<h3>Nx curves</h3><p class="note">Hover to read Nx of every assembly at that point; click a legend entry to '
             'hide or show it. Solid hap1, dashed hap2, dotted primary.</p><div id="nx" class="nx"></div></section>')

# ---------------- CiFi contact analysis (bin/cifi_stats.py, after McGinty et al. 2026, Dennis lab)
DATA["cifi"] = []
for run in RUNS:
    run["cifi_data"] = []
    for js in run["cifi"]:
        try: d = json.load(open(js))
        except ValueError: continue
        if not d.get("summary"): continue
        m = re.search(r"\.(hap[12])\.", os.path.basename(js))
        lab = f'{run["label"]}_{m[1]}' if m else f'{run["label"]}_{os.path.basename(js).split(".cifi_stats")[0]}'
        d["label"], d["hap"], d["png"] = lab, m[1] if m else "", js.replace(".cifi_stats.json", ".cifi_stats.png")
        run["cifi_data"].append(d)
        DATA["cifi"].append({k: d[k] for k in ("label", "hap", "ps", "cardinality", "summary")})
for d in DATA["cifi"]:      # same hue per run as the assembly plots
    d["color"] = PALETTE[(run_of(d["label"]) or 0) % len(PALETTE)]
CIFI_COLS = [("Haplotype", "label", None), ("CiFi reads", "cifi_reads", "int"), ("Segments", "segments", "int"),
             ("% MAPQ ≥ 1", "pct_segments_mapq_ge1", "pct"), ("% MAPQ ≥ 30", "pct_segments_mapq_ge30", "pct"),
             ("Median usable seg./read", "median_usable_segments", "int"), ("% multi-way reads", "pct_multiway_reads", "pct"),
             ("Contacts/read", "contacts_per_read", "pct"), ("Contacts", "contacts", "int"),
             ("% cis", "pct_cis", "pct"), ("% cis ≥ 10 kb", "pct_cis_ge10kb", "pct"), ("% cis ≥ 1 Mb", "pct_cis_ge1Mb", "pct"),
             ("% trans adjacent", "pct_trans_adjacent", "pct"), ("% trans non-adjacent", "pct_trans_non_adjacent", "pct"),
             ("% reads > 1 chromosome", "pct_reads_multi_chrom", "pct"), ("P(s) slope", "ps_slope_50kb_5Mb", "pct")]
def cfmt(v, kind):
    if v is None: return "NA"
    return f"{v:,.0f}" if kind == "int" else f"{v:,.2f}" if isinstance(v, float) else str(v)
MARKDUP = []      # pbmarkdup: run, library, reads, unique, duplicates, % duplicates (TOTAL row, or the only library)
for run in RUNS:
    for f in run["markdup"]:
        rs = [l.rstrip("\n").split("\t") for l in open(f)][1:]
        rs = [r for r in rs if len(r) >= 5]
        tot = [r for r in rs if r[0] == "TOTAL"] or rs[-1:]
        for r in tot: MARKDUP.append([run["label"]] + r)
CIFI_SUMMARY = ""    # cross-run part of the CiFi section; per-haplotype parts are added below
if DATA["cifi"]:
    body = [[(d["label"], "")] + [(cfmt(d["summary"].get(k), kind), "", d["summary"].get(k)) for _, k, kind in CIFI_COLS[1:]]
            for d in DATA["cifi"]]
    CIFI_SUMMARY = (
        '<h3>Contacts</h3>'
        '<p class="note">CiFi segments mapped to each haplotype&#39;s contigs and lifted to the YaHS scaffolds; every pair of '
        'usable segments (MAPQ ≥ the threshold) in a read is a contact, as in <code>cifi contacts</code> and McGinty et al. '
        '(Nat Commun 2026). <b>Adjacent</b> segments are neighbours along the CiFi read (direct ligation junctions); '
        '<b>non-adjacent</b> pairs are the extra, multi-way contacts that long reads add. cis / trans count only scaffolds '
        'that are chromosome-scale; trans between two large scaffolds can also mean two pieces of one chromosome that were '
        'not joined. A P(s) slope near −1 is typical of chromatin; a flat curve suggests random ligation or a misassembly.</p>'
        + table([c[0] for c in CIFI_COLS], body)
        + '<h3>Contact distance decay P(s)</h3><p class="note">Contact probability per bp of distance (log–log); hover for '
          'values, switch between all pairs and adjacent / non-adjacent pairs.</p>'
          '<div class="seg"><button data-k="all" class="on">all pairs</button><button data-k="adjacent">adjacent</button>'
          '<button data-k="non_adjacent">non-adjacent</button></div><div id="cifi-ps" class="nx"></div>'
          '<h3>Usable segments per CiFi read</h3><p class="note">Share of CiFi reads with k segments mapped at the '
          'MAPQ threshold (k ≥ 2 give contacts, k ≥ 3 are multi-way).</p><div id="cifi-card" class="nx"></div>')

# ---------------- ONT reads mapped to hap1 / hap2 (--map_ont): samtools coverage + flagstat
ONT, ONT_SEQ, ONT_GAP = [], [], []   # per assembly; per sequence >= 1 Mb; joins without ONT support
for run in RUNS:
    for cf in run.get("ontcov", []):
        base = os.path.basename(cf)[:-len(".ont.coverage.tsv")]
        m = re.search(r"\.(hap[12])\.", base)
        lab = f'{run["label"]}_{m[1]}' if m else f'{run["label"]}_{base}'
        seqs = []
        for l in open(cf):
            if l.startswith("#"): continue
            f = l.rstrip("\n").split("\t")
            if len(f) < 9: continue
            L = int(f[2]) - int(f[1]) + 1
            seqs.append((f[0], L, int(f[3]), int(f[4]), float(f[6]), float(f[8])))
        tot = sum(x[1] for x in seqs)
        if not tot: continue
        fs = {}
        ff = cf[:-len(".coverage.tsv")] + ".flagstat.tsv"
        if os.path.exists(ff):
            for l in open(ff):
                f = l.rstrip("\n").split("\t")
                if len(f) >= 3: fs[f[2]] = f[0]
        # ONT support of the joins (N gaps), from ont_gaps.sh: supported / conflict / no_evidence
        gp = cf[:-len(".coverage.tsv")] + ".gaps.tsv"
        gaps_ = [l.rstrip("\n").split("\t") for l in open(gp)][1:] if os.path.exists(gp) else None
        gst = None if gaps_ is None else {k: sum(1 for g in gaps_ if len(g) >= 6 and g[5] == k)
                                         for k in ("supported", "conflict", "no_evidence")}
        if gaps_: ONT_GAP += [[lab, g[0], int(g[1]), int(g[2]), int(g[3]), int(g[4]), g[5]]
                              for g in gaps_ if len(g) >= 6 and g[5] != "supported"]
        depth = sum(x[1] * x[4] for x in seqs) / tot
        big = [x for x in seqs if x[1] >= 1_000_000]
        ONT.append([lab, num(fs.get("primary")), num(fs.get("primary mapped %")), depth,
                    100 * sum(x[3] for x in seqs) / tot, len(seqs),
                    sum(1 for x in big if x[4] < 0.5 * depth or x[4] > 2 * depth), base,
                    None if gst is None else len(gaps_), gst])
        ONT_SEQ += [[lab, x[0], x[1], x[4], 100 * x[3] / x[1], x[5], x[4] / depth if depth else None] for x in big]
ONT_HTML = ""
if ONT:
    f2 = lambda v, d=2: "NA" if v is None else f"{v:,.{d}f}"
    ONT_HTML = (
        '<section id="ont"><h2>ONT read coverage</h2>'
        '<p class="note">ONT reads mapped with minimap2 to each haplotype (scaffolds with CiFi, contigs otherwise); '
        'primary alignments, from <code>samtools flagstat</code> and <code>samtools coverage</code>. Depth is the '
        'length-weighted mean over all sequences. A sequence ≥ 1 Mb at less than half or more than twice that depth '
        'is flagged: possible collapsed repeat, duplication, or contamination. BAMs (with CSI index) are in '
        '<code>ont_mapping/</code> for IGV. <b>Joins</b>: N gaps between contigs (the YaHS joins). A join is '
        '<b>supported</b> when ≥ 2 reads align ≥ 1 kb into both sides on the same strand, a <b>conflict</b> when no '
        'read crosses it and ≥ 2 reads stop at it with ≥ 1 kb clipped (the rest of the read is elsewhere: possible '
        'misjoin); otherwise there is no evidence either way (e.g. a repeat longer than the reads).</p>'
        + table(["Haplotype", "Reads", "% mapped", "Mean depth (×)", "% bases covered", "Sequences", "Flagged ≥ 1 Mb",
                 "Joins", "ONT-supported", "Conflict", "No evidence"],
                [[(o[0], ""), (f2(o[1], 0), "", o[1]), (f2(o[2]), "", o[2]), (f2(o[3], 1), "", o[3]),
                  (f2(o[4]), "", o[4]), (f"{o[5]:,}", "", o[5]), (str(o[6]), "", o[6])]
                 + ([("NA", "", None)] * 4 if o[9] is None else
                    [(f"{o[8]:,}", "", o[8]),
                     (f'{o[9]["supported"]:,} ({100 * o[9]["supported"] / o[8]:.0f}%)' if o[8] else "0", "", o[9]["supported"]),
                     (f'{o[9]["conflict"]:,}', "worse" if o[9]["conflict"] else "", o[9]["conflict"]),
                     (f'{o[9]["no_evidence"]:,}', "", o[9]["no_evidence"])]) for o in ONT])
        + ('' if not ONT_GAP else '<details><summary>Joins not supported by ONT reads '
           f'({sum(1 for g in ONT_GAP if g[6] == "conflict")} conflict, {sum(1 for g in ONT_GAP if g[6] != "conflict")} no evidence)'
           '</summary><p class="note">Conflicts first: check them in the IGV session (<code>IGV session (<code>igv/</code>, track lt;hapIGV session (<code>igv/</code>, track gt;.igv.xml</code>, track '
           '&quot;ONT support of joins&quot;) and the contact map in Juicebox.</p>' +
           table(["Haplotype", "Scaffold", "Gap start", "Gap end", "Spanning reads", "Broken reads", "Status"],
                 [[(g[0], ""), (g[1], ""), (f"{g[2]:,}", "", g[2]), (f"{g[3]:,}", "", g[3]), (str(g[4]), "", g[4]),
                   (str(g[5]), "", g[5]), (g[6].replace("_", " "), "worse" if g[6] == "conflict" else "")]
                  for g in sorted(ONT_GAP, key=lambda g: (g[6] != "conflict", g[0], g[1], g[2]))],
                 numeric_from=2) + '</details>')
        + ('' if not ONT_SEQ else '<h3>Sequences ≥ 1 Mb</h3><p class="note">Click a header to sort; flagged rows '
           'are outside 0.5–2× the haplotype&#39;s mean depth.</p>' +
           table(["Haplotype", "Sequence", "Length (Mb)", "Mean depth (×)", "% covered", "Mean MAPQ", "Depth / mean"],
                 [[(r[0], ""), (r[1], ""), (f2(r[2] / 1e6), "", r[2]), (f2(r[3], 1), "", r[3]), (f2(r[4]), "", r[4]),
                   (f2(r[5], 1), "", r[5]),
                   (f2(r[6]), "worse" if r[6] is not None and (r[6] < 0.5 or r[6] > 2) else "", r[6])] for r in ONT_SEQ],
                 numeric_from=2))
        + '</section>')

# ---------------- one section per topic; with several runs, one subsection per run
multi = len(RUNS) > 1
def run_head(run, level="h3"):
    return f'<{level}>{esc(run["label"])}</{level}>' if multi else ""
def links(files, base):
    return " · ".join(f'<a href="{esc(os.path.relpath(f, outdir))}">{esc(os.path.relpath(f, base))}</a>' for f in files)
def hap_of(path):
    m = re.search(r"\.(hap[12])\.", os.path.basename(path))
    return m[1] if m else ""

spectra = []
for run in RUNS:
    pngs = [p for p in run["spectra"] if os.path.getsize(p)]
    if pngs:
        spectra.append(run_head(run) + '<div class="grid">' + "".join(fig(p, esc(os.path.basename(p))) for p in pngs) + "</div>")
if spectra:
    parts.append('<section id="spectra"><h2>Merqury k-mer spectra</h2>'
                 '<p class="note">Read k-mers by copy number in the assembly (spectra-cn) and by haplotype (spectra-asm).</p>'
                 + "".join(spectra) + "</section>")

# scaffold layout per haplotype: YaHS AGP (contigs, orientation, gaps) + telomere ends of the scaffolds
def read_agp(fn):
    rows, idx = [], {}
    for line in open(fn):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9: continue
        s, b, e = f[0], int(f[1]) - 1, int(f[2])
        if s not in idx: idx[s] = len(rows); rows.append([s, 0, [], [], 0, 0])   # name, length, contigs, gaps, telo
        r = rows[idx[s]]; r[1] = max(r[1], e)
        if f[4] in ("N", "U"): r[3].append([b, e])
        else: r[2].append([f[5], b, e, f[8], int(f[6]) - 1, int(f[7])])        # contig, start, end, orient, cstart, cend
    return rows, idx
def agp_of(run, label):
    """contigs of a scaffolded assembly by scaffold: {scaffold: [[start, end, contig, orient], ...]}"""
    if not label.endswith(".scaffolds"): return None
    fn = [f for f in run["agp"] if os.path.basename(f) == label[:-len(".scaffolds")] + ".yahs_scaffolds_final.agp"]
    if not fn: return None
    return {r[0]: [[c[1], c[2], c[0], c[3]] for c in r[2]] for r in read_agp(fn[0])[0]}

dots = []
for run in RUNS:
    if not run["dotplots"]: continue
    sec = [run_head(run)]
    rs = [r for _, r, _ in run["dotplots"] if len(r) >= 10]
    if rs:
        sec.append(table(["Query", "Target", "Query aligned (%)", "Target aligned (%)", "Inverted (Mb)", "Alignments"],
                         [[(r[0], ""), (r[1], ""), (r[5], "", num(r[5])), (r[7], "", num(r[7])),
                           (mb(r[8]), "", num(r[8])), (intf(r[9]), "", num(r[9]))] for r in rs], 2, run=run["label"]))
    figs = []
    for png, r, js in run["dotplots"]:
        cap = esc(" vs ".join(r[:2])) if r else esc(os.path.basename(png))
        d = None
        if js:
            try: d = json.load(open(js))
            except ValueError: d = None
        if d and d.get("tseq") is not None:
            key = f"dp{len(DATA['dp'])}"
            for ax, seqs in (("tctg", d["tseq"]), ("qctg", d["qseq"])):   # contig strips along the axes
                a = agp_of(run, d.get(ax[0], ""))
                if a: d[ax] = [a.get(q[0], []) for q in seqs]
            DATA["dp"][key] = d
            figs.append(f'<figure class="dpfig"><div class="dp" data-key="{key}"></div><figcaption>{cap}</figcaption></figure>')
        elif os.path.getsize(png):
            figs.append(fig(png, cap))
    sec.append('<div class="grid dpgrid">' + "".join(figs) + "</div>")
    dots.append("".join(sec))
if dots:
    parts.append('<section id="dotplots"><h2>Dot plots</h2><p class="note">Hover to highlight the target and query sequences '
                 'under the cursor; near an alignment, all alignments of that query–target pair and a summary of them; drag to zoom, double-click to reset. Blue same orientation, red inverted. Query sequences are ordered along the '
                 'target and flipped (−) where most of them align reversed. Scaffold plots have a contig strip on each axis '
                 '(YaHS AGP; alternating shades, a break = contig join); hover it to highlight that contig.</p>' + "".join(dots) + "</section>")

def pair_haps(run, n1, r1, n2, r2):
    """hap2 scaffolds matched to hap1 scaffolds: from the hap2-vs-hap1 scaffold dot plot (most aligned bp; reversed
    when most of it aligns on the minus strand), else by size rank. Returns (groups, method):
    groups = [[hap1 index or -1, [[hap2 index, flipped], ...]], ...] in hap1 length order, unmatched hap2 last."""
    i1 = {r[0]: i for i, r in enumerate(r1)}; i2 = {r[0]: i for i, r in enumerate(r2)}
    d = None
    for _png, _row, js in run["dotplots"]:
        try: j = json.load(open(js)) if js else None
        except ValueError: j = None
        if j and j.get("t") == n1 + ".scaffolds" and j.get("q") == n2 + ".scaffolds" and j.get("segs"): d = j; break
    m = {}                                                # hap2 index -> [hap1 index, flipped, mean hap1 position]
    if d:
        acc = {}
        for ti, qi, ts, te, qs, qe, plus in d["segs"]:
            a = acc.setdefault(qi, {}).setdefault(ti, [0, 0, 0.0])        # bp, plus bp, position * bp
            w = qe - qs; a[0] += w; a[1] += w * plus; a[2] += (ts + te) / 2 * w
        for qi, ts_ in acc.items():
            ti, a = max(ts_.items(), key=lambda kv: kv[1][0])
            q, t = d["qseq"][qi][0], d["tseq"][ti][0]
            if q in i2 and t in i1: m[i2[q]] = [i1[t], a[1] * 2 < a[0], a[2] / a[0]]
        method = "alignment"
    else:                                                 # no dot plot: pair the n-th longest of each haplotype
        o1 = sorted(range(len(r1)), key=lambda i: -r1[i][1]); o2 = sorted(range(len(r2)), key=lambda i: -r2[i][1])
        for a, b in zip(o1, o2):
            if r1[a][1] >= 1e6 and r2[b][1] >= 1e6: m[b] = [a, False, 0]
        method = "size"
    groups = []
    for a in sorted(range(len(r1)), key=lambda i: -r1[i][1]):
        hs = sorted((b for b, v in m.items() if v[0] == a), key=lambda b: m[b][2])
        groups.append([a, [[b, int(m[b][1])] for b in hs]])
    groups += [[-1, [[b, 0]]] for b in sorted((b for b in range(len(r2)) if b not in m), key=lambda b: -r2[b][1])]
    return groups, method

scafs = []
for run in RUNS:
    figs, byhap = [], {}
    for fn in run["agp"]:
        name = os.path.basename(fn)[:-len(".yahs_scaffolds_final.agp")]
        srows, idx = read_agp(fn)
        tt = os.path.join(run["path"], "eval", "telomeres", f"{name}.scaffolds.telo.tsv")
        if os.path.exists(tt):
            for line in list(open(tt))[1:]:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 6 and f[0] in idx: srows[idx[f[0]]][4:6] = [int(f[4]), int(f[5])]
        if not srows: continue
        key = f"scaf{len(DATA['scaf'])}"
        h = hap_of(fn)
        # asm: the Assemblies-table label of these scaffolds (colour + filter in the scaffold scatter plot)
        DATA["scaf"][key] = {"rows": srows, "hap": h, "run": run["label"], "asm": f'{run["label"]}Scaf_{h}' if h else ""}
        if h: byhap[h] = (name, srows)
        lab = f'{run["label"]}_{h}' if h else f'{run["label"]}_{name}'
        figs.append(f'<figure class="scfig"><div class="scaf" data-key="{key}"></div><figcaption>{esc(lab)}: '
                    f'{sum(1 for r in srows if r[1] >= 1e6)} scaffolds ≥ 1 Mb of {len(srows)}</figcaption></figure>')
    if "hap1" in byhap and "hap2" in byhap:              # hap1 / hap2 side by side, per matched chromosome
        (n1, r1), (n2, r2) = byhap["hap1"], byhap["hap2"]
        groups, method = pair_haps(run, n1, r1, n2, r2)
        key = f"pair{len(DATA['pair'])}"
        DATA["pair"][key] = {"h1": r1, "h2": r2, "groups": groups, "by": method}
        how = ("hap2 scaffolds matched to the hap1 scaffold they align to best (hap2 vs hap1 scaffold dot plot), "
               "reverse-complemented (−) where they align reversed" if method == "alignment" else
               "no hap2 vs hap1 scaffold dot plot: paired by size rank only, which can mismatch chromosomes")
        # replaces the per-haplotype figures (same scaffolds, drawn per chromosome pair); kept for single-haplotype runs
        figs = [f'<figure class="scfig pairfig"><div class="scpair" data-key="{key}"></div><figcaption>'
                f'{esc(run["label"])}: hap1 and hap2 paired; {how}.</figcaption></figure>']
    if figs: scafs.append(run_head(run) + '<div class="grid scgrid">' + "".join(figs) + "</div>")
if scafs:
    parts.append('<section id="scaffolds"><h2>Scaffolds</h2>'
                 '<p class="note">Every YaHS scaffold ≥ 1 Mb: length against the number of contigs it is built from, '
                 'all runs in one plot (contig-only assemblies have no scaffolds and are not shown). Well-assembled chromosomes sit '
                 'low (few contigs) as filled circles (telomere at both ends). Follows the table filter; click a run in the legend '
                 'to hide or show it, hover a point for details.</p><div id="scsc" class="nx"></div>'
                 '<h3>Layout</h3><p class="note">YaHS scaffolds, hap1 and hap2 paired by '
                 'chromosome when both are scaffolded (one figure per haplotype otherwise): contigs '
                 'in alternating colours (arrow: orientation), gaps as dark marks, '
                 '<span style="color:#d62728">■</span> telomere at that end; length and number of contigs at the '
                 'right. Hover a contig or scaffold name for details. In the paired view, a hollow '
                 '<span style="color:#d62728">□</span> marks a telomere that the other haplotype has and this one '
                 'lacks, and a contig count in red differs from the other haplotype (1:1 pairs only).</p>' + "".join(scafs) + "</section>")

# CiFi: library QC -> contact summary across runs -> per run and haplotype: contact map, analysis, scaffolds
lib = []
for run in RUNS:
    if run["cifi_html"]:
        lib.append(f'<p>{esc(run["label"]) + ": " if multi else ""}cifi reports: '
                   + links(run["cifi_html"], os.path.join(run["path"], "cifi")) + "</p>")
if MARKDUP:
    lib.append('<p class="note">PCR duplicates (pbmarkdup): CiFi reads sequenced from the same molecule, before '
               'digestion (the paper reports 1.8% for its GM12878 DpnII library). Rerun with <code>--cifi_rmdup</code> '
               'to remove them.</p>' + table(["Run", "Reads", "Unique molecules", "Duplicate reads", "% duplicates"],
               [[(m[0], ""), (f"{int(m[2]):,}", "", int(m[2])), (f"{int(m[3]):,}", "", int(m[3])),
                 (f"{int(m[4]):,}", "", int(m[4])), (m[5], "", num(m[5]))] for m in MARKDUP]))
per_hap = []
for run in RUNS:
    maps = [p for p in run["contact"] if os.path.getsize(p)]
    data = {d["hap"]: d for d in run["cifi_data"]}
    haps = sorted(set(data) | {hap_of(p) for p in maps})
    if not haps: continue
    sec = []
    for h in haps:   # labels are <run>_<hap>, so no extra run heading
        d = data.get(h)
        name = d["label"] if d else f'{run["label"]}_{h}' if h else run["label"]
        sec.append(f'<h3>{esc(name)}</h3>')
        hm = [p for p in maps if hap_of(p) == h]
        hic = d.get("hic") if d else None
        if hic and hic.get("levels"):
            key = f"hic{len(DATA['hic'])}"
            DATA["hic"][key] = hic
            sec.append(f'<figure class="hicfig"><div class="hic" data-key="{key}"></div><figcaption>{esc(name)}: '
                       f'CiFi contacts (segment MAPQ ≥ {hic.get("min_mapq", 0)}) on the YaHS scaffolds, '
                       f'{hic["pairs"]:,} contacts</figcaption></figure>')
            if hm: sec.append('<details><summary>PretextSnapshot image</summary><div class="grid">' + "".join(
                fig(p, esc(os.path.basename(p))) for p in hm) + "</div></details>")
        elif hm: sec.append('<div class="grid">' + "".join(
            fig(p, "contact map (PretextSnapshot): " + esc(os.path.basename(p))) for p in hm) + "</div>")
        if d and os.path.exists(d["png"]) and os.path.getsize(d["png"]):   # multi-panel figure: full width
            sec.append(fig(d["png"], f'{esc(d["label"])}: segments per read, P(s), MAPQ and contact classes, coverage'))
        rows_ = d.get("scaffolds", []) if d else []
        if rows_:
            cis = sorted(r[6] for r in rows_ if r[6] is not None)
            med = cis[len(cis) // 2] if cis else None
            body = [[(r[0], ""), (f"{r[1] / 1e6:,.2f}", "", r[1]), (f"{r[2]:,}", "", r[2]), (f"{r[3]:,.1f}", "", r[3]),
                     (f"{r[4]:,}", "", r[4]), (f"{r[5]:,}", "", r[5]),
                     ("NA" if r[6] is None else f"{r[6]:.2f}", "low" if med and r[6] is not None and r[6] < 0.8 * med else "", r[6])]
                    for r in rows_]
            sec.append(f'<p class="note">Chromosome-scale scaffolds (≥ {d.get("min_scaf", 0) / 1e6:g} Mb). '
                       '<span class="lowkey">Highlighted</span>: % cis below 80% of the median, worth a look in '
                       'Juicebox (possible misjoin or a scaffold that belongs to another chromosome).</p>'
                       + table(["Scaffold", "Length (Mb)", "Segments", "Segments/Mb", "cis contacts",
                                "trans contacts", "% cis"], body, run=run["label"], hap=h or None, kind="scaffolds"))
    per_hap.append("".join(sec))
HAS_CIFI = bool(lib or CIFI_SUMMARY or per_hap)
if HAS_CIFI:
    parts.append('<section id="cifi"><h2>CiFi</h2><p class="note">Library QC, contact statistics of all haplotypes, '
                 'then per haplotype: contact map of the YaHS scaffolds (solid lines scaffolds, dashed lines contigs; scroll to zoom, drag to pan, double-click '
                 'to reset; finer bins load as you zoom), contact analysis and per-scaffold contacts.</p>'
                 + ('<h3>Library</h3>' + "".join(lib) if lib else "")
                 + CIFI_SUMMARY
                 + "".join(per_hap)
                 + "</section>")

parts.append(ONT_HTML)

# reads: per read type (HiFi, ONT, CiFi) from faster --table (reads/reads.stats.tsv) + faster-report links.
# Coverage = bases / size of the run's hap1 contigs (a haploid genome size estimate).
KINDS = {"hifi": "HiFi", "ont": "ONT", "cifi": "CiFi"}
asm_size = {r[0]: num(r[1]) for r in rows}
rtot, rfile = [], []
for run in RUNS:
    if not os.path.exists(run["readstats"]): continue
    gsize = asm_size.get(f'{run["label"]}_hap1')
    by_kind = {}
    for l in list(open(run["readstats"]))[1:]:
        f = l.rstrip("\n").split("\t")
        if len(f) < 13: continue
        fn = os.path.basename(f[0]); k = fn.split(".")[0]
        kind = KINDS.get(k, k)
        n, bp, mean, n50, q20, q30 = int(f[1]), int(f[2]), float(f[6]), int(f[10]), float(f[11]), float(f[12])
        rfile.append([run["label"], kind, fn[len(k) + 1:] if k in KINDS else fn, n, bp, mean, n50, q30])
        t = by_kind.setdefault(kind, [0, 0, 0, 0, 0.0, 0.0])          # files, reads, bases, max N50, Q20 bp, Q30 bp
        t[0] += 1; t[1] += n; t[2] += bp; t[3] = max(t[3], n50); t[4] += q20 * bp; t[5] += q30 * bp
    for kind, t in by_kind.items():
        rtot.append([run["label"], kind, t[0], t[1], t[2], t[2] / t[1] if t[1] else None, t[3],
                     t[4] / t[2] if t[2] else None, t[5] / t[2] if t[2] else None, t[2] / gsize if gsize else None])
rlinks = [(f'<p>{esc(run["label"]) + ": " if multi else ""}faster-report: '
           + links(run["reads"], os.path.join(run["path"], "reads")) + "</p>") for run in RUNS if run["reads"]]
HAS_READS = bool(rtot or rlinks)
if HAS_READS:
    g = lambda v, d=1: "NA" if v is None else f"{v:,.{d}f}"
    lead = [("Run", 0)] if multi else []
    sec = ['<section id="reads"><h2>Reads</h2>']
    if rtot:
        sec.append('<p class="note">Starting reads (before filters), from <code>faster --table</code>. N50 is the largest '
                   'of the files; Q20 / Q30 are base-weighted over the files. Coverage = bases / hap1 contig size of the '
                   'run (a haploid genome size estimate): per haplotype it is half of that for HiFi and ONT.</p>'
                   + table([c for c, _ in lead] + ["Type", "Files", "Reads (n)", "Bases (Gb)", "Mean length (kb)",
                            "N50 (kb)", "Q20 (%)", "Q30 (%)", "Coverage (×)"],
                           [[(r[0], "")] * bool(multi) + [(r[1], ""), (str(r[2]), "", r[2]), (f"{r[3]:,}", "", r[3]),
                             (g(r[4] / 1e9, 2), "", r[4]), (g(r[5] and r[5] / 1e3), "", r[5]), (g(r[6] / 1e3), "", r[6]),
                             (g(r[7], 2), "", r[7]), (g(r[8], 2), "", r[8]), (g(r[9]), "", r[9])] for r in rtot],
                           1 + bool(multi)))
        sec.append('<details><summary>Per file</summary>' + table(
            [c for c, _ in lead] + ["Type", "File", "Reads (n)", "Bases (Gb)", "Mean length (kb)", "N50 (kb)", "Q30 (%)"],
            [[(r[0], "")] * bool(multi) + [(r[1], ""), (r[2], ""), (f"{r[3]:,}", "", r[3]), (g(r[4] / 1e9, 2), "", r[4]),
              (g(r[5] / 1e3), "", r[5]), (g(r[6] / 1e3), "", r[6]), (g(r[7], 2), "", r[7])] for r in rfile],
            2 + bool(multi)) + "</details>")
    sec += rlinks
    parts.insert(0, "".join(sec) + "</section>")

# commands: the nextflow command of each run (pipeline_info/command.txt) and this script's
cmds = [(f'nxf-cifi run {run["label"]}', open(run["command"]).read().strip()) for run in RUNS if os.path.exists(run["command"])]
cmds.append(("report", os.environ.get("REPORT_CMD", "compare_assemblies.sh").strip()))
CMD_HTML = ('<details class="cmd"><summary>Command used</summary>'
            + "".join(f'<p class="note">{esc(t)}</p><pre>{esc(c)}</pre>' for t, c in cmds) + "</details>")

CSS = r"""
:root { --bg:#f7f7f5; --fg:#1d1f21; --muted:#676c73; --line:#e1e2de; --card:#ffffff; --accent:#1f5fa8; --best:#e3edf9; --worse:#fbe5e1; --hover:#f1f4f8; --band:#ecede8; }
* { box-sizing:border-box }
html { color-scheme: light }
body { margin:0; background:var(--bg); color:var(--fg); font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif }
header, main { max-width:1400px; margin:0 auto; padding:0 16px }
header { padding-top:28px; position:sticky; top:-80px; background:var(--bg); z-index:5 }
h1 { font-size:24px; margin:0 0 4px } h2 { font-size:19px; margin:28px 0 8px; padding-top:8px }
section + section h2 { border-top:1px solid var(--line); padding-top:24px }
h3 { font-size:15px; margin:22px 0 6px }
.note, figcaption, .meta { color:var(--muted); font-size:13px }
.note { margin:0 0 10px }
nav { padding:10px 0; border-bottom:1px solid var(--line); overflow-x:auto; white-space:nowrap }
nav a { color:var(--accent); text-decoration:none; margin-right:18px }
a { color:var(--accent) }
code { font-size:12px; word-break:break-all }
details { margin:10px 0 }
details > summary { cursor:pointer; color:var(--accent) }
details.cmd { margin:16px 0 0 }
pre { font:12px/1.45 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; background:var(--card); border:1px solid var(--line);
      border-radius:6px; padding:10px 12px; margin:4px 0 10px; white-space:pre-wrap; word-break:break-all }
.scroll { overflow-x:auto; width:fit-content; max-width:100%; border:1px solid var(--line); border-radius:6px; background:var(--card) }
table { border-collapse:collapse; font-size:13px }
th, td { padding:6px 10px; border-bottom:1px solid var(--line); white-space:nowrap; text-align:left }
th { font-weight:600; background:#fafafa; cursor:pointer; user-select:none; position:relative; padding-right:20px }
th:hover { background:var(--hover) }
/* column groups (assemblies table): one colour per topic, header and cells */
th.grp { cursor:default; text-align:center; font-size:12px; letter-spacing:.02em; padding:4px 10px; color:#3b3f45 }
th.grp::after { content:none }
th.g-scaf, th.g-scaf:hover { background:#dbe7f6 } td.g-scaf { background:#f5f9fd }
th.g-ctg, th.g-ctg:hover { background:#dbf0df } td.g-ctg { background:#f5fbf6 }
th.g-telo, th.g-telo:hover { background:#f4ead6 } td.g-telo { background:#fcf9f3 }
th.g-busco, th.g-busco:hover { background:#e9e1f5 } td.g-busco { background:#faf8fd }
th.g-merq, th.g-merq:hover { background:#f7dfdb } td.g-merq { background:#fdf7f6 }
th.cg:hover { filter:brightness(.95) }
th.cg0, td.cg0 { border-left:1px solid #d3d6da }
th::after { content:"↕"; position:absolute; right:6px; color:#b5b9bf; font-size:11px }
th[data-dir="asc"]::after { content:"▲"; color:var(--accent) } th[data-dir="desc"]::after { content:"▼"; color:var(--accent) }
th.num, td.num { text-align:right; font-variant-numeric:tabular-nums }
td.best { font-weight:700; background:var(--best) }
td.worse { font-weight:700; background:var(--worse) }
#flt { display:flex; flex-wrap:wrap; gap:6px 14px; align-items:center; padding:8px 0; border-bottom:1px solid var(--line); font-size:13px; color:var(--muted) }
#flt:empty { display:none }
#flt details.ms { position:relative }
#flt details.ms summary { list-style:none; cursor:pointer; font-size:12px; padding:3px 22px 3px 8px; border:1px solid var(--line); border-radius:6px; background:var(--card); color:var(--fg); position:relative; white-space:nowrap }
#flt details.ms summary::-webkit-details-marker { display:none }
#flt details.ms summary::after { content:'▾'; position:absolute; right:7px; color:var(--muted) }
#flt details.ms.on summary { border-color:var(--accent); color:var(--accent) }
#flt .msbox { position:absolute; top:calc(100% + 3px); left:0; z-index:20; min-width:170px; max-height:300px; overflow:auto; background:var(--card); border:1px solid var(--line); border-radius:6px; box-shadow:0 4px 14px rgba(0,0,0,.12); padding:6px 0 }
#flt .msbox label { display:flex; gap:6px; align-items:center; padding:3px 10px; color:var(--fg); font-size:12px; cursor:pointer; white-space:nowrap }
#flt .msbox label:hover { background:var(--hover) }
#flt .msact { display:flex; gap:10px; padding:0 10px 5px; margin-bottom:4px; border-bottom:1px solid var(--line) }
#flt .msact button { font:inherit; font-size:12px; border:0; background:none; color:var(--accent); cursor:pointer; padding:0 }
#flt input[type=search] { font:inherit; font-size:12px; padding:3px 8px; border:1px solid var(--line); border-radius:6px; width:170px; background:var(--card); color:var(--fg) }
#flt .cnt { font-size:12px } #flt .clr { font:inherit; font-size:12px; border:0; background:none; color:var(--accent); cursor:pointer; padding:0 }
tbody tr.fx { display:none }
tbody tr:hover td { background:var(--hover) } tbody tr:hover td.best { background:var(--best) }
tbody tr:last-child td { border-bottom:none }
table.runs tbody tr.rb td:not(.best):not(.worse):not(.low):not(.cg) { background:var(--band) }
table.runs tbody tr.rb td.cg:not(.best):not(.worse):not(.low) { box-shadow:inset 0 0 0 999px rgba(70,70,50,.06) }   /* every other run */
figure { margin:0 0 18px; background:var(--card); border:1px solid var(--line); border-radius:6px; padding:10px }
figure img { width:100%; height:auto; display:block }
figcaption { padding-top:6px }
.grid { display:grid; grid-template-columns:repeat(auto-fill, minmax(420px, 1fr)); gap:16px }
.dpgrid { grid-template-columns:repeat(auto-fill, minmax(min(100%, 560px), 1fr)) }
@media (max-width:600px) { .grid { grid-template-columns:1fr } }
.bars { display:grid; grid-template-columns:repeat(auto-fill, minmax(min(100%, 400px), 1fr)); gap:14px }
.nx { display:grid; grid-template-columns:repeat(auto-fill, minmax(min(100%, 520px), 1fr)); gap:14px }
.panel { background:var(--card); border:1px solid var(--line); border-radius:6px; padding:10px }
.panel svg { width:100%; height:auto; display:block; overflow:visible }
svg text { font-family:inherit; fill:#3b3f45 }
svg .ptitle { font-size:13px; font-weight:600; fill:#1d1f21 }
svg .tick { font-size:10px; fill:#80858c } svg .lab { font-size:11px } svg .val { font-size:10px; fill:#555 }
svg .na { font-size:10px; fill:#a0a4aa; font-style:italic }
svg .grid { stroke:#eceded } svg .sep { stroke:#c4c7cc; stroke-width:1 } svg .axis { stroke:#b9bdc2 }
svg .hit { fill:transparent } svg g.row { cursor:default }
.bars.hl g.row:not(.on) { opacity:.3 } .bars g.row.on .hit { fill:var(--hover) }
svg .cross { stroke:#8d9298; stroke-dasharray:3 3 }
svg .n50 { stroke:#c9ccd0; stroke-dasharray:4 4 }
.legend { display:flex; flex-wrap:wrap; gap:4px 14px; margin-top:10px; font-size:12px; grid-column:1/-1 }
.legend span { cursor:pointer; display:inline-flex; align-items:center; gap:6px; user-select:none }
.legend span.off { opacity:.35; text-decoration:line-through }
.legend.static span { cursor:default }
.legend svg { width:26px; height:10px }
.dp { position:relative; width:100% }
.hicfig { max-width:820px }
.scgrid { grid-template-columns:repeat(auto-fill, minmax(min(100%, 560px), 1fr)) }
.scaf svg { display:block; width:100%; height:auto; overflow:visible }
.scaf .row:hover .lab, .scpair .row:hover .lab { font-weight:700 }
.scgrid .pairfig { grid-column:1/-1 }
.scpair .more { font:12px inherit; margin-bottom:6px; padding:2px 8px; border:1px solid var(--line); background:#fff; border-radius:4px; cursor:pointer }
.scaf .more { font:12px inherit; margin-bottom:6px; padding:2px 8px; border:1px solid var(--line); background:#fff; border-radius:4px; cursor:pointer }
.hic { position:relative; width:100% }
.hic canvas { display:block; cursor:grab; max-width:100%; touch-action:none } .hic canvas.drag { cursor:grabbing }
.hic .reset { position:absolute; top:4px; left:4px; font:12px inherit; padding:2px 8px; border:1px solid var(--line); background:#fff; border-radius:4px; cursor:pointer; display:none }
.hic .res { position:absolute; top:6px; right:8px; font-size:11px; color:var(--muted) }
.dp canvas { display:block; cursor:crosshair; max-width:100% }
.dp .sel { position:absolute; border:1px solid var(--accent); background:rgba(31,95,168,.10); display:none; pointer-events:none }
.dp .reset { position:absolute; top:4px; left:4px; font:12px inherit; padding:2px 8px; border:1px solid var(--line); background:#fff; border-radius:4px; cursor:pointer; display:none }
#tip { position:fixed; z-index:20; display:none; pointer-events:none; background:#fff; color:var(--fg); border:1px solid #d5d8dc;
       border-radius:6px; padding:8px 10px; font-size:12px; line-height:1.45; box-shadow:0 4px 14px rgba(0,0,0,.12); max-width:380px }
#tip table { font-size:12px } #tip td { padding:1px 6px 1px 0; border:none; background:none }
.sw { display:inline-block; width:10px; height:10px; border-radius:2px; margin-right:6px; vertical-align:-1px }
h4 { font-size:14px; margin:18px 0 6px }
td.low { background:#fde8e6; font-weight:700 } .lowkey { background:#fde8e6; padding:0 4px; border-radius:3px }
.seg { display:inline-flex; border:1px solid var(--line); border-radius:6px; overflow:hidden; margin:0 0 10px }
.seg button { font:inherit; font-size:12px; padding:4px 12px; border:0; background:#fff; color:var(--fg); cursor:pointer }
.seg button + button { border-left:1px solid var(--line) }
.seg button.on { background:var(--accent); color:#fff }
footer { max-width:1400px; margin:40px auto 30px; padding:0 16px; color:var(--muted); font-size:12px }
"""

JS = r"""
(function () {
const D = JSON.parse(document.getElementById('report-data').textContent);
const tip = document.getElementById('tip');
function showTip(h, ev) {
  tip.innerHTML = h; tip.style.display = 'block';
  const r = tip.getBoundingClientRect(), p = 14;
  let x = ev.clientX + p, y = ev.clientY + p;
  if (x + r.width > innerWidth - 8) x = ev.clientX - r.width - p;
  if (y + r.height > innerHeight - 8) y = ev.clientY - r.height - p;
  tip.style.left = Math.max(4, x) + 'px'; tip.style.top = Math.max(4, y) + 'px';
}
function hideTip() { tip.style.display = 'none'; }
const fmt = (v, d = 2) => v == null ? 'NA' : v.toLocaleString(undefined, {minimumFractionDigits: d, maximumFractionDigits: d});
const fint = v => v == null ? 'NA' : Math.round(v).toLocaleString();
const flen = bp => bp == null ? 'NA' : bp >= 1e6 ? fmt(bp / 1e6, 2) + ' Mb' : bp >= 1e3 ? fmt(bp / 1e3, 1) + ' kb' : Math.round(bp) + ' bp';
const sw = c => `<span class="sw" style="background:${c}"></span>`;
const escH = s => String(s).replace(/[&<>"]/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c]));
const NS = 'http://www.w3.org/2000/svg';
function el(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) if (attrs[k] !== '' && attrs[k] != null) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function svgPoint(svg, ev) {
  const p = svg.createSVGPoint(); p.x = ev.clientX; p.y = ev.clientY;
  return p.matrixTransform(svg.getScreenCTM().inverse());
}
function niceTicks(lo, hi, n) {
  const span = hi - lo; if (!(span > 0)) return [lo];
  let step = Math.pow(10, Math.floor(Math.log10(span / n)));
  const err = span / n / step;
  if (err >= 7.5) step *= 10; else if (err >= 3.5) step *= 5; else if (err >= 1.5) step *= 2;
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) out.push(+t.toFixed(10));
  return out;
}
const short = v => Number.isInteger(v) ? v.toLocaleString() : Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2);
const dash = h => h === 'hap2' ? '5 3' : h === 'primary' ? '1.5 2.5' : '';

/* ---------------------------------------------------------------- sortable tables */
document.querySelectorAll('table.sortable').forEach(t => {
  const ths = [...t.tHead.rows[t.tHead.rows.length - 1].cells];
  ths.forEach((th, i) => th.addEventListener('click', () => {
    const asc = th.dataset.dir !== 'asc';
    ths.forEach(h => delete h.dataset.dir);
    th.dataset.dir = asc ? 'asc' : 'desc';
    const key = r => { const c = r.cells[i], v = c.dataset.v; return v === undefined ? c.textContent : v === '' ? null : +v; };
    const tb = t.tBodies[0];
    [...tb.rows].sort((a, b) => {
      const x = key(a), y = key(b);
      if (x === null || y === null) return x === y ? 0 : x === null ? 1 : -1;    // NA always last
      const c = typeof x === 'number' ? x - y : String(x).localeCompare(String(y), undefined, {numeric: true});
      return asc ? c : -c;
    }).forEach(r => tb.appendChild(r));
  }));
});

/* ---------------------------------------------------------------- table filter (header bar) */
// FLT: assembly labels hidden by the filter (rows of the Assemblies table); plots redraw on the 'flt' event
const FLT = new Set();
(() => {
  const bar = document.getElementById('flt'), rows = [...document.querySelectorAll('table.sortable tbody tr[data-lab]')];
  if (!rows.length) return;
  const vals = k => [...new Set(rows.map(r => r.dataset[k]).filter(Boolean))];
  const order = {hap1: 0, hap2: 1, primary: 2, contigs: 0, scaffolds: 1};
  const dims = [['run', 'Run', vals('run')], ['hap', 'Haplotype', vals('hap').sort((a, b) => order[a] - order[b])],
                ['kind', 'Assembly', vals('kind').sort((a, b) => order[a] - order[b])]].filter(d => d[2].length > 1);
  const off = {run: new Set(), hap: new Set(), kind: new Set()};   // unticked values per dimension
  const q = document.createElement('input'), cnt = document.createElement('span'), clr = document.createElement('button');
  q.type = 'search'; q.placeholder = 'filter rows (text)…'; cnt.className = 'cnt'; clr.className = 'clr'; clr.textContent = 'show all';
  const menus = [];
  function apply() {
    const t = q.value.trim().toLowerCase(); let n = 0;
    for (const r of rows) {
      const d = r.dataset, txt = [d.run, d.hap, d.kind, d.lab].join(' ').toLowerCase();
      const hide = Object.keys(off).some(k => d[k] && off[k].has(d[k])) || (t && !txt.includes(t));
      r.classList.toggle('fx', !!hide); n += !!hide;
    }
    menus.forEach(m => m());
    FLT.clear();
    document.querySelectorAll('#assemblies tbody tr.fx[data-lab]').forEach(r => FLT.add(r.dataset.lab));
    document.dispatchEvent(new Event('flt'));
    cnt.textContent = n ? `${n} of ${rows.length} rows hidden` : ''; clr.style.display = n ? '' : 'none';
  }
  bar.append(Object.assign(document.createElement('span'), {textContent: 'Filter tables:'}));
  for (const [k, lab, vs] of dims) {
    // multi-select dropdown: <details> with a checkbox per value, plus all / none
    const dd = document.createElement('details'), sm = document.createElement('summary'), box = document.createElement('div');
    dd.className = 'ms'; box.className = 'msbox';
    const cbs = vs.map(v => {
      const l = document.createElement('label'), c = document.createElement('input');
      c.type = 'checkbox'; c.checked = true;
      c.addEventListener('change', () => { c.checked ? off[k].delete(v) : off[k].add(v); apply(); });
      l.append(c, ' ' + v); box.append(l); return [c, v];
    });
    const act = document.createElement('div'); act.className = 'msact';
    for (const [txt, f] of [['all', () => off[k].clear()], ['none', () => vs.forEach(v => off[k].add(v))]]) {
      const b = document.createElement('button'); b.textContent = txt; b.type = 'button';
      b.addEventListener('click', () => { f(); apply(); }); act.append(b);
    }
    box.prepend(act);
    menus.push(() => {
      cbs.forEach(([c, v]) => c.checked = !off[k].has(v));
      const on = vs.filter(v => !off[k].has(v));
      sm.textContent = `${lab}: ${on.length === vs.length ? 'all' : !on.length ? 'none' : on.length <= 2 ? on.join(', ') : on.length + ' of ' + vs.length}`;
      dd.classList.toggle('on', on.length < vs.length);
    });
    dd.append(sm, box); bar.append(dd);
  }
  document.addEventListener('click', ev => bar.querySelectorAll('details.ms[open]').forEach(d => { if (!d.contains(ev.target)) d.open = false; }));
  q.addEventListener('input', apply);
  clr.addEventListener('click', () => { Object.values(off).forEach(o => o.clear()); q.value = ''; apply(); });
  bar.append(q, cnt, clr); apply();
})();

/* ---------------------------------------------------------------- metric bars */
const V = (r, k) => r.v[k];
const METRICS = [
  {t: 'Scaffold N50 (Mb)', get: r => V(r, 'scaffold_N50') == null ? null : V(r, 'scaffold_N50') / 1e6,
   tip: r => `Scaffold N50 <b>${flen(V(r, 'scaffold_N50'))}</b><br>${fint(V(r, 'scaffolds'))} scaffolds · longest ${flen(V(r, 'longest'))}<br>assembly size ${flen(V(r, 'size_bp'))}`},
  {t: 'Contig N50 (Mb)', get: r => V(r, 'contig_N50') == null ? null : V(r, 'contig_N50') / 1e6,
   tip: r => `Contig N50 <b>${flen(V(r, 'contig_N50'))}</b><br>${fint(V(r, 'contigs'))} contigs · ${fint(V(r, 'gaps'))} gaps`},
  {t: '% bp in scaffolds ≥ 1 Mb', get: r => V(r, 'pct_bp_in_ge_1Mb'),
   tip: r => `<b>${fmt(V(r, 'pct_bp_in_ge_1Mb'))}%</b> of ${flen(V(r, 'size_bp'))}<br>in ${fint(V(r, 'scaf_ge_1Mb'))} scaffolds ≥ 1 Mb`},
  {t: 'T2T sequences', get: r => V(r, 'T2T'),
   tip: r => `<b>${fint(V(r, 'T2T'))}</b> telomere-to-telomere<br>${fint(V(r, 'one_telo'))} with a telomere at one end`},
  {t: 'BUSCO complete (S+D, %)', zoom: true, max: 100,
   get: r => V(r, 'BUSCO_S') == null ? null : V(r, 'BUSCO_S') + (V(r, 'BUSCO_D') || 0),
   tip: r => `Complete <b>${fmt(V(r, 'BUSCO_S') + (V(r, 'BUSCO_D') || 0))}%</b><br>single ${fmt(V(r, 'BUSCO_S'))}% · duplicated ${fmt(V(r, 'BUSCO_D'))}%<br>fragmented ${fmt(V(r, 'BUSCO_F'))}% · missing ${fmt(V(r, 'BUSCO_M'))}%`},
  {t: 'Merqury QV', zoom: true, get: r => V(r, 'QV'),
   tip: r => `QV <b>${fmt(V(r, 'QV'), 1)}</b> (≈ 1 error per ${fint(Math.pow(10, V(r, 'QV') / 10))} bp)<br>k-mer completeness ${fmt(V(r, 'kmer_compl_pct'))}%`},
];
function bars() {
  const host = document.getElementById('bars'); if (!host || !D.rows.length) return;
  host.innerHTML = '';
  const rows = D.rows.filter(r => !FLT.has(r.label));
  if (!rows.length) { host.innerHTML = '<p class="note">All assemblies hidden by the filter.</p>'; return; }
  const n = rows.length, rh = 20, top = 28, lw = 124, W = 420, bw = W - lw - 44, H = top + n * rh + 20;
  const hl = i => {
    host.classList.toggle('hl', i >= 0);
    host.querySelectorAll('g.row').forEach(g => g.classList.toggle('on', +g.dataset.i === i));
  };
  METRICS.forEach(m => {
    const vals = rows.map(m.get), known = vals.filter(v => v != null);
    let lo = 0, hi = known.length ? Math.max(...known) : 1;
    if (m.zoom && known.length) {
      lo = Math.max(0, Math.floor(Math.min(...known) - 3));
      hi = m.max ? m.max : Math.ceil(Math.max(...known) + 1);
    } else hi = hi > 0 ? hi * 1.05 : 1;
    const sx = v => lw + (Math.min(Math.max(v, lo), hi) - lo) / (hi - lo) * bw;
    const box = document.createElement('div'); box.className = 'panel'; host.appendChild(box);
    const svg = el('svg', {viewBox: `0 0 ${W} ${H}`}, box);
    el('text', {x: 0, y: 14, class: 'ptitle'}, svg).textContent = m.t;
    niceTicks(lo, hi, 4).forEach(t => {
      const x = sx(t);
      el('line', {x1: x, x2: x, y1: top - 4, y2: top + n * rh, class: 'grid'}, svg);
      el('text', {x, y: top + n * rh + 14, class: 'tick', 'text-anchor': 'middle'}, svg).textContent = short(t);
    });
    el('line', {x1: lw, x2: lw, y1: top - 4, y2: top + n * rh, class: 'axis'}, svg);
    rows.forEach((r, i) => {
      const y = top + i * rh, v = vals[i];
      if (i && r.group !== rows[i - 1].group) el('line', {x1: 0, x2: W, y1: y, y2: y, class: 'sep'}, svg);   // next assembly
      const g = el('g', {class: 'row', 'data-i': i}, svg);
      el('rect', {x: 0, y, width: W, height: rh, class: 'hit'}, g);
      el('text', {x: lw - 6, y: y + rh / 2 + 4, 'text-anchor': 'end', class: 'lab'}, g).textContent = r.label;
      if (v != null) {
        el('rect', {x: lw, y: y + 3, width: Math.max(sx(v) - lw, 1), height: rh - 6, rx: 2, fill: r.color,
                    'fill-opacity': r.hap === 'primary' ? 0.4 : 0.95}, g);
        el('text', {x: sx(v) + 4, y: y + rh / 2 + 4, class: 'val'}, g).textContent = short(v);
      } else el('text', {x: lw + 4, y: y + rh / 2 + 4, class: 'na'}, g).textContent = 'NA';
      g.addEventListener('mousemove', ev => { hl(i); showTip(`${sw(r.color)}<b>${escH(r.label)}</b><br>${v == null ? 'not computed for this assembly' : m.tip(r)}`, ev); });
      g.addEventListener('mouseleave', () => { hl(-1); hideTip(); });
    });
  });
}

/* ---------------------------------------------------------------- Nx curves */
const nxHidden = new Set();                         // legend toggles, kept across filter redraws
function nxCurves() {
  const host = document.getElementById('nx');
  if (!host || !D.rows.some(r => D.nx[r.label])) return;
  host.innerHTML = '';
  const rows = D.rows.filter(r => D.nx[r.label] && !FLT.has(r.label));
  if (!rows.length) { host.innerHTML = '<p class="note">All assemblies hidden by the filter.</p>'; return; }
  const all = [];
  rows.forEach(r => ['scaf', 'ctg'].forEach(k => D.nx[r.label][k].forEach(p => all.push(p[1]))));
  if (!all.length) return;
  const l0 = Math.floor(Math.log10(Math.max(1, Math.min(...all)))), l1 = Math.max(l0 + 1, Math.ceil(Math.log10(Math.max(...all))));
  const W = 560, H = 340, ml = 58, mr = 14, mt = 28, mb = 40, pw = W - ml - mr, ph = H - mt - mb;
  const sx = x => ml + x / 100 * pw, sy = y => mt + ph - (Math.log10(Math.max(y, 1)) - l0) / (l1 - l0) * ph;
  const hidden = nxHidden, paths = [];
  const at = (pts, x) => { let a = 0, b = pts.length - 1; while (a < b) { const c = (a + b + 1) >> 1; if (pts[c][0] <= x) a = c; else b = c - 1; } return pts[a]; };
  [['scaf', 'Scaffold Nx'], ['ctg', 'Contig Nx (split at N gaps)']].forEach(([k, title]) => {
    const box = document.createElement('div'); box.className = 'panel'; host.appendChild(box);
    const svg = el('svg', {viewBox: `0 0 ${W} ${H}`}, box);
    el('text', {x: 0, y: 14, class: 'ptitle'}, svg).textContent = title;
    for (let e = l0; e <= l1; e++) {
      for (const m of [1, 2, 5]) {
        const v = m * Math.pow(10, e); if (v < Math.pow(10, l0) || v > Math.pow(10, l1)) continue;
        el('line', {x1: ml, x2: ml + pw, y1: sy(v), y2: sy(v), class: 'grid'}, svg);
        if (m === 1) el('text', {x: ml - 6, y: sy(v) + 3, class: 'tick', 'text-anchor': 'end'}, svg).textContent = flen(v).replace('.00', '').replace('.0 ', ' ');
      }
    }
    for (let x = 0; x <= 100; x += 10) {
      el('line', {x1: sx(x), x2: sx(x), y1: mt, y2: mt + ph, class: 'grid'}, svg);
      el('text', {x: sx(x), y: mt + ph + 14, class: 'tick', 'text-anchor': 'middle'}, svg).textContent = x;
    }
    el('line', {x1: sx(50), x2: sx(50), y1: mt, y2: mt + ph, class: 'n50'}, svg);
    el('rect', {x: ml, y: mt, width: pw, height: ph, fill: 'none', class: 'axis'}, svg);
    el('text', {x: ml + pw / 2, y: H - 6, class: 'tick', 'text-anchor': 'middle'}, svg).textContent = 'x (% of assembly)';
    const yl = el('text', {x: 12, y: mt + ph / 2, class: 'tick', 'text-anchor': 'middle', transform: `rotate(-90 12 ${mt + ph / 2})`}, svg);
    yl.textContent = 'sequence length';
    const lines = rows.map(r => {
      const pts = D.nx[r.label][k];
      let d = '';
      pts.forEach((p, i) => { d += i ? `H${sx(p[0]).toFixed(1)}V${sy(p[1]).toFixed(1)}` : `M${sx(p[0]).toFixed(1)},${sy(p[1]).toFixed(1)}`; });
      const path = el('path', {d, fill: 'none', stroke: r.color, 'stroke-width': 2, 'stroke-dasharray': dash(r.hap)}, svg);
      if (hidden.has(r.label)) path.style.display = 'none';
      paths.push({label: r.label, path});
      return {r, pts, dot: el('circle', {r: 3.2, fill: r.color, stroke: '#fff', 'stroke-width': 1, visibility: 'hidden'}, svg)};
    });
    const cross = el('line', {y1: mt, y2: mt + ph, class: 'cross', visibility: 'hidden'}, svg);
    const hit = el('rect', {x: ml, y: mt, width: pw, height: ph, fill: 'transparent'}, svg);
    hit.addEventListener('mousemove', ev => {
      const p = svgPoint(svg, ev), x = Math.min(100, Math.max(0, (p.x - ml) / pw * 100));
      cross.setAttribute('x1', sx(x)); cross.setAttribute('x2', sx(x)); cross.setAttribute('visibility', 'visible');
      const vals = [];
      lines.forEach(l => {
        if (hidden.has(l.r.label) || !l.pts.length) { l.dot.setAttribute('visibility', 'hidden'); return; }
        const q = at(l.pts, x);
        l.dot.setAttribute('cx', sx(x)); l.dot.setAttribute('cy', sy(q[1])); l.dot.setAttribute('visibility', 'visible');
        vals.push([l.r, q]);
      });
      vals.sort((a, b) => b[1][1] - a[1][1]);
      showTip(`<b>N${x.toFixed(0)}</b> — ${k === 'scaf' ? 'scaffolds' : 'contigs'}<table>` + vals.map(([r, q]) =>
        `<tr><td>${sw(r.color)}${escH(r.label)}</td><td class="num"><b>${flen(q[1])}</b></td><td class="num" style="color:#888">#${q[2].toLocaleString()}</td></tr>`).join('') + '</table>', ev);
    });
    hit.addEventListener('mouseleave', () => { hideTip(); cross.setAttribute('visibility', 'hidden'); lines.forEach(l => l.dot.setAttribute('visibility', 'hidden')); });
  });
  const leg = document.createElement('div'); leg.className = 'legend'; host.appendChild(leg);
  rows.forEach(r => {
    const s = document.createElement('span');
    s.classList.toggle('off', hidden.has(r.label));
    s.innerHTML = `<svg viewBox="0 0 26 10"><line x1="1" y1="5" x2="25" y2="5" stroke="${r.color}" stroke-width="2.2" stroke-dasharray="${dash(r.hap)}"/></svg>${escH(r.label)}`;
    s.addEventListener('click', () => {
      const off = !hidden.has(r.label);
      off ? hidden.add(r.label) : hidden.delete(r.label);
      s.classList.toggle('off', off);
      paths.filter(p => p.label === r.label).forEach(p => p.path.style.display = off ? 'none' : '');
    });
    leg.appendChild(s);
  });
}

/* ---------------------------------------------------------------- dot plots (canvas, zoomable) */
function dotplot(host, d) {
  const T = d.tseq, Q = d.qseq;
  const end = a => a.length ? a[a.length - 1][1] + a[a.length - 1][2] : 1;
  const tTot = end(T), qTot = end(Q);
  const S = d.segs.map(s => {
    const [ti, qi, ts, te, qs, qe, plus] = s, t = T[ti], q = Q[qi], L = q[2], flip = q[3] === 1;
    const p = (plus ? [[ts, qs], [te, qe]] : [[ts, qe], [te, qs]]).map(([x, y]) => [t[1] + x, q[1] + (flip ? L - y : y)]);
    return {x0: p[0][0], y0: p[0][1], x1: p[1][0], y1: p[1][1], fwd: (plus === 1) !== flip, s};
  });
  const cv = document.createElement('canvas'), sel = document.createElement('div'), btn = document.createElement('button');
  sel.className = 'sel'; btn.className = 'reset'; btn.textContent = 'reset zoom';
  host.append(cv, sel, btn);
  const ML = 62, MR = 96, MT = 76, MB = 44;
  let full = {x0: 0, x1: tTot, y0: 0, y1: qTot}, view = {...full}, size = 0, dpr = 1, hover = null, drag = null;
  let hl = null;                                     // {t, q}: target / query sequence under the cursor
  // contigs of scaffolded assemblies (YaHS AGP), in plot coordinates: strips along the top and right edges
  const SW = 6, ctgs = (seqs, cs, isQ) => !cs ? null : seqs.flatMap((q, i) => (cs[i] || []).map((c, k) => {
    const [b, e] = isQ && q[3] === 1 ? [q[2] - c[1], q[2] - c[0]] : [c[0], c[1]];
    return {a: q[1] + b, b: q[1] + e, name: c[2], o: c[3], k, n: cs[i].length, seq: q};
  }));
  const TC = ctgs(T, d.tctg, false), QC = ctgs(Q, d.qctg, true);
  let hc = null;                                     // contig under the cursor on a strip
  const TL = TC ? SW + 4 : 0, QL = QC ? SW + 4 : 0;  // label offsets past the strips
  const pw = () => size - ML - MR, ph = () => size - MT - MB;
  const X = x => ML + (x - view.x0) / (view.x1 - view.x0) * pw(), Y = y => MT + ph() - (y - view.y0) / (view.y1 - view.y0) * ph();
  const iX = px => view.x0 + (px - ML) / pw() * (view.x1 - view.x0), iY = py => view.y0 + (MT + ph() - py) / ph() * (view.y1 - view.y0);
  const trunc = (s, n) => s.length > n ? s.slice(0, n - 1) + '…' : s;
  const find = (a, v) => { let lo = 0, hi = a.length - 1; while (lo < hi) { const m = (lo + hi + 1) >> 1; if (a[m][1] <= v) lo = m; else hi = m - 1; } return a[lo]; };
  function draw() {
    const c = cv.getContext('2d');
    c.setTransform(dpr, 0, 0, dpr, 0, 0); c.clearRect(0, 0, size, size);
    c.fillStyle = '#fff'; c.fillRect(0, 0, size, size);
    c.save(); c.beginPath(); c.rect(ML, MT, pw(), ph()); c.clip();
    if (hl) {                                          // highlight the sequences under the cursor
      c.fillStyle = 'rgba(240,160,30,.13)';
      c.fillRect(X(hl.t[1]), MT, X(hl.t[1] + hl.t[2]) - X(hl.t[1]), ph());
      c.fillRect(ML, Y(hl.q[1] + hl.q[2]), pw(), Y(hl.q[1]) - Y(hl.q[1] + hl.q[2]));
    }
    c.strokeStyle = '#e6e8eb'; c.lineWidth = 1; c.beginPath();
    T.forEach(t => { const x = Math.round(X(t[1])) + .5; c.moveTo(x, MT); c.lineTo(x, MT + ph()); });
    Q.forEach(q => { const y = Math.round(Y(q[1])) + .5; c.moveTo(ML, y); c.lineTo(ML + pw(), y); });
    c.stroke();
    if (hc) {                                          // contig hovered on a strip: its span across the plot
      c.fillStyle = 'rgba(63,120,181,.12)';
      if (hc.ax === 't') c.fillRect(X(hc.a), MT, X(hc.b) - X(hc.a), ph());
      else c.fillRect(ML, Y(hc.b), pw(), Y(hc.a) - Y(hc.b));
    }
    c.lineWidth = 1.5; c.lineCap = 'round';
    for (const fwd of [true, false]) {
      c.strokeStyle = fwd ? '#1f5fa8' : '#c8413b'; c.beginPath();
      for (const s of S) if (s.fwd === fwd) { c.moveTo(X(s.x0), Y(s.y0)); c.lineTo(X(s.x1), Y(s.y1)); }
      c.stroke();
    }
    if (hover) {                                       // all alignments of the hovered query-target pair
      c.strokeStyle = 'rgba(240,160,30,.9)'; c.lineWidth = 5; c.beginPath();
      for (const s of S) if (s.s[0] === hover.s[0] && s.s[1] === hover.s[1]) { c.moveTo(X(s.x0), Y(s.y0)); c.lineTo(X(s.x1), Y(s.y1)); }
      c.stroke();
    }
    c.restore();
    c.strokeStyle = '#b9bdc2'; c.lineWidth = 1; c.strokeRect(ML + .5, MT + .5, pw(), ph());
    // contig strips: neighbouring contigs of a scaffold alternate dark / light blue, hovered one orange
    const strip = (cs, ax) => {
      if (!cs) return;
      for (const g of cs) {
        c.fillStyle = hc && hc.g === g ? '#e8a13a' : g.k % 2 ? '#9cc3e6' : '#3f78b5';
        if (ax === 't') {
          const a = Math.max(X(g.a), ML), b = Math.min(X(g.b), ML + pw());
          if (b > a) c.fillRect(a, MT - SW - 2, Math.max(b - a - (b - a > 3 ? 1 : 0), .5), SW);   // 1 px gap = join
        } else {
          const a = Math.min(Y(g.a), MT + ph()), b = Math.max(Y(g.b), MT);
          if (a > b) c.fillRect(ML + pw() + 2, b, SW, Math.max(a - b - (a - b > 3 ? 1 : 0), .5));
        }
      }
    };
    strip(TC, 't'); strip(QC, 'q');
    c.fillStyle = '#6b7078'; c.font = '11px -apple-system, Segoe UI, Roboto, sans-serif';
    const tk = (a, b) => { const st = niceTicks(a / 1e6, b / 1e6, 5); return st.map(t => [t, short(t)]); };
    c.textAlign = 'center';
    tk(view.x0, view.x1).forEach(([t, s]) => { const x = X(t * 1e6); if (x >= ML - 1 && x <= ML + pw() + 1) { c.fillText(s, x, MT + ph() + 15); c.fillRect(x, MT + ph(), 1, 4); } });
    c.textAlign = 'right';
    tk(view.y0, view.y1).forEach(([t, s]) => { const y = Y(t * 1e6); if (y >= MT - 1 && y <= MT + ph() + 1) { c.fillText(s, ML - 7, y + 4); c.fillRect(ML - 4, y, 4, 1); } });
    c.fillStyle = '#3b3f45'; c.font = '12px -apple-system, Segoe UI, Roboto, sans-serif'; c.textAlign = 'center';
    c.fillText(`${d.t} (Mb)`, ML + pw() / 2, size - 8);
    c.save(); c.translate(14, MT + ph() / 2); c.rotate(-Math.PI / 2); c.fillText(`${d.q} (Mb)`, 0, 0); c.restore();
    c.font = '10px -apple-system, Segoe UI, Roboto, sans-serif'; c.fillStyle = '#4a4f56';
    T.forEach(t => {
      const a = Math.max(X(t[1]), ML), b = Math.min(X(t[1] + t[2]), ML + pw());
      if (b - a < 11 && !(hl && hl.t === t)) return;
      const on = hl && hl.t === t;
      c.save(); c.translate((a + b) / 2 + 3, MT - 5 - TL); c.rotate(-Math.PI / 2); c.textAlign = 'left';
      c.font = (on ? 'bold ' : '') + '10px -apple-system, Segoe UI, Roboto, sans-serif'; c.fillStyle = on ? '#9a5a00' : '#4a4f56';
      c.fillText(trunc(t[0], 12), 0, 0); c.restore();
    });
    c.textAlign = 'left';
    Q.forEach(q => {
      const a = Math.min(Y(q[1]), MT + ph()), b = Math.max(Y(q[1] + q[2]), MT);
      if (a - b < 11 && !(hl && hl.q === q)) return;
      const on = hl && hl.q === q;
      c.font = (on ? 'bold ' : '') + '10px -apple-system, Segoe UI, Roboto, sans-serif'; c.fillStyle = on ? '#9a5a00' : '#4a4f56';
      c.fillText(trunc(q[0], 13) + (q[3] ? ' (−)' : ''), ML + pw() + 5 + QL, (a + b) / 2 + 3);
    });
  }
  function layout() {
    const w = Math.max(320, Math.min(host.clientWidth, 760));
    if (w === size) return;
    size = w; dpr = window.devicePixelRatio || 1;
    cv.width = w * dpr; cv.height = w * dpr; cv.style.width = w + 'px'; cv.style.height = w + 'px';
    draw();
  }
  const inside = (x, y) => x >= ML && x <= ML + pw() && y >= MT && y <= MT + ph();
  const pos = ev => { const r = cv.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
  function nearest(px, py) {
    let best = null, bd = 36;                          // within 6 px
    for (const s of S) {
      const ax = X(s.x0), ay = Y(s.y0), bx = X(s.x1), by = Y(s.y1), dx = bx - ax, dy = by - ay;
      const l2 = dx * dx + dy * dy, t = l2 ? Math.max(0, Math.min(1, ((px - ax) * dx + (py - ay) * dy) / l2)) : 0;
      const ex = ax + t * dx - px, ey = ay + t * dy - py, dd = ex * ex + ey * ey;
      if (dd < bd) { bd = dd; best = s; }
    }
    return best;
  }
  cv.addEventListener('mousemove', ev => {
    const [px, py] = pos(ev);
    if (drag) {
      const x = Math.min(drag[0], px), y = Math.min(drag[1], py);
      Object.assign(sel.style, {display: 'block', left: x + 'px', top: y + 'px', width: Math.abs(px - drag[0]) + 'px', height: Math.abs(py - drag[1]) + 'px'});
      hideTip(); return;
    }
    // on a contig strip: that contig
    const onT = TC && px >= ML && px <= ML + pw() && py >= MT - SW - 3 && py < MT;
    const onQ = QC && py >= MT && py <= MT + ph() && px > ML + pw() && px <= ML + pw() + SW + 3;
    if (onT || onQ) {
      const v = onT ? iX(px) : iY(py), g = (onT ? TC : QC).find(g => g.a <= v && v < g.b);
      const nh = g ? {g, ax: onT ? 't' : 'q', a: g.a, b: g.b} : null;
      if ((nh && nh.g) !== (hc && hc.g) || hover || hl) { hc = nh; hover = hl = null; draw(); }
      if (g) showTip(`<b>${escH(g.name)}</b> (${g.o === '-' ? 'reversed' : 'forward'}) ${flen(g.b - g.a)}<br>` +
                     `${escH(g.seq[0])}: contig ${g.k + 1} of ${g.n}`, ev);
      else hideTip();
      return;
    }
    if (hc) { hc = null; draw(); }
    if (!inside(px, py)) { if (hover || hl) { hover = hl = null; draw(); } hideTip(); return; }
    const h = nearest(px, py);
    const nt = h ? T[h.s[0]] : find(T, iX(px)), nq = h ? Q[h.s[1]] : find(Q, iY(py));
    if (h !== hover || !hl || hl.t !== nt || hl.q !== nq) { hover = h; hl = {t: nt, q: nq}; draw(); }
    if (h) {
      // summary of the query-target pair: single alignments are split at the 10 Mb mapping chunks
      const ti = h.s[0], qi = h.s[1], t = T[ti], q = Q[qi];
      let n = 0, bp = 0, inv = 0, q0 = Infinity, q1 = 0, t0 = Infinity, t1 = 0;
      for (const s of S) if (s.s[0] === ti && s.s[1] === qi) {
        const [, , ts, te, qs, qe] = s.s, L = qe - qs;
        n++; bp += L; if (!s.fwd) inv += L;
        q0 = Math.min(q0, qs); q1 = Math.max(q1, qe); t0 = Math.min(t0, ts); t1 = Math.max(t1, te);
      }
      showTip(`<b>${escH(q[0])}</b>${q[3] ? ' (−)' : ''} ${flen(q0)}–${flen(q1)} of ${flen(q[2])}<br>` +
              `<b>${escH(t[0])}</b> ${flen(t0)}–${flen(t1)} of ${flen(t[2])}<br>` +
              `${flen(bp)} aligned (${fmt(100 * bp / q[2], 1)}% of query) · ${n} alignment${n > 1 ? 's' : ''}<br>` +
              `<span style="color:#1f5fa8">collinear ${flen(bp - inv)}</span> · <span style="color:#c8413b">inverted ${flen(inv)}</span>`, ev);
    } else {
      const x = iX(px), y = iY(py), t = find(T, x), q = find(Q, y);
      const qp = q[3] ? q[2] - (y - q[1]) : y - q[1];
      const ct = TC && TC.find(g => g.a <= x && x < g.b), cq = QC && QC.find(g => g.a <= y && y < g.b);
      showTip(`<span style="color:#888">${escH(d.t)}</span> ${escH(t[0])}: ${flen(x - t[1])}${ct ? ` · ${escH(ct.name)}` : ''}<br>` +
              `<span style="color:#888">${escH(d.q)}</span> ${escH(q[0])}${q[3] ? ' (−)' : ''}: ${flen(qp)}${cq ? ` · ${escH(cq.name)}` : ''}`, ev);
    }
  });
  cv.addEventListener('mouseleave', () => { hideTip(); if (hover || hl || hc) { hover = hl = hc = null; draw(); } });
  cv.addEventListener('mousedown', ev => { const p = pos(ev); if (inside(...p)) { drag = p; ev.preventDefault(); } });
  window.addEventListener('mouseup', ev => {
    if (!drag) return;
    const [px, py] = pos(ev), [qx, qy] = drag; drag = null; sel.style.display = 'none';
    if (Math.abs(px - qx) < 6 || Math.abs(py - qy) < 6) return;
    const cx = v => Math.min(ML + pw(), Math.max(ML, v)), cy = v => Math.min(MT + ph(), Math.max(MT, v));
    view = {x0: iX(cx(Math.min(px, qx))), x1: iX(cx(Math.max(px, qx))), y0: iY(cy(Math.max(py, qy))), y1: iY(cy(Math.min(py, qy)))};
    btn.style.display = 'block'; hover = hl = null; draw();
  });
  const reset = () => { view = {...full}; btn.style.display = 'none'; draw(); };
  cv.addEventListener('dblclick', reset); btn.addEventListener('click', reset);
  if (!S.length) host.insertAdjacentHTML('beforeend', '<p class="note">no alignments to show</p>');
  layout();
  if (window.ResizeObserver) new ResizeObserver(layout).observe(host); else window.addEventListener('resize', layout);
  const cap = host.parentElement.querySelector('figcaption');
  if (cap) cap.insertAdjacentHTML('beforeend', ` — ${d.qpct}% of query / ${d.tpct}% of target in alignments ≥ ${flen(d.min_aln)}` +
    (d.shown < d.total ? `; ${d.shown.toLocaleString()} longest of ${d.total.toLocaleString()} alignments drawn` : '') +
    (d.min_seq ? `; sequences < ${flen(d.min_seq)} not drawn` : ''));
}

/* ---------------------------------------------------------------- interactive contact map (HiGlass-like) */
// h: {genome, pairs, scaffolds: [[name, offset, len]], levels: [{n, bp, vmax, data}]} from cifi_stats.py;
// data = gzip + base64 of the upper triangle, log-scaled counts 0-255. Levels are decoded when first needed.
function hicLevel(L) {
  if (!L.p) L.p = (async () => {
    const raw = Uint8Array.from(atob(L.data), c => c.charCodeAt(0));
    const t = new Uint8Array(await new Response(new Blob([raw]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());
    const n = L.n, sq = new Uint8Array(n * n), hist = new Uint32Array(256);
    for (let i = 0, k = 0; i < n; i++) for (let j = i; j < n; j++, k++) { const v = t[k]; sq[i * n + j] = v; sq[j * n + i] = v; hist[v]++; }
    let nz = t.length - hist[0], c = 0, hi = 255;                    // colour scale: 99.5th percentile of non-zero bins
    for (let v = 255; v > 0; v--) { c += hist[v]; if (c >= nz * 0.005) { hi = v; break; } }
    L.hi = Math.max(hi, 1); L.sq = sq;
  })();
  return L.p;
}
function hicMap(host, h) {
  const G = h.genome, SC = h.scaffolds, LV = h.levels;               // LV: coarse -> fine
  const CT = (h.contigs || []).slice().sort((a, b) => a[1] - b[1]);  // [name, genome offset, length, orient]
  const cv = document.createElement('canvas'), btn = document.createElement('button'), res = document.createElement('div');
  btn.className = 'reset'; btn.textContent = 'reset zoom'; res.className = 'res';
  host.append(cv, btn, res);
  if (typeof DecompressionStream === 'undefined') { host.innerHTML = '<p class="note">this browser cannot decompress the contact map; see the PretextSnapshot image</p>'; return; }
  const ML = 14, MR = 104, MT = 84, MB = 14;
  let view = {x0: 0, y0: 0, s: G}, size = 0, dpr = 1, drag = null, moved = false;
  const pw = () => size - ML - MR;
  const X = x => ML + (x - view.x0) / view.s * pw(), Y = y => MT + (y - view.y0) / view.s * pw();
  const iX = px => view.x0 + (px - ML) / pw() * view.s, iY = py => view.y0 + (py - MT) / pw() * view.s;
  const at = (A, g) => { let lo = 0, hi = A.length - 1; while (lo < hi) { const m = (lo + hi + 1) >> 1; if (A[m][1] <= g) lo = m; else hi = m - 1; } return A[lo]; };
  const scafAt = g => at(SC, g);
  const ctgAt = g => { if (!CT.length) return null; const c = at(CT, g); return g < c[1] + c[2] ? c : null; };   // null in a gap
  // coarsest level with bins no larger than ~1.5 px (LV is coarse -> fine); the finest when zoomed in further
  const level = () => { const want = view.s / pw(); return LV.find(l => l.bp <= want * 1.5) || LV[LV.length - 1]; };
  const off = document.createElement('canvas');
  function draw() {
    const L = level(), have = L.sq ? L : [...LV].reverse().find(l => l.bp >= L.bp && l.sq);
    if (!L.sq) hicLevel(L).then(draw);
    const c = cv.getContext('2d'), P = Math.round(pw() * dpr);
    c.setTransform(1, 0, 0, 1, 0, 0); c.fillStyle = '#fff'; c.fillRect(0, 0, cv.width, cv.height);
    if (have) {
      off.width = off.height = P;
      const oc = off.getContext('2d'), img = oc.createImageData(P, P), px = img.data, n = have.n, sq = have.sq, hi = have.hi;
      const bins = new Int32Array(P);
      for (let k = 0; k < P; k++) { const g = view.x0 + (k + .5) / P * view.s; bins[k] = g < 0 || g >= G ? -1 : Math.min(n - 1, Math.floor(g / have.bp)); }
      const rows = new Int32Array(P);
      for (let k = 0; k < P; k++) { const g = view.y0 + (k + .5) / P * view.s; rows[k] = g < 0 || g >= G ? -1 : Math.min(n - 1, Math.floor(g / have.bp)); }
      for (let r = 0, o = 0; r < P; r++) {
        const bi = rows[r];
        for (let k = 0; k < P; k++, o += 4) {
          const bj = bins[k];
          if (bi < 0 || bj < 0) { px[o] = px[o + 1] = px[o + 2] = 247; px[o + 3] = 255; continue; }
          const f = Math.min(1, sq[bi * n + bj] / hi), w = Math.round(255 * (1 - f));
          px[o] = 255 - Math.round(f * 40); px[o + 1] = w; px[o + 2] = w; px[o + 3] = 255;    // white -> red
        }
      }
      oc.putImageData(img, 0, 0);
      c.imageSmoothingEnabled = false; c.drawImage(off, ML * dpr, MT * dpr);
    }
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.save(); c.beginPath(); c.rect(ML, MT, pw(), pw()); c.clip();
    // contig boundaries inside scaffolds (dashed), drawn where neighbouring boundaries are >= 3 px apart
    const sStart = new Set(SC.map(x => x[1]));
    c.strokeStyle = 'rgba(30,60,110,.30)'; c.lineWidth = 1; c.setLineDash([3, 3]); c.beginPath();
    let last = -1e9;
    for (const [, o] of CT) {
      if (sStart.has(o)) continue;
      const x = X(o); if (x < ML || x > ML + pw() || x - last < 3) continue; last = x;
      const xr = Math.round(x) + .5, yr = Math.round(Y(o)) + .5;
      c.moveTo(xr, MT); c.lineTo(xr, MT + pw()); c.moveTo(ML, yr); c.lineTo(ML + pw(), yr);
    }
    c.stroke(); c.setLineDash([]);
    c.strokeStyle = 'rgba(30,60,110,.45)'; c.lineWidth = 1; c.beginPath();
    for (const [, o] of SC) { const x = Math.round(X(o)) + .5, y = Math.round(Y(o)) + .5; c.moveTo(x, MT); c.lineTo(x, MT + pw()); c.moveTo(ML, y); c.lineTo(ML + pw(), y); }
    c.stroke(); c.restore();
    c.strokeStyle = '#b9bdc2'; c.strokeRect(ML + .5, MT + .5, pw(), pw());
    c.font = '10px -apple-system, Segoe UI, Roboto, sans-serif'; c.fillStyle = '#4a4f56';
    for (const [nm, o, len] of SC) {
      const a = Math.max(X(o), ML), b = Math.min(X(o + len), ML + pw());
      if (b - a >= 11) { c.save(); c.translate((a + b) / 2 + 3, MT - 5); c.rotate(-Math.PI / 2); c.textAlign = 'left'; c.fillText(nm.replace(/^scaffold_/, 's'), 0, 0); c.restore(); }
      const u = Math.max(Y(o), MT), v = Math.min(Y(o + len), MT + pw());
      if (v - u >= 11) { c.textAlign = 'left'; c.fillText(nm.replace(/^scaffold_/, 's'), ML + pw() + 5, (u + v) / 2 + 3); }
    }
    res.textContent = `${flen((have || L).bp)} bins` + (L.sq ? '' : ' · loading finer…');
  }
  function layout() {
    const w = Math.max(320, Math.min(host.clientWidth, 820));
    if (w === size) return;
    size = w; dpr = window.devicePixelRatio || 1;
    cv.width = cv.height = w * dpr; cv.style.width = cv.style.height = w + 'px';
    draw();
  }
  const pos = ev => { const r = cv.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
  const inside = (x, y) => x >= ML && x <= ML + pw() && y >= MT && y <= MT + pw();
  const clamp = () => { view.s = Math.min(G, Math.max(view.s, LV[LV.length - 1].bp * 20));
    view.x0 = Math.min(Math.max(view.x0, 0), G - view.s); view.y0 = Math.min(Math.max(view.y0, 0), G - view.s); };
  const zoomed = () => { btn.style.display = view.s < G ? 'block' : 'none'; };
  cv.addEventListener('wheel', ev => {
    const [px, py] = pos(ev); if (!inside(px, py)) return;
    ev.preventDefault(); hideTip();
    const gx = iX(px), gy = iY(py), f = Math.exp(ev.deltaY * (ev.deltaMode ? 0.05 : 0.0015));
    view.s *= f; view.x0 = gx - (px - ML) / pw() * view.s; view.y0 = gy - (py - MT) / pw() * view.s;
    clamp(); zoomed(); draw();
  }, {passive: false});
  cv.addEventListener('mousedown', ev => { const p = pos(ev); if (inside(...p)) { drag = {p, x0: view.x0, y0: view.y0}; moved = false; cv.classList.add('drag'); ev.preventDefault(); } });
  window.addEventListener('mouseup', () => { drag = null; cv.classList.remove('drag'); });
  cv.addEventListener('mousemove', ev => {
    const [px, py] = pos(ev);
    if (drag) {
      const k = view.s / pw(); view.x0 = drag.x0 - (px - drag.p[0]) * k; view.y0 = drag.y0 - (py - drag.p[1]) * k;
      moved = true; clamp(); hideTip(); draw(); return;
    }
    if (!inside(px, py)) { hideTip(); return; }
    const gx = iX(px), gy = iY(py), L = [...LV].reverse().find(l => l.bp >= level().bp && l.sq);
    const sx = scafAt(gx), sy = scafAt(gy);
    let cnt = '';
    if (L) {
      const v = L.sq[Math.min(L.n - 1, Math.floor(gy / L.bp)) * L.n + Math.min(L.n - 1, Math.floor(gx / L.bp))];
      cnt = `<br>contacts ≈ <b>${fint(Math.expm1(v * Math.log1p(L.vmax) / 255))}</b> in ${flen(L.bp)} × ${flen(L.bp)}`;
    }
    const cx = ctgAt(gx), cy = ctgAt(gy);
    const ctg = (c, g) => !CT.length ? '' : c ? ` <span style="color:#888">· ${escH(c[0])} (${c[3]}) ${flen(c[3] === '-' ? c[1] + c[2] - g : g - c[1])}</span>` : ' <span style="color:#888">· gap</span>';
    showTip(`<b>${escH(sx[0])}</b> ${flen(gx - sx[1])}${ctg(cx, gx)}<br><b>${escH(sy[0])}</b> ${flen(gy - sy[1])}${ctg(cy, gy)}${cnt}` +
            (sx === sy ? '' : '<br><span style="color:#888">between scaffolds</span>'), ev);
  });
  cv.addEventListener('mouseleave', hideTip);
  const reset = () => { view = {x0: 0, y0: 0, s: G}; zoomed(); draw(); };
  cv.addEventListener('dblclick', reset); btn.addEventListener('click', reset);
  hicLevel(LV[0]).then(() => { layout(); draw(); });
  if (window.ResizeObserver) new ResizeObserver(layout).observe(host); else window.addEventListener('resize', layout);
}

/* ---------------------------------------------------------------- scaffold layout (YaHS AGP) */
// d.rows: [name, length, contigs [[ctg, start, end, orient, cstart, cend]], gaps [[start, end]], telo_start, telo_end]
const SC = {W: 560, LW: 118, RW: 96, RH: 16, BH: 10, C: ['#3f78b5', '#9cc3e6'], TEL: '#d62728'};
const teloTxt = (a, b) => a && b ? 'both ends (T2T)' : a ? 'start' : b ? 'end' : 'none';
// one scaffold row at y, sc = px per bp. o: label, labFill, flip (draw reverse-complemented), tip (extra tooltip line),
// miss [left, right] (the other haplotype has a telomere at that drawn end, this one not), diff (contig count differs)
function scafRow(svg, r, y, sc, o = {}) {
  const {RH, BH, C, TEL} = SC, W = o.W || SC.W, LW = o.LW || SC.LW, [nm, L, ctgs, gaps, ts, te] = r, f = !!o.flip;
  const X = p => LW + (f ? L - p : p) * sc;                // scaffold position -> x, mirrored when flipped
  const g = el('g', {class: 'row'}, svg);
  el('rect', {x: 0, y, width: W, height: RH, fill: 'transparent'}, g);
  const lab = el('text', {x: LW - 8, y: y + RH / 2 + 4, 'text-anchor': 'end', class: 'lab', 'font-size': 11, fill: o.labFill || '#3b3f45'}, g);
  lab.textContent = o.label || nm.replace(/^scaffold_/, 's');
  const len = el('text', {x: LW + L * sc + 9, y: y + RH / 2 + 4, 'font-size': 10, fill: '#6b7078'}, g);
  len.textContent = `${flen(L)} · `;
  el('tspan', o.diff ? {fill: TEL, 'font-weight': 700} : {}, len).textContent = `${ctgs.length} ctg`;
  ctgs.forEach((c, k) => {
    const a = X(c[1]), b = X(c[2]), x = Math.min(a, b), w = Math.max(Math.abs(b - a), .6);
    const rc = el('rect', {x, y: y + (RH - BH) / 2, width: w, height: BH, fill: C[k % 2]}, g);
    if (w > 12) {                                         // orientation arrow
      const cx = x + w / 2, cy = y + RH / 2, d = (c[3] === '-') !== f ? -1 : 1;
      el('path', {d: `M${cx - 3 * d},${cy - 3} L${cx + 3 * d},${cy} L${cx - 3 * d},${cy + 3}`, fill: 'none', stroke: '#fff', 'stroke-width': 1.3}, g);
    }
    rc.addEventListener('mousemove', ev => showTip(
      `<b>${escH(c[0])}</b> (${c[3]}) · ${flen(c[2] - c[1])}<br>contig ${fint(c[4] + 1)}–${fint(c[5])}` +
      (c[5] - c[4] < c[2] - c[1] || c[4] > 0 ? ' (part)' : '') + `<br>${escH(nm)} ${flen(c[1])}–${flen(c[2])} · contig ${k + 1} of ${ctgs.length}` +
      (f ? '<br>drawn reverse-complemented' : ''), ev));
    rc.addEventListener('mouseleave', hideTip);
  });
  gaps.forEach(([a, b]) => el('rect', {x: X((a + b) / 2) - .6, y: y + 1, width: 1.2, height: RH - 2, fill: '#2b2f35'}, g));
  const [tl, tr] = f ? [te, ts] : [ts, te], miss = o.miss || [];
  [[tl, LW - 5, miss[0]], [tr, LW + L * sc, miss[1]]].forEach(([t, x, m]) => {
    if (t) el('rect', {x, y: y + 1, width: 5, height: RH - 2, fill: TEL}, g);
    else if (m) el('rect', {x: x + .6, y: y + 1.6, width: 3.8, height: RH - 3.2, fill: 'none', stroke: TEL, 'stroke-width': 1.2}, g);
  });
  lab.addEventListener('mousemove', ev => showTip(`<b>${escH(nm)}</b> · ${flen(L)}<br>${ctgs.length} contig${ctgs.length > 1 ? 's' : ''}, ` +
    `${gaps.length} gap${gaps.length === 1 ? '' : 's'}<br>telomere: ${teloTxt(ts, te)}` + (o.tip ? `<br>${o.tip}` : ''), ev));
  lab.addEventListener('mouseleave', hideTip);
}
function moreButton(host, showAll, rest, txt, onclick) {   // show all / only large toggle, above a figure
  const b = document.createElement('button'); b.className = 'more';
  b.textContent = showAll ? txt[0] : txt[1]; b.onclick = onclick; host.appendChild(b);
}
function scafLayout(host, d) {
  const all = d.rows, big = all.filter(r => r[1] >= 1e6), rest = all.length - big.length;
  const restBp = all.filter(r => r[1] < 1e6).reduce((a, r) => a + r[1], 0);
  let showAll = !big.length;
  const {W, LW, RW, RH} = SC, bw = W - LW - RW;
  function draw() {
    host.innerHTML = '';
    if (rest && big.length) moreButton(host, showAll, rest, ['show only scaffolds ≥ 1 Mb',
      `show all (+${rest.toLocaleString()} scaffolds < 1 Mb, ${flen(restBp)})`], () => { showAll = !showAll; draw(); });
    const rows = showAll ? all : big, maxL = Math.max(...rows.map(r => r[1]), 1);
    const svg = el('svg', {viewBox: `0 0 ${W} ${rows.length * RH + 4}`}, host);
    rows.forEach((r, i) => scafRow(svg, r, i * RH + 2, bw / maxL));
  }
  draw();
}
// hap1 and hap2 per matched chromosome. d: {h1, h2: rows, groups: [[hap1 index or -1, [[hap2 index, flipped], ...]]], by}
function pairLayout(host, d) {
  const {RW, RH} = SC, LW = 150, GAP = 10;               // LW: chromosome label + hap / scaffold labels
  let W = SC.W, bw = W - LW - RW;
  const ends = (r, f) => f ? [r[5], r[4]] : [r[4], r[5]];            // telomeres at the drawn left / right end
  const sn = nm => nm.replace(/^scaffold_/, 's');
  const isBig = ([a, hs]) => (a >= 0 && d.h1[a][1] >= 1e6) || hs.some(([b]) => d.h2[b][1] >= 1e6);
  const big = d.groups.filter(isBig), rest = d.groups.length - big.length;
  let showAll = !big.length;
  function draw() {
    W = Math.max(SC.W, Math.floor(host.clientWidth || SC.W)); bw = W - LW - RW;
    host.innerHTML = '';
    if (rest && big.length) moreButton(host, showAll, rest, ['show only chromosomes with a scaffold ≥ 1 Mb',
      `show all (+${rest.toLocaleString()} smaller scaffolds)`], () => { showAll = !showAll; draw(); });
    const items = [], bands = [];
    let y = 2, chr = 0;
    (showAll ? d.groups : big).forEach(([a, hs], gi) => {
      const y0 = y, r1 = a >= 0 ? d.h1[a] : null, one = !!r1 && hs.length === 1;
      const r2 = one ? d.h2[hs[0][0]] : null, e1 = r1 && ends(r1, 0), e2 = one && ends(r2, hs[0][1]);
      const nc = one && r1[2].length !== r2[2].length;
      if (r1) {
        items.push({r: r1, y, o: {label: 'hap1 ' + sn(r1[0]), miss: one ? [!e1[0] && e2[0], !e1[1] && e2[1]] : null, diff: nc,
                    tip: hs.length ? `hap2: ${hs.map(([b, f]) => escH(d.h2[b][0]) + (f ? ' (−)' : '')).join(', ')}` : 'no hap2 scaffold matched'}});
        y += RH;
      }
      hs.forEach(([b, f]) => {
        items.push({r: d.h2[b], y, o: {label: 'hap2 ' + sn(d.h2[b][0]) + (f ? ' −' : ''), labFill: '#5d6168', flip: f, diff: nc,
                    miss: one ? [!e2[0] && e1[0], !e2[1] && e1[1]] : null,
                    tip: (r1 ? `matched to hap1 ${escH(r1[0])}` : 'no hap1 match') + (f ? ', drawn reverse-complemented' : '')}});
        y += RH;
      });
      bands.push([y0 - 2, y - y0 + 4, r1 ? `chr ${++chr}` : 'unplaced', gi % 2 === 0]);
      y += GAP;
    });
    const maxL = Math.max(1, ...items.map(it => it.r[1]));
    const svg = el('svg', {viewBox: `0 0 ${W} ${y}`}, host);
    bands.forEach(([by, bh, lab, shade]) => {           // one band per chromosome: label at the left
      if (shade) el('rect', {x: 0, y: by, width: W, height: bh, rx: 3, style: 'fill:var(--band)', opacity: .55}, svg);
      el('line', {x1: 0, x2: W, y1: by + bh + GAP / 2, y2: by + bh + GAP / 2, stroke: '#e1e2de'}, svg);
      el('text', {x: 6, y: by + bh / 2 + 4, 'font-size': 12, 'font-weight': 700, fill: lab === 'unplaced' ? '#8a8f96' : '#1d1f21'}, svg).textContent = lab;
    });
    items.forEach(it => scafRow(svg, it.r, it.y, bw / maxL, Object.assign({W, LW}, it.o)));
  }
  draw();
  let lastW = W;                                          // redraw when the page width changes
  if (window.ResizeObserver) new ResizeObserver(() => {
    const w = Math.max(SC.W, Math.floor(host.clientWidth || SC.W));
    if (w !== lastW) { lastW = w; draw(); }
  }).observe(host);
}

/* ---------------------------------------------------------------- scaffold scatter: length vs contigs */
// all runs in one plot: colour = run, fill = telomeres (filled both ends, half one, hollow none); hap1 and
// hap2 alike (haplotype in the tooltip). Legend entries (one per run) hide / show a run.
const scHidden = new Set();                                // legend toggles, kept across filter redraws
// colour-blind safe (Okabe-Ito), distinct at marker size; run colours of the other plots are too dark / similar here
const SCCOL = ['#0072b2', '#e69f00', '#009e73', '#cc79a7', '#56b4e9', '#d55e00', '#7f7f7f', '#000000'];
function scMark(svg, x, y, tel, col, size = 1) {
  // circle in the run colour: filled = telomeres at both ends, half-filled = one end, hollow = none
  const r = 4.8 * size, g = el('g', {}, svg);
  // semi-transparent fill and outline so that overlapping points stay visible
  el('circle', {cx: x, cy: y, r, fill: tel === 2 ? col : 'transparent', 'fill-opacity': .5, stroke: col,
                'stroke-opacity': .85, 'stroke-width': 1.6 * size}, g);
  if (tel === 1) el('path', {d: `M${x},${y - r} A${r},${r} 0 0 0 ${x},${y + r} Z`, fill: col, 'fill-opacity': .5}, g);   // left half
  return g;
}
function scafScatter() {
  const host = document.getElementById('scsc'); if (!host) return;
  host.innerHTML = '';
  const pts = [], sets = [];
  Object.values(D.scaf).forEach(s => {
    const r = D.rows.find(r => r.label === s.asm);
    if (!r || FLT.has(s.asm)) return;
    sets.push({s, r});
    s.rows.forEach(q => { if (q[1] >= 1e6) pts.push({s, r, q, n: q[2].length, tel: !!q[4] + !!q[5]}); });
  });
  if (!pts.length) { host.innerHTML = '<p class="note">No scaffolds ≥ 1 Mb in the assemblies left by the filter.</p>'; return; }
  const W = 640, H = 380, ml = 54, mr = 14, mt = 28, mb = 40, pw = W - ml - mr, ph = H - mt - mb;
  const lx = pts.map(p => Math.log10(p.q[1])), x0 = Math.floor(Math.min(...lx) * 2) / 2, x1 = Math.max(x0 + .5, Math.ceil(Math.max(...lx) * 2) / 2);
  const ySteps = [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000], nMax = Math.max(...pts.map(p => p.n));
  const yTop = ySteps.find(v => v >= nMax * 1.15) || nMax * 1.15;
  const sx = L => ml + (Math.log10(L) - x0) / (x1 - x0) * pw;
  const sy = n => mt + ph - 10 - Math.log(n) / Math.log(yTop) * (ph - 10);   // 10 px above the axis for 1 contig
  const box = document.createElement('div'); box.className = 'panel'; host.appendChild(box);
  const svg = el('svg', {viewBox: `0 0 ${W} ${H}`}, box);
  el('text', {x: 0, y: 14, class: 'ptitle'}, svg).textContent = 'Scaffolds ≥ 1 Mb: length vs contigs';
  for (let e = Math.floor(x0); e <= Math.ceil(x1); e++) for (const m of [1, 2, 5]) {
    const v = m * Math.pow(10, e), lv = Math.log10(v); if (lv < x0 - 1e-9 || lv > x1 + 1e-9) continue;
    el('line', {x1: sx(v), x2: sx(v), y1: mt, y2: mt + ph, class: 'grid'}, svg);
    el('text', {x: sx(v), y: mt + ph + 14, class: 'tick', 'text-anchor': 'middle'}, svg).textContent = flen(v).replace('.00', '');
  }
  [1, ...ySteps].filter(v => v <= yTop).forEach(v => {
    el('line', {x1: ml, x2: ml + pw, y1: sy(v), y2: sy(v), class: 'grid'}, svg);
    el('text', {x: ml - 6, y: sy(v) + 3, class: 'tick', 'text-anchor': 'end'}, svg).textContent = v;
  });
  el('rect', {x: ml, y: mt, width: pw, height: ph, fill: 'none', class: 'axis'}, svg);
  el('text', {x: ml + pw / 2, y: H - 6, class: 'tick', 'text-anchor': 'middle'}, svg).textContent = 'scaffold length (log)';
  el('text', {x: 12, y: mt + ph / 2, class: 'tick', 'text-anchor': 'middle', transform: `rotate(-90 12 ${mt + ph / 2})`}, svg).textContent = 'contigs per scaffold (log)';
  const runs = [...new Set(Object.values(D.scaf).map(s => s.run))];   // colour by run order, stable under the filter
  const rcol = Object.fromEntries(runs.map((r, i) => [r, SCCOL[i % SCCOL.length]]));
  const marks = [];
  pts.sort((a, b) => a.tel - b.tel).forEach((p, i) => {   // T2T drawn last (on top)
    const x = sx(p.q[1]), y = sy(p.n) + ((i * 7919) % 9 - 4) * .7;    // small fixed jitter: many scaffolds share a count
    const m = scMark(svg, x, y, p.tel, rcol[p.s.run]);
    if (scHidden.has(p.s.run)) m.style.display = 'none';
    marks.push([p.s.run, m]);
    m.addEventListener('mousemove', ev => showTip(`${sw(rcol[p.s.run])}<b>${escH(p.s.run)}</b> ${p.s.hap} · ${escH(p.q[0])}<br>${flen(p.q[1])} · ` +
      `${p.n} contig${p.n > 1 ? 's' : ''}, ${p.q[3].length} gap${p.q[3].length === 1 ? '' : 's'}<br>telomere: ${teloTxt(p.q[4], p.q[5])}`, ev));
    m.addEventListener('mouseleave', hideTip);
  });
  // legend: one entry per run (click to hide / show) with its T2T counts, then the shape / fill key
  const leg = document.createElement('div'); leg.className = 'legend'; host.appendChild(leg);
  const ico = (tel, col) => { const s = document.createElementNS(NS, 'svg'); s.setAttribute('viewBox', '0 0 26 12'); scMark(s, 13, 6, tel, col, .85); return s; };
  [...new Set(sets.map(x => x.s.run))].forEach(run => {
    const rp = pts.filter(p => p.s.run === run), col = rcol[run];
    const t2t = h => `${rp.filter(p => p.s.hap === h && p.tel === 2).length}/${rp.filter(p => p.s.hap === h).length}`;
    const sp = document.createElement('span'); sp.classList.toggle('off', scHidden.has(run));
    sp.append(ico(2, col), `${run} (T2T hap1 ${t2t('hap1')}, hap2 ${t2t('hap2')})`);
    sp.addEventListener('click', () => {
      const off = !scHidden.has(run); off ? scHidden.add(run) : scHidden.delete(run);
      sp.classList.toggle('off', off); marks.filter(([r]) => r === run).forEach(([, m]) => m.style.display = off ? 'none' : '');
    });
    leg.appendChild(sp);
  });
  const key = document.createElement('div'); key.className = 'legend static'; host.appendChild(key);
  [[2, 'telomeres: both ends (T2T)'], [1, 'one end'], [0, 'none']].forEach(([t, txt]) => {
    const sp = document.createElement('span'); sp.append(ico(t, '#6b7078'), txt); key.appendChild(sp);
  });
}

/* ---------------------------------------------------------------- generic line chart (CiFi) */
function lineChart(host, series, o) {
  // series: [{label, color, dash, pts: [[x, y], ...]}]; o: {xlog, ylog, xfmt, yfmt, xlabel, ylabel, title, tipHead}
  host.innerHTML = '';
  const all = series.flatMap(s => s.pts).filter(p => (!o.xlog || p[0] > 0) && (!o.ylog || p[1] > 0));
  if (!all.length) { host.innerHTML = '<p class="note">no data</p>'; return; }
  const W = 560, H = 330, ml = 62, mr = 14, mt = 26, mb = 40, pw = W - ml - mr, ph = H - mt - mb;
  const tx = v => o.xlog ? Math.log10(v) : v, ty = v => o.ylog ? Math.log10(v) : v;
  let x0 = Math.min(...all.map(p => tx(p[0]))), x1 = Math.max(...all.map(p => tx(p[0])));
  let y0 = Math.min(...all.map(p => ty(p[1]))), y1 = Math.max(...all.map(p => ty(p[1])));
  if (o.xlog) { x0 = Math.floor(x0); x1 = Math.ceil(x1); } if (o.ylog) { y0 = Math.floor(y0); y1 = Math.ceil(y1); }
  if (!o.ylog) y0 = Math.min(0, y0);
  if (x1 === x0) x1 = x0 + 1; if (y1 === y0) y1 = y0 + 1;
  const sx = v => ml + (tx(v) - x0) / (x1 - x0) * pw, sy = v => mt + ph - (ty(v) - y0) / (y1 - y0) * ph;
  const box = document.createElement('div'); box.className = 'panel'; host.appendChild(box);
  const svg = el('svg', {viewBox: `0 0 ${W} ${H}`}, box);
  if (o.title) el('text', {x: 0, y: 14, class: 'ptitle'}, svg).textContent = o.title;
  const xt = o.xlog ? Array.from({length: x1 - x0 + 1}, (_, i) => Math.pow(10, x0 + i)) : niceTicks(x0, x1, 6);
  const yt = o.ylog ? Array.from({length: y1 - y0 + 1}, (_, i) => Math.pow(10, y0 + i)) : niceTicks(y0, y1, 5);
  xt.forEach(t => { el('line', {x1: sx(t), x2: sx(t), y1: mt, y2: mt + ph, class: 'grid'}, svg);
    el('text', {x: sx(t), y: mt + ph + 14, class: 'tick', 'text-anchor': 'middle'}, svg).textContent = o.xfmt(t); });
  yt.forEach(t => { el('line', {x1: ml, x2: ml + pw, y1: sy(t), y2: sy(t), class: 'grid'}, svg);
    el('text', {x: ml - 6, y: sy(t) + 3, class: 'tick', 'text-anchor': 'end'}, svg).textContent = o.yfmt(t); });
  el('rect', {x: ml, y: mt, width: pw, height: ph, fill: 'none', class: 'axis'}, svg);
  el('text', {x: ml + pw / 2, y: H - 6, class: 'tick', 'text-anchor': 'middle'}, svg).textContent = o.xlabel;
  el('text', {x: 12, y: mt + ph / 2, class: 'tick', 'text-anchor': 'middle', transform: `rotate(-90 12 ${mt + ph / 2})`}, svg).textContent = o.ylabel;
  const hidden = new Set();
  const lines = series.map(s => {
    const pts = s.pts.filter(p => (!o.xlog || p[0] > 0) && (!o.ylog || p[1] > 0));
    const d = pts.map((p, i) => `${i ? 'L' : 'M'}${sx(p[0]).toFixed(1)},${sy(p[1]).toFixed(1)}`).join('');
    const path = el('path', {d, fill: 'none', stroke: s.color, 'stroke-width': 1.7, 'stroke-dasharray': s.dash || ''}, svg);
    return {s, pts, path, dot: el('circle', {r: 3.2, fill: s.color, stroke: '#fff', 'stroke-width': 1, visibility: 'hidden'}, svg)};
  });
  const cross = el('line', {y1: mt, y2: mt + ph, class: 'cross', visibility: 'hidden'}, svg);
  const hit = el('rect', {x: ml, y: mt, width: pw, height: ph, fill: 'transparent'}, svg);
  hit.addEventListener('mousemove', ev => {
    const p = svgPoint(svg, ev), xv = x0 + (p.x - ml) / pw * (x1 - x0);
    const rows = [];
    let xs = null;
    lines.forEach(l => {
      if (hidden.has(l.s.label) || !l.pts.length) { l.dot.setAttribute('visibility', 'hidden'); return; }
      let best = l.pts[0];
      l.pts.forEach(q => { if (Math.abs(tx(q[0]) - xv) < Math.abs(tx(best[0]) - xv)) best = q; });
      xs = xs === null ? best[0] : xs;
      l.dot.setAttribute('cx', sx(best[0])); l.dot.setAttribute('cy', sy(best[1])); l.dot.setAttribute('visibility', 'visible');
      rows.push([l.s, best]);
    });
    if (xs === null) return;
    cross.setAttribute('x1', sx(xs)); cross.setAttribute('x2', sx(xs)); cross.setAttribute('visibility', 'visible');
    rows.sort((a, b) => b[1][1] - a[1][1]);
    showTip(`<b>${o.tipHead(xs)}</b><table>` + rows.map(([s, q]) =>
      `<tr><td>${sw(s.color)}${escH(s.label)}</td><td class="num"><b>${o.tipVal(q[1])}</b></td></tr>`).join('') + '</table>', ev);
  });
  hit.addEventListener('mouseleave', () => { hideTip(); cross.setAttribute('visibility', 'hidden'); lines.forEach(l => l.dot.setAttribute('visibility', 'hidden')); });
  const leg = document.createElement('div'); leg.className = 'legend'; box.appendChild(leg);
  series.forEach(s => {
    const it = document.createElement('span');
    it.innerHTML = `<svg viewBox="0 0 26 10"><line x1="1" y1="5" x2="25" y2="5" stroke="${s.color}" stroke-width="2.2" stroke-dasharray="${s.dash || ''}"/></svg>${escH(s.label)}`;
    it.addEventListener('click', () => {
      const off = !hidden.has(s.label); off ? hidden.add(s.label) : hidden.delete(s.label);
      it.classList.toggle('off', off);
      lines.filter(l => l.s.label === s.label).forEach(l => l.path.style.display = off ? 'none' : '');
    });
    leg.appendChild(it);
  });
}
const sci = v => { const e = Math.round(Math.log10(v)); return e === 0 ? '1' : `1e${e}`; };
function cifiCharts() {
  const C = D.cifi || []; if (!C.length) return;
  const psHost = document.getElementById('cifi-ps'), cardHost = document.getElementById('cifi-card');
  const drawPs = k => lineChart(psHost, C.map(c => ({label: c.label, color: c.color, dash: dash(c.hap), pts: c.ps[k] || []})),
    {xlog: true, ylog: true, xfmt: v => flen(v).replace('.00', '').replace('.0 ', ' '), yfmt: sci,
     xlabel: 'genomic distance s', ylabel: 'contact probability per bp',
     title: {all: 'P(s), all pairs', adjacent: 'P(s), adjacent segments', non_adjacent: 'P(s), non-adjacent segments'}[k],
     tipHead: x => `s ≈ ${flen(x)}`, tipVal: y => y.toExponential(2)});
  drawPs('all');
  document.querySelectorAll('.seg button').forEach(b => b.addEventListener('click', () => {
    document.querySelectorAll('.seg button').forEach(x => x.classList.toggle('on', x === b)); drawPs(b.dataset.k);
  }));
  lineChart(cardHost, C.map(c => {
    const tot = c.cardinality.reduce((a, p) => a + p[1], 0) || 1;
    return {label: c.label, color: c.color, dash: dash(c.hap), pts: c.cardinality.filter(p => p[1] > 0).map(p => [p[0], 100 * p[1] / tot])};
  }), {xlog: false, ylog: true, xfmt: v => v, yfmt: v => v >= 1 ? v + '%' : v + '%', xlabel: 'usable segments per read (k)',
       ylabel: '% of reads', title: 'Segments per CiFi read', tipHead: x => `k = ${x}`, tipVal: y => y.toFixed(2) + '%'});
}

bars();
nxCurves();
scafScatter();
document.addEventListener('flt', () => { hideTip(); bars(); nxCurves(); scafScatter(); });
cifiCharts();
document.querySelectorAll('.dp').forEach(h => dotplot(h, D.dp[h.dataset.key]));
document.querySelectorAll('.hic').forEach(h => hicMap(h, D.hic[h.dataset.key]));
document.querySelectorAll('.scaf').forEach(h => scafLayout(h, D.scaf[h.dataset.key]));
document.querySelectorAll('.scpair').forEach(h => pairLayout(h, D.pair[h.dataset.key]));
})();
"""

nav = "".join(f'<a href="#{i}">{t}</a>' for i, t, show in [("reads", "Reads", HAS_READS),
    ("assemblies", "Assemblies", True), ("plots", "Plots", True), ("spectra", "k-mer spectra", bool(spectra)),
    ("dotplots", "Dot plots", bool(dots)), ("scaffolds", "Scaffolds", bool(scafs)), ("cifi", "CiFi", HAS_CIFI), ("ont", "ONT coverage", bool(ONT)),
    ] if show)
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
data_json = json.dumps(DATA, separators=(",", ":")).replace("</", "<\\/")
page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light">
<title>{esc(title)}</title>
<style>{CSS}</style></head><body>
<header><h1>{esc(title)}</h1>
<p class="meta">{len(rows)} assemblies · {len(RUNS)} run(s) · generated {now}</p>
<nav>{nav}</nav><div id="flt"></div></header>
<main>{CMD_HTML}{"".join(parts)}</main>
<footer>compare_assemblies.sh (nxf-cifi) · tables: comparison.tsv, hap_pairs.tsv · static plots: comparison.png, nx_curves.png</footer>
<div id="tip"></div>
<script type="application/json" id="report-data">{data_json}</script>
<script>{JS}</script>
</body></html>
"""
report = os.path.join(outdir, "report.html")
with open(report, "w") as fh: fh.write(page)
print(f"Report: {report}")
PY
