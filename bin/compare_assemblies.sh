#!/usr/bin/env bash
#
# Side-by-side comparison of nxf-cifi runs, e.g. HiFi+ONT (no CiFi) vs HiFi+ONT+CiFi
# (--dual-scaf) vs YaHS scaffolds. Needs only bash + python3: it collects the results the
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
# HTML report (report.html, self-contained): tables + plots above, and per run the dot plots
#   (dotplots/), merqury spectra, CiFi contact maps (scaffolds/) and links to the reads QC reports
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
        spectra=sorted(f for f in glob.glob(os.path.join(ev, "merqury", "*.png"))
                       if re.search(r"^(haps|primary)\.spectra-(cn|asm)\.fl\.png$", os.path.basename(f))),
        contact=sorted(glob.glob(os.path.join(path, "scaffolds", "**", "*.png"), recursive=True)),
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

PALETTE = ["#2f6fb3", "#e07b24", "#2e9a52", "#c8413b", "#7d5bb5", "#8a5a44", "#d0619f", "#6f747b", "#a7a426", "#1fa3b5"]
grp = lambda l: re.sub(r"_(hap[12]|primary)$", "", l)
groups = list(dict.fromkeys(grp(r[0]) for r in rows))
gcol = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(groups)}

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

DATA = {"rows": [], "nx": {}, "dp": {}}
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
        ("Scaffolds", lambda r: intf(r[col["scaffolds"]]), lambda r: num(r[col["scaffolds"]]), None),
        ("Scaffold N50 (Mb)", lambda r: mb(r[col["scaffold_N50"]]), lambda r: num(r[col["scaffold_N50"]]), max),
        ("Longest (Mb)", lambda r: mb(r[col["longest"]]), lambda r: num(r[col["longest"]]), max),
        ("≥ 1 Mb (n)", lambda r: intf(r[col["scaf_ge_1Mb"]]), lambda r: num(r[col["scaf_ge_1Mb"]]), None),
        ("% bp in ≥ 1 Mb", lambda r: r[col["pct_bp_in_ge_1Mb"]], lambda r: num(r[col["pct_bp_in_ge_1Mb"]]), max),
        ("Contigs", lambda r: intf(r[col["contigs"]]), lambda r: num(r[col["contigs"]]), None),
        ("Contig N50 (Mb)", lambda r: mb(r[col["contig_N50"]]), lambda r: num(r[col["contig_N50"]]), max),
        ("Gaps", lambda r: intf(r[col["gaps"]]), lambda r: num(r[col["gaps"]]), None),
        ("T2T", lambda r: r[col["T2T"]], lambda r: num(r[col["T2T"]]), max),
        ("One telomere", lambda r: r[col["one_telo"]], lambda r: num(r[col["one_telo"]]), None),
        ("BUSCO S", lambda r: r[col["BUSCO_S"]], lambda r: num(r[col["BUSCO_S"]]), max),
        ("BUSCO D", lambda r: r[col["BUSCO_D"]], lambda r: num(r[col["BUSCO_D"]]), None),
        ("BUSCO F", lambda r: r[col["BUSCO_F"]], lambda r: num(r[col["BUSCO_F"]]), None),
        ("BUSCO M", lambda r: r[col["BUSCO_M"]], lambda r: num(r[col["BUSCO_M"]]), min),
        ("QV", lambda r: r[col["QV"]], lambda r: num(r[col["QV"]]), max),
        ("k-mer compl. (%)", lambda r: r[col["kmer_compl_pct"]], lambda r: num(r[col["kmer_compl_pct"]]), max)]

def table(head, body, numeric_from=1):
    """body: rows of (text, css class, sort value or None); click a header to sort"""
    h = "".join(f'<th class="{"num" if i >= numeric_from else ""}">{esc(c)}</th>' for i, c in enumerate(head))
    def td(i, cell):
        text, cls, v = (list(cell) + ["", None])[:3]
        cls = ("num " if i >= numeric_from else "") + cls
        dv = "" if i < numeric_from else f' data-v="{"" if v is None else v}"'
        return f'<td class="{cls.strip()}"{dv}>{esc(str(text))}</td>'
    b = "".join("<tr>" + "".join(td(i, c) for i, c in enumerate(r)) + "</tr>" for r in body)
    return f'<div class="scroll"><table class="sortable"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'

best = {}
for name, _, key, better in COLS:
    if key and better and len(rows) > 1:
        vals = [key(r) for r in rows if key(r) is not None]
        if len(set(vals)) > 1: best[name] = better(vals)
body = [[(f(r), "best" if name in best and key(r) == best[name] else "", key(r) if key else None)
         for name, f, key, _ in COLS] for r in rows]

parts = [f'<section id="assemblies"><h2>Assemblies</h2>'
         f'<p class="note">Click a column header to sort. Best value per column highlighted (when the assemblies differ). '
         f'Scaffolds inherit BUSCO and QV from their contigs. All columns: <code>comparison.tsv</code>.</p>'
         f'{table([c[0] for c in COLS], body)}']
if pairs:
    parts.append('<h3>Haplotype pairs (hap1 + hap2)</h3>' +
                 table(["Run", "Combined QV", "Combined k-mer completeness (%)"],
                       [[(p[0], ""), (p[1], "", num(p[1])), (p[2], "", num(p[2]))] for p in pairs]))
parts.append("</section>")
parts.append('<section id="plots"><h2>Comparison plots</h2>'
             '<p class="note">Hover a bar for details; the same assembly is highlighted in every panel.</p>'
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
cgroups = list(dict.fromkeys(re.sub(r"_hap[12]$", "", d["label"]) for d in DATA["cifi"]))
for d in DATA["cifi"]:
    d["color"] = PALETTE[cgroups.index(re.sub(r"_hap[12]$", "", d["label"])) % len(PALETTE)]
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
if DATA["cifi"] or MARKDUP:
    body = [[(d["label"], "")] + [(cfmt(d["summary"].get(k), kind), "", d["summary"].get(k)) for _, k, kind in CIFI_COLS[1:]]
            for d in DATA["cifi"]]
    parts.append(
        '<section id="cifi"><h2>CiFi contacts</h2>'
        '<p class="note">CiFi segments mapped to each haplotype&#39;s contigs and lifted to the YaHS scaffolds; every pair of '
        'usable segments (MAPQ ≥ the threshold) in a read is a contact, as in <code>cifi contacts</code> and McGinty et al. '
        '(Nat Commun 2026). <b>Adjacent</b> segments are neighbours along the CiFi read (direct ligation junctions); '
        '<b>non-adjacent</b> pairs are the extra, multi-way contacts that long reads add. cis / trans count only scaffolds '
        'that are chromosome-scale; trans between two large scaffolds can also mean two pieces of one chromosome that were '
        'not joined. A P(s) slope near −1 is typical of chromatin; a flat curve suggests random ligation or a misassembly.</p>'
        + ((table([c[0] for c in CIFI_COLS], body)) if body else "")
        + (('<h3>PCR duplicates (pbmarkdup)</h3><p class="note">CiFi reads sequenced from the same molecule, before '
            'digestion (the paper reports 1.8% for its GM12878 DpnII library). Rerun with <code>--cifi_rmdup</code> to '
            'remove them.</p>' + table(["Run", "Reads", "Unique molecules", "Duplicate reads", "% duplicates"],
              [[(m[0], ""), (f"{int(m[2]):,}", "", int(m[2])), (f"{int(m[3]):,}", "", int(m[3])),
                (f"{int(m[4]):,}", "", int(m[4])), (m[5], "", num(m[5]))] for m in MARKDUP])) if MARKDUP else "")
        + ('' if not DATA["cifi"] else '<h3>Contact distance decay P(s)</h3><p class="note">Contact probability per bp of distance (log–log); hover for '
          'values, switch between all pairs and adjacent / non-adjacent pairs.</p>'
          '<div class="seg"><button data-k="all" class="on">all pairs</button><button data-k="adjacent">adjacent</button>'
          '<button data-k="non_adjacent">non-adjacent</button></div><div id="cifi-ps" class="nx"></div>'
          '<h3>Usable segments per CiFi read</h3><p class="note">Share of CiFi reads with k segments mapped at the '
          'MAPQ threshold (k ≥ 2 give contacts, k ≥ 3 are multi-way).</p><div id="cifi-card" class="nx"></div>') + '</section>')

# ---------------- ONT reads mapped to hap1 / hap2 (--map_ont): samtools coverage + flagstat
ONT, ONT_SEQ = [], []   # per assembly; per sequence >= 1 Mb
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
        depth = sum(x[1] * x[4] for x in seqs) / tot
        big = [x for x in seqs if x[1] >= 1_000_000]
        ONT.append([lab, num(fs.get("primary")), num(fs.get("primary mapped %")), depth,
                    100 * sum(x[3] for x in seqs) / tot, len(seqs),
                    sum(1 for x in big if x[4] < 0.5 * depth or x[4] > 2 * depth), base])
        ONT_SEQ += [[lab, x[0], x[1], x[4], 100 * x[3] / x[1], x[5], x[4] / depth if depth else None] for x in big]
if ONT:
    f2 = lambda v, d=2: "NA" if v is None else f"{v:,.{d}f}"
    parts.append(
        '<section id="ont"><h2>ONT read coverage</h2>'
        '<p class="note">ONT reads mapped with minimap2 to each haplotype (scaffolds with CiFi, contigs otherwise); '
        'primary alignments, from <code>samtools flagstat</code> and <code>samtools coverage</code>. Depth is the '
        'length-weighted mean over all sequences. A sequence ≥ 1 Mb at less than half or more than twice that depth '
        'is flagged: possible collapsed repeat, duplication, or contamination. BAMs (with CSI index) are in '
        '<code>ont_mapping/</code> for IGV.</p>'
        + table(["Haplotype", "Reads", "% mapped", "Mean depth (×)", "% bases covered", "Sequences", "Flagged ≥ 1 Mb"],
                [[(o[0], ""), (f2(o[1], 0), "", o[1]), (f2(o[2]), "", o[2]), (f2(o[3], 1), "", o[3]),
                  (f2(o[4]), "", o[4]), (f"{o[5]:,}", "", o[5]), (str(o[6]), "", o[6])] for o in ONT])
        + ('' if not ONT_SEQ else '<h3>Sequences ≥ 1 Mb</h3><p class="note">Click a header to sort; flagged rows '
           'are outside 0.5–2× the haplotype&#39;s mean depth.</p>' +
           table(["Haplotype", "Sequence", "Length (Mb)", "Mean depth (×)", "% covered", "Mean MAPQ", "Depth / mean"],
                 [[(r[0], ""), (r[1], ""), (f2(r[2] / 1e6), "", r[2]), (f2(r[3], 1), "", r[3]), (f2(r[4]), "", r[4]),
                   (f2(r[5], 1), "", r[5]),
                   (f2(r[6]), "worse" if r[6] is not None and (r[6] < 0.5 or r[6] > 2) else "", r[6])] for r in ONT_SEQ],
                 numeric_from=2))
        + '</section>')

for run in RUNS:
    lab = run["label"]
    sec = [f'<section id="run-{esc(lab)}"><h2>Run {esc(lab)}</h2><p class="note"><code>{esc(run["path"])}</code></p>']
    if run["reads"]:
        sec.append("<p>Reads QC: " + " · ".join(
            f'<a href="{esc(os.path.relpath(f, outdir))}">{esc(os.path.basename(f))}</a>' for f in run["reads"]) + "</p>")
    if run["dotplots"]:
        rs = [r for _, r, _ in run["dotplots"] if len(r) >= 10]
        sec.append('<h3>Dot plots</h3><p class="note">Hover an alignment for coordinates; drag to zoom, double-click to '
                   'reset. Blue same orientation, red inverted. Query sequences are ordered along the target and '
                   'flipped (−) where most of them align reversed.</p>')
        if rs:
            sec.append(table(["Query", "Target", "Query aligned (%)", "Target aligned (%)", "Inverted (Mb)", "Alignments"],
                             [[(r[0], ""), (r[1], ""), (r[5], "", num(r[5])), (r[7], "", num(r[7])),
                               (mb(r[8]), "", num(r[8])), (intf(r[9]), "", num(r[9]))] for r in rs], 2))
        figs = []
        for png, r, js in run["dotplots"]:
            cap = esc(" vs ".join(r[:2])) if r else esc(os.path.basename(png))
            d = None
            if js:
                try: d = json.load(open(js))
                except ValueError: d = None
            if d and d.get("tseq") is not None:
                key = f"dp{len(DATA['dp'])}"
                DATA["dp"][key] = d
                figs.append(f'<figure class="dpfig"><div class="dp" data-key="{key}"></div><figcaption>{cap}</figcaption></figure>')
            elif os.path.getsize(png):
                figs.append(fig(png, cap))
        sec.append('<div class="grid dpgrid">' + "".join(figs) + "</div>")
    if run["spectra"]:
        sec.append('<h3>Merqury k-mer spectra</h3><div class="grid">' +
                   "".join(fig(p, esc(os.path.basename(p))) for p in run["spectra"] if os.path.getsize(p)) + "</div>")
    if run["contact"]:
        sec.append('<h3>CiFi contact maps (scaffolds)</h3><div class="grid">' +
                   "".join(fig(p, esc(os.path.basename(p))) for p in run["contact"] if os.path.getsize(p)) + "</div>")
    if run["cifi_data"] or run["cifi_html"]:
        sec.append("<h3>CiFi contact analysis</h3>")
        if run["cifi_html"]:
            sec.append('<p>cifi reports: ' + " · ".join(
                f'<a href="{esc(os.path.relpath(f, outdir))}">{esc(os.path.relpath(f, os.path.join(run["path"], "cifi")))}</a>'
                for f in run["cifi_html"]) + "</p>")
        for d in run["cifi_data"]:
            rows_ = [r for r in d.get("scaffolds", [])]
            cis = sorted(r[6] for r in rows_ if r[6] is not None)
            med = cis[len(cis) // 2] if cis else None
            body = [[(r[0], ""), (f"{r[1] / 1e6:,.2f}", "", r[1]), (f"{r[2]:,}", "", r[2]), (f"{r[3]:,.1f}", "", r[3]),
                     (f"{r[4]:,}", "", r[4]), (f"{r[5]:,}", "", r[5]),
                     ("NA" if r[6] is None else f"{r[6]:.2f}", "low" if med and r[6] is not None and r[6] < 0.8 * med else "", r[6])]
                    for r in rows_]
            sec.append(f'<h4>{esc(d["label"])}</h4>')
            if os.path.exists(d["png"]) and os.path.getsize(d["png"]):
                sec.append(fig(d["png"], f'{esc(d["label"])}: segments per read, P(s), MAPQ and contact classes, coverage'))
            if body:
                sec.append(f'<p class="note">Chromosome-scale scaffolds (≥ {d.get("min_scaf", 0) / 1e6:g} Mb). '
                           '<span class="lowkey">Highlighted</span>: % cis below 80% of the median, worth a look in '
                           'Juicebox (possible misjoin or a scaffold that belongs to another chromosome).</p>'
                           + table(["Scaffold", "Length (Mb)", "Segments", "Segments/Mb", "cis contacts",
                                    "trans contacts", "% cis"], body))
    if len(sec) == 1 + bool(run["reads"]):
        sec.append('<p class="note">No dot plots, spectra, contact maps or CiFi analysis found for this run.</p>')
    parts.append("".join(sec) + "</section>")

CSS = r"""
:root { --bg:#f7f7f5; --fg:#1d1f21; --muted:#676c73; --line:#e1e2de; --card:#ffffff; --accent:#1f5fa8; --best:#e3edf9; --worse:#fbe5e1; --hover:#f1f4f8; }
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
.scroll { overflow-x:auto; width:fit-content; max-width:100%; border:1px solid var(--line); border-radius:6px; background:var(--card) }
table { border-collapse:collapse; font-size:13px }
th, td { padding:6px 10px; border-bottom:1px solid var(--line); white-space:nowrap; text-align:left }
th { font-weight:600; background:#fafafa; cursor:pointer; user-select:none; position:relative; padding-right:20px }
th:hover { background:var(--hover) }
th::after { content:"↕"; position:absolute; right:6px; color:#b5b9bf; font-size:11px }
th[data-dir="asc"]::after { content:"▲"; color:var(--accent) } th[data-dir="desc"]::after { content:"▼"; color:var(--accent) }
th.num, td.num { text-align:right; font-variant-numeric:tabular-nums }
td.best { font-weight:700; background:var(--best) }
td.worse { font-weight:700; background:var(--worse) }
tbody tr:hover td { background:var(--hover) } tbody tr:hover td.best { background:var(--best) }
tbody tr:last-child td { border-bottom:none }
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
svg .grid { stroke:#eceded } svg .axis { stroke:#b9bdc2 }
svg .hit { fill:transparent } svg g.row { cursor:default }
.bars.hl g.row:not(.on) { opacity:.3 } .bars g.row.on .hit { fill:var(--hover) }
svg .cross { stroke:#8d9298; stroke-dasharray:3 3 }
svg .n50 { stroke:#c9ccd0; stroke-dasharray:4 4 }
.legend { display:flex; flex-wrap:wrap; gap:4px 14px; margin-top:10px; font-size:12px; grid-column:1/-1 }
.legend span { cursor:pointer; display:inline-flex; align-items:center; gap:6px; user-select:none }
.legend span.off { opacity:.35; text-decoration:line-through }
.legend svg { width:26px; height:10px }
.dp { position:relative; width:100% }
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
  const ths = [...t.tHead.rows[0].cells];
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
  const n = D.rows.length, rh = 20, top = 28, lw = 124, W = 420, bw = W - lw - 44, H = top + n * rh + 20;
  const hl = i => {
    host.classList.toggle('hl', i >= 0);
    host.querySelectorAll('g.row').forEach(g => g.classList.toggle('on', +g.dataset.i === i));
  };
  METRICS.forEach(m => {
    const vals = D.rows.map(m.get), known = vals.filter(v => v != null);
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
    D.rows.forEach((r, i) => {
      const y = top + i * rh, v = vals[i];
      const g = el('g', {class: 'row', 'data-i': i}, svg);
      el('rect', {x: 0, y, width: W, height: rh, class: 'hit'}, g);
      el('text', {x: lw - 6, y: y + rh / 2 + 4, 'text-anchor': 'end', class: 'lab'}, g).textContent = r.label;
      if (v != null) {
        el('rect', {x: lw, y: y + 3, width: Math.max(sx(v) - lw, 1), height: rh - 6, rx: 2, fill: r.color,
                    'fill-opacity': r.hap === 'hap2' ? 0.62 : r.hap === 'primary' ? 0.4 : 0.95}, g);
        el('text', {x: sx(v) + 4, y: y + rh / 2 + 4, class: 'val'}, g).textContent = short(v);
      } else el('text', {x: lw + 4, y: y + rh / 2 + 4, class: 'na'}, g).textContent = 'NA';
      g.addEventListener('mousemove', ev => { hl(i); showTip(`${sw(r.color)}<b>${escH(r.label)}</b><br>${v == null ? 'not computed for this assembly' : m.tip(r)}`, ev); });
      g.addEventListener('mouseleave', () => { hl(-1); hideTip(); });
    });
  });
}

/* ---------------------------------------------------------------- Nx curves */
function nxCurves() {
  const host = document.getElementById('nx');
  const rows = D.rows.filter(r => D.nx[r.label]);
  if (!host || !rows.length) return;
  const all = [];
  rows.forEach(r => ['scaf', 'ctg'].forEach(k => D.nx[r.label][k].forEach(p => all.push(p[1]))));
  if (!all.length) return;
  const l0 = Math.floor(Math.log10(Math.max(1, Math.min(...all)))), l1 = Math.max(l0 + 1, Math.ceil(Math.log10(Math.max(...all))));
  const W = 560, H = 340, ml = 58, mr = 14, mt = 28, mb = 40, pw = W - ml - mr, ph = H - mt - mb;
  const sx = x => ml + x / 100 * pw, sy = y => mt + ph - (Math.log10(Math.max(y, 1)) - l0) / (l1 - l0) * ph;
  const hidden = new Set(), paths = [];
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
      const path = el('path', {d, fill: 'none', stroke: r.color, 'stroke-width': 1.7, 'stroke-dasharray': dash(r.hap)}, svg);
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
    c.strokeStyle = '#e6e8eb'; c.lineWidth = 1; c.beginPath();
    T.forEach(t => { const x = Math.round(X(t[1])) + .5; c.moveTo(x, MT); c.lineTo(x, MT + ph()); });
    Q.forEach(q => { const y = Math.round(Y(q[1])) + .5; c.moveTo(ML, y); c.lineTo(ML + pw(), y); });
    c.stroke();
    c.lineWidth = 1.5; c.lineCap = 'round';
    for (const fwd of [true, false]) {
      c.strokeStyle = fwd ? '#1f5fa8' : '#c8413b'; c.beginPath();
      for (const s of S) if (s.fwd === fwd) { c.moveTo(X(s.x0), Y(s.y0)); c.lineTo(X(s.x1), Y(s.y1)); }
      c.stroke();
    }
    if (hover) {
      c.strokeStyle = 'rgba(240,160,30,.9)'; c.lineWidth = 5; c.beginPath();
      c.moveTo(X(hover.x0), Y(hover.y0)); c.lineTo(X(hover.x1), Y(hover.y1)); c.stroke();
    }
    c.restore();
    c.strokeStyle = '#b9bdc2'; c.lineWidth = 1; c.strokeRect(ML + .5, MT + .5, pw(), ph());
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
      if (b - a < 11) return;
      c.save(); c.translate((a + b) / 2 + 3, MT - 5); c.rotate(-Math.PI / 2); c.textAlign = 'left'; c.fillText(trunc(t[0], 12), 0, 0); c.restore();
    });
    c.textAlign = 'left';
    Q.forEach(q => {
      const a = Math.min(Y(q[1]), MT + ph()), b = Math.max(Y(q[1] + q[2]), MT);
      if (a - b < 11) return;
      c.fillText(trunc(q[0], 13) + (q[3] ? ' (−)' : ''), ML + pw() + 5, (a + b) / 2 + 3);
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
    if (!inside(px, py)) { if (hover) { hover = null; draw(); } hideTip(); return; }
    const h = nearest(px, py);
    if (h !== hover) { hover = h; draw(); }
    if (h) {
      const [ti, qi, ts, te, qs, qe, plus] = h.s, t = T[ti], q = Q[qi];
      showTip(`<b>${escH(q[0])}</b> ${fint(qs)}–${fint(qe)}<br><b>${escH(t[0])}</b> ${fint(ts)}–${fint(te)}<br>` +
              `${flen(qe - qs)} · strand ${plus ? '+' : '−'} · <span style="color:${h.fwd ? '#1f5fa8' : '#c8413b'}">${h.fwd ? 'collinear' : 'inverted'}</span>`, ev);
    } else {
      const x = iX(px), y = iY(py), t = find(T, x), q = find(Q, y);
      const qp = q[3] ? q[2] - (y - q[1]) : y - q[1];
      showTip(`<span style="color:#888">${escH(d.t)}</span> ${escH(t[0])}: ${flen(x - t[1])}<br><span style="color:#888">${escH(d.q)}</span> ${escH(q[0])}${q[3] ? ' (−)' : ''}: ${flen(qp)}`, ev);
    }
  });
  cv.addEventListener('mouseleave', () => { hideTip(); if (hover) { hover = null; draw(); } });
  cv.addEventListener('mousedown', ev => { const p = pos(ev); if (inside(...p)) { drag = p; ev.preventDefault(); } });
  window.addEventListener('mouseup', ev => {
    if (!drag) return;
    const [px, py] = pos(ev), [qx, qy] = drag; drag = null; sel.style.display = 'none';
    if (Math.abs(px - qx) < 6 || Math.abs(py - qy) < 6) return;
    const cx = v => Math.min(ML + pw(), Math.max(ML, v)), cy = v => Math.min(MT + ph(), Math.max(MT, v));
    view = {x0: iX(cx(Math.min(px, qx))), x1: iX(cx(Math.max(px, qx))), y0: iY(cy(Math.max(py, qy))), y1: iY(cy(Math.min(py, qy)))};
    btn.style.display = 'block'; hover = null; draw();
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
cifiCharts();
document.querySelectorAll('.dp').forEach(h => dotplot(h, D.dp[h.dataset.key]));
})();
"""

nav = '<a href="#assemblies">Assemblies</a><a href="#plots">Plots</a>' + ('<a href="#cifi">CiFi contacts</a>' if DATA["cifi"] or MARKDUP else "") + ('<a href="#ont">ONT coverage</a>' if ONT else "") + "".join(
    f'<a href="#run-{esc(r["label"])}">{esc(r["label"])}</a>' for r in RUNS)
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
data_json = json.dumps(DATA, separators=(",", ":")).replace("</", "<\\/")
page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light">
<title>{esc(title)}</title>
<style>{CSS}</style></head><body>
<header><h1>{esc(title)}</h1>
<p class="meta">{len(rows)} assemblies · {len(RUNS)} run(s) · generated {now}</p>
<nav>{nav}</nav></header>
<main>{"".join(parts)}</main>
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
