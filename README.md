# nxf-cifi

Genome assembly with hifiasm from HiFi and ONT ultra-long reads. With CiFi reads, it also phases the haplotypes and scaffolds them to chromosome scale. Every process runs in a container (Docker by default).

```
without --cifi : reads QC -> hifiasm --ul -> contigs (primary, hap1, hap2) -> evaluation
with --cifi    : reads QC -> pbmarkdup (PCR duplicates) -> cifi qc -> cifi digest (segments + Hi-C-like R1/R2)
                 -> hifiasm --ul --dual-scaf --h1 R1 --h2 R2 -> phased contigs (hap1, hap2)
                 -> per haplotype: minimap2 map-hifi -> cifi contacts -> yahs scaffolds
                    -> Juicebox files (juicer pre -a [+ juicer_tools .hic]) -> PretextMap + PretextSnapshot
                    -> CiFi contact analysis (MAPQ, segments/read, cis/trans, distance decay P(s), coverage)
                 -> evaluation of contigs and scaffolds
evaluation     : N50 / contigs / gaps, telomeres (TTTAGGG)n, compleasm (BUSCO genes), merqury (QV, completeness)
dot plots      : minimap2 hap2 vs hap1 (contigs and scaffolds); with --ref every assembly vs the reference
reads QC       : faster-report (github.com/angelovangel/faster-report), one report for all input reads (HiFi, ONT, CiFi)
report         : report/report.html, a single self-contained page with all of the above
```

## Run

```bash
# HiFi + ONT UL
nextflow run nxf-cifi --hifi hifi.bam --ont ont_ul.fastq.gz --prefix myplant --outdir asm_noCiFi

# HiFi + ONT UL + CiFi (phasing + scaffolding)
nextflow run nxf-cifi --hifi hifi.bam --ont ont_ul.fastq.gz --cifi cifi.bam --enzyme HindIII \
    --prefix myplant --outdir asm_CiFi --meryl_db asm_noCiFi/eval/merqury/hifi.k21.meryl
```

With `--cifi`, each scaffolded haplotype also gets a CiFi contact analysis, modelled on McGinty et al. (Dennis lab,
Nat Commun 2026; github.com/mydennislab/CiFi):
- % of segments at MAPQ ≥ 1/10/20/30/60 (the paper's mappability comparison)
- usable segments per CiFi read, % multi-way reads (≥ 3 segments), contacts per read
- cis vs trans contacts between chromosome-scale scaffolds, cis by distance (< 1 kb … ≥ 1 Mb), reads touching > 1 chromosome
- contact distance decay P(s) (log-log, slope 50 kb–5 Mb) for all pairs, adjacent segments (direct ligations)
  and non-adjacent segments (the multi-way contacts long reads add)
- segment coverage along the scaffolds and a per-scaffold cis/trans table (a scaffold with low % cis is worth a look in Juicebox)

Contacts are every pair of segments of a read (k segments → k(k-1)/2), as in `cifi contacts` and the paper.
`--cifi_stats_mapq` (1), `--cifi_stats_min_scaf` (1 Mb), `--cifi_stats_bin` (100 kb) and `--cifi_stats_max_pairs`
(20 M; segment statistics always use all reads) tune it; `--skip_cifi_stats` turns it off.
The CiFi reads also go through `pbmarkdup` (github.com/PacificBiosciences/pbmarkdup) for the PCR-duplicate rate, as in
the paper (1.8% for its GM12878 DpnII library). By default the rate is only reported; `--cifi_rmdup` removes the duplicates
before digestion, and `--skip_cifi_markdup` skips the step. pbmarkdup expects PacBio HiFi BAM/FASTQ; for CiFi given as FASTQ
the reads are merged into one BAM first.
Not included from the paper: TAD calling (TopDom) and Hi-C correlation, which need a reference annotation or Illumina Hi-C.

Add `--ref ref.fa` to dot plot every assembly against a reference. `--dotplot_min_aln` (10 kb) and
`--dotplot_min_seq` (500 kb) drop short alignments and sequences from the plots; `--skip_dotplot` turns them off.

Add `--map_ont` to map the ONT reads to hap1 and hap2 (the scaffolds with `--cifi`, the contigs otherwise) with minimap2
`-x lr:hq` (Q20+ ONT reads). `ont_mapping/` gets a sorted BAM with a
CSI index (chromosomes can exceed the BAI limit) for IGV, `samtools coverage` per sequence and `samtools flagstat`.
The report shows the mapping rate, mean depth and breadth per haplotype, and flags sequences ≥ 1 Mb at less than half or
more than twice the haplotype's mean depth (collapsed repeats, duplications, contamination). With ultra-long data this
step takes hours.

### CiFi on an existing assembly

`--assembly hap1.fa[,hap2.fa]` skips hifiasm and runs only the CiFi part on contigs you already have:
duplicate rate, digestion, segment mapping, YaHS scaffolding, Juicebox/Pretext files and the CiFi contact analysis,
followed by the usual evaluation (stats, telomeres, BUSCO, dot plots) and the report.

```bash
nextflow run nxf-cifi --assembly asm_noCiFi/fasta/myplant.bp.hap1.p_ctg.fa,asm_noCiFi/fasta/myplant.bp.hap2.p_ctg.fa \
    --cifi cifi.bam --enzyme HindIII --prefix myplant --outdir asm_CiFiScaf \
    --meryl_db asm_noCiFi/eval/merqury/hifi.k21.meryl
```

- One or two FASTAs (.gz ok), in hap1, hap2 order. They are renamed `<prefix>.asm.hap1.p_ctg` / `hap2`, so outputs,
  dot plots and `compare_assemblies.sh` work as for a hifiasm run (`-d CiFiScaf=asm_CiFiScaf`).
- `--cifi` is required. `--hifi` or `--meryl_db` turns on merqury QV; without either, merqury is skipped.
  `--ont` is optional (reads QC only).
- The CiFi reads are used for scaffolding only. Phasing with CiFi (`hifiasm --dual-scaf`) needs a full run from reads.

Inputs can be comma-separated lists or quoted globs (`--hifi 'reads/*.bam'`).
Set `--max_cpus` / `--max_memory` to the machine size; hifiasm gets the full amount and the other steps are capped at it.
Add `-resume` to continue a run, `-profile singularity,slurm` on a cluster.
Run `nextflow run nxf-cifi --help` for all options (`--help <param>` for one, `--showHidden` for the advanced ones).
Options are defined in `nextflow_schema.json` and checked at start by the nf-schema plugin (Nextflow ≥ 25.04);
`--monochrome_logs` turns colours off. The defaults match `run_hifiasm.sh`.

## Output

```
report/report.html                tables, Nx + metric plots, CiFi contact analysis, dot plots, k-mer spectra, contact maps, links
summary.tsv                       per assembly: contiguity, telomeres, BUSCO, QV, k-mer completeness
reads/                            reads.faster-report.html: one report for all reads (+ merged CiFi BAM)
cifi/qc  cifi/digest              cifi qc report, digest report, segments, R1/R2
cifi/markdup                      pbmarkdup duplicate summary (<prefix>.cifi.markdup.tsv) + log
cifi/contacts                     cifi contacts report + stats per haplotype
cifi/stats                        <hap>.cifi_stats.{tsv,json,png}, <hap>.cifi_scaffolds.tsv (contact analysis)
hifiasm/                          contig GFAs + hifiasm log
fasta/                            contigs (+ .fai) and <hap>.scaffolds.fa
scaffolds/<hap>/                  yahs output, <hap>_JBAT.* (Juicebox), .hic, .pretext + snapshot PNGs
eval/{stats,telomeres,compleasm,merqury}
dotplots/                         <query>_vs_<target>.dotplot.{png,json,tsv} (+ .paf.gz); tsv = % of each side aligned, inverted bp
pipeline_info/                    timeline, resource report, trace
```

After curation in Juicebox (load `<hap>_JBAT.hic`, import `<hap>_JBAT.assembly`, export `review.assembly`):
`juicer post -o <hap>_JBAT <hap>_JBAT.review.assembly <hap>_JBAT.liftover.agp fasta/<hap>.fa`

## Compare runs

`bin/compare_assemblies.sh` compares finished runs side by side. It runs no assembly or evaluation tools. It only collects what the pipeline already computed (`eval/stats`, `eval/compleasm`, `eval/merqury`, `fasta/*.fai`), so it needs just bash, python3 and matplotlib:

```bash
bash nxf-cifi/bin/compare_assemblies.sh -d noCiFi=asm_noCiFi -d CiFi=asm_CiFi -d CiFiL0=asm_CiFi_l0 -o compare_3runs
```

Each `-d label=outdir` adds `<label>_hap1/_hap2` and, for scaffolded runs, `<label>Scaf_hap1/_hap2` (`-p` also adds the primary assembly).
`-a label=fasta` adds any other assembly. For those, contiguity and telomeres are computed from the FASTA, and BUSCO/QV are NA.

Output: `comparison.tsv` (per assembly), `hap_pairs.tsv` (combined hap1+hap2 QV / completeness + spectra-cn plot path),
`nx_curves.png` (scaffold + contig Nx), `comparison.png` (N50s, % bp in ≥ 1 Mb, T2T, BUSCO complete, QV) and
`report.html`: one self-contained page (light theme, works offline) with sortable tables, interactive metric bars and Nx curves
(hover for details, click the legend to hide assemblies), each run's interactive dot plots (hover an alignment for its
coordinates, drag to zoom, double-click to reset), merqury spectra, contact maps and links to its reads QC reports. Use `-P` to skip the plots and `-t` to set the report title.
The pipeline runs the same script on each run to make `report/report.html`.

## Containers

| Step | Image |
| :-- | :-- |
| BAM → FASTQ (HiFi; CiFi for reads QC), CiFi BAM | `quay.io/biocontainers/samtools:1.21--h96c455f_1` |
| ONT length filter, GFA → FASTA | `quay.io/biocontainers/seqkit:2.10.1--he881be0_0` |
| reads QC | `docker.io/aangeloo/faster-report:latest` |
| CiFi PCR duplicates | `quay.io/biocontainers/pbmarkdup:1.2.0--h9ee0642_0` |
| cifi qc / digest / contacts | `quay.io/biocontainers/cifi:1.1.0--py312h0c6b66a_0` |
| hifiasm | `quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0` |
| segment mapping, dot plot alignment, ONT mapping | `community.wave.seqera.io/library/minimap2_samtools_procps-ng:a6a2be08b29bce9f` |
| yahs + juicer | `quay.io/biocontainers/yahs:1.2.2--h577a1d6_1` |
| juicer_tools (`--juicer_jar`) | `eclipse-temurin:21-jre` |
| PretextMap / PretextSnapshot | `quay.io/biocontainers/pretextmap:0.2.4--hd63eeec_0`, `pretextsnapshot:0.1.1--hc52dbad_0` |
| stats, telomeres, summary, dot plots, CiFi contact analysis, report | `community.wave.seqera.io/library/python_matplotlib_procps-ng:e896924bb279b00c` |
| compleasm | `quay.io/biocontainers/compleasm:0.2.9--pyhdfd78af_0` |
| meryl / merqury | `quay.io/biocontainers/merqury:1.3--hdfd78af_4` |
