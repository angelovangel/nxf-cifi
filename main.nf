#!/usr/bin/env nextflow
/*
 * nxf-cifi: plant genome assembly with hifiasm (HiFi + ONT ultra-long), optionally with CiFi
 * (PacBio long-read Hi-C) for phasing and chromosome-scale scaffolding, followed by evaluation.
 * Nextflow port of run_hifiasm.sh; every process runs in a container (docker by default).
 *
 * Without --cifi : reads QC -> hifiasm --ul -> contigs (primary, hap1, hap2) -> evaluation
 * With --cifi    : reads QC -> cifi qc -> cifi digest (segments + Hi-C-like R1/R2)
 *                  -> hifiasm --ul --dual-scaf --h1 R1 --h2 R2 -> phased contigs (hap1, hap2)
 *                  -> per haplotype: minimap2 map-hifi (segments) -> cifi contacts -> yahs scaffolds
 *                     -> Juicebox files (juicer pre -a [+ juicer_tools .hic]) -> PretextMap contact map
 *                  -> CiFi contact analysis per haplotype (MAPQ, segments/read, cis/trans, P(s), coverage)
 *                  -> evaluation of contigs and scaffolds
 * Evaluation     : contiguity, telomeres (TTTAGGG)n, compleasm (BUSCO genes), merqury (QV, completeness),
 *                  minimap2 dot plots (hap2 vs hap1, optionally every assembly vs --ref)
 * Report         : report/report.html (bin/compare_assemblies.sh on this run)
 * Reads QC       : faster-report (github.com/angelovangel/faster-report) for HiFi, ONT and CiFi reads
 */

nextflow.enable.dsl = 2

// --help, parameter validation and the run summary come from nextflow_schema.json (nf-schema plugin)
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

def readFiles(p) {
    def pats = p instanceof List ? p : p.toString().tokenize(',')
    return Channel.fromPath(pats.collect { it.toString().trim() }, checkIfExists: true)
}

def readType(p) {
    return p.toString().contains('.bam') ? 'bam' : 'fastq'
}

/* =========================================================================== reads */
process BAM2FASTQ {
    tag "$kind: $bam.name"
    label 'mid'
    input:  tuple val(kind), path(bam)
    output: tuple val(kind), path("${bam.baseName}.fastq.gz")
    script:
    "samtools fastq -@ ${task.cpus} $bam | bgzip -@ ${task.cpus} > ${bam.baseName}.fastq.gz"
    stub:
    """
    touch ${bam.baseName}.fastq.gz
    """
}

process ONT_FILTER {
    tag "$fq.name"
    label 'mid'
    input:  path fq
    output: path "*.min${params.ont_minlen}.fastq.gz"
    script:
    def base = fq.name.replaceAll(/\.f(ast)?[aq](\.gz)?$/, '')
    "seqkit seq -j ${task.cpus} -m ${params.ont_minlen} $fq -o ${base}.min${params.ont_minlen}.fastq.gz"
    stub:
    """
    touch ${fq.simpleName}.min${params.ont_minlen}.fastq.gz
    """
}

process FASTER_REPORT {
    tag "$kind"
    label 'mid'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/reads", mode: params.publish_mode
    input:  tuple val(kind), path('reads/*'), val(filetype)
    output: path "${kind}.faster-report.html"
    script:
    def type = kind == 'ont' ? 'ont' : 'pacbio'
    """
    # as in angelovangel/faster-report main.nf: assets + run info from the first read header
    cp /temp/custom.css ./custom.css 2>/dev/null || true
    mkdir -p www && cp -r /temp/www/*.png ./www/ 2>/dev/null || true
    IFS=, read -r d_platform d_flowcell d_rundate d_basecall d_mods d_ft < <(get-header-data.sh reads) || true
    opts=()
    [[ -n "\${d_rundate:-}"  && "\$d_rundate"  != NA ]] && opts+=(-d "\$d_rundate")
    [[ -n "\${d_flowcell:-}" && "\$d_flowcell" != NA ]] && opts+=(-f "\$d_flowcell")
    [[ -n "\${d_basecall:-}" && "\$d_basecall" != NA ]] && opts+=(-b "\$d_basecall")
    [[ -n "\${d_mods:-}"     && "\$d_mods"     != - ]]  && opts+=(-m "\$d_mods")
    /temp/faster-report.R -p reads -r '\\.(fastq|fq|fasta)(\\.gz)?\$' -t $type -y $filetype \\
        "\${opts[@]}" -x ${params.report_subsample} -o ${kind}.faster-report.html
    """
    stub:
    """
    touch ${kind}.faster-report.html
    """
}

/* =========================================================================== CiFi */
process CIFI_MARKDUP {
    label 'mid'
    errorStrategy { params.cifi_rmdup ? 'terminate' : 'ignore' }   // only a QC number unless duplicates are removed
    publishDir "${params.outdir}/cifi/markdup", mode: params.publish_mode, pattern: '*.{tsv,log}'
    input:  path bam
    output:
        path "${params.prefix}.cifi.markdup.tsv", emit: tsv
        path "${params.prefix}.cifi.pbmarkdup.log"
        path "${params.prefix}.cifi.dedup.bam", emit: bam, optional: true
    script:
    def out = params.cifi_rmdup ? "${params.prefix}.cifi.dedup.bam" : 'marked.tmp.bam'
    def log = "${params.prefix}.cifi.pbmarkdup.log"
    """
    pbmarkdup -j ${task.cpus} --log-level INFO ${params.cifi_rmdup ? '--rmdup' : ''} $bam $out > $log 2>&1 \\
        || { tail -n 20 $log >&2; exit 1; }
    rm -f marked.tmp.bam*
    # summary table: LIBRARY  READS  UNIQUE MOLECULES (%)  DUPLICATE READS (%)
    awk -v OFS='\\t' 'BEGIN { print "library", "reads", "unique_molecules", "duplicate_reads", "pct_duplicates" }
        /[0-9]+ +[0-9]+ [(][0-9.]+%[)] +[0-9]+ [(][0-9.]+%[)]/ { n = NF; p = \$n; gsub(/[()%]/, "", p); print \$(n-5), \$(n-4), \$(n-3), \$(n-1), p }' \\
        $log > ${params.prefix}.cifi.markdup.tsv
    """
    stub:
    """
    printf "library\\treads\\tunique_molecules\\tduplicate_reads\\tpct_duplicates\\nTOTAL\\t100\\t98\\t2\\t2.0\\n" > ${params.prefix}.cifi.markdup.tsv
    touch ${params.prefix}.cifi.pbmarkdup.log ${params.cifi_rmdup ? "${params.prefix}.cifi.dedup.bam" : ''}
    """
}

process CIFI_BAM {
    label 'mid'
    publishDir "${params.outdir}/reads", mode: params.publish_mode
    input:  path(reads, stageAs: 'in/*')
    output: path "${params.prefix}.cifi.bam"
    script: "cifi_bam.sh ${params.prefix}.cifi.bam ${task.cpus} ${reads}"
    stub:
    """
    touch ${params.prefix}.cifi.bam
    """
}

process CIFI_QC {
    label 'mid'
    publishDir "${params.outdir}/cifi", mode: params.publish_mode
    input:  path bam
    output: path 'qc'
    script: "cifi qc $bam -o qc -e ${params.enzyme} -n 0 --min-sites 0"
    stub:
    """
    mkdir qc; touch qc/qc.json
    """
}

process CIFI_DIGEST {
    label 'mid'
    publishDir "${params.outdir}/cifi", mode: params.publish_mode
    input:  path bam
    output:
        path "digest/${params.prefix}.segments.fastq.gz", emit: segments
        tuple path("digest/${params.prefix}_R1.fastq*"), path("digest/${params.prefix}_R2.fastq*"), emit: hic
        path 'digest/*'
    script:
    def p = "digest/${params.prefix}"
    """
    mkdir -p digest
    cifi digest $bam -e ${params.enzyme} -o $p -m ${params.cifi_min_seg} -l ${params.cifi_min_seg_len} \\
        --segments-out ${p}.segments.fastq.gz --no-strip-overhang --gzip
    """
    stub:
    """
    mkdir digest; touch digest/${params.prefix}.segments.fastq.gz digest/${params.prefix}_R1.fastq.gz digest/${params.prefix}_R2.fastq.gz
    """
}

/* =========================================================================== assembly */
process HIFIASM {
    label 'hifiasm'
    publishDir "${params.outdir}/hifiasm", mode: params.publish_mode
    input:
        path(hifi, stageAs: 'hifi/*')
        path(ont, stageAs: 'ont/*')
        path(hic, stageAs: 'hic/*')       // [] or [R1, R2]
    output:
        path "${params.prefix}.*ctg.gfa", emit: gfa
        path "${params.prefix}.hifiasm.log"
    script:
    def hicFiles = [hic].flatten()
    def args = ["-o ${params.prefix}", "-t ${task.cpus}", "--n-hap ${params.ploidy}", "--telo-m ${params.telo_motif_r}"]
    if (params.purge != null) args << "-l ${params.purge}"
    args << "--ul ${[ont].flatten().join(',')}"
    if (hicFiles) args << "--dual-scaf --h1 ${hicFiles[0]} --h2 ${hicFiles[1]}"
    if (params.hifiasm_args) args << params.hifiasm_args
    """
    hifiasm --version > hifiasm.version
    hifiasm ${args.join(' ')} ${[hifi].flatten().join(' ')} 2> ${params.prefix}.hifiasm.log \\
        || { tail -n 30 ${params.prefix}.hifiasm.log >&2; exit 1; }
    """
    stub:
    """
    for h in "" hap1. hap2.; do printf "S\\tctg1\\tACGTACGT\\n" > ${params.prefix}.${hic ? "hic" : "bp"}.\${h}p_ctg.gfa; done; touch ${params.prefix}.hifiasm.log
    """
}

process ASM_IMPORT {
    tag "$name"
    publishDir "${params.outdir}/fasta", mode: params.publish_mode
    input:  tuple val(name), path(fa, stageAs: 'in/*')
    output: tuple val(name), path("${name}.fa"), path("${name}.fa.fai")
    script:
    """
    seqkit seq -w 80 $fa > ${name}.fa
    seqkit faidx ${name}.fa
    """
    stub:
    """
    printf ">ctg1\\nACGTACGT\\n" > ${name}.fa; printf "ctg1\\t8\\t6\\t8\\t9\\n" > ${name}.fa.fai
    """
}

process GFA2FASTA {
    tag "$gfa.baseName"
    publishDir "${params.outdir}/fasta", mode: params.publish_mode
    input:  path gfa
    output: tuple val(gfa.baseName), path("${gfa.baseName}.fa"), path("${gfa.baseName}.fa.fai")
    script:
    """
    awk '/^S/{print ">"\$2; print \$3}' $gfa | seqkit seq -w 80 > ${gfa.baseName}.fa
    seqkit faidx ${gfa.baseName}.fa
    """
}

/* =========================================================================== scaffolding */
process MAP_SEGMENTS {
    tag "$name"
    label 'big'
    input:  tuple val(name), path(fa), path(fai), path(segments)
    output:
        tuple val(name), path("${name}.segments.ns.bam"), emit: bam
        tuple val(name), path("${name}.segments.tsv.gz"), emit: table
    script:
    def mm2t = task.cpus > 5 ? task.cpus - 4 : 1
    def st = task.cpus > 4 ? 4 : 1
    """
    minimap2 -t $mm2t -ax map-hifi --secondary=no --no-hash-name $fa $segments 2> minimap2.log \\
        | samtools sort -n -@ $st -m 2G -T tmp.sort -o ${name}.segments.ns.bam
    # qname flag contig pos mapq ref_span of primary alignments, for CIFI_STATS
    samtools view -@ $st -F 2304 ${name}.segments.ns.bam \\
        | awk -v OFS='\\t' '{ span = 0; c = \$6
                              while (match(c, /[0-9]+[MDN=X]/)) { span += substr(c, RSTART, RLENGTH - 1); c = substr(c, RSTART + RLENGTH) }
                              print \$1, \$2, \$3, \$4, \$5, span }' \\
        | gzip -c > ${name}.segments.tsv.gz
    """
    stub:
    """
    touch ${name}.segments.ns.bam; echo "" | gzip -c > ${name}.segments.tsv.gz
    """
}

process CIFI_CONTACTS {
    tag "$name"
    label 'mid'
    publishDir "${params.outdir}/cifi/contacts", mode: params.publish_mode, pattern: '*_contacts_*'
    input:  tuple val(name), path(bam)
    output:
        tuple val(name), path("${name}.contacts.bed"), emit: bed
        path "${name}.contacts_contacts_*", optional: true
    script: "cifi contacts $bam --format bed -q ${params.cifi_mapq} -o ${name}.contacts.bed -t ${task.cpus}"
    stub:
    """
    touch ${name}.contacts.bed ${name}.contacts_contacts_stats.json
    """
}

process YAHS {
    tag "$name"
    label 'mid'
    publishDir "${params.outdir}/scaffolds/${name}", mode: params.publish_mode, pattern: '{*.yahs*,*_JBAT*,*.log}'
    publishDir "${params.outdir}/fasta", mode: params.publish_mode, pattern: '*.scaffolds.fa'
    input:  tuple val(name), path(fa), path(fai), path(bed)
    output:
        tuple val(name), path("${name}.scaffolds.fa"), emit: scaffolds
        tuple val(name), path("${name}_JBAT.txt"), path("${name}_JBAT.chrom.sizes"), emit: jbat
        tuple val(name), path("${name}.scaffolds.pairs.gz"), emit: pairs, optional: true
        tuple val(name), path("${name}.yahs_scaffolds_final.agp"), emit: agp
        path "{*.yahs*,*_JBAT*,*.log}"
    script:
    "scaffold_yahs.sh $name $fa $fai $bed ${params.cifi_mapq} ${params.contig_ec ? 1 : 0} ${params.skip_pretext ? 0 : 1}"
    stub:
    """
    printf ">s1\\nACGT\\n" > ${name}.scaffolds.fa; printf "s1\\t1\\t4\\t1\\tW\\tctg1\\t1\\t4\\t+\\n" > ${name}.yahs_scaffolds_final.agp; touch ${name}_JBAT.txt ${name}_JBAT.chrom.sizes ${name}.scaffolds.pairs.gz ${name}.yahs.bin
    """
}

process CIFI_STATS {
    tag "$name"
    label 'mid'
    publishDir "${params.outdir}/cifi/stats", mode: params.publish_mode
    input:  tuple val(name), path(table), path(agp)
    output: path "${name}.cifi_{stats.tsv,stats.json,stats.png,scaffolds.tsv}"
    script:
    """
    cifi_stats.py ${name} ${table} ${agp} ${params.cifi_stats_mapq} ${params.cifi_stats_max_pairs} \\
        ${params.cifi_stats_min_scaf} ${params.cifi_stats_bin}
    """
    stub:
    """
    touch ${name}.cifi_stats.tsv ${name}.cifi_stats.png ${name}.cifi_scaffolds.tsv; echo "{}" > ${name}.cifi_stats.json
    """
}

process JUICER_TOOLS {
    tag "$name"
    label 'big'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/scaffolds/${name}", mode: params.publish_mode
    input:
        tuple val(name), path(txt), path(sizes)
        path jar
    output: path "${name}_JBAT.hic"
    script:
    "java -Xmx${Math.max(task.memory.toGiga() - 2, 1)}g -jar $jar pre $txt ${name}_JBAT.hic $sizes"
    stub:
    """
    touch ${name}_JBAT.hic
    """
}

process PRETEXTMAP {
    tag "$name"
    label 'mid'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/scaffolds/${name}", mode: params.publish_mode
    input:  tuple val(name), path(pairs)
    output: tuple val(name), path("${name}.scaffolds.pretext")
    script: "zcat $pairs | PretextMap -o ${name}.scaffolds.pretext"
    stub:
    """
    touch ${name}.scaffolds.pretext
    """
}

process PRETEXTSNAPSHOT {
    tag "$name"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/scaffolds/${name}", mode: params.publish_mode
    input:  tuple val(name), path(pretext)
    output: path '*.png'
    script: "PretextSnapshot -m $pretext -f png -r 1080 --sequences '=full' -o . --prefix ${name}."
    stub:
    """
    touch ${name}.FullMap.png
    """
}

/* =========================================================================== evaluation */
process ASM_STATS {
    tag "$name"
    publishDir "${params.outdir}/eval/stats",     mode: params.publish_mode, pattern: '*.{stats.tsv,gaps.bed}'
    publishDir "${params.outdir}/eval/telomeres", mode: params.publish_mode, pattern: '*.telo.tsv'
    input:  tuple val(name), path(fa)
    output:
        path "${name}.stats.tsv", emit: stats
        path "${name}.{telo.tsv,gaps.bed}"
    script:
    "asm_stats.py $name $fa ${params.telo_motif_f} ${params.telo_motif_r} ${params.telo_window} ${params.telo_min}"
}

process COMPLEASM_DOWNLOAD {
    output: path 'mb_downloads'
    script: "compleasm download ${params.lineage} -L mb_downloads"
    stub:
    """
    mkdir mb_downloads
    """
}

process COMPLEASM {
    tag "$name"
    label 'big'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/eval/compleasm", mode: params.publish_mode
    input:
        tuple val(name), path(fa)
        path lineage_dir
    output: path "${name}.compleasm.txt", emit: summary
    script:
    """
    compleasm run -a $fa -o $name -l ${params.lineage} -L $lineage_dir -t ${task.cpus}
    cp ${name}/summary.txt ${name}.compleasm.txt
    """
    stub:
    """
    printf "S:90.00%%, 9\\nD:1.00%%, 1\\nF:1.00%%, 1\\nM:8.00%%, 8\\n" > ${name}.compleasm.txt
    """
}

process MERYL_COUNT {
    label 'big'
    publishDir "${params.outdir}/eval/merqury", mode: params.publish_mode
    input:  path reads
    output: path "hifi.k${params.kmer}.meryl"
    script:
    """
    meryl count k=${params.kmer} threads=${task.cpus} memory=${task.memory.toGiga()} \\
        $reads output hifi.k${params.kmer}.meryl
    """
    stub:
    """
    mkdir hifi.k${params.kmer}.meryl
    """
}

process MERQURY {
    tag "$set"
    label 'big'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/eval/merqury", mode: params.publish_mode, saveAs: { it.endsWith('.meryl') ? null : it }
    input:
        tuple val(set), path(fas)
        path meryl_db
    output:
        path "${set}.{qv,completeness.stats}", emit: tables
        path "*.{qv,stats,png,hist,bed,wig,log}", emit: files
    script:
    """
    if [ -f /usr/local/env-activate.sh ]; then set +u; . /usr/local/env-activate.sh; set -u; fi
    export OMP_NUM_THREADS=${task.cpus}
    merqury.sh $meryl_db ${fas} $set
    """
    stub:
    """
    for f in ${fas}; do b=\$(basename \$f .fa); printf "%s\\t1\\t2\\t55.5\\t1e-6\\n" \$b >> ${set}.qv; printf "%s\\tall\\t1\\t2\\t98.7\\n" \$b >> ${set}.completeness.stats; done; printf "Both\\t1\\t2\\t56.6\\t1e-6\\n" >> ${set}.qv; touch ${set}.log
    """
}

process DOTPLOT_ALIGN {
    tag "$name"
    label 'big'
    publishDir "${params.outdir}/dotplots", mode: params.publish_mode, pattern: '*.paf.gz'
    input:  tuple val(name), val(qlabel), path(qfa, stageAs: 'q/*'), val(tlabel), path(tfa, stageAs: 't/*')
    output: tuple val(name), val(qlabel), val(tlabel), path("${name}.paf.gz"), path('q.len'), path('t.len')
    script:
    """
    minimap2 -x ${params.dotplot_preset} -t ${task.cpus} --secondary=no $tfa $qfa 2> minimap2.log | gzip -c > ${name}.paf.gz
    for x in q t; do
        zcat -f \$x/* | awk '/^>/{if(n!="")print n"\\t"l; n=substr(\$1,2); l=0; next}{l+=length(\$0)} END{if(n!="")print n"\\t"l}' > \$x.len
    done
    """
    stub:
    """
    echo "" | gzip -c > ${name}.paf.gz; printf "ctg1\\t8\\n" > q.len; cp q.len t.len
    """
}

process DOTPLOT {
    tag "$name"
    publishDir "${params.outdir}/dotplots", mode: params.publish_mode
    input:  tuple val(name), val(qlabel), val(tlabel), path(paf), path(qlen), path(tlen)
    output: path "${name}.dotplot.{png,tsv,json}"
    script:
    "dotplot.py $name $paf $qlen $tlen ${params.dotplot_min_aln} ${params.dotplot_min_seq} $qlabel $tlabel"
    stub:
    """
    touch ${name}.dotplot.png; echo "{}" > ${name}.dotplot.json; printf "query\ttarget\n${qlabel}\t${tlabel}\n" > ${name}.dotplot.tsv
    """
}

process SUMMARY {
    publishDir params.outdir, mode: params.publish_mode
    input:
        path stats
        path busco
        path merq
    output: path 'summary.tsv'
    script:
    """
    summarize_run.py
    column -t summary.tsv
    """
}

process REPORT {
    publishDir params.outdir, mode: params.publish_mode
    input:  // staged in the layout of the output dir, so compare_assemblies.sh reads it as a finished run
        path(fasta,     stageAs: 'fasta/*')
        path(stats,     stageAs: 'eval/stats/*')
        path(telo,      stageAs: 'eval/telomeres/*')
        path(busco,     stageAs: 'eval/compleasm/*')
        path(merq,      stageAs: 'eval/merqury/*')
        path(dotplots,  stageAs: 'dotplots/*')
        path(contact,   stageAs: 'scaffolds/*')
        path(readsqc,   stageAs: 'reads/*')
        path(cifistats, stageAs: 'cifi/stats/*')
        path(markdup,   stageAs: 'cifi/markdup/*')
    output: path 'report'
    script:
    def label = params.prefix.toString().replaceAll(/[^A-Za-z0-9_]/, '_')
    "compare_assemblies.sh -d ${label}=. -p -o report -t 'nxf-cifi: ${params.prefix}'"
    stub:
    """
    mkdir report; touch report/report.html
    """
}

/* =========================================================================== workflow */
workflow {
    validateParameters()
    log.info paramsSummaryLog(workflow)
    def asmFiles = params.assembly ? params.assembly.toString().tokenize(',').collect { it.trim() } : []
    if (params.assembly) {
        if (!params.cifi) error "--assembly needs --cifi (CiFi scaffolding of an existing assembly); see --help"
        if (asmFiles.size() > 2) error "--assembly takes one or two FASTAs (hap1[,hap2]), got ${asmFiles.size()}"
    } else if (!params.hifi || !params.ont) error "--hifi and --ont are required (or --assembly + --cifi); see --help"
    // name of the assembly set: hifiasm bp (no CiFi) / hic (CiFi phasing) or asm (--assembly)
    def mode = params.assembly ? 'asm' : params.cifi ? 'hic' : 'bp'

    // ---------------------------------------------------------------- reads
    // BAM -> FASTQ: HiFi always (hifiasm input), CiFi only for the reads QC (cifi itself reads the BAM)
    ch_cifi_qc_in = params.cifi && !params.skip_reads_qc ? readFiles(params.cifi).map { ['cifi', it] } : Channel.empty()
    ch_conv = (params.hifi ? readFiles(params.hifi).map { ['hifi', it] } : Channel.empty()).mix(ch_cifi_qc_in)
        .branch { _k, f -> bam: f.name.endsWith('.bam'); other: true }
    ch_fq = BAM2FASTQ(ch_conv.bam).mix(ch_conv.other)
    ch_hifi = ch_fq.filter { k, _f -> k == 'hifi' }.map { _k, f -> f }.collect()
    ch_ont_raw = params.ont ? readFiles(params.ont) : Channel.empty()
    ch_ont = (params.ont_minlen as int > 0 ? ONT_FILTER(ch_ont_raw) : ch_ont_raw).collect()

    // reads QC of all starting reads: HiFi, ONT and CiFi (raw, before any filter)
    ch_readsqc = Channel.empty()
    if (!params.skip_reads_qc) {
        ch_qc = ch_hifi.map { ['hifi', it, readType(params.hifi)] }
            .mix(ch_ont_raw.collect().map { ['ont', it, 'fastq'] })          // HiFi / ONT are optional with --assembly
        if (params.cifi)
            ch_qc = ch_qc.mix(ch_fq.filter { k, _f -> k == 'cifi' }.map { _k, f -> f }.collect().map { ['cifi', it, readType(params.cifi)] })
        ch_readsqc = FASTER_REPORT(ch_qc)
    }

    // ---------------------------------------------------------------- CiFi
    ch_hic = Channel.value([])
    ch_segments = Channel.empty()
    ch_markdup = Channel.empty()
    if (params.cifi) {
        def cifiFiles = readFiles(params.cifi).collect()
        // a single BAM is used as is; FASTQ is imported, several files are concatenated
        ch_cifi_bam = cifiFiles.branch { single: it.size() == 1 && it[0].name.endsWith('.bam'); other: true }
        ch_bam = ch_cifi_bam.single.map { it[0] }.mix(CIFI_BAM(ch_cifi_bam.other))
        // PCR duplicates (pbmarkdup, as in the CiFi paper); --cifi_rmdup removes them before digestion
        if (!params.skip_cifi_markdup || params.cifi_rmdup) {
            CIFI_MARKDUP(ch_bam)
            ch_markdup = CIFI_MARKDUP.out.tsv
            if (params.cifi_rmdup) ch_bam = CIFI_MARKDUP.out.bam
        }
        CIFI_QC(ch_bam)
        CIFI_DIGEST(ch_bam)
        ch_hic = CIFI_DIGEST.out.hic.map { r1, r2 -> [r1, r2] }
        ch_segments = CIFI_DIGEST.out.segments
    }

    // ---------------------------------------------------------------- hifiasm, or an existing assembly
    if (params.assembly) {
        // hap1[,hap2] in the given order -> <prefix>.asm.hap1.p_ctg / .hap2.p_ctg (same names as a hifiasm run)
        def asm = asmFiles.withIndex().collect { f, i -> ["${params.prefix}.asm.hap${i + 1}.p_ctg", file(f, checkIfExists: true)] }
        asm.each { n, f -> log.info "  ${f} -> fasta/${n}.fa" }
        ch_ctg = ASM_IMPORT(Channel.fromList(asm))
    } else {
        HIFIASM(ch_hifi, ch_ont, ch_hic)
        ch_ctg = GFA2FASTA(HIFIASM.out.gfa.flatten())    // [name, fa, fai]
    }

    // ---------------------------------------------------------------- scaffolding
    ch_scaf = Channel.empty()
    ch_contact = Channel.empty()
    ch_cifi_stats = Channel.empty()
    if (params.cifi) {
        ch_haps = ch_ctg.filter { n, _fa, _fai -> n.endsWith(".${mode}.hap1.p_ctg") || n.endsWith(".${mode}.hap2.p_ctg") }
        MAP_SEGMENTS(ch_haps.combine(ch_segments))
        CIFI_CONTACTS(MAP_SEGMENTS.out.bam)
        YAHS(ch_haps.join(CIFI_CONTACTS.out.bed))
        if (!params.skip_cifi_stats) ch_cifi_stats = CIFI_STATS(MAP_SEGMENTS.out.table.join(YAHS.out.agp))
        ch_scaf = YAHS.out.scaffolds.map { n, fa -> [n + '.scaffolds', fa] }
        if (params.juicer_jar) JUICER_TOOLS(YAHS.out.jbat, Channel.value(file(params.juicer_jar, checkIfExists: true)))
        if (!params.skip_pretext) ch_contact = PRETEXTSNAPSHOT(PRETEXTMAP(YAHS.out.pairs))
    }

    // ---------------------------------------------------------------- evaluation
    ch_contigs = ch_ctg.map { n, fa, _fai -> [n, fa] }
    ASM_STATS(ch_contigs.mix(ch_scaf))
    ch_telo = ASM_STATS.out[1].flatten().filter { it.name.endsWith('.telo.tsv') }
    ch_gaps = ASM_STATS.out[1].flatten().filter { it.name.endsWith('.gaps.bed') }

    // compleasm + merqury on contigs only: scaffolds carry the same sequence
    ch_busco = Channel.empty()
    if (!params.skip_busco) {
        ch_lin = params.lineage_dir ? Channel.value(file(params.lineage_dir, checkIfExists: true)) : COMPLEASM_DOWNLOAD()
        ch_busco = COMPLEASM(ch_contigs.filter { n, _fa -> !n.endsWith('a_ctg') }, ch_lin).summary
    }

    ch_merq = Channel.empty()
    ch_merq_files = Channel.empty()
    if (!params.skip_merqury && (params.meryl_db || params.hifi)) {     // merqury needs HiFi reads or a meryl DB
        ch_db = params.meryl_db ? Channel.value(file(params.meryl_db, checkIfExists: true)) : MERYL_COUNT(ch_hifi)
        ch_sets = ch_contigs
            .filter { n, _fa -> n.endsWith(".${mode}.hap1.p_ctg") || n.endsWith(".${mode}.hap2.p_ctg") }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .filter { it.size() == 2 }
            .map { l -> ['haps', l.collect { it[1] }] }
            .mix(ch_contigs.filter { n, _fa -> n.endsWith(".${mode}.p_ctg") }.map { _n, fa -> ['primary', [fa]] })
        MERQURY(ch_sets, ch_db)
        ch_merq = MERQURY.out.tables
        ch_merq_files = MERQURY.out.files
    }

    // dot plots: hap2 vs hap1 (contigs, scaffolds) and, with --ref, every assembly vs the reference
    ch_dot = Channel.empty()
    if (!params.skip_dotplot) {
        def ch_asm = ch_contigs.mix(ch_scaf)
        def h1 = ch_asm.filter { n, _fa -> n.contains(".${mode}.hap1.") }.map { n, fa -> [n.replace('.hap1.', '.hapN.'), n, fa] }
        def h2 = ch_asm.filter { n, _fa -> n.contains(".${mode}.hap2.") }.map { n, fa -> [n.replace('.hap2.', '.hapN.'), n, fa] }
        ch_pairs = h2.join(h1).map { _k, n2, fa2, n1, fa1 -> [n2.replace('.hap2.', '.hap2_vs_hap1.'), n2, fa2, n1, fa1] }
        if (params.ref) {
            def ref = file(params.ref, checkIfExists: true)
            def refName = ref.name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
            ch_pairs = ch_pairs.mix(ch_asm.filter { n, _fa -> !n.endsWith('a_ctg') }
                                          .map { n, fa -> ["${n}_vs_${refName}", n, fa, refName, ref] })
        }
        ch_dot = DOTPLOT(DOTPLOT_ALIGN(ch_pairs))
    }

    SUMMARY(
        ASM_STATS.out.stats.collect(),
        ch_busco.collect().ifEmpty([]),
        ch_merq.flatten().collect().ifEmpty([])
    )

    REPORT(
        ch_ctg.flatMap { _n, fa, fai -> [fa, fai] }.mix(ch_scaf.map { _n, fa -> fa }).collect(),
        ASM_STATS.out.stats.mix(ch_gaps).collect(),
        ch_telo.collect(),
        ch_busco.collect().ifEmpty([]),
        ch_merq_files.flatten().collect().ifEmpty([]),
        ch_dot.flatten().collect().ifEmpty([]),
        ch_contact.flatten().collect().ifEmpty([]),
        ch_readsqc.collect().ifEmpty([]),
        ch_cifi_stats.flatten().collect().ifEmpty([]),
        ch_markdup.collect().ifEmpty([])
    )
}
