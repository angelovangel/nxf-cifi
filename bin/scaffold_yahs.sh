#!/usr/bin/env bash
# YaHS scaffolding of one haplotype from CiFi contacts + Juicebox (JBAT) files + pairs for PretextMap.
# usage: scaffold_yahs.sh NAME CONTIGS.fa CONTIGS.fa.fai CONTACTS.bed MAPQ CONTIG_EC(0/1) PRETEXT(0/1)
# writes NAME.scaffolds.fa, NAME.yahs*, NAME_JBAT.{txt,assembly,liftover.agp,chrom.sizes,log},
#        NAME.scaffolds.pairs.gz (if PRETEXT=1)
set -euo pipefail
name=$1; fa=$2; fai=$3; bed=$4; mapq=$5; ec=$6; pretext=$7

args=(-q "$mapq"); (( ec == 0 )) && args+=(--no-contig-ec)
yp="$name.yahs"
yahs "${args[@]}" -o "$yp" -v 1 "$fa" "$bed" > yahs.log 2>&1 || { tail -n 30 yahs.log >&2; exit 1; }
cp "${yp}_scaffolds_final.fa" "$name.scaffolds.fa"

# Juicebox files for manual curation
jb="${name}_JBAT"
juicer pre -a -q "$mapq" -o "$jb" "$yp.bin" "${yp}_scaffolds_final.agp" "$fai" > "$jb.log" 2>&1
grep PRE_C_SIZE "$jb.log" | awk '{print $2" "$3}' > "$jb.chrom.sizes"

# scaffold-level pairs (juicer pre short format: str1 chr1 pos1 frag1 str2 chr2 pos2 frag2)
if (( pretext )); then
  { echo "## pairs format v1.0"
    awk '/^>/{if(n!="")printf "#chromsize: %s %d\n",n,l; n=substr($1,2); l=0; next}
         {l+=length($0)} END{if(n!="")printf "#chromsize: %s %d\n",n,l}' "$name.scaffolds.fa"
    echo "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2"
    juicer pre -q "$mapq" "$yp.bin" "${yp}_scaffolds_final.agp" "$fai" 2> juicer_pre_scaf.log \
      | awk -v OFS='\t' '{print ".",$2,$3,$6,$7,($1==0?"+":"-"),($5==0?"+":"-")}'
  } | gzip -c > "$name.scaffolds.pairs.gz"
fi
