#!/usr/bin/env bash
set -euo pipefail
D=${1:?usage: $0 <datadir>}
mkdir -p "$D/chip" "$D/ref"
G=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/GSE128072/suppl
for m in H3K27ac H3K4me3; do
  f="$D/chip/GSE128072_Pr1_${m}_FeatureCounts.txt.gz"
  [[ -s $f ]] || curl -fsSL -o "$f" "$G/GSE128072_Pr1_${m}_FeatureCounts.txt.gz"
  gunzip -c "$f" |
    awk 'BEGIN{OFS="\t"} /^#/ || $1=="Geneid" {next}
         { n=split($2,c,";"); split($3,s,";"); split($4,e,";");
           for (i=1;i<=n;i++) print c[i], s[i]-1, e[i], $1 }' > "$D/chip/${m}.bed"
  echo "[info] $m: $(wc -l < "$D/chip/${m}.bed") regions -> $D/chip/${m}.bed"
done
c="$D/ref/hg19ToHg38.over.chain.gz"
[[ -s $c ]] || curl -fsSL -o "$c" https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
echo "[done] $D/chip and $c"
