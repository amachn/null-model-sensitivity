#!/usr/bin/env bash
# Re-run of Figure S7 / Table S3 (Wenz et al. 2026) from the authors' Zenodo ATAC peaks.
# Histone peaks: ENCODE (Roadmap) human liver Histone ChIP-seq, GRCh38 "pseudoreplicated peaks"
# (2021 re-processing), all donors pooled per mark -> reproduces the observed counts of Table S3 exactly.
# Observed = number of rows of `bedtools intersect -a ATAC -b marks` (ATAC-mark overlap pairs).
# Null = 1000 x `bedtools shuffle -chrom` (size- and chromosome-matched random regions).
# Usage (from P03/): bash code/figS7_roadmap_rerun.sh [n_shuffle=1000] [n_parallel=4]
set -euo pipefail
N=${1:-1000}; P=${2:-4}
OUT=results/S7_rerun; mkdir -p $OUT/marks $OUT/shuffles
GENOME=data/hg38.autosomes.genome
ATAC=$OUT/atac_peaks.chr.sorted.bed
[ -s $ATAC ] || gzip -dc data/genrichAllPeaks_m10_g50_9.15.21.noBLnarrowpeak.gz | \
  awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$4}' | sort -k1,1 -k2,2n > $ATAC

# pooled mark files
MARKS=$(cut -f5 data/encode_liver_histone_used/files_used.tsv | sort -u)
for m in $MARKS; do
  for a in $(awk -F'\t' -v m=$m '$5==m{print $2}' data/encode_liver_histone_used/files_used.tsv); do
    gzcat data/encode_liver_histone_used/$a.bed.gz | cut -f1-3
  done | sort -k1,1 -k2,2n > $OUT/marks/$m.bed
done

count_all() {  # $1 = region bed ; prints one count per mark (sorted mark order)
  for m in $MARKS; do bedtools intersect -sorted -a "$1" -b $OUT/marks/$m.bed | wc -l | tr -d ' '; done | paste -sd, -
}
echo "mark,$(echo $MARKS | tr ' ' ',')" > $OUT/observed.csv
echo "observed,$(count_all $ATAC)" >> $OUT/observed.csv
cat $OUT/observed.csv

export OUT GENOME ATAC MARKS
export -f count_all
one() {
  i=$1; f=$OUT/shuffles/$i.csv
  [ -s $f ] && return
  t=$(mktemp); bedtools shuffle -chrom -seed $i -i $ATAC -g $GENOME | sort -k1,1 -k2,2n > $t
  echo "$i,$(count_all $t)" > $f; rm -f $t
}
export -f one
seq 1 $N | xargs -P $P -I{} bash -c 'one {}'
( echo "iter,$(echo $MARKS | tr ' ' ',')"; cat $OUT/shuffles/*.csv | sort -t, -k1,1n ) > $OUT/shuffled_counts.csv
echo "done: $(($(wc -l < $OUT/shuffled_counts.csv)-1)) shuffles"
