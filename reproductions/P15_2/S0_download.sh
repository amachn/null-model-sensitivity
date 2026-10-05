#!/usr/bin/env bash
set -uo pipefail
D=${1:?usage: $0 <datadir>}
mkdir -p "$D"/{zenodo,ref}
dl() { [[ -s $2 ]] && return 0; wget -c -q --show-progress -O "$2" "$1" || { echo "[warn] failed: $1"; rm -f "$2"; }; }

Z=https://zenodo.org/records/18329519/files
dl "$Z/genrichAllPeaks_m10_g50_9.15.21.noBLnarrowpeak.gz?download=1" "$D/zenodo/genrichAllPeaks_m10_g50_9.15.21.noBLnarrowpeak.gz"
dl "$Z/rasqual_liver_175samples_fdr5_caqtls_4.1.24.txt.gz?download=1" "$D/zenodo/rasqual_liver_175samples_fdr5_caqtls_4.1.24.txt.gz"
dl "$Z/Files_column_description.rtf?download=1" "$D/zenodo/Files_column_description.rtf"

dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes "$D/ref/hg38.chrom.sizes"
dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz "$D/ref/hg38.gap.txt.gz"
dl https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz "$D/ref/hg38-blacklist.v2.bed.gz"
dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz "$D/ref/hg38.fa.gz"

echo "[done] $D  (next: Rscript S1_prepare_query.R ...)"
