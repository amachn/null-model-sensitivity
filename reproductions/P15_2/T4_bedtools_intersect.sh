#!/usr/bin/env bash
set -euo pipefail

if [[ ${1:-} == --worker ]]; then
  i=$2; OUT=$3; G=$4; EXCL=$5
  q=$OUT/tmp/q.$i.bed
  if [[ $i -eq 0 ]]; then
    cp "$OUT/atac.bed" "$q"
  else
    args=(-i "$OUT/atac.bed" -g "$G" -chrom -seed "$i")
    [[ $EXCL != none ]] && args+=(-excl "$EXCL")
    bedtools shuffle "${args[@]}" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' > "$q"
  fi
  for c in "$OUT"/chip/*.bed; do
    lab=$(basename "$c" .bed)
    bedtools intersect -a "$q" -b "$c" -wo |
      awk -v i="$i" -v lab="$lab" '{if(!($4 in s)){s[$4]=1; n++} p++; bp+=$NF}
           END{printf "%d\t%s\t%d\t%d\t%d\n", i, lab, n+0, p+0, bp+0}'
  done > "$OUT/tmp/res.$i.tsv"
  rm -f "$q"
  exit 0
fi

ATAC=""; CHIPS=(); CHAIN=""; NOLIFT=0; G=""; N=1000; T=4; OUT="results_T4"; EXCL=none; MERGE=1
while [[ $# -gt 0 ]]; do
  case $1 in
    --atac) ATAC=$2; shift 2;;       --chip) CHIPS+=("$2"); shift 2;;
    --chain) CHAIN=$2; shift 2;;     --no-lift) NOLIFT=1; shift;;
    --genome) G=$2; shift 2;;        --n) N=$2; shift 2;;
    --threads) T=$2; shift 2;;       --out) OUT=$2; shift 2;;
    --excl) EXCL=$2; shift 2;;       --no-merge) MERGE=0; shift;;
    *) echo "unknown option $1"; exit 1;;
  esac
done
[[ -n $ATAC && ${#CHIPS[@]} -gt 0 && -n $G ]] || { echo "need --atac, --chip (1+), --genome"; exit 1; }
[[ $NOLIFT -eq 1 || -n $CHAIN ]] || { echo "need --chain (hg19->hg38) or --no-lift"; exit 1; }
command -v bedtools >/dev/null || { echo "bedtools not found"; exit 1; }
[[ $NOLIFT -eq 1 ]] || command -v CrossMap >/dev/null || { echo "CrossMap not found (conda install -c bioconda crossmap)"; exit 1; }
mkdir -p "$OUT"/chip "$OUT"/tmp
AUTO='^chr([1-9]|1[0-9]|2[0-2])$'
cat_any() { if [[ $1 == *.gz ]]; then gunzip -c "$1"; else cat "$1"; fi; }

cut -f1-3 "$ATAC" | awk -v re="$AUTO" '$1 ~ re' | LC_ALL=C sort -k1,1 -k2,2n |
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"p"NR}' > "$OUT/atac.bed"
echo "[info] ATAC peaks: $(wc -l < "$OUT/atac.bed") (paper: 2,518,633)"

for f in "${CHIPS[@]}"; do
  lab=$(basename "$f"); lab=${lab%.gz}; lab=${lab%.bed}; lab=${lab%.narrowPeak}; lab=${lab%.broadPeak}
  raw=$OUT/tmp/$lab.raw.bed
  cat_any "$f" | grep -vE '^(track|browser|#)' | cut -f1-3 |
    awk 'BEGIN{OFS="\t"} {c=$1; if (c !~ /^chr/) c="chr"c; print c,$2,$3}' > "$raw"
  n0=$(wc -l < "$raw")
  if [[ $NOLIFT -eq 0 ]]; then
    CrossMap bed "$CHAIN" "$raw" "$raw.hg38" > /dev/null 2>&1; mv "$raw.hg38" "$raw"
  fi
  awk -v re="$AUTO" '$1 ~ re' "$raw" | LC_ALL=C sort -k1,1 -k2,2n > "$raw.s"
  if [[ $MERGE -eq 1 ]]; then bedtools merge -i "$raw.s" > "$OUT/chip/$lab.bed"; else cp "$raw.s" "$OUT/chip/$lab.bed"; fi
  echo "[info] $lab: $n0 peaks -> $(wc -l < "$raw") after liftover -> $(wc -l < "$OUT/chip/$lab.bed") autosomal$([[ $MERGE -eq 1 ]] && echo ", merged")"
  rm -f "$raw" "$raw.s"
done

cat > "$OUT/run_settings.txt" <<EOF
atac=$ATAC
chip=${CHIPS[*]}
liftover=$([[ $NOLIFT -eq 1 ]] && echo none || echo "CrossMap $CHAIN")
null=bedtools shuffle -chrom (size and chromosome matched)
exclusion=$EXCL
iterations=$N
chip_merged=$MERGE
count_definition=distinct ATAC peaks with >=1 bp overlap (primary); overlap pairs also reported
bp_definition=sum of overlapping bp from bedtools intersect -wo
p_value=proportion of shuffles >= observed (paper); enrichment=observed/median(shuffles)
EOF

echo "[info] observed + $N shuffles on $T threads"
start=$(date +%s)
bash "$0" --worker 0 "$OUT" "$G" "$EXCL"
bash "$0" --worker 1 "$OUT" "$G" "$EXCL"
per=$(( $(date +%s) - start ))
echo "[info] ~$(( per * N / 2 / T / 60 + 1 )) min for all shuffles (rough estimate)"
seq 2 "$N" | xargs -P "$T" -I{} bash "$0" --worker {} "$OUT" "$G" "$EXCL"
{ printf "iter\tmark\tn_atac_overlapping\tn_pairs\tbp_overlap\n"; cat "$OUT"/tmp/res.*.tsv | sort -k1,1n -k2,2; } > "$OUT/iterations.tsv"
rm -rf "$OUT/tmp"

python3 - "$OUT" <<'PY'
import csv, statistics, sys, os
out = sys.argv[1]
rows = list(csv.DictReader(open(os.path.join(out, "iterations.tsv")), delimiter="\t"))
pub = {"H3K27ac": 1.86, "H3K4me3": 1.83}
res = []
for mark in sorted({r["mark"] for r in rows}):
    for metric in ("n_atac_overlapping", "bp_overlap", "n_pairs"):
        obs = [float(r[metric]) for r in rows if r["mark"] == mark and r["iter"] == "0"][0]
        null = [float(r[metric]) for r in rows if r["mark"] == mark and r["iter"] != "0"]
        med = statistics.median(null)
        res.append(dict(mark=mark, metric=metric, observed=int(obs), null_median=med,
                        fold=round(obs / med, 3) if med else float("nan"),
                        p_paper=round(sum(x >= obs for x in null) / len(null), 4),
                        n_shuffles=len(null),
                        published_fold=next((v for k, v in pub.items() if k.lower() in mark.lower()), "")))
with open(os.path.join(out, "TableS4_reproduction.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(res[0]), delimiter="\t"); w.writeheader(); w.writerows(res)
print(f"\n{'mark':<14}{'metric':<22}{'observed':>11}{'null_median':>13}{'fold':>8}{'p_paper':>9}{'published':>11}")
for r in res:
    print(f"{r['mark']:<14}{r['metric']:<22}{r['observed']:>11}{r['null_median']:>13.1f}{r['fold']:>8}{r['p_paper']:>9}{str(r['published_fold']):>11}")
print(f"\n[done] {out}/TableS4_reproduction.tsv  (p_paper = 0 means p < 1/{res[0]['n_shuffles']})")
PY
