# P11 – Reproduction: Iotchkova et al. 2019 (GARFIELD), Nat Genet 51:343-353, doi:10.1038/s41588-018-0322-6

Assignment step covered: reproduce a published genomic-enrichment result with the authors' own
code and data (the later null-model re-test with EGRE / GAT / regioneR is not part of this
deliverable).

**Note:** the paper never uses regioneR or GAT. It compares five enrichment tools — GARFIELD,
GREGOR, LDSC, fgwas, GoShifter — across 21 GWAS traits and DHS / H3K27ac / H3K4me3 annotations
(Fig. 4). GARFIELD is the paper's own method (Fig. 3), and is the primary reproduction target here.

Every file below is either (a) required to run the reproduction, or (b) direct evidence of its
result. Files that are easy to regenerate from the two (prep files kept only as intermediate build
products of `make`, per-panel image variants, etc.) were removed — see "Reproducing this."

```
P11/
├── code/
│   ├── garfield_v2/                          authors' released GARFIELD source (only the scripts actually run)
│   ├── authors_Figure4_ST4_ST5_SF3.R         authors' original Fig. 4 script (reference only)
│   └── fig4_recreate.py                      my Fig. 4 recreation code
├── data/
│   └── SuppTable5.xlsx                       authors' published data
└── results/
    ├── fig3a_garfield_rerun/                 PRIMARY: authors' code, rerun by me
    └── fig4_recreation/                      SECONDARY: re-plot of authors' own published table
```

## Primary reproduction: `results/fig3a_garfield_rerun/` — Fig. 3a (GARFIELD, height, DHS hotspots)

Compiled and ran the authors' unmodified GARFIELD v2 source (`code/garfield_v2/`, only the
`DATADIR`/`INPUTNAME`/`OUTDIR` paths edited in the `garfield` wrapper) on the authors' publicly
released reference data, for the height GWAS shipped as the package's own example. Full pipeline:
`garfield-prep-chr` (LD pruning/annotation) → `garfield-Meff-Padj.R` (multiple-testing correction)
→ `garfield-test.R` (enrichment test) → `garfield-plot.R` (the authors' own plotting code, which
sources `garfield-plot-function.R`). `code/garfield_v2/` here holds only the four scripts the
wrapper actually calls, plus `makefile`/`garfield-prep-chr.cpp` to build the pruning tool and the
authors' `LICENSE`/`README`; their other released utilities (custom-GWAS formatting, custom
annotation building, a post-hoc variant-extraction script) were never invoked and are not included.

| file | what it is |
|---|---|
| `paper_Fig3a_height_published.png` | published Fig. 3a, cropped from the paper |
| `garfield_rerun_HGT_wheel_plot.png` | same plot, drawn by the authors' `garfield-plot.R` from my rerun |
| `garfield_rerun_vs_published_scatter.png` | quantitative check: rerun vs. published −log10 P, all 424 DHS annotations, T<1e-8 |
| `garfield_rerun_vs_ST5_HGT.csv` / `..._summary.csv` | the numbers behind that scatter plot |
| `raw_output/garfield.test.GIANT_HEIGHT.out` | raw rerun output (all annotations, all thresholds) — ground truth the scatter/CSVs above were derived from |
| `raw_output/garfield.Meff.GIANT_HEIGHT.out` | effective-annotation count (Meff = 492.56) |
| `raw_output/garfield.test.GIANT_HEIGHT.out.Hotspots.pdf` | source PDF for the rerun wheel plot |

**Result:** Pearson r = 0.86 (rerun vs. published −log10 P); 363/364 published significant DHS
enrichments recovered (99.7%); median offset +1.94 −log10 P (rerun slightly more significant).

**Correction (kept for transparency):** an earlier pass tried a second binning (`m1,n15,t5`) copied
from `garfield_run_dhs.sh` in the authors' `manuscript_custom_code`, on the assumption it was the
real Fig. 4 setting. On closer reading, that script line actually calls a different, unpublished
tool (`garfield-perm-feat.R`) and writes a file named `...timing.out` — it was the authors' CPU-timing
benchmark, not the real analysis, and `m1` disables MAF binning entirely, which inflated every
P-value by ~13 −log10 units. That run and its outputs were discarded. The result above uses only the
standard, documented `garfield-test.R` default settings (`-b m5,n5,t5`).

A shipped Crohn's example (`cd-meta`) was also run as a sanity check, but is not the paper's actual
CD GWAS, so it wasn't a like-for-like comparison and was dropped from this folder rather than kept
as inconclusive filler.

## Secondary: `results/fig4_recreation/` — five-tool comparison, from the authors' published table

Not an independent computation — a re-plot of the authors' own summarized results, since the raw
per-method files behind Fig. 4 were never released. Included to show *why* null-model choice
matters (the theme of the assignment): the five tools disagree substantially on the same data.

| file | what it is |
|---|---|
| `Fig4_recreated.png` | recreated figure, all 4 panels (from `data/SuppTable5.xlsx` via `code/fig4_recreate.py`) |
| `validation_vs_paper.csv` | median/max/mean significant DHS enrichments per method — match the paper exactly |

Per-panel image variants (a/b/c/d individually) and the intermediate per-trait CSVs are not stored
here since they're one command away — see "Reproducing this."

## Deviations / limits
1. Fig. 4 recreation uses Supp. Table 5 (authors' summarised output), not the raw per-method results (unpublished); box size in panels b–d is approximated as scaled −log10 P (no raw enrichment statistic is published).
2. Fig. 4a bars at "1 method" differ slightly from the printed figure (GARFIELD ≈0.13 vs. ≈0.10 printed; LDSC ≈0.05 vs. ≈0.17 printed); bars for 2–5 agree.
3. GARFIELD re-run (Fig. 3a) is close but not exact: the paper used stricter r²≥0.01 LD-pruning tags (`tags/r001`, unreleased); the public archive ships only r²≥0.1 tags, which likely explains the small residual +1.94 median offset and 41 extra calls.
4. Only GARFIELD was independently re-run; GREGOR, LDSC, fgwas, GoShifter were not (their inputs were never published).
5. The rerun's quantitative check (r=0.86, 363/364) covers only the T<1e-8 layer of Fig. 3a, since that's the one threshold Supp. Table 5 publishes; the wheel-plot image itself covers all 8 thresholds, straight from the rerun's own output.

## Reproducing this
```
cd code/garfield_v2 && make                     # compiles garfield-prep-chr
# download & extract the authors' garfield-data.tar.gz from ebi.ac.uk/birney-srv/GARFIELD/
DATADIR=/path/to/garfield-data INPUTNAME=GIANT_HEIGHT OUTDIR=../../results/fig3a_garfield_rerun/raw_output ./garfield
cd ../.. && python3 code/fig4_recreate.py       # regenerates Fig4_recreated.png/pdf, all 4 individual panels, and all intermediate CSVs
```
