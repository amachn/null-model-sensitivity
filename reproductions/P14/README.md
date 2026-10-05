P14 - Pearl et al. 2025
Reproduction Status: awaiting QC

## Target

* Target result: Figure 3E ~ enrichment of EZH2 and histone-mark (H3K27ac, H3K27me3, H3K4me3, H3K9me3)
  ChIP-seq peaks in HTT ChIP-seq peak regions
* Method: custom Python re-implementation of GAT (Genomic Association Tester) — nucleotide-overlap
  permutation test (no GAT/analysis code was released; see Issues)
* Permutations: 10,000 (original study used 100,000)
* Genome build: mm10

### Code

* No original analysis code was released for this result. The following scripts were written
  from scratch, based on the paper's Methods, to reproduce it:
  * `enrichlib.py` -> `scripts/` — interval algebra, permutation null, overlap statistics (the GAT re-implementation)
  * `pearl2025.py` -> `scripts/` — loads the paper's peak sets and runs the reproduction
  * `workspace_fit.py`, `make_workspace_bed.py`, `diagnose_workspace.py` -> `scripts/` — used to investigate the
    unspecified GAT workspace

### Data

* HTT ChIP-seq peak coordinates (`all_reproducible`, n = 9,624) and EZH2/H3K27ac/H3K27me3/H3K4me3/H3K9me3
  ChIP-seq peak coordinates were obtained from the authors' GitHub repository (`seth-ament/ament-carroll-collab`)
  and stored in `data/`
* ENCODE mm10 blacklist (v2) and Umap k36 mappability track were obtained externally and stored in `data/ref/`

## Environment

* Operating System: macOS Sequoia 15.7.4
* Hardware: MacBook Pro, Apple M2 Pro, 16 GB RAM
* Editor: VS Code
* Python Version: 3.11.4
* Key packages: numpy, pandas, openpyxl (see `requirements.txt`)

### Environment Setup

The environment is stored in a `requirements.txt` file; to recreate the environment used in this reproduction, run:

```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Results

The reproduced results qualitatively reproduced the overall pattern in Figure 3E: enrichment of HTT peaks at
EZH2, H3K27ac and H3K4me3, and depletion at H3K27me3 (17 of 20 peak-set x mark comparisons matched the
published enrichment/depletion call). Effect sizes showed consistent quantitative deviation from the original
(log2FC differing by a median of 0.3-0.6, depending on workspace). H3K9me3 did not reproduce: the original
reports no significant overlap, while our result is significant and its direction (enriched vs depleted)
depends on which workspace is used. This is most likely explained by the unspecified genome mask (below)
rather than a difference in analysis.

## Issues / Notes

* Genome mask (workspace) used in the permutation test
  * The original Methods state permutations were run "within the mappable genome" but do not specify or provide
    this mask, so it could not be reconstructed exactly.
  * Results were generated for two candidate workspaces: whole genome minus the ENCODE blacklist, and a strict
    Umap k36 multi-read mappability mask minus the blacklist.
  * The mappability mask's deterministic overlap counts were closest to the published values (median error
    0.6% vs 5.7% for blacklist-only), so it appears closest to what was originally used, though this is not
    confirmed.
* No analysis code was available
  * The Data and resource availability statement links a GitHub repository, but it contains only peak-calling
    and caQTL-mapping scripts, not the enrichment/permutation analysis.
* Extra HTT peak files
  * The peak set used in Figure 3E (9,624 peaks, "all reproducible") is not included in the paper's
    supplementary tables; Table S1 contains only a smaller, differently filtered subset (4,900 peaks).
  * The 9,624-peak set was instead obtained from the authors' GitHub repository.
* Ambiguous peak-set labeling in Table S7
  * Table S7 lists three genotype-specific peak subsets (WT-specific, mHTT-specific, shared) per mark with no
    column identifying which row is which, and the row order is not consistent across marks.
  * Rows were assigned by matching each row's expected overlap to the peak-set size it implies; this
    assignment is inferred, not confirmed by the authors.