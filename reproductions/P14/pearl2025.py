from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from enrichlib import Genome, IntervalSet, run_enrichment

HERE = Path(__file__).resolve().parent
BLACKLIST_URL = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz"

B2020 = "Htt_ChIPSeq_Files/09_B2020_WT_QQ_data_analysis"
_STEM = "Htt_ChIPSeq.downsample_30000000.filtered_p1e3.merged_peaks.any_2_samples.final"
HTT_ALL = f"{B2020}/{_STEM}.non_blacklisted.bed"
HTT_COND = f"{B2020}/{_STEM}.conditionally_reproducible.txt"
MARK_DIR = "Htt_ChIPSeq_Files/06_ActiveMotif_Histone_Marks"
MARK_FILE = "{m}_Striata_intervals.merged.any_3_samples.final.annotated_with_description.txt"
MARKS = ["EZH2", "H3K27Ac", "H3K27me3", "H3K4me3", "H3K9me3"]
S7_COND_BY_EXPECTED_RANK = ["mHTT-specific", "WT-mHTT-shared", "WT-specific"]
S7_COND_ORDER = S7_COND_BY_EXPECTED_RANK

SEED_OFFSET = {"all_reproducible": 0, "cond:mHTT-specific": 1, "cond:WT-mHTT-shared": 2, "cond:WT-specific": 3}

def load_marks(repo: Path, genome: Genome, start_offset: int) -> dict[str, IntervalSet]:
    out = {}
    for m in MARKS:
        path = repo / MARK_DIR / MARK_FILE.format(m=m)
        df = pd.read_csv(path, sep="\t", usecols=["Chr", "Start", "End"], dtype={"Chr": str})
        out[m] = IntervalSet.from_frame(genome, df, chrom="Chr", start="Start", end="End", start_offset=start_offset)
    return out

def load_peaksets(repo: Path, genome: Genome) -> dict[str, IntervalSet]:
    sets = {"all_reproducible": IntervalSet.from_bed(genome, repo / HTT_ALL)}
    cond = pd.read_csv(repo / HTT_COND, sep="\t", header=None, usecols=[0, 1, 2, 4],
                       names=["chrom", "start", "end", "label"], dtype={0: str})
    for label in S7_COND_ORDER:
        sets[f"cond:{label}"] = IntervalSet.from_frame(genome, cond[cond.label == label])
    return sets

def fetch_blacklist(dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        try:
            urllib.request.urlretrieve(BLACKLIST_URL, dest)
        except Exception as exc:
            sys.exit(f"Could not download the mm10 blacklist ({exc}).\n"
                     f"Download {BLACKLIST_URL} yourself and pass --blacklist PATH.")
    return dest

def build_workspace(spec: str, genome: Genome, blacklist: Path | None, out: Path) -> IntervalSet:
    whole = genome.whole()
    if spec == "genome":
        return whole
    if spec == "genome-minus-blacklist":
        bl = IntervalSet.from_bed(genome, blacklist or fetch_blacklist(out / "ref" / "mm10-blacklist.v2.bed.gz"))
        return whole.subtract(bl)
    return whole.intersect(IntervalSet.from_bed(genome, spec))

def load_s7(path: Path) -> pd.DataFrame:
    s7 = pd.read_excel(path)
    for c in ["Pvalue", "Qvalue", "Adjusted.Pvalue", "Expected"]:
        s7[c] = pd.to_numeric(s7[c], errors="coerce")
    s7 = s7[s7["Annotation"].isin(MARKS)].copy()
    s7["kind"] = s7["Peakset"].str.rsplit(".", n=1).str[-1]
    out = []
    for _, r in s7[s7.kind == "all_reproducible"].iterrows():
        out.append((r, "all_reproducible", "exact"))
    for ann, g in s7[s7.kind == "conditionally_reproducible"].groupby("Annotation", sort=False):
        ranks = g["Expected"].rank(method="first").astype(int)
        for (_, r), k in zip(g.iterrows(), ranks):
            out.append((r, f"cond:{S7_COND_BY_EXPECTED_RANK[k - 1]}", "by expected-rank"))
    return pd.DataFrame([dict(peakset=key, annotation=r["Annotation"], s7_row_assignment=how,
                              s7_observed=r["Observed"], s7_expected=r["Expected"], s7_l2fold=r["Log2FC"],
                              s7_pvalue=r["Pvalue"], s7_adj_pvalue=r["Adjusted.Pvalue"]) for r, key, how in out])

def attach_s7(results: pd.DataFrame, s7_path: Path) -> pd.DataFrame:
    results = results.drop(columns=[c for c in results.columns if c.startswith("s7_") or c.endswith("_vs_s7")
                                    or c in ("l2fold_minus_s7", "implied_workspace_bp_upper")])
    results = results.merge(load_s7(s7_path), on=["peakset", "annotation"], how="left")
    results["obs_ratio_vs_s7"] = results["observed"] / results["s7_observed"]
    results["exp_ratio_vs_s7"] = results["expected"] / results["s7_expected"]
    results["l2fold_minus_s7"] = results["l2fold"] - results["s7_l2fold"]
    results["implied_workspace_bp_upper"] = results["workspace_bp"] * results["expected"] / results["s7_expected"]
    return results

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", type=Path, help="path to a clone of seth-ament/ament-carroll-collab (not needed with --from-results)")
    ap.add_argument("--from-results", type=Path,
                    help="skip the permutation run: reload an earlier fig3e_reproduction.tsv and redo the S7 comparison")
    ap.add_argument("--s7", type=Path, help="Table S7 xlsx; adds comparison columns")
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--chrom-sizes", type=Path, default=HERE / "data" / "mm10.chrom.sizes")
    ap.add_argument("--workspace", default="genome,genome-minus-blacklist",
                    help="comma-separated: genome | genome-minus-blacklist | BED path")
    ap.add_argument("--blacklist", type=Path, help="local mm10 blacklist BED(.gz); otherwise downloaded")
    ap.add_argument("--n-perm", type=int, default=10_000, help="paper used 100,000")
    ap.add_argument("--counter", default="nucleotide-overlap", choices=["nucleotide-overlap", "segment-overlap"])
    ap.add_argument("--peaksets", default="all,cond", help="'all' (9,624 peaks) and/or 'cond' (3 genotype groups)")
    ap.add_argument("--marks-start-offset", type=int, default=0,
                    help="added to mark-peak starts; use -1 if the annotated tables are 1-based (effect <0.5%%)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv)

    a.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if a.from_results:
        results = pd.read_csv(a.from_results, sep="\t")
    else:
        if not a.repo:
            ap.error("--repo is required unless --from-results is given")
        results = run_all(a)
    if a.s7:
        results = attach_s7(results, a.s7)
    return report(a, results, t0)

def run_all(a) -> pd.DataFrame:
    genome = Genome.from_file(a.chrom_sizes)
    marks = load_marks(a.repo, genome, a.marks_start_offset)
    peaksets = load_peaksets(a.repo, genome)
    want = set(a.peaksets.split(","))
    peaksets = {k: v for k, v in peaksets.items()
                if (k == "all_reproducible" and "all" in want) or (k.startswith("cond:") and "cond" in want)}

    print(f"genome: {len(genome.chroms)} chromosomes, {genome.total:,} bp (global axis)")
    for k, v in {**peaksets, **marks}.items():
        note = f"  ({v.dropped} rows on unknown chromosomes dropped)" if v.dropped else ""
        print(f"  {k:<26s} {len(v):>7,d} merged intervals  {v.bp:>13,d} bp{note}")

    frames = []
    for wspec in a.workspace.split(","):
        ws = build_workspace(wspec, genome, a.blacklist, a.out)
        wname = Path(wspec).name if wspec not in ("genome", "genome-minus-blacklist") else wspec
        print(f"\nworkspace = {wname}: {len(ws):,} intervals, {ws.bp:,} bp")
        for pname, seg in peaksets.items():
            t = time.time()
            res = run_enrichment(seg, marks, ws, n_perm=a.n_perm, counter=a.counter, seed=a.seed + SEED_OFFSET[pname])
            res.insert(0, "peakset", pname)
            res.insert(0, "workspace", wname)
            res["workspace_bp"] = ws.bp
            frames.append(res)
            print(f"  {pname:<26s} {len(seg):>6,d} segments  {time.time() - t:5.1f}s")
    return pd.concat(frames, ignore_index=True)

def report(a, results: pd.DataFrame, t0: float) -> int:
    tsv = a.out / "fig3e_reproduction.tsv"
    results.to_csv(tsv, sep="\t", index=False)
    if not a.from_results:
      (a.out / "run_info.json").write_text(json.dumps(
          dict(argv=sys.argv, n_perm=a.n_perm, seed=a.seed, counter=a.counter, workspaces=a.workspace,
               numpy=np.__version__, pandas=pd.__version__, seconds=round(time.time() - t0, 1)), indent=2))

    show = ["workspace", "peakset", "annotation", "observed", "expected", "l2fold", "pvalue"]
    if a.s7:
        show += ["s7_observed", "obs_ratio_vs_s7", "s7_l2fold", "l2fold_minus_s7", "exp_ratio_vs_s7"]
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.float_format", "{:,.3f}".format):
        print("\n" + results[show].to_string(index=False))
    if a.s7:
        d = results[(results.peakset == "all_reproducible") & (results.annotation != "H3K9me3")]
        for w, g in d.groupby("workspace"):
            print(f"[diagnostic] workspace={w}: implied S7 workspace <= {g.implied_workspace_bp_upper.min() / 1e9:.2f}-"
                  f"{g.implied_workspace_bp_upper.max() / 1e9:.2f} Gb (EZH2/H3K27ac/H3K27me3/H3K4me3), "
                  f"vs {g.workspace_bp.iloc[0] / 1e9:.2f} Gb used here")
    print(f"\nwrote {tsv}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
