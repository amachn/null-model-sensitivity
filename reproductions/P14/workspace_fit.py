from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

from enrichlib import Genome, IntervalSet
from pearl2025 import MARKS, load_marks, load_peaksets, load_s7

HERE = Path(__file__).resolve().parent

def build_candidate(spec: str, genome: Genome) -> IntervalSet:
    ws = genome.whole()
    for tok in (t.strip() for t in spec.split(",")):
        if tok in ("", "genome"):
            continue
        sign, path = tok[0], tok[1:]
        if sign not in "+-":
            raise ValueError(f"bad token {tok!r}: expected +PATH, -PATH or 'genome'")
        bed = IntervalSet.from_bed(genome, path)
        ws = ws.intersect(bed) if sign == "+" else ws.subtract(bed)
    return ws

def evaluate(ws: IntervalSet, peaksets: Dict[str, IntervalSet], marks: Dict[str, IntervalSet]) -> pd.DataFrame:
    marks_w = {m: a.intersect(ws) for m, a in marks.items()}
    rows = []
    for pname, seg in peaksets.items():
        seg_w = seg.intersect(ws)
        for m, a in marks_w.items():
            obs = int(a.overlap_bp(seg_w.starts, seg_w.ends).sum())
            rows.append(dict(peakset=pname, annotation=m, observed=obs, segment_bp=seg_w.bp, annotation_bp=a.bp,
                             expected_approx=seg_w.bp * a.bp / ws.bp))
    return pd.DataFrame(rows)

def score(detail: pd.DataFrame, target: pd.DataFrame) -> dict:
    d = detail.merge(target, on=["peakset", "annotation"], how="inner")
    d["obs_err"] = np.abs(np.log(d.observed / d.s7_observed))
    d["exp_err"] = np.abs(np.log(d.expected_approx / d.s7_expected))
    core = d[d.annotation != "H3K9me3"]
    k9 = d[d.annotation == "H3K9me3"]
    return dict(
        n_rows=len(d),
        exact_observed_rows=int((d.observed == d.s7_observed).sum()),
        obs_err_median_4marks=float(core.obs_err.median()), obs_err_max_4marks=float(core.obs_err.max()),
        obs_err_median_H3K9me3=float(k9.obs_err.median()) if len(k9) else np.nan,
        exp_err_median_4marks=float(core.exp_err.median()), exp_err_median_H3K9me3=float(k9.exp_err.median()) if len(k9) else np.nan,
    ), d

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", required=True, type=Path)
    ap.add_argument("--s7", required=True, type=Path)
    ap.add_argument("--candidate", action="append", required=True, metavar="NAME:SPEC")
    ap.add_argument("--chrom-sizes", type=Path, default=HERE / "data" / "mm10.chrom.sizes")
    ap.add_argument("--marks-start-offset", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("workspace_fit"))
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)

    genome = Genome.from_file(a.chrom_sizes)
    marks = load_marks(a.repo, genome, a.marks_start_offset)
    peaksets = load_peaksets(a.repo, genome)
    target = load_s7(a.s7)

    summ, details = [], []
    for item in a.candidate:
        name, spec = item.split(":", 1)
        ws = build_candidate(spec, genome)
        if ws.bp == 0 or len(ws) == 0:
            print(f"skipping candidate {name!r} ({spec}): workspace is empty (0 bp) -- threshold too strict?")
            continue
        s, d = score(evaluate(ws, peaksets, marks), target)
        summ.append(dict(candidate=name, spec=spec, workspace_bp=ws.bp, frac_of_genome=ws.bp / genome.whole().bp,
                         n_intervals=len(ws), mean_interval_bp=ws.bp / len(ws), **s))
        d.insert(0, "candidate", name)
        details.append(d)
    summary = pd.DataFrame(summ).sort_values("obs_err_median_4marks").reset_index(drop=True)
    summary.to_csv(a.out / "workspace_fit_summary.tsv", sep="\t", index=False)
    pd.concat(details).to_csv(a.out / "workspace_fit_detail.tsv", sep="\t", index=False)

    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.float_format", "{:.4f}".format):
        print(summary.drop(columns="spec").to_string(index=False))
    best = summary.iloc[0]
    print(f"\nbest by Observed: {best.candidate}  (median |ln ratio| {best.obs_err_median_4marks:.4f}; "
          f"{int(best.exact_observed_rows)}/{int(best.n_rows)} rows exactly equal to S7)")
    print("An exact match on Observed (0 error, all rows equal) would identify the workspace; anything else is only closer.")
    frag = summary[summary.mean_interval_bp < 5_000]
    if len(frag):
        print(f"\ncaution: {', '.join(frag.candidate)} {'is' if len(frag) == 1 else 'are'} fragmented (mean interval < 5 kb), so the "
              f"closed-form exp_err_* columns are unreliable there (see the docstring). Confirm with pearl2025.py --n-perm 500.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
