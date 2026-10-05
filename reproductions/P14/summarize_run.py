from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ALPHA = 0.05

def call(l2fold: float, p: float) -> str:
    return "n.s." if p >= ALPHA else ("enriched" if l2fold > 0 else "depleted")

def load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    need = {"s7_l2fold", "s7_pvalue", "s7_observed", "s7_expected"}
    if not need <= set(df.columns):
        raise SystemExit(f"{path} has no S7 comparison columns; rerun pearl2025.py with --s7 (or use --from-results).")
    df = df.dropna(subset=["s7_l2fold"]).copy()
    df["ours"] = [call(l, p) for l, p in zip(df.l2fold, df.pvalue)]
    df["s7_call"] = [call(l, p) for l, p in zip(df.s7_l2fold, df.s7_pvalue)]
    df["agree"] = df.ours == df.s7_call
    return df

def summarize(df: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for w, g in df.groupby("workspace", sort=False):
        core = g[g.annotation != "H3K9me3"]
        rows.append({
            "run": label, "workspace": w, "n_perm": int(g.n_perm.iloc[0]), "rows": len(g),
            "calls agree": f"{int(g.agree.sum())}/{len(g)}",
            "median dL2FC": g.l2fold_minus_s7.median(),
            "median |dL2FC|": g.l2fold_minus_s7.abs().median(),
            "Spearman": g.l2fold.rank().corr(g.s7_l2fold.rank()),
            "median Obs/S7": g.obs_ratio_vs_s7.median(),
            "median Exp/S7 (4 marks)": core.exp_ratio_vs_s7.median(),
            "median Exp/S7 (H3K9me3)": g[g.annotation == "H3K9me3"].exp_ratio_vs_s7.median(),
        })
    return pd.DataFrame(rows)

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", type=Path, help="fig3e_reproduction.tsv of the run to summarize")
    ap.add_argument("--compare", type=Path, action="append", default=[], help="earlier run(s) to put alongside")
    a = ap.parse_args(argv)

    runs = [(a.results.parent.name or "run", load(a.results))]
    for c in a.compare:
        runs.append((c.parent.name or "compare", load(c)))
    seen: dict[str, int] = {}
    labelled = []
    for lab, df in runs:
        seen[lab] = seen.get(lab, 0) + 1
        labelled.append((lab if seen[lab] == 1 else f"{lab}#{seen[lab]}", df))

    table = pd.concat([summarize(df, lab) for lab, df in labelled], ignore_index=True)
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.float_format", "{:.3f}".format):
        print("closeness to Table S7 (calls at p < 0.05; dL2FC = our log2FC minus S7's)\n")
        print(table.to_string(index=False))

        main_label, main_df = labelled[0]
        cols = ["workspace", "peakset", "annotation", "l2fold", "pvalue", "ours", "s7_l2fold", "s7_pvalue", "s7_call"]
        bad = main_df[~main_df.agree][cols]
        print(f"\ndisagreeing rows in {main_label} ({len(bad)}):")
        print(bad.round(4).to_string(index=False) if len(bad) else "  none")

        k9 = main_df[(main_df.peakset == "all_reproducible") & (main_df.annotation == "H3K9me3")]
        print(f"\nall_reproducible x H3K9me3 in {main_label}   (S7: log2FC {k9.s7_l2fold.iloc[0]:.3f}, p = {k9.s7_pvalue.iloc[0]:.3g}, {k9.s7_call.iloc[0]})")
        print(k9[["workspace", "observed", "s7_observed", "expected", "s7_expected", "l2fold", "pvalue", "ours"]].round(4).to_string(index=False))
    print(f"\nnote: our p-values floor at 1/(n_perm+1) ({1 / (int(main_df.n_perm.iloc[0]) + 1):.1e}); S7's floor at 1e-05.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
