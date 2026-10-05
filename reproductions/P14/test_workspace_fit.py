import numpy as np
import pandas as pd

from enrichlib import Genome, IntervalSet
from workspace_fit import build_candidate, evaluate, score

G = Genome({"c1": 400_000, "c2": 300_000, "c3": 200_000})

def rand_set(rng, n, lo, hi):
    rows = []
    for _ in range(n):
        c = rng.choice(list(G.sizes))
        s = int(rng.integers(0, G.sizes[c] - hi))
        rows.append((c, s, s + int(rng.integers(lo, hi))))
    return IntervalSet.from_frame(G, pd.DataFrame(rows, columns=["chrom", "start", "end"]))

def make_world(seed=0):
    rng = np.random.default_rng(seed)
    peaks = {"all_reproducible": rand_set(rng, 600, 100, 2000), "cond:WT-specific": rand_set(rng, 200, 100, 1500)}
    marks = {m: rand_set(rng, 900, 200, 6000) for m in ["EZH2", "H3K27Ac", "H3K27me3", "H3K4me3", "H3K9me3"]}
    hidden_excl = rand_set(rng, 40, 5_000, 25_000)
    return peaks, marks, hidden_excl

def target_from(ws, peaks, marks):
    d = evaluate(ws, peaks, marks)
    return d[["peakset", "annotation"]].assign(s7_observed=d.observed, s7_expected=d.expected_approx)

def test_hidden_workspace_is_recovered_exactly_and_decoys_are_worse():
    peaks, marks, excl = make_world()
    hidden = G.whole().subtract(excl)
    target = target_from(hidden, peaks, marks)

    s_true, _ = score(evaluate(hidden, peaks, marks), target)
    assert s_true["exact_observed_rows"] == s_true["n_rows"] == 10
    assert s_true["obs_err_median_4marks"] == 0.0 and s_true["exp_err_median_4marks"] == 0.0

    rng = np.random.default_rng(99)
    decoys = {
        "whole genome": G.whole(),
        "different mask, same size": G.whole().subtract(rand_set(rng, 40, 5_000, 25_000)),
        "hidden mask + extra hole": hidden.subtract(rand_set(rng, 10, 5_000, 25_000)),
    }
    for name, ws in decoys.items():
        s, _ = score(evaluate(ws, peaks, marks), target)
        assert s["exact_observed_rows"] < s["n_rows"], name
        assert s["obs_err_median_4marks"] > 0 or s["obs_err_max_4marks"] > 0, name

def test_build_candidate_grammar(tmp_path):
    a = tmp_path / "keep.bed"
    b = tmp_path / "drop.bed"
    a.write_text("c1\t0\t100000\nc2\t0\t100000\n")
    b.write_text("c1\t40000\t60000\n")
    ws = build_candidate(f"+{a},-{b}", G)
    assert ws.bp == 200_000 - 20_000
    assert build_candidate("genome", G).bp == G.whole().bp
