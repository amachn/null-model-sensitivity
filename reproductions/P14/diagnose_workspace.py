from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from enrichlib import Genome, IntervalSet
from pearl2025 import load_peaksets

HERE = Path(__file__).resolve().parent

def valid_starts(ws: IntervalSet, piece_lengths: np.ndarray) -> np.ndarray:
    lens = ws.ends - ws.starts
    desc = np.sort(lens)[::-1]
    prefix = np.concatenate([[0], np.cumsum(desc)]).astype(np.int64)
    n_fit = np.searchsorted(-desc, -piece_lengths, side="right")
    return prefix[n_fit] - n_fit * (piece_lengths - 1)

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", required=True, type=Path)
    ap.add_argument("--workspace", required=True, type=Path, help="workspace BED(.gz), e.g. from make_workspace_bed.py")
    ap.add_argument("--chrom-sizes", type=Path, default=HERE / "data" / "mm10.chrom.sizes")
    ap.add_argument("--n-perm", type=int, default=10_000, help="permutations to assume when estimating failures")
    ap.add_argument("--rounds", type=int, default=10_000, help="rejection rounds of the original sampler")
    ap.add_argument("--top", type=int, default=8, help="how many hardest pieces to list")
    a = ap.parse_args(argv)

    genome = Genome.from_file(a.chrom_sizes)
    ws = IntervalSet.from_bed(genome, a.workspace)
    wl = ws.ends - ws.starts
    print(f"workspace: {len(ws):,} intervals, {ws.bp:,} bp ({ws.bp / genome.whole().bp:.1%} of the genome)")
    q = np.percentile(wl, [50, 90, 99, 99.9])
    print(f"  interval length: mean {wl.mean():,.0f} bp | median {q[0]:,.0f} | 90th pct {q[1]:,.0f} | 99th {q[2]:,.0f} | "
          f"99.9th {q[3]:,.0f} | longest {wl.max():,}")
    print(f"  intervals >= 10 kb: {(wl >= 10_000).sum():,} | >= 50 kb: {(wl >= 50_000).sum():,} | >= 100 kb: {(wl >= 100_000).sum():,}")

    for name, peaks in load_peaksets(a.repo, genome).items():
        pieces = peaks.intersect(ws)
        pl = pieces.ends - pieces.starts
        if len(pieces) == 0:
            print(f"\n{name}: no piece survives clipping to this workspace")
            continue
        valid = valid_starts(ws, pl)
        p = valid / ws.bp
        p_fail = np.exp(-p * a.rounds)
        expected_unplaced = float(p_fail.sum() * a.n_perm)
        order = np.argsort(valid)[: a.top]
        print(f"\n{name}: {len(peaks):,} peaks -> {len(pieces):,} pieces after clipping to the workspace "
              f"(longest piece {pl.max():,} bp; {int((pl >= 5_000).sum()):,} pieces >= 5 kb)")
        print(f"  hardest pieces (fewest valid start positions):")
        print(f"    {'length_bp':>10} {'valid_starts':>13} {'chance/draw':>12} {'P(unplaced after ' + str(a.rounds) + ' draws)':>30}")
        for i in order:
            print(f"    {pl[i]:>10,} {valid[i]:>13,} {p[i]:>12.1e} {p_fail[i]:>30.3g}")
        verdict = ("the ORIGINAL rejection-only sampler would fail" if expected_unplaced >= 1
                   else "the original rejection-only sampler would probably have coped")
        print(f"  expected cells still unplaced after {a.rounds:,} rounds over {a.n_perm:,} permutations: "
              f"{expected_unplaced:,.1f}  -> {verdict}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
