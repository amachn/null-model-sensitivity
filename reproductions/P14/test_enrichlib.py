import numpy as np
import pandas as pd
import pytest

from enrichlib import (Genome, IntervalSet, UniformWorkspaceSampler, bh_adjust,
                       run_enrichment)

G = Genome({"c1": 1000, "c2": 700})

def bitmap(iset: IntervalSet) -> np.ndarray:
    m = np.zeros(G.total + 5, dtype=bool)
    for s, e in zip(iset.starts, iset.ends):
        m[s:e] = True
    return m

def random_set(rng, n, maxlen, chroms=("c1", "c2")):
    rows = []
    for _ in range(n):
        c = rng.choice(chroms)
        s = int(rng.integers(0, G.sizes[c] - 1))
        e = min(G.sizes[c], s + int(rng.integers(1, maxlen)))
        rows.append((c, s, e))
    return IntervalSet.from_frame(G, pd.DataFrame(rows, columns=["chrom", "start", "end"]))

def test_merge_is_disjoint_sorted_and_conserves_coverage():
    rng = np.random.default_rng(1)
    for _ in range(20):
        a = random_set(rng, 40, 120)
        assert np.all(a.starts[1:] > a.ends[:-1])
        assert a.bp == bitmap(a).sum()

def test_chromosome_boundary_does_not_merge():
    df = pd.DataFrame({"chrom": ["c1", "c2"], "start": [900, 0], "end": [1000, 100]})
    assert len(IntervalSet.from_frame(G, df)) == 2

def test_intersect_subtract_match_bitmaps():
    rng = np.random.default_rng(2)
    for _ in range(30):
        a, b = random_set(rng, 30, 150), random_set(rng, 30, 150)
        assert np.array_equal(bitmap(a.intersect(b)), bitmap(a) & bitmap(b))
        assert np.array_equal(bitmap(a.subtract(b)) & bitmap(G.whole()), bitmap(a) & ~bitmap(b))

def test_overlap_bp_matches_bitmap():
    rng = np.random.default_rng(3)
    for _ in range(30):
        ann = random_set(rng, 25, 200)
        seg = random_set(rng, 25, 100)
        bm = bitmap(ann)
        expect = np.array([bm[s:e].sum() for s, e in zip(seg.starts, seg.ends)])
        assert np.array_equal(ann.overlap_bp(seg.starts, seg.ends), expect)

def test_sampler_falls_back_when_workspace_is_too_fragmented_for_rejection_sampling():
    rng = np.random.default_rng(11)
    n_frag = 20_000
    starts = np.sort(rng.integers(0, 9_000_000, size=n_frag))
    ends = starts + rng.integers(5, 40, size=n_frag)
    starts = np.append(starts, [5_000_000]); ends = np.append(ends, [5_004_000])
    ws = IntervalSet(G, starts, ends)
    lengths = np.array([3800] + [10] * 40)
    out = UniformWorkspaceSampler(ws).place(np.random.default_rng(3), lengths, n_perm=500)
    ends_out = out + lengths
    for col, L in enumerate(lengths):
        assert (ws.overlap_bp(out[:, [col]], ends_out[:, [col]]).ravel() == L).all()

def test_exact_sampler_raises_a_clear_error_when_truly_impossible():
    ws = IntervalSet(G, [0], [100])
    with pytest.raises(RuntimeError, match="no single workspace interval"):
        UniformWorkspaceSampler(ws)._place_exact(np.random.default_rng(0), np.array([500] * 10))

def test_exact_sampler_is_uniform_over_every_valid_start():
    ws = IntervalSet(G, [0, 300, 700, 1200, 1500], [120, 330, 740, 1500, 1520])
    s = UniformWorkspaceSampler(ws)
    for L in (10, 30, 45, 120, 300):
        valid = np.concatenate([np.arange(a, b - L + 1) for a, b in zip(ws.starts, ws.ends) if b - a >= L])
        n = 40_000 * len(valid) if len(valid) < 10 else 200_000
        got = s._place_exact(np.random.default_rng(L), np.full(n, L))
        assert np.isin(got, valid).all(), f"L={L}: placed a segment somewhere it does not fit"
        counts = np.bincount(np.searchsorted(valid, got), minlength=len(valid))
        assert counts.min() > 0, f"L={L}: some valid start positions are never chosen"
        expect = n / len(valid)
        assert np.abs(counts - expect).max() < 6 * np.sqrt(expect), f"L={L}: not uniform"

def test_exact_and_rejection_stages_agree_in_distribution():
    ws = IntervalSet(G, [0, 300, 700], [200, 330, 900])
    lengths = np.array([20, 60, 150])
    a = UniformWorkspaceSampler(ws)
    b = UniformWorkspaceSampler(ws)
    b.REJECTION_ROUNDS = 0
    xa = a.place(np.random.default_rng(1), lengths, 60_000)
    xb = b.place(np.random.default_rng(2), lengths, 60_000)
    for col in range(3):
        vals = np.union1d(xa[:, col], xb[:, col])
        ca = np.bincount(np.searchsorted(vals, xa[:, col]), minlength=len(vals))
        cb = np.bincount(np.searchsorted(vals, xb[:, col]), minlength=len(vals))
        keep = (ca + cb) >= 10
        chi2 = ((ca[keep] - cb[keep]) ** 2 / (ca[keep] + cb[keep])).sum()
        dof = keep.sum() - 1
        assert chi2 < dof + 5 * np.sqrt(2 * dof), f"column {col}: chi2={chi2:.0f} on {dof} dof"

def test_ordinary_workspace_never_leaves_the_rejection_stage_but_fragmented_one_switches_early():
    lengths = np.array([100] * 200 + [200] * 100)
    ordinary = UniformWorkspaceSampler(IntervalSet(G, [0, 500], [450, 1000]))
    ordinary.place(np.random.default_rng(0), lengths, 200)
    assert ordinary.last_exact_cells == 0

    rng = np.random.default_rng(5)
    starts = np.sort(rng.integers(0, 900_000, size=40_000))
    frag = IntervalSet(G, np.append(starts, [500_000]), np.append(starts + rng.integers(5, 40, size=40_000), [504_000]))
    s = UniformWorkspaceSampler(frag)
    lengths = np.array([3000] * 40 + [10] * 200)
    out = s.place(np.random.default_rng(1), lengths, 600)
    assert s.last_exact_cells > 0 and s.last_rounds < 20
    assert (frag.overlap_bp(out, out + lengths) == lengths).all()

def test_sampler_preserves_length_and_stays_inside_workspace():
    ws = G.whole().subtract(IntervalSet(G, [300, G.offsets[1] + 100], [400, G.offsets[1] + 200]))
    lengths = np.array([5, 20, 60, 60, 99, 100])
    starts = UniformWorkspaceSampler(ws).place(np.random.default_rng(0), lengths, 400)
    ends = starts + lengths
    bm = bitmap(ws)
    for r in range(starts.shape[0]):
        for s, e in zip(starts[r], ends[r]):
            assert bm[s:e].all()

def test_sampler_is_uniform_over_valid_starts():
    ws = IntervalSet(G, [0, 500], [200, 600])
    starts = UniformWorkspaceSampler(ws).place(np.random.default_rng(0), np.array([50]), 60_000).ravel()
    in_first = (starts < 200).mean()
    assert in_first == pytest.approx(151 / (151 + 51), abs=0.01)

def test_expected_overlap_matches_analytic_value():
    ann = IntervalSet(G, [0, 500], [300, 800])
    ws = G.whole()
    seg = IntervalSet(G, [10, 900, 950], [40, 930, 990])
    res = run_enrichment(seg, {"a": ann}, ws, n_perm=4000, seed=1).iloc[0]
    lengths = np.array([30, 30, 40])
    exp = 0.0
    bm = bitmap(ann)
    for L in lengths:
        vals = []
        for off, size in [(G.offsets[0], 1000), (G.offsets[1], 700)]:
            for s in range(off, off + size - L + 1):
                vals.append(bm[s:s + L].sum())
        exp += np.mean(vals)
    assert res["expected"] == pytest.approx(exp, rel=0.05)

def test_pvalue_direction_and_floor():
    ann = IntervalSet(G, [100], [220])
    ws = G.whole()
    on = IntervalSet(G, [100, 130, 160, 190], [125, 155, 185, 215])
    res_hi = run_enrichment(on, {"a": ann}, ws, n_perm=999, seed=0).iloc[0]
    assert res_hi["l2fold"] > 1 and res_hi["pvalue"] == pytest.approx(1 / 1000)
    off = IntervalSet(G, [500, 550, 600], [540, 590, 640])
    big = IntervalSet(G, [0], [700])
    res_lo = run_enrichment(off, {"b": big.intersect(G.whole().subtract(IntervalSet(G, [480], [660])))},
                            ws, n_perm=999, seed=0).iloc[0]
    assert res_lo["observed"] == 0 and res_lo["l2fold"] < 0

def test_segment_counter_and_clipping_to_workspace():
    ann = IntervalSet(G, [100], [200])
    seg = IntervalSet(G, [90, 150, 500], [110, 160, 520])
    ws = IntervalSet(G, [0, 400], [180, 1000])
    res = run_enrichment(seg, {"a": ann}, ws, n_perm=200, counter="segment-overlap", seed=0).iloc[0]
    assert res["observed"] == 2
    res_bp = run_enrichment(seg, {"a": ann}, ws, n_perm=200, seed=0).iloc[0]
    assert res_bp["observed"] == 20

def test_bh_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.005])
    assert np.allclose(bh_adjust(p), [0.02, 0.04, 0.04, 0.02])

def test_reproducible_with_seed():
    ann = IntervalSet(G, [100], [300])
    seg = IntervalSet(G, [10, 400], [50, 460])
    a = run_enrichment(seg, {"a": ann}, G.whole(), n_perm=300, seed=7)
    b = run_enrichment(seg, {"a": ann}, G.whole(), n_perm=300, seed=7)
    assert a.equals(b)
