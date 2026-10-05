from __future__ import annotations

import gzip
from pathlib import Path
from typing import Callable, Dict, Optional, Sequence

import numpy as np
import pandas as pd

GAP = 1_000

class Genome:

    def __init__(self, sizes: Dict[str, int]):
        self.chroms = list(sizes)
        self.sizes = {c: int(v) for c, v in sizes.items()}
        lengths = np.array([self.sizes[c] for c in self.chroms], dtype=np.int64)
        self.offsets = np.concatenate([[0], np.cumsum(lengths + GAP)[:-1]]).astype(np.int64)
        self._off = dict(zip(self.chroms, self.offsets.tolist()))
        self.total = int(self.offsets[-1] + lengths[-1])

    @classmethod
    def from_file(cls, path: str | Path, keep: Optional[Sequence[str]] = None) -> "Genome":
        sizes = {}
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                c, n = line.split()[:2]
                if keep is None or c in keep:
                    sizes[c] = int(n)
        return cls(sizes)

    def whole(self) -> "IntervalSet":
        starts = self.offsets.copy()
        ends = starts + np.array([self.sizes[c] for c in self.chroms], dtype=np.int64)
        return IntervalSet(self, starts, ends, merged=True)

    def to_global(self, chrom: Sequence[str], pos: np.ndarray) -> np.ndarray:
        off = pd.Series(chrom).map(self._off)
        out = off.to_numpy(dtype="float64")
        ok = ~np.isnan(out)
        res = np.full(len(out), -1, dtype=np.int64)
        res[ok] = out[ok].astype(np.int64) + np.asarray(pos, dtype=np.int64)[ok]
        return res

    def to_local(self, gpos: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        idx = np.searchsorted(self.offsets, gpos, side="right") - 1
        chroms = np.array(self.chroms, dtype=object)[idx]
        return chroms, gpos - self.offsets[idx]

def _merge(starts: np.ndarray, ends: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    keep = ends > starts
    starts, ends = starts[keep], ends[keep]
    if len(starts) == 0:
        return starts.astype(np.int64), ends.astype(np.int64)
    order = np.argsort(starts, kind="stable")
    s, e = starts[order], ends[order]
    cummax = np.maximum.accumulate(e)
    new = np.empty(len(s), dtype=bool)
    new[0] = True
    new[1:] = s[1:] > cummax[:-1]
    idx = np.flatnonzero(new)
    ms = s[idx]
    me = np.append(cummax[idx[1:] - 1], cummax[-1])
    return ms.astype(np.int64), me.astype(np.int64)

class IntervalSet:

    def __init__(self, genome: Genome, starts, ends, merged: bool = False):
        starts = np.asarray(starts, dtype=np.int64)
        ends = np.asarray(ends, dtype=np.int64)
        if not merged:
            starts, ends = _merge(starts, ends)
        self.genome, self.starts, self.ends = genome, starts, ends
        self._cum = np.concatenate([[0], np.cumsum(ends - starts)]).astype(np.int64)
        self.dropped = 0

    @classmethod
    def from_frame(cls, genome: Genome, df: pd.DataFrame, chrom="chrom", start="start", end="end",
                   start_offset: int = 0, end_offset: int = 0) -> "IntervalSet":
        g_s = genome.to_global(df[chrom].astype(str).to_numpy(), df[start].to_numpy() + start_offset)
        g_e = genome.to_global(df[chrom].astype(str).to_numpy(), df[end].to_numpy() + end_offset)
        ok = (g_s >= 0) & (g_e >= 0)
        obj = cls(genome, g_s[ok], g_e[ok])
        obj.dropped = int((~ok).sum())
        return obj

    @classmethod
    def from_bed(cls, genome: Genome, path: str | Path, start_offset: int = 0, end_offset: int = 0) -> "IntervalSet":
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            df = pd.read_csv(fh, sep="\t", header=None, comment="#", usecols=[0, 1, 2],
                             names=["chrom", "start", "end"], dtype={0: str})
        return cls.from_frame(genome, df, start_offset=start_offset, end_offset=end_offset)

    def __len__(self) -> int:
        return len(self.starts)

    @property
    def bp(self) -> int:
        return int(self._cum[-1])

    def to_frame(self) -> pd.DataFrame:
        chrom, s = self.genome.to_local(self.starts)
        _, e = self.genome.to_local(self.ends - 1)
        return pd.DataFrame({"chrom": chrom, "start": s, "end": e + 1})

    def intersect(self, other: "IntervalSet") -> "IntervalSet":
        if len(self) == 0 or len(other) == 0:
            return IntervalSet(self.genome, [], [], merged=True)
        i0 = np.searchsorted(other.ends, self.starts, side="right")
        i1 = np.searchsorted(other.starts, self.ends, side="left")
        cnt = np.maximum(i1 - i0, 0)
        tot = int(cnt.sum())
        if tot == 0:
            return IntervalSet(self.genome, [], [], merged=True)
        a_idx = np.repeat(np.arange(len(self)), cnt)
        first = np.cumsum(cnt) - cnt
        b_idx = np.arange(tot) - np.repeat(first, cnt) + np.repeat(i0, cnt)
        s = np.maximum(self.starts[a_idx], other.starts[b_idx])
        e = np.minimum(self.ends[a_idx], other.ends[b_idx])
        keep = e > s
        return IntervalSet(self.genome, s[keep], e[keep], merged=True)

    def complement(self) -> "IntervalSet":
        s = np.concatenate([[0], self.ends])
        e = np.concatenate([self.starts, [self.genome.total]])
        keep = e > s
        return IntervalSet(self.genome, s[keep], e[keep], merged=True)

    def subtract(self, other: "IntervalSet") -> "IntervalSet":
        return self.intersect(other.complement())

    def cover(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=np.int64)
        k = np.searchsorted(self.starts, x, side="right")
        out = self._cum[k].copy()
        has = k > 0
        out[has] -= np.maximum(self.ends[k[has] - 1] - x[has], 0)
        return out

    def overlap_bp(self, seg_starts: np.ndarray, seg_ends: np.ndarray) -> np.ndarray:
        shp = np.shape(seg_starts)
        s, e = np.ravel(seg_starts), np.ravel(seg_ends)
        return (self.cover(e) - self.cover(s)).reshape(shp)

def _count_nucleotide(overlap: np.ndarray) -> np.ndarray:
    return overlap.sum(axis=-1)

def _count_segment(overlap: np.ndarray) -> np.ndarray:
    return (overlap > 0).sum(axis=-1)

COUNTERS: Dict[str, Callable[[np.ndarray], np.ndarray]] = {
    "nucleotide-overlap": _count_nucleotide,
    "segment-overlap": _count_segment,
}

class UniformWorkspaceSampler:

    name = "uniform-workspace"
    REJECTION_ROUNDS = 100
    SWITCH_ACCEPT_RATE = 0.25
    MIN_CELLS_FOR_SWITCH = 1000

    def __init__(self, workspace: IntervalSet):
        self.ws = workspace
        lens = workspace.ends - workspace.starts
        self._cs = np.cumsum(lens) - lens
        self._total = int(lens.sum())
        self._maxlen = int(lens.max()) if len(lens) else 0
        self._exact_ready = False
        self.last_rounds = 0
        self.last_exact_cells = 0

    def _prepare_exact(self) -> None:
        if self._exact_ready:
            return
        lens = self.ws.ends - self.ws.starts
        order = np.argsort(-lens, kind="stable")
        self._ord_starts = self.ws.starts[order]
        self._neg_lens = -lens[order]
        self._prefix = np.concatenate([[0], np.cumsum(lens[order])]).astype(np.int64)
        self._exact_ready = True

    def _place_exact(self, rng: np.random.Generator, seg_len: np.ndarray) -> np.ndarray:
        self._prepare_exact()
        seg_len = np.asarray(seg_len, dtype=np.int64)
        n_fit = np.searchsorted(self._neg_lens, -seg_len, side="right")
        if (n_fit == 0).any():
            bad = int(seg_len[n_fit == 0].max())
            raise RuntimeError(f"no single workspace interval is at least {bad} bp long, so a segment of that "
                               f"length cannot be placed inside it (longest interval available: {self._maxlen} bp)")
        shrink = seg_len - 1
        total = self._prefix[n_fit] - n_fit * shrink
        u = rng.integers(0, total)
        lo, hi = np.zeros_like(n_fit), n_fit - 1
        for _ in range(int(np.ceil(np.log2(max(len(self._neg_lens), 2)))) + 1):
            mid = (lo + hi + 1) >> 1
            ok = (self._prefix[mid] - mid * shrink) <= u
            lo = np.where(ok, mid, lo)
            hi = np.where(ok, hi, mid - 1)
        return self._ord_starts[lo] + (u - (self._prefix[lo] - lo * shrink))

    def place(self, rng: np.random.Generator, lengths: np.ndarray, n_perm: int) -> np.ndarray:
        lengths = np.asarray(lengths, dtype=np.int64)
        n = len(lengths)
        if n and lengths.max() > self._maxlen:
            raise ValueError(f"longest segment ({lengths.max()} bp) exceeds the longest workspace interval ({self._maxlen} bp)")
        L = np.broadcast_to(lengths, (n_perm, n))
        starts = np.empty((n_perm, n), dtype=np.int64)
        todo = np.ones((n_perm, n), dtype=bool)
        flat_starts, flat_todo, flat_L = starts.reshape(-1), todo.reshape(-1), np.ascontiguousarray(L).reshape(-1)
        self.last_rounds, self.last_exact_cells = 0, 0
        for _ in range(self.REJECTION_ROUNDS):
            idx = np.flatnonzero(flat_todo)
            if idx.size == 0:
                return starts
            u = rng.integers(0, self._total, size=idx.size)
            j = np.searchsorted(self._cs, u, side="right") - 1
            pos = self.ws.starts[j] + (u - self._cs[j])
            ok = pos + flat_L[idx] <= self.ws.ends[j]
            flat_starts[idx[ok]] = pos[ok]
            flat_todo[idx[ok]] = False
            self.last_rounds += 1
            if idx.size >= self.MIN_CELLS_FOR_SWITCH and ok.sum() < self.SWITCH_ACCEPT_RATE * idx.size:
                break
        idx = np.flatnonzero(flat_todo)
        if idx.size:
            self.last_exact_cells = int(idx.size)
            flat_starts[idx] = self._place_exact(rng, flat_L[idx])
        return starts

def bh_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out

def run_enrichment(
    segments: IntervalSet,
    annotations: Dict[str, IntervalSet],
    workspace: IntervalSet,
    n_perm: int = 10_000,
    counter: str = "nucleotide-overlap",
    sampler=None,
    seed: int = 0,
    batch_elements: int = 3_000_000,
    progress: Optional[Callable[[int, int], None]] = None,
) -> pd.DataFrame:
    if counter not in COUNTERS:
        raise ValueError(f"unknown counter {counter!r}; choose from {list(COUNTERS)}")
    count = COUNTERS[counter]

    seg = segments.intersect(workspace)
    ann = {k: a.intersect(workspace) for k, a in annotations.items()}
    lengths = seg.ends - seg.starts
    n = len(seg)
    if n == 0:
        raise ValueError("no segments remain inside the workspace")

    observed = {k: int(count(a.overlap_bp(seg.starts, seg.ends))) for k, a in ann.items()}

    sampler = sampler or UniformWorkspaceSampler(workspace)
    rng = np.random.default_rng(seed)
    sims = {k: np.empty(n_perm, dtype=np.int64) for k in ann}
    batch = max(1, min(n_perm, batch_elements // max(n, 1)))
    done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        starts = sampler.place(rng, lengths, b)
        ends = starts + lengths
        for k, a in ann.items():
            sims[k][done:done + b] = count(a.overlap_bp(starts, ends))
        done += b
        if progress:
            progress(done, n_perm)

    rows = []
    for k, a in ann.items():
        s, obs = sims[k], observed[k]
        exp = float(s.mean())
        hits = int((s >= obs).sum()) if obs >= exp else int((s <= obs).sum())
        p = min(1.0, (hits + 1) / (n_perm + 1))
        fold = (obs + 1.0) / (exp + 1.0)
        rows.append(dict(
            annotation=k, observed=obs, expected=exp,
            CI95low=float(np.percentile(s, 2.5)), CI95high=float(np.percentile(s, 97.5)),
            stddev=float(s.std(ddof=0)), fold=fold, l2fold=float(np.log2(fold)),
            pvalue=p, track_nsegments=n, track_size=seg.bp,
            annotation_size=a.bp, annotation_nsegments=len(a),
            n_perm=n_perm, counter=counter, null=getattr(sampler, "name", type(sampler).__name__),
        ))
    out = pd.DataFrame(rows)
    out["qvalue"] = bh_adjust(out["pvalue"].to_numpy())
    return out
