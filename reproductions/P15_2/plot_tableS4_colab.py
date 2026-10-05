import io, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from google.colab import files
    up = files.upload()
    it = pd.read_csv(io.BytesIO(next(iter(up.values()))), sep="\t")
    IN_COLAB = True
except ImportError:
    it = pd.read_csv(os.environ.get("T4_ITER", "results_T4/iterations.tsv"), sep="\t")
    IN_COLAB = False

PUBLISHED = {"H3K27ac": 1.86, "H3K4me3": 1.83}
METRICS = {"bp_overlap": "Overlapping base pairs",
           "n_atac_overlapping": "ATAC peaks overlapping (≥1 bp)",
           "n_pairs": "ATAC–ChIP overlap pairs"}
COL = {"null": "#D63F86", "obs": "#76A33B"}
marks = sorted(it.mark.unique())
obs = it[it.iter == 0].set_index("mark")
null = it[it.iter != 0]
N = null.iter.nunique()
saved = []

def stats(mark, metric):
    o = obs.loc[mark, metric]; x = null[null.mark == mark][metric].to_numpy()
    med = np.median(x)
    return dict(obs=o, med=med, fold=o / med, p_enr=(x >= o).mean(), p_dep=(x <= o).mean(),
                lo=np.percentile(x, 2.5), hi=np.percentile(x, 97.5), x=x)

def ptxt(p):
    return f"p < {1 / N:.3g}" if p == 0 else f"p = {p:.3g}"

def save(fig, name):
    fig.savefig(name, dpi=200, bbox_inches="tight"); saved.append(name); plt.show()

fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
for ax, (m, lab) in zip(axes, METRICS.items()):
    x = np.arange(len(marks)); w = 0.38
    S = [stats(k, m) for k in marks]
    scale = 1e6 if m == "bp_overlap" else 1e3
    ax.bar(x - w / 2, [s["med"] / scale for s in S], w, color=COL["null"], label="Random region median",
           yerr=[[(s["med"] - s["lo"]) / scale for s in S], [(s["hi"] - s["med"]) / scale for s in S]], capsize=3)
    ax.bar(x + w / 2, [s["obs"] / scale for s in S], w, color=COL["obs"], label="ATAC peak regions")
    for i, s in enumerate(S):
        ax.text(x[i], max(s["obs"], s["hi"]) / scale * 1.03, f"{s['fold']:.2f}×", ha="center", fontsize=10,
                fontweight="bold", color="black" if s["fold"] >= 1 else "firebrick")
    ax.set_xticks(x); ax.set_xticklabels(marks); ax.set_title(lab, fontsize=11)
    ax.set_ylabel("Mb" if m == "bp_overlap" else "thousands")
    ax.set_ylim(0, max(max(s["obs"], s["hi"]) for s in S) / scale * 1.18); ax.grid(axis="y", alpha=0.3)
axes[0].legend(frameon=False, loc="upper right", fontsize=8)
fig.suptitle(f"Liver ATAC peaks vs GSE128072 histone marks: three overlap definitions "
             f"({N} size- and chromosome-matched shuffles; error bars = 95% of shuffles)", fontweight="bold")
fig.tight_layout(); save(fig, "TableS4_paper_style_bars.png")

fig, axes = plt.subplots(len(marks), 3, figsize=(15, 3.6 * len(marks)))
axes = np.atleast_2d(axes)
for r, k in enumerate(marks):
    for c, (m, lab) in enumerate(METRICS.items()):
        ax = axes[r, c]; s = stats(k, m)
        scale = 1e6 if m == "bp_overlap" else 1e3
        ax.hist(s["x"] / scale, bins=40, color=COL["null"], alpha=0.8, label=f"{N} shuffles")
        ax.axvline(s["obs"] / scale, color=COL["obs"], lw=3, label="observed")
        ax.axvline(s["med"] / scale, color="k", lw=1, ls="--", label="shuffle median")
        direction = "enriched" if s["fold"] > 1 else "depleted"
        p = s["p_enr"] if s["fold"] > 1 else s["p_dep"]
        txt = f"{s['fold']:.2f}× ({direction})\n{ptxt(p)}"
        if m == "bp_overlap" and k in PUBLISHED:
            txt += f"\npublished: {PUBLISHED[k]}×"
        ax.text(0.02, 0.95, txt, transform=ax.transAxes, va="top", fontsize=9,
                bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.9))
        ax.set_title(f"{k} — {lab}", fontsize=10)
        ax.set_xlabel("Mb" if m == "bp_overlap" else "thousands")
        lo = min(s["x"].min(), s["obs"]) / scale; hi = max(s["x"].max(), s["obs"]) / scale
        pad = (hi - lo) * 0.08; ax.set_xlim(lo - pad, hi + pad)
axes[0, 0].legend(frameon=False, fontsize=8, loc="center left")
fig.suptitle("Observed overlap vs the shuffle null: same data, different overlap definitions",
             fontweight="bold", y=1.0)
fig.tight_layout(); save(fig, "TableS4_null_distributions.png")

fig, ax = plt.subplots(figsize=(7, 3.2))
ypos = np.arange(len(marks))
style = {"bp_overlap": ("#1f77b4", "o"), "n_atac_overlapping": ("#ff7f0e", "s"), "n_pairs": ("#2ca02c", "^")}
for m, lab in METRICS.items():
    ax.scatter([stats(k, m)["fold"] for k in marks], ypos, s=70, color=style[m][0], marker=style[m][1], label=lab, zorder=3)
ax.scatter([PUBLISHED.get(k, np.nan) for k in marks], ypos, s=160, facecolors="none", edgecolors="k", linewidths=2,
           marker="D", label="Published (Table S4)", zorder=4)
ax.axvline(1, color="grey", lw=1)
ax.set_yticks(ypos); ax.set_yticklabels(marks); ax.set_xlabel("Fold enrichment (observed / median of shuffles)")
ax.set_xlim(0.8, 2.05); ax.grid(axis="x", alpha=0.3)
ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
ax.set_title("Which overlap definition reproduces the published folds?", fontweight="bold")
fig.tight_layout(); save(fig, "TableS4_fold_by_metric.png")

rows = []
for k in marks:
    for m in METRICS:
        s = stats(k, m)
        rows.append(dict(mark=k, metric=m, observed=s["obs"], null_median=s["med"], fold=round(s["fold"], 3),
                         p_enrichment=s["p_enr"], p_depletion=s["p_dep"],
                         published_fold=PUBLISHED.get(k) if m == "bp_overlap" else None))
summary = pd.DataFrame(rows); print(summary.to_string(index=False))
summary.to_csv("TableS4_summary.csv", index=False); saved.append("TableS4_summary.csv")
if IN_COLAB:
    import shutil
    os.makedirs("figs_T4", exist_ok=True)
    for f in saved: shutil.copy(f, "figs_T4")
    shutil.make_archive("TableS4_figures", "zip", "figs_T4"); files.download("TableS4_figures.zip")
else:
    print("Saved:", ", ".join(saved))
