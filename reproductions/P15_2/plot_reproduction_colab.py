import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MAIN = """figure\tannotation\tpub_obs\tpub_null\tpub_p\tA_obs\tA_null\tA_p\tC1_obs\tC1_null\tC1_p\tC2_obs\tC2_null\tC2_p\tB_obs\tB_null\tB_p
S10\tenhancers_fantom\t415\t322\t0.000999\t599\t397\t0.002\t599\t397\t0.002\t415\t210\t0.002\t415\t224\t0.020
S10\tgenes_1to5kb\t2635\t2356\t0.000999\t3407\t2331\t0.002\t3407\t2333\t0.002\t2635\t1822\t0.002\t2635\t1921\t0.020
S10\tgenes_3UTRs\t288\t387\t0.000999\t519\t564\t0.044\t519\t563\t0.038\t288\t284\t0.819\t288\t304\t0.436
S10\tgenes_5UTRs\t188\t209\t0.063\t543\t302\t0.002\t543\t302\t0.002\t188\t83\t0.002\t188\t89\t0.020
S10\tgenes_cds\t403\t378\t0.098\t1061\t905\t0.002\t1061\t906\t0.002\t403\t197\t0.002\t403\t212\t0.020
S10\tgenes_exons\t1579\t1406\t0.000999\t2691\t1882\t0.002\t2691\t1880\t0.002\t1579\t932\t0.002\t1579\t985\t0.020
S10\tgenes_firstexons\t895\t713\t0.000999\t1537\t798\t0.002\t1537\t797\t0.002\t895\t398\t0.002\t895\t416\t0.020
S10\tgenes_intergenic\t3776\t3428\t0.000999\t4836\t5208\t0.002\t4836\t5208\t0.002\t3776\t4139\t0.002\t3776\t3712\t0.158
S10\tgenes_intronexonboundaries\t1639\t1596\t0.111\t2233\t1613\t0.002\t2233\t1612\t0.002\t1639\t1122\t0.002\t1639\t1186\t0.020
S10\tgenes_introns\t6420\t6921\t0.000999\t8075\t8258\t0.002\t8075\t8257\t0.004\t6420\t6523\t0.054\t6420\t6923\t0.020
S10\tgenes_promoters\t1469\t1280\t0.000999\t1890\t1088\t0.002\t1890\t1085\t0.002\t1469\t808\t0.002\t1469\t853\t0.020
S11\tH3K27ac\t4283\t4093\t0.057\t933\t597\t0.002\t933\t599\t0.002\t870\t445\t0.002\t870\t479\t0.020
S11\tH3K27me3\t619\t960\t0.001\t100\t104\t0.733\t100\t104\t0.661\t77\t66\t0.188\t77\t71\t0.495
S11\tH3K36me3\t2825\t4631.5\t0.001\t395\t726\t0.002\t395\t725\t0.002\t306\t505\t0.002\t306\t542\t0.020
S11\tH3K4me1\t7132\t6675.5\t0.01\t1201\t985\t0.002\t1201\t984\t0.002\t1085\t708\t0.002\t1085\t768\t0.020
S11\tH3K4me3\t3492\t3586\t0.244\t457\t306\t0.002\t457\t306\t0.002\t436\t222\t0.002\t436\t236\t0.020
S11\tH3K9ac\t1567\t1435\t0.032\t393\t243\t0.002\t393\t243\t0.002\t373\t178\t0.002\t373\t192\t0.020
S11\tH3K9me3\t1257\t1036\t0.004\t177\t196\t0.172\t177\t196\t0.172\t137\t137\t1.000\t137\t140\t0.772
"""
VARIANTS = """mark\tpub_obs\tpub_null\tpseudorep_obs\tpseudorep_null\trelaxed_obs\trelaxed_null\tunion_grch38_obs\tunion_grch38_null\tunion_all_builds_obs\tunion_all_builds_null
H3K27ac\t4283\t4093\t933\t595\t1081\t777.5\t1134\t791\t1804\t1418
H3K27me3\t619\t960\t100\t105\t555\t735\t564\t737\t1196\t1318
H3K36me3\t2825\t4631.5\t395\t735\t617\t1111.5\t644\t1130\t1656\t1974
H3K4me1\t7132\t6675.5\t1201\t992.5\t1438\t1309\t1483\t1337\t2419\t2272
H3K4me3\t3492\t3586\t457\t306\t959\t822.5\t977\t828\t1612\t1444
H3K9ac\t1567\t1435\t393\t240\t659\t559.5\t675\t560.5\t1159\t1025
H3K9me3\t1257\t1036\t177\t203\t617\t790.5\t640\t800.5\t1254\t1416
"""
df = pd.read_csv(io.StringIO(MAIN), sep="\t")
var = pd.read_csv(io.StringIO(VARIANTS), sep="\t")

CONDS = {"pub": "Published (Tables S6/S7)",
         "A": "A: annotatr null, 1 bp",
         "C1": "C1: paper-literal shuffle, 1 bp",
         "C2": "C2: paper shuffle, 200 bp",
         "B": "B: shuffle excl. blacklist, 200 bp"}
STYLE = {"pub": ("k", "D", 75), "A": ("#1f77b4", "o", 50), "C1": ("#17becf", "^", 50),
         "C2": ("#2ca02c", "v", 50), "B": ("#ff7f0e", "s", 50)}
REF = "C2"
COLORS = {"null": "#D63F86", "obs": "#76A33B"}
for c in CONDS:
    df[f"{c}_fold"] = df[f"{c}_obs"] / df[f"{c}_null"]
    df[f"{c}_dir"] = np.select([df[f"{c}_fold"] > 1, df[f"{c}_fold"] < 1], ["enriched", "depleted"], "none")
saved = []

def save(fig, name):
    fig.savefig(name, dpi=200, bbox_inches="tight")
    saved.append(name)
    plt.show()

def paper_style(fig_id, title):
    d = df[df.figure == fig_id].reset_index(drop=True)
    x = np.arange(len(d)); w = 0.4
    ymax = d[[f"{c}_{k}" for c in CONDS for k in ("obs", "null")]].max().max() * 1.08
    fig, axes = plt.subplots(1, len(CONDS), figsize=(5.2 * len(CONDS), 5.5), sharey=True)
    for ax, (c, label) in zip(axes, CONDS.items()):
        ax.bar(x - w / 2, d[f"{c}_null"], w, color=COLORS["null"], label="Random region median")
        ax.bar(x + w / 2, d[f"{c}_obs"], w, color=COLORS["obs"], label="caQTL peaks")
        if True:
            for i, p in enumerate(d[f"{c}_p"]):
                if p < 0.05:
                    ax.text(x[i], max(d[f"{c}_obs"][i], d[f"{c}_null"][i]) + ymax * 0.01, "*", ha="center", fontsize=13)
        ax.set_title(label, fontsize=11); ax.set_xticks(x)
        ax.set_xticklabels(d.annotation, rotation=55, ha="right", fontsize=9)
        ax.set_ylim(0, ymax); ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("Number of peaks")
    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle(f"{title}   (* empirical p < 0.05)", fontsize=13, fontweight="bold")
    fig.tight_layout()
    save(fig, f"{fig_id}_paper_style_comparison.png")

paper_style("S10", "Fig S10: caQTL peaks vs genomic annotations")
paper_style("S11", "Fig S11: caQTL peaks vs liver histone marks")

def fold_dotplot(fig_id, title):
    d = df[df.figure == fig_id].reset_index(drop=True)
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(9, 0.5 * len(d) + 2))
    for c, label in CONDS.items():
        col, m, s = STYLE[c]
        ax.scatter(np.log2(d[f"{c}_fold"]), y, color=col, marker=m, s=s, label=label, zorder=3)
    for i in range(len(d)):
        flips = [c for c in CONDS if c != "pub" and d[f"{c}_dir"][i] != d["pub_dir"][i]]
        if flips:
            ax.text(1.2, y[i], "flip: " + ",".join(flips), va="center",
                    fontsize=8, color="firebrick")
    ax.axvline(0, color="grey", lw=1)
    ax.set_yticks(y); ax.set_yticklabels(d.annotation)
    ax.set_xlabel("log2(observed / random median)   <- depleted | enriched ->")
    ax.set_xlim(-1.1, 1.6); ax.grid(axis="x", alpha=0.3)
    ax.set_title(title, fontweight="bold")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=8, ncol=3)
    fig.tight_layout()
    save(fig, f"{fig_id}_fold_comparison.png")

fold_dotplot("S10", "Fig S10: fold enrichment, published vs reproduced")
fold_dotplot("S11", "Fig S11: fold enrichment, published vs reproduced")

d = df.copy()
d["obs_pct"] = 100 * (d[f"{REF}_obs"] - d.pub_obs) / d.pub_obs
d["null_pct"] = 100 * (d[f"{REF}_null"] - d.pub_null) / d.pub_null
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw={"width_ratios": [11, 7]})
for ax, fid in zip(axes, ("S10", "S11")):
    s = d[d.figure == fid].reset_index(drop=True); x = np.arange(len(s)); w = 0.4
    ax.bar(x - w / 2, s.obs_pct, w, color=COLORS["obs"], label="observed count")
    ax.bar(x + w / 2, s.null_pct, w, color=COLORS["null"], label="random-set median")
    ax.axhline(0, color="k", lw=1); ax.axhspan(-5, 5, color="grey", alpha=0.15, label="within ±5%")
    ax.set_xticks(x); ax.set_xticklabels(s.annotation, rotation=55, ha="right", fontsize=9)
    ax.set_title(f"Fig {fid}"); ax.grid(axis="y", alpha=0.3)
    if (s.obs_pct.abs() < 0.5).all():
        ax.set_title(f"Fig {fid}\nobserved counts: exact match to published (0% difference)", color="#4d7a1f")
axes[0].set_ylabel(f"% difference from published ({REF})")
axes[0].legend(frameon=False)
fig.suptitle(f"How far off is the reproduction? ({CONDS[REF]})", fontweight="bold")
fig.tight_layout()
save(fig, "pct_difference_from_published.png")

inputs = ["pseudorep", "relaxed", "union_grch38", "union_all_builds"]
v = var.copy()
v["pub_fold"] = v.pub_obs / v.pub_null
for k in inputs:
    v[f"{k}_fold"] = v[f"{k}_obs"] / v[f"{k}_null"]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
y = np.arange(len(v))
ax1.scatter(np.log2(v.pub_fold), y, color="k", marker="D", s=80, label="Published", zorder=4)
cols = ["#1f77b4", "#2ca02c", "#9467bd", "#d62728"]
for k, c in zip(inputs, cols):
    ax1.scatter(np.log2(v[f"{k}_fold"]), y, color=c, s=45, label=k, zorder=3)
ax1.axvline(0, color="grey"); ax1.set_yticks(y); ax1.set_yticklabels(v.mark)
ax1.set_xlabel("log2 fold (observed / random median)"); ax1.grid(axis="x", alpha=0.3)
ax1.set_title("Fold enrichment by mark-file input"); ax1.legend(frameon=False, fontsize=8)
x = np.arange(len(v)); w = 0.16
ax2.bar(x - 2 * w, v.pub_obs, w, color="k", label="Published")
for j, (k, c) in enumerate(zip(inputs, cols)):
    ax2.bar(x + (j - 1) * w, v[f"{k}_obs"], w, color=c, label=k)
ax2.set_xticks(x); ax2.set_xticklabels(v.mark, rotation=45, ha="right")
ax2.set_ylabel("caQTL peaks overlapping mark"); ax2.set_title("Observed counts by mark-file input")
ax2.legend(frameon=False, fontsize=8); ax2.grid(axis="y", alpha=0.3)
fig.suptitle("Fig S11 input sensitivity (annotatr null, 1-bp overlap, 10 random sets)", fontweight="bold")
fig.tight_layout()
save(fig, "S11_input_variants.png")

DISP = """figure\tannotation\tpublished\treproduced
S10\tenhancers_fantom\t0.93\t0.94
S10\tgenes_1to5kb\t0.96\t0.98
S10\tgenes_3UTRs\t0.97\t0.99
S10\tgenes_5UTRs\t1.00\t0.99
S10\tgenes_cds\t0.94\t0.95
S10\tgenes_exons\t0.87\t0.99
S10\tgenes_firstexons\t0.90\t0.98
S10\tgenes_intergenic\t0.90\t0.94
S10\tgenes_intronexonboundaries\t0.90\t0.92
S10\tgenes_introns\t0.84\t0.86
S10\tgenes_promoters\t0.94\t1.00
S11\tH3K27ac\t2.25\t0.94
S11\tH3K27me3\t2.06\t1.01
S11\tH3K36me3\t3.43\t0.99
S11\tH3K4me1\t3.07\t0.98
S11\tH3K4me3\t2.46\t0.99
S11\tH3K9ac\t1.82\t0.97
S11\tH3K9me3\t2.43\t1.00
"""
dp = pd.read_csv(io.StringIO(DISP), sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw={"width_ratios": [11, 7]}, sharey=True)
for ax, fid, lab in zip(axes, ("S10", "S11"), ("Fig S10 / Table S6: genomic categories", "Fig S11 / Table S7: histone marks")):
    q = dp[dp.figure == fid].reset_index(drop=True); x = np.arange(len(q)); w = 0.38
    ax.bar(x - w / 2, q.published, w, color="#444444", label="Published null")
    ax.bar(x + w / 2, q.reproduced, w, color="#76A33B", label="This reproduction (C2 for S10, C1 for S11)")
    ax.axhline(1, color="firebrick", ls="--", lw=1.2, label="independent placement (ratio = 1)")
    ax.set_xticks(x); ax.set_xticklabels(q.annotation, rotation=55, ha="right", fontsize=9)
    ax.set_title(lab, fontsize=11); ax.grid(axis="y", alpha=0.3)
axes[0].set_ylabel("SD of random-set counts / binomial SD")
axes[0].legend(frameon=False, fontsize=8, loc="upper left")
fig.suptitle("Null spread: S6 matches independent, chromosome-matched placement; S7 is 1.8–3.4× wider",
             fontweight="bold")
fig.tight_layout(); save(fig, "null_spread_S6_vs_S7.png")

cols = ["figure", "annotation"] + [f"{c}_fold" for c in CONDS] + [f"{c}_dir" for c in CONDS]
summary = df[cols].round(2)
for c in CONDS:
    if c != "pub":
        summary[f"{c}_matches_pub"] = summary[f"{c}_dir"] == summary.pub_dir
print(summary.to_string(index=False))
print("\nDirection matches with published bars:")
for fid in ("S10", "S11"):
    s_ = summary[summary.figure == fid]
    print(f"  {fid}: " + ", ".join(f"{c} {int(s_[f'{c}_matches_pub'].sum())}/{len(s_)}" for c in CONDS if c != "pub"))
summary.to_csv("reproduction_summary.csv", index=False); saved.append("reproduction_summary.csv")
try:
    import shutil, os
    from google.colab import files
    os.makedirs("figs", exist_ok=True)
    for f in saved:
        shutil.copy(f, "figs")
    shutil.make_archive("reproduction_figures", "zip", "figs")
    files.download("reproduction_figures.zip")
except ImportError:
    print("Saved:", ", ".join(saved))
