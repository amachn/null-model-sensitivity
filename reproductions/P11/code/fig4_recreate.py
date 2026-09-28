"""Recreate Fig. 4 of Iotchkova et al. 2019 (Nat Genet 51:343-353, GARFIELD).

Input : data/SuppTable5.xlsx  (sheets Hotspots = ST5a, H3K27ac = ST5b, H3K4me3 = ST5c;
        columns Annotation, Tissue, Trait, GOSHIFTER, GREGOR, FGWAS, GARFIELD, LDSC
        holding -log10 P per method)
Logic : follows the authors' Figure4_ST4_ST5_SF3.R (manuscript_custom_code.tar.gz).
Output: results/Fig4_recreated.png/pdf, results/Fig4{a,b,c,d}_*.png, results/*.csv

Deviation: ST5 stores only -log10 P, not the raw enrichment statistic the authors used
for box size in panels b-d. Box size here is -log10 P scaled 0-1 per method and trait.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "SuppTable5.xlsx")
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

# order and colours from the authors' script
METHODS = ["GARFIELD", "GREGOR", "LDSC", "FGWAS", "GOSHIFTER"]
COL = {"GARFIELD": "#FF0000", "GREGOR": "#FFA500", "LDSC": "#EEDC82",
       "FGWAS": "#43CD80", "GOSHIFTER": "#4F94CD"}
BAR_ORDER = ["GREGOR", "GARFIELD", "LDSC", "FGWAS", "GOSHIFTER"]
BAR_COL = {"GREGOR": "orange", "GARFIELD": "red", "LDSC": "#EEDC82",
           "FGWAS": "lightgreen", "GOSHIFTER": "#4F94CD"}
TRAITS = ["HGT", "BMI", "HDL", "TG", "LDL", "TC", "MCH", "MCV", "MCHC", "RBC", "HGB",
          "PCV", "PLT", "MPV", "CD", "IBD", "UC", "FPI", "FG", "HbA1C", "SCZ"]
# significance thresholds used by the authors (P values)
THR = {"Hotspots": 0.000257, "H3K27ac": 0.00047, "H3K4me3": 0.00047}


def load(sheet):
    d = pd.read_excel(DATA, sheet_name=sheet, header=2, usecols=range(8))
    d = d.rename(columns={"GOSHIFTER": "GOSHIFTER"})
    d = d[d["Trait"].isin(TRAITS)].copy()
    d["Tissue"] = d["Tissue"].astype(str).str.replace("_", " ")
    for m in METHODS:
        d[m] = pd.to_numeric(d[m], errors="coerce")
    return d


def long(d, sheet):
    l = d.melt(id_vars=["Annotation", "Tissue", "Trait"], value_vars=METHODS,
               var_name="Method", value_name="nlogP")
    l["sig"] = l["nlogP"] > -np.log10(THR[sheet])
    return l


def fig4a(ax, l):
    s = l[l["sig"]].copy()
    s["Numbers"] = s.groupby(["Annotation", "Trait"])["Method"].transform("count")
    cnt = s.groupby(["Numbers", "Method"]).size().rename("All").reset_index()
    cnt["Norm"] = cnt.groupby("Numbers")["All"].transform("sum")
    cnt["Prop"] = cnt["All"] / cnt["Norm"]
    cnt["Height"] = cnt["Prop"] * cnt["Numbers"]   # authors plot Prop*Numbers
    w = 0.16
    for i, m in enumerate(BAR_ORDER):
        sub = cnt[cnt["Method"] == m].set_index("Numbers")
        xs = np.arange(1, 6)
        ys = [sub["Height"].get(n, 0) for n in xs]
        ax.bar(xs + (i - 2) * w, ys, w, color=BAR_COL[m], edgecolor="black",
               linewidth=0.5, label=m)
    ax.set_xlabel("Number of methods declaring enrichment")
    ax.set_ylabel("Proportion of enrichments\nattributable to method")
    ax.set_ylim(0, 1.05)
    ax.legend(title="Method", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return cnt


def grid(ax, l, sheet, title=None, fs=5):
    """Tissue x (trait, method) grid; box = best annotation per tissue."""
    l = l.copy()
    # scale -log10P to 0-1 per trait and method across all annotations
    g = l.groupby(["Trait", "Method"])["nlogP"]
    l["size"] = (l["nlogP"] - g.transform("min")) / (g.transform("max") - g.transform("min"))
    best = (l.sort_values("nlogP", ascending=False)
             .groupby(["Trait", "Method", "Tissue"], as_index=False).first())
    tissues = sorted(best["Tissue"].unique(), key=str.lower)
    ty = {t: i for i, t in enumerate(tissues)}
    nm = len(METHODS)
    for ti, tr in enumerate(TRAITS):
        for mi, m in enumerate(METHODS):
            x = ti * (nm + 1) + mi
            ax.add_patch(plt.Rectangle((x - .5, -.5), 1, len(tissues), color="#E5E5E5", lw=0))
    for _, r in best.iterrows():
        x = TRAITS.index(r["Trait"]) * (nm + 1) + METHODS.index(r["Method"])
        col = COL[r["Method"]] if r["sig"] else "#B3B3B3"
        area = 0.04 + 0.96 * (r["size"] if pd.notna(r["size"]) else 0)
        ax.scatter(x, ty[r["Tissue"]], s=fs * 8 * area, marker="s", color=col, lw=0)
    ax.set_xlim(-1, len(TRAITS) * (nm + 1))
    ax.set_ylim(len(tissues) - .5, -.5)
    ax.set_yticks(range(len(tissues)))
    ax.set_yticklabels(tissues, fontsize=fs + 1)
    ax.set_xticks([ti * (nm + 1) + (nm - 1) / 2 for ti in range(len(TRAITS))])
    ax.set_xticklabels(TRAITS, fontsize=7)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    ax.set_ylabel("Tissue")
    for sp in ax.spines.values():
        sp.set_visible(False)
    if title:
        ax.set_title(title, loc="left", fontsize=9, pad=14)
    return best


def summary(l, sheet):
    s = (l[l["sig"]].groupby(["Method", "Trait"]).size().unstack(fill_value=0)
         .reindex(index=METHODS, columns=TRAITS, fill_value=0))
    return s


def main():
    sheets = {"Hotspots": "DHS hotspots", "H3K27ac": "H3K27ac", "H3K4me3": "H3K4me3"}
    L = {s: long(load(s), s) for s in sheets}

    # ---- stand-alone panels
    fig, ax = plt.subplots(figsize=(6, 3.4))
    cnt = fig4a(ax, L["Hotspots"])
    fig.tight_layout(); fig.savefig(f"{OUT}/Fig4a_recreated.png", dpi=200); plt.close(fig)
    cnt.to_csv(f"{OUT}/Fig4a_counts.csv", index=False)

    for key, name, h in (("Hotspots", "b", 11), ("H3K27ac", "c", 4.5), ("H3K4me3", "d", 4.5)):
        fig, ax = plt.subplots(figsize=(15, h))
        best = grid(ax, L[key], key)
        fig.tight_layout(); fig.savefig(f"{OUT}/Fig4{name}_recreated.png", dpi=200); plt.close(fig)
        best.to_csv(f"{OUT}/Fig4{name}_best_per_tissue.csv", index=False)

    # ---- composite figure, layout like the paper
    fig = plt.figure(figsize=(15, 19))
    gs = fig.add_gridspec(4, 1, height_ratios=[2, 8, 3.3, 3.3], hspace=0.3)
    axa = fig.add_subplot(gs[0])
    fig4a(axa, L["Hotspots"]); axa.set_title("a", loc="left", fontweight="bold", x=-0.12)
    for i, (key, lab) in enumerate((("Hotspots", "b  DHS hotspots"), ("H3K27ac", "c  H3K27ac"),
                                    ("H3K4me3", "d  H3K4me3"))):
        ax = fig.add_subplot(gs[i + 1])
        grid(ax, L[key], key, title=lab)
    axa.text(1.01, -0.45, "Coloured box: significant\nGrey box: not significant\nBox size: -log10 P scaled\n0-1 per method and trait\n(panels b-d)",
             transform=axa.transAxes, fontsize=8, va="bottom")
    fig.suptitle("Recreated Fig. 4 (GARFIELD, Iotchkova et al. 2019) from Supplementary Table 5",
                 y=0.93, fontsize=11)
    fig.savefig(f"{OUT}/Fig4_recreated.png", dpi=170, bbox_inches="tight")
    fig.savefig(f"{OUT}/Fig4_recreated.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---- validation against numbers stated in the paper / authors' script
    s = summary(L["Hotspots"], "Hotspots")
    v = pd.DataFrame({"median": s.median(axis=1), "max": s.max(axis=1), "mean": s.mean(axis=1).round(2)})
    v["paper_median"] = pd.Series({"GARFIELD": 10, "GREGOR": 24, "LDSC": 5, "FGWAS": 5, "GOSHIFTER": 0})
    v["paper_max"] = pd.Series({"GARFIELD": 364, "GREGOR": 398, "LDSC": 144, "FGWAS": 327, "GOSHIFTER": 5})
    v.to_csv(f"{OUT}/validation_vs_paper.csv")
    s.to_csv(f"{OUT}/significant_DHS_per_trait_method.csv")
    print(v.to_string())
    print("\nFig 4a counts (bar height = Prop*Numbers):")
    print(cnt.pivot(index="Numbers", columns="Method", values="Height").round(2).to_string())


if __name__ == "__main__":
    main()
