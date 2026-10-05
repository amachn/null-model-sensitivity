#!/usr/bin/env python3
import csv, io, re, sys
import numpy as np
import pandas as pd

N_QUERY = 14076
HEAD_KEYS = {"annotation", "epigenomicmark"}

def tables_in(path):
    if path.lower().endswith((".csv", ".txt", ".tsv")):
        with open(path, newline="", encoding="utf-8-sig") as fh:
            text = fh.read()
        delim = "\t" if text.count("\t") > text.count(",") else ","
        rows = list(csv.reader(io.StringIO(text), delimiter=delim))
        width = max(len(r) for r in rows)
        raw = pd.DataFrame([r + [""] * (width - len(r)) for r in rows]).replace("", np.nan)
        sheets = {"(file)": raw}
    else:
        sheets = pd.read_excel(path, sheet_name=None, header=None, dtype=str)
    for name, raw in sheets.items():
        raw = raw.dropna(how="all")
        for i, row in raw.iterrows():
            first = str(row.iloc[0]).strip().lower()
            if first in HEAD_KEYS:
                tab = raw.loc[i:].copy()
                tab.columns = [str(c).strip() for c in tab.iloc[0]]
                tab = tab.iloc[1:]
                tab = tab[tab.iloc[:, 0].notna() & ~tab.iloc[:, 0].astype(str).str.lower().isin(HEAD_KEYS)]
                title = str(raw.iloc[0, 0])[:90] if str(raw.iloc[0, 0]).lower() not in HEAD_KEYS else name
                yield name, title, tab
                break

for path in sys.argv[1:]:
    for sheet, title, tab in tables_in(path):
        sh = [c for c in tab.columns if re.match(r"^Shuffled_\d+(_Count)?$", c)]
        if not sh:
            continue
        m = tab[sh].apply(pd.to_numeric, errors="coerce").to_numpy()
        obs_col = tab.columns[1]
        out = pd.DataFrame({
            "annotation": tab.iloc[:, 0].str.replace("hg38_", "", regex=False).values,
            "observed": pd.to_numeric(tab[obs_col], errors="coerce").values,
            "n_shuffles": np.sum(~np.isnan(m), axis=1),
            "median": np.nanmedian(m, axis=1),
            "sd": np.round(np.nanstd(m, axis=1, ddof=1), 1),
        })
        p = out["median"] / N_QUERY
        out["binomial_sd"] = np.round(np.sqrt(N_QUERY * p * (1 - p)), 1)
        out["sd_ratio"] = np.round(out["sd"] / out["binomial_sd"], 2)
        out["z_published_null"] = np.round((out["observed"] - out["median"]) / out["sd"], 1)
        out["z_binomial_null"] = np.round((out["observed"] - out["median"]) / out["binomial_sd"], 1)
        print(f"\n== {path} | sheet '{sheet}' | {title}")
        print(out.to_string(index=False))
print("\nsd_ratio ~1: shuffles behave like independent random placements; >1: extra variance (different randomization).")
