#!/usr/bin/env python3
import argparse, csv, gzip, os, sys, urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import S2_fetch_encode_liver_marks as s2

VARIANTS = {
    "pseudorep_one":    lambda f: f["assembly"] == "GRCh38" and f["output_type"] == "pseudoreplicated peaks",
    "relaxed":          lambda f: f["assembly"] == "GRCh38" and f["output_type"] == "peaks",
    "union_grch38":     lambda f: f["assembly"] == "GRCh38",
    "union_all_builds": lambda f: True,
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="marks_manifest.tsv from S2")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    raw = os.path.join(a.out, "_raw")
    os.makedirs(raw, exist_ok=True)
    log = open(os.path.join(a.out, "variant_files.tsv"), "w")
    log.write("variant\tlabel\texperiment\tfile\tassembly\toutput_type\tn_lines\n")

    for r in csv.DictReader(open(a.manifest), delimiter="\t"):
        files = [f for f in s2.search(type="File", dataset=f"/experiments/{r['experiment']}/")
                 if f.get("file_format") == "bed"]
        files.sort(key=lambda f: f["accession"])
        local = {}
        for f in files:
            dest = os.path.join(raw, f"{f['accession']}.bed.gz")
            if not os.path.exists(dest):
                urllib.request.urlretrieve(s2.BASE + f["href"], dest)
            local[f["accession"]] = dest
        print(f"[info] {r['label']}: {len(files)} bed files "
              f"({sum(f['assembly'] == 'GRCh38' for f in files)} GRCh38, {sum(f['assembly'] == 'hg19' for f in files)} hg19)")

        for v, keep in VARIANTS.items():
            sel = [f for f in files if keep(f)]
            if v == "pseudorep_one":
                sel = [f for f in sel if f["accession"] == r["file"]] or sel[:1]
            if not sel:
                print(f"[warn] {v}/{r['label']}: no matching files"); continue
            vdir = os.path.join(a.out, v)
            os.makedirs(vdir, exist_ok=True)
            with gzip.open(os.path.join(vdir, f"{r['label']}.bed.gz"), "wt") as out:
                for f in sel:
                    n = 0
                    with gzip.open(local[f["accession"]], "rt") as fh:
                        for line in fh:
                            if line.startswith(("track", "browser", "#")):
                                continue
                            out.write(line); n += 1
                    log.write(f"{v}\t{r['label']}\t{r['experiment']}\t{f['accession']}\t{f['assembly']}\t{f['output_type']}\t{n}\n")
    log.close()
    print(f"[done] variants in {a.out}/ : {', '.join(VARIANTS)}  (file list: variant_files.tsv)")

if __name__ == "__main__":
    main()
