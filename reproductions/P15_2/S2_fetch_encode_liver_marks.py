#!/usr/bin/env python3
import argparse, collections, json, os, sys, urllib.error, urllib.parse, urllib.request

BASE = "https://www.encodeproject.org"
MARKS = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3"]
TERMS = ["liver", "right lobe of liver"]
OUTPUT_PREF = ["replicated peaks", "pseudoreplicated peaks", "pseudo-replicated peaks", "stable peaks", "peaks"]

def get_json(path, **params):
    params.setdefault("format", "json")
    url = BASE + path + "?" + urllib.parse.urlencode(params, doseq=True)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)

def search(**params):
    try:
        return get_json("/search/", limit="all", frame="object", status="released", **params).get("@graph", [])
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return []
        raise

def target_label(t):
    return t.rstrip("/").split("/")[-1].replace("-human", "") if t else None

def list_experiments():
    exps = []
    for term in TERMS:
        exps += search(**{"type": "Experiment", "assay_title": "Histone ChIP-seq",
                          "biosample_ontology.term_name": term,
                          "replicates.library.biosample.donor.organism.scientific_name": "Homo sapiens"})
    roadmap = {e["accession"] for term in TERMS for e in search(**{
        "type": "Experiment", "assay_title": "Histone ChIP-seq", "biosample_ontology.term_name": term,
        "award.project": "Roadmap"})}
    rows = []
    for e in exps:
        refs = [s.rstrip("/").split("/")[-1] for s in e.get("related_series", []) if "reference-epigenomes" in s]
        rows.append(dict(accession=e["accession"], mark=target_label(e.get("target")),
                         group=refs[0] if refs else "biosample:" + e.get("biosample_summary", "?"),
                         project="Roadmap" if e["accession"] in roadmap else "ENCODE",
                         biosample=e.get("biosample_summary", "")))
    return rows

def rank_groups(rows):
    g = collections.defaultdict(list)
    for r in rows:
        g[r["group"]].append(r)
    out = []
    for name, rs in g.items():
        marks = {r["mark"] for r in rs}
        out.append(dict(group=name, n_s11_marks=len(marks & set(MARKS)), missing=",".join(sorted(set(MARKS) - marks)),
                        extra=len(marks - set(MARKS)), roadmap=any(r["project"] == "Roadmap" for r in rs),
                        biosample=rs[0]["biosample"], rows=rs))
    out.sort(key=lambda x: (-x["n_s11_marks"], not x["roadmap"], x["extra"]))
    return out

def pick_file(exp_acc):
    fs = search(**{"type": "File", "dataset": f"/experiments/{exp_acc}/", "assembly": "GRCh38", "file_format": "bed"})
    if not fs:
        return None
    def score(f):
        ot = f.get("output_type", "")
        return (OUTPUT_PREF.index(ot) if ot in OUTPUT_PREF else len(OUTPUT_PREF), not f.get("preferred_default", False))
    return sorted(fs, key=score)[0]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--group", help="reference epigenome accession or 'biosample:...' label to use")
    ap.add_argument("--list-only", action="store_true")
    ap.add_argument("--donor", help="only use experiments whose biosample contains this text, "
                    "e.g. 'female adult (25 years)'; groups such as ENCSR000FYQ mix several donors")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    rows = list_experiments()
    groups = rank_groups(rows)
    with open(os.path.join(a.out, "encode_liver_candidates.tsv"), "w") as fh:
        fh.write("group\tn_s11_marks\tmissing\troadmap\tbiosample\n")
        for g in groups:
            fh.write(f"{g['group']}\t{g['n_s11_marks']}\t{g['missing']}\t{g['roadmap']}\t{g['biosample']}\n")
    print(f"{'group':<28}{'S11 marks':>10}  {'roadmap':<8} missing / biosample")
    for g in groups[:15]:
        print(f"{g['group']:<28}{g['n_s11_marks']:>7}/7   {str(g['roadmap']):<8} {g['missing'] or '-'} / {g['biosample'][:60]}")
    if a.list_only:
        return

    chosen = next((g for g in groups if g["group"] == a.group), None) if a.group else groups[0]
    if chosen is None:
        sys.exit(f"group {a.group} not found; see encode_liver_candidates.tsv")
    print(f"\n[info] using {chosen['group']} ({chosen['n_s11_marks']}/7 marks; {chosen['biosample']})")
    donors = collections.defaultdict(set)
    for r in chosen["rows"]:
        donors[r["biosample"]].add(r["mark"])
    for d, ms in sorted(donors.items(), key=lambda x: -len(x[1] & set(MARKS))):
        print(f"       donor: {d}  ({len(ms & set(MARKS))}/7 S11 marks)")
    if a.donor:
        print(f"[info] restricting to donor matching '{a.donor}'")
    elif len(donors) > 1:
        print("[warn] this group mixes donors; without --donor the first experiment per mark is used")
    man = open(os.path.join(a.out, "marks_manifest.tsv"), "w")
    man.write("label\tmark\texperiment\tfile\toutput_type\tgroup\tbiosample\n")
    for mark in MARKS:
        exps = [r for r in chosen["rows"] if r["mark"] == mark and (not a.donor or a.donor in r["biosample"])]
        if not exps:
            print(f"[warn] {mark}: not in this group"); continue
        f = None
        for e in exps:
            f = pick_file(e["accession"])
            if f: break
        if not f:
            print(f"[warn] {mark}: no GRCh38 peak BED"); continue
        label = f"{mark}-human"
        dest = os.path.join(a.out, f"{label}.bed.gz")
        urllib.request.urlretrieve(BASE + f["href"], dest)
        man.write(f"{label}\t{mark}\t{e['accession']}\t{f['accession']}\t{f.get('output_type','')}\t{chosen['group']}\t{e['biosample']}\n")
        print(f"[ok] {label}: {e['accession']} {f['accession']} ({f.get('output_type','')})")
    man.close()
    print(f"[done] {a.out}/marks_manifest.tsv")

if __name__ == "__main__":
    main()
