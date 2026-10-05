from __future__ import annotations

import argparse
import gzip
from pathlib import Path

from enrichlib import Genome
from workspace_fit import build_candidate

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec", required=True, help="same +PATH,-PATH,... syntax as workspace_fit.py --candidate. "
                    "If the spec starts with '-' (e.g. a pure subtraction), write it as --spec=-path,... "
                    "(with the '='), or argparse will mistake the leading '-' for a flag.")
    ap.add_argument("--chrom-sizes", type=Path, default=Path("data/mm10.chrom.sizes"))
    ap.add_argument("--out", type=Path, required=True, help="output BED(.gz) path")
    a = ap.parse_args(argv)

    genome = Genome.from_file(a.chrom_sizes)
    ws = build_candidate(a.spec, genome)
    df = ws.to_frame()

    opener = gzip.open if str(a.out).endswith(".gz") else open
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with opener(a.out, "wt") as fh:
        df.to_csv(fh, sep="\t", header=False, index=False)

    print(f"wrote {a.out}: {len(df):,} intervals, {ws.bp:,} bp ({ws.bp / genome.whole().bp:.2%} of the genome)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
