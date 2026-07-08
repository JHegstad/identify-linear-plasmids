#!/usr/bin/env python3
"""
Combine per-sample identify_linear_plasmids.py outputs (written to
SAMPLE_DIR/linear_plasmid/*.tsv[.json]) into one batch-wide TSV/JSON.

Each per-sample run is invoked in single-sample mode (either --hybracter-dir
discovering one sample, or plain -i), so its TSV has no 'sample' column
(that column is only added by the script's own multi-sample batch loop).
This script adds it back (from the report's own filename stem) and
concatenates everything.

Usage:
    combine_batch_results.py SRC_DIR [-o OUT_PREFIX]
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("src_dir", help="Root dir containing SAMPLE_DIR/linear_plasmid/*.tsv")
    ap.add_argument("-o", "--output", default="linear_plasmid_batch_summary",
                     help="Output prefix (default: linear_plasmid_batch_summary)")
    args = ap.parse_args()

    src = Path(args.src_dir)
    # Exactly one main report TSV per linear_plasmid/ dir — not byproduct
    # files written to the same dir, e.g. <sample>.amrfinder.tsv. Matched by
    # content (not filename == parent dir name) since different wrapper
    # scripts may name the report differently from the sample directory.
    tsvs = []
    for lp_dir in sorted(src.glob("*/linear_plasmid")):
        candidates = [p for p in lp_dir.glob("*.tsv") if not p.name.endswith(".amrfinder.tsv")]
        if len(candidates) != 1:
            print(f"[WARN] {lp_dir}: expected exactly 1 main report TSV, found "
                  f"{len(candidates)}; skipping", file=sys.stderr)
            continue
        tsvs.append(candidates[0])
    if not tsvs:
        print(f"[ERROR] No per-sample TSVs found under {src}/*/linear_plasmid/", file=sys.stderr)
        sys.exit(1)

    frames = []
    combined_json = []
    for tsv_path in tsvs:
        sample = tsv_path.stem
        try:
            df = pd.read_csv(tsv_path, sep="\t")
        except pd.errors.EmptyDataError:
            # Every contig in this sample was gated to 0 / filtered below
            # --min-score, so identify_linear_plasmids.py wrote a columnless
            # file (pd.DataFrame([]).to_csv() writes just a newline, not 0
            # bytes). Nothing to report for this sample.
            print(f"[INFO] {sample}: no contigs reported (empty TSV) — skipping")
            continue
        if "sample" not in df.columns:
            df.insert(0, "sample", sample)
        frames.append(df)

        json_path = tsv_path.with_suffix(".json")
        if json_path.exists():
            with open(json_path) as fh:
                records = json.load(fh)
            for r in records:
                r.setdefault("sample", sample)
            combined_json.extend(records)

    if not frames:
        print("[INFO] No sample had any reported contigs — nothing to combine.")
        sys.exit(0)

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.sort_values("score", ascending=False)

    out_tsv = f"{args.output}.tsv"
    combined.to_csv(out_tsv, sep="\t", index=False)
    print(f"[OUTPUT] Combined TSV  → {out_tsv}  ({len(combined)} contigs, {len(tsvs)} samples)")

    if combined_json:
        out_json = f"{args.output}.json"
        with open(out_json, "w") as fh:
            json.dump(combined_json, fh, indent=2)
        print(f"[OUTPUT] Combined JSON → {out_json}")

    print("\nConfidence summary:")
    print(combined.groupby(["sample", "confidence"]).size().unstack(fill_value=0)
          .reindex(columns=["HIGH", "MEDIUM", "LOW"], fill_value=0))


if __name__ == "__main__":
    main()
