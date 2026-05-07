# SKANI plasmid database

Pre-sketched database for use with `--skani-db skani_db/sketches`.

## Current references (`refs/`)

| File | Source | Description |
|------|--------|-------------|
| `pELF1_LC495616.1.fasta` | LC495616/LC495616.1.fasta | *E. faecium* AA708 pELF1, complete genome (NCBI LC495616.1) |
| `51525510_plasmid00002.fasta` | 51525510/51525510_plasmid00002.fasta | Sample 51525510 plasmid assembly |

## Rebuilding the sketches

```bash
skani sketch skani_db/refs/*.fasta -o skani_db/sketches
```

The output directory must not already exist — delete `skani_db/sketches/` first if rebuilding.

## Adding new references

1. Copy or symlink the FASTA into `skani_db/refs/` with a descriptive name that includes the
   plasmid family (e.g. `pELF2_accession.fasta`, `pBSSB1_accession.fasta`). The filename stem
   is used as a fallback for keyword matching if the sequence header is absent.
2. Delete `skani_db/sketches/` and re-run the sketch command above.

## Keyword matching

`interpret_blast_hits()` is reused for SKANI results. It matches the full sequence header
(`Ref_name` column from skani ≥0.3) against two-tier keyword lists:

- **Tier 2 (20 pts):** `\bpelf`, `\bpbssb` — named linear plasmid families
- **Tier 1 (15 pts):** `linear plasmid`, `linear chromosome`, `telomere`, …

Ensure reference FASTA headers or filenames contain these terms.
