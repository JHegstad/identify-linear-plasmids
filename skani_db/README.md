# SKANI plasmid database

Pre-sketched database for use with `--skani-db skani_db/sketches`.

## Current references (`refs/`)

31 pELF-family *E. faecium* linear plasmids from NCBI (source: `ncbi_pelf.fasta.gz`,
48 records deduplicated down to 31 unique sequences — GenBank/RefSeq duplicate
accessions of the same plasmid were dropped, keeping the first-encountered record).

| File | Accession | Description |
|------|-----------|--------------|
| `pELF1_LC495616.1.fasta` | LC495616.1 | *E. faecium* AA708 plasmid pELF1, complete genome |
| `pELF_2001564_CP086645.1.fasta` | CP086645.1 | *E. faecium* strain 2001564 plasmid pELF_2001564 |
| `pELF_2008300_CP115972.1.fasta` | CP115972.1 | *E. faecium* strain 2008300 plasmid pELF_2008300 |
| `pELF_AA242_AP026630.1.fasta` | AP026630.1 | *E. faecium* AA242 plasmid pELF_AA242 |
| `pELF_AA290_AP026634.1.fasta` | AP026634.1 | *E. faecium* AA290 plasmid pELF_AA290 |
| `pELF_AA55_AP026606.1.fasta` | AP026606.1 | *E. faecium* AA55 plasmid pELF_AA55 |
| `pELF_AA610_AP026643.1.fasta` | AP026643.1 | *E. faecium* AA610 plasmid pELF_AA610 |
| `pELF_AA818_AP026652.1.fasta` | AP026652.1 | *E. faecium* AA818 plasmid pELF_AA818 |
| `pELF_AA94_AP026614.1.fasta` | AP026614.1 | *E. faecium* AA94 plasmid pELF_AA94 |
| `pELF_AA96_AP026619.1.fasta` | AP026619.1 | *E. faecium* AA96 plasmid pELF_AA96 |
| `pELF_BSI_SJ40_CP076332.1.fasta` | CP076332.1 | *E. faecium* strain BSI_SJ40 plasmid pELF_BSI_SJ40 |
| `pELF_GK923_AP026657.1.fasta` | AP026657.1 | *E. faecium* GK923 plasmid pELF_GK923 (identical seq to GK941/GK961) |
| `pELF_JHP35_AP026581.1.fasta` | AP026581.1 | *E. faecium* JHP35 plasmid pELF_JHP35 |
| `pELF_JHP36_AP026588.1.fasta` | AP026588.1 | *E. faecium* JHP36 plasmid pELF_JHP36 (identical seq to JHP38) |
| `pELF_JHP38_AP026595.1.fasta` | AP026595.1 | *E. faecium* JHP38 plasmid pELF_JHP38 |
| `pELF_JHP80_AP026602.1.fasta` | AP026602.1 | *E. faecium* JHP80 plasmid pELF_JHP80 |
| `pELF_JHP9_AP026568.1.fasta` | AP026568.1 | *E. faecium* JHP9 plasmid pELF_JHP9 (identical seq to JHP10) |
| `pELF_mdr_AP026773.1.fasta` | AP026773.1 | *E. faecium* NUITM-VRE1 plasmid pELF_mdr |
| `pELF_N21-03014_CP197454.1.fasta` | CP197454.1 | *E. faecium* strain N21-03014 plasmid pELF_N21-03014 |
| `pELF_NZ_PV113235.1.fasta` | NZ_PV113235.1 | *E. faecium* strain VR2 plasmid pELF |
| `pELF_NZ_PV113236.1.fasta` | NZ_PV113236.1 | *E. faecium* strain VR3 plasmid pELF |
| `pELF_NZ_PV113237.1.fasta` | NZ_PV113237.1 | *E. faecium* strain VR15 plasmid pELF |
| `pELF_NZ_PV113238.1.fasta` | NZ_PV113238.1 | *E. faecium* strain VR27 plasmid pELF |
| `pELF_NZ_PV113239.1.fasta` | NZ_PV113239.1 | *E. faecium* strain VR40 plasmid pELF |
| `pELF_SJ10_CP076469.1.fasta` | CP076469.1 | *E. faecium* strain SJ10 plasmid pELF_SJ10 |
| `pELF_USZ_NZ_OU015710.1.fasta` | NZ_OU015710.1 | *E. faecium* isolate USZ_VRE32_P32 plasmid pELF_USZ |
| `pELF_USZ_NZ_OU016038.1.fasta` | NZ_OU016038.1 | *E. faecium* isolate USZ_VRE53_P46 plasmid pELF_USZ |
| `pELF_V12_AP031245.1.fasta` | AP031245.1 | *E. faecium* NUITM-VRE12 plasmid pELF_V12 |
| `pELF_V19_AP031257.1.fasta` | AP031257.1 | *E. faecium* NUITM-VRE19 plasmid pELF_V19 |
| `pELF_V3_AP031223.1.fasta` | AP031223.1 | *E. faecium* NUITM-VRE3 plasmid pELF_V3 |
| `pELF_V6_AP031234.1.fasta` | AP031234.1 | *E. faecium* NUITM-VRE6 plasmid pELF_V6 |

Previously this db also carried `51525510_plasmid00002.fasta` (an unrelated sample
plasmid, not pELF-family) — dropped when the db was rebuilt from `ncbi_pelf.fasta.gz`
(2026-08-25) so all references are pELF-family plasmids.

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
