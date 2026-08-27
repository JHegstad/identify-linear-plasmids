# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment Setup

```bash
conda env create -f environment.yml
conda activate linear_plasmids_env
```

Key dependencies: `biopython`, `numpy`, `pandas`, `pysam` (Python); `minimap2`, `samtools`, `blast` (external tools via bioconda).

## Running the Script

Basic usage (sequence-only analysis):
```bash
python identify_linear_plasmids.py -i assembly.fasta
```

With long-read auto-mapping (recommended for full evidence):
```bash
python identify_linear_plasmids.py -i assembly.fasta \
    --longread-fastq reads.fastq.gz --longread-preset map-ont --longread-threads 8
```

With short-read mapping (interleaved PE):
```bash
python identify_linear_plasmids.py -i assembly.fasta \
    --shortread-fastq reads.interleaved.fastq.gz --shortread-threads 8
```

With short-read mapping (separate R1/R2):
```bash
python identify_linear_plasmids.py -i assembly.fasta \
    --shortread-r1 R1.fastq.gz --shortread-r2 R2.fastq.gz
```

If both `--shortread-*` and `--longread-fastq` are supplied (and no explicit `--bam`), short reads are mapped first and take priority for the single shared BAM (`coverage_drop` + `copy_number`) — long reads are only auto-mapped as a fallback if short-read mapping isn't attempted or fails. This holds even when the assembly itself came from a long-read-only assembler (Flye/Autocycler/Hybracter-long): the coverage-drop check's wide-window signal is actually cleaner on Illumina data (adapter ligation is blocked outright by the hairpin fold) than on ONT (mappability is only reduced, not blocked, across the palindrome).

With all optional inputs:
```bash
python identify_linear_plasmids.py -i assembly.fasta \
    --bam reads.bam --gfa assembly.gfa \
    --annot prokka.gff --blast-db plsdb.fasta \
    --skani-db skani_db/sketches \
    --amrfinder --amrfinder-organism Enterococcus_faecium \
    --ref-gc 37.5 --chromosome-contigs chr1 \
    --json
```

Test data is in `LC495616/`: `LC495616.1.fasta` (assembly), `DRR194273.fastq.gz` (long reads), `PROKKA_04152026/PROKKA_04152026.gff` (annotation), `linear_plasmid_report.sorted.bam` (pre-mapped BAM).

Example run against test data:
```bash
python identify_linear_plasmids.py \
    -i LC495616/LC495616.1.fasta \
    --bam LC495616/linear_plasmid_report.sorted.bam \
    --annot LC495616/PROKKA_04152026/PROKKA_04152026.gff \
    --skani-db skani_db/sketches \
    --json -o LC495616/linear_plasmid_report
```

With a Hybracter output directory (batch mode, one sample per assembly — auto-discovers each sample's FINAL_OUTPUT FASTA, per_contig_stats.tsv, Flye/Plassembler GFA, and Hybracter's own QC'd long reads for auto-BAM-mapping; mutually exclusive with `-i`, and `--bam`/`--longread-fastq`/`--shortread-fastq`/`--shortread-r1`/`--shortread-r2` are not supported alongside it):
```bash
python identify_linear_plasmids.py --hybracter-dir hybracter_out/ \
    --annot prokka.gff --skani-db skani_db/sketches --json -o batch_report
```
Produces one combined `<prefix>.tsv`/`.json` across all samples, with a `sample` column.

To supply matching short reads per sample in batch mode (e.g. Illumina data for samples assembled long-read-only via Flye/Autocycler/Hybracter-long), use `--shortread-dir`: a directory of FASTQs auto-matched to each sample by filename (`{sample}_R1.fastq.gz`/`{sample}_R2.fastq.gz`, `{sample}_1.fastq.gz`/`{sample}_2.fastq.gz`, or interleaved `{sample}.fastq.gz`; `.fq`/`.fq.gz` also accepted). A matched sample's short reads are mapped and take priority over Hybracter's own long reads for that sample's BAM; samples with no match fall back to long-read auto-mapping as before.
```bash
python identify_linear_plasmids.py --hybracter-dir hybracter_out/ \
    --shortread-dir illumina_reads/ \
    --annot prokka.gff --skani-db skani_db/sketches --json -o batch_report
```

## Architecture

Single-file script (`identify_linear_plasmids.py`) with a modular pipeline:

| Module | Function(s) | Purpose |
|--------|-------------|---------|
| 0 | `parse_fasta_header`, `assess_header_metadata` | Extract circular/copy-number metadata from assembler FASTA headers (Unicycler, Plassembler, Hybracter) |
| 0b | `discover_hybracter_samples`, `parse_hybracter_contig_stats`, `find_shortreads_for_sample` | `--hybracter-dir` batch mode: per-sample FASTA/GFA/reads discovery under a Hybracter output root; `find_shortreads_for_sample` matches `--shortread-dir` files to a sample by filename for short-read-priority BAM mapping |
| 2 | `detect_terminal_hairpins` | Localized per-end hairpin/TIR fold-back detector (ported from linear-plasmid-hairpin-tools' `find_hairpins.py`, 2026-07) — replaces the old whole-window `detect_self_complement_ends`, which tested a symmetric-hairpin hypothesis that never fires on genuinely asymmetric ends (pELF1-type) |
| 2c | `detect_coverage_drop_ends` | Detect coverage drop at contig ends via BAM (hairpin inaccessibility artefact) |
| 2d | `detect_boundary_clip_signature` | Detect soft/hard-clip pile-up exactly at contig ends via BAM (circular-seam artefact — the counterpart to 2c's depth-drop check). **Not scored** (removed 2026-08, same day as added): once genuine Illumina data became available for the 6 circular negative controls, 4/6 came back falsely "consistent with linear" and the ratio ranges overlap the positive controls — no threshold cleanly separates them. Evidence still computed and left in JSON (`evidence["boundary_clip"]`) for diagnostic purposes only; see SCORING_WEIGHTS comment block |
| 3 | `gc_content`, `assess_gc` | GC% deviation from reference |
| 4 | `classify_size` | Size range check against known linear plasmid families |
| 5 | `parse_annotation`, `screen_genes` | Screen GFF3/Prokka TSV for linear-plasmid-associated gene keywords |
| 5a | `detect_tra_operon` | Screen annotation for the pELF2 conjugation (`tra`) operon, 12 genes (orf1-orf12) per Table 1 of Kurushima et al. 2026, PLOS Pathogens (10.1371/journal.ppat.1013937); only genes with a reported domain (traC/traD/traG plus orf1/orf7/orf11) are keyword-matchable, the rest have no recognisable domain per the paper itself. Reported in TSV (`tra_operon_hits`/`tra_operon_n`/`tra_operon_essential_n`) and JSON (`evidence["tra_operon"]`) — informational only, not scored |
| 5b | `run_amrfinder`, `interpret_amr_hits` | `--amrfinder`: AMR gene screening via NCBI AMRFinderPlus, run once per assembly and attributed to each contig by its `Contig id`. Reported in TSV/JSON (`amr_hit_count`/`amr_genes`/`amr_classes`) — informational only, not scored (same precedent as `resistance_genes`, commit 863d4b2) |
| 6 | `estimate_copy_number` | BAM-based copy number relative to chromosome |
| 7 | `parse_gfa_topology`, `is_linear_in_gfa` | Strand-aware GFA topology classification (ported from linear-plasmid-hairpin-tools' `autocycler_dotplot_classify.py`, 2026-07) — walks the graph per connected component (circular/linear/fragmented) instead of the old link-degree count, which could misclassify a same-strand circular self-loop as linear evidence |
| 7b | `annotate_gfa_hairpins` | `--annotate-gfa-hairpins`: writes Autocycler-style hairpin links into a GFA copy for Bandage visualization (diagnostic only, not scored) |
| 8 | `run_blast`, `interpret_blast_hits` | BLAST against plasmid DB (PLSDB); two-tier scoring |
| 9 | `compute_score` | Composite weighted scoring → confidence call (HIGH/MEDIUM/LOW/NONE) |

**Scoring system**: `SCORING_WEIGHTS` dict maps evidence categories to integer weights, rescaled (2026-07) so the achievable ceiling is 100 (`compute_score`'s `max_possible` accounts for the mutually-exclusive blast/skani tiers). `CONFIDENCE_THRESHOLDS` maps HIGH ≥35, MEDIUM ≥20, LOW ≥8. Two categories (`circular_flag_absent`, `self_complement_end`) were removed from scoring during revalidation — they were unreachable/non-discriminating on known-positive linear plasmids (see comments above `SCORING_WEIGHTS`). `CONTIG_SPECIFIC_SCORES` lists gene-based scores that are suppressed when annotation cannot be matched per-contig. Contigs with `circular=true` in header have gene-based scoring blocked (`CIRCULAR_DISQUALIFIES_GENE_SCORING`).

**Annotation handling**: `screen_genes` tries to match the contig ID against the annotation's `contig` column. For multi-contig annotations with no match (Prokka TSV prefix ≠ FASTA header), gene scoring is skipped — use GFF3 (`--annot`) for reliable per-contig mapping.

**Outputs**: `<prefix>.tsv` always; `<prefix>.json` with `--json`; `<prefix>.sorted.bam` when `--longread-fastq` triggers auto-mapping.
