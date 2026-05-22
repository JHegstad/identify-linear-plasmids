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

With all optional inputs:
```bash
python identify_linear_plasmids.py -i assembly.fasta \
    --bam reads.bam --gfa assembly.gfa \
    --annot prokka.gff --blast-db plsdb.fasta \
    --skani-db skani_db/sketches \
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

## Architecture

Single-file script (`identify_linear_plasmids.py`) with a modular pipeline:

| Module | Function(s) | Purpose |
|--------|-------------|---------|
| 0 | `parse_fasta_header`, `assess_header_metadata` | Extract circular/copy-number metadata from assembler FASTA headers (Unicycler, Plassembler, Hybracter) |
| 2 | `detect_self_complement_ends` | Detect hairpin/palindromic ends via end-window reverse-complement identity |
| 2c | `detect_coverage_drop_ends` | Detect coverage drop at contig ends via BAM (hairpin inaccessibility artefact) |
| 3 | `gc_content`, `assess_gc` | GC% deviation from reference |
| 4 | `classify_size` | Size range check against known linear plasmid families |
| 5 | `parse_annotation`, `screen_genes` | Screen GFF3/Prokka TSV for linear-plasmid-associated gene keywords |
| 6 | `estimate_copy_number` | BAM-based copy number relative to chromosome |
| 7 | `parse_gfa_topology`, `is_linear_in_gfa` | GFA assembly graph linear topology detection |
| 8 | `run_blast`, `interpret_blast_hits` | BLAST against plasmid DB (PLSDB); two-tier scoring |
| 9 | `compute_score` | Composite weighted scoring → confidence call (HIGH/MEDIUM/LOW/NONE) |

**Scoring system**: `SCORING_WEIGHTS` dict maps evidence categories to integer weights. `CONFIDENCE_THRESHOLDS` maps HIGH ≥70, MEDIUM ≥40, LOW ≥15. `CONTIG_SPECIFIC_SCORES` lists gene-based scores that are suppressed when annotation cannot be matched per-contig. Contigs with `circular=true` in header have gene-based scoring blocked (`CIRCULAR_DISQUALIFIES_GENE_SCORING`).

**Annotation handling**: `screen_genes` tries to match the contig ID against the annotation's `contig` column. For multi-contig annotations with no match (Prokka TSV prefix ≠ FASTA header), gene scoring is skipped — use GFF3 (`--annot`) for reliable per-contig mapping.

**Outputs**: `<prefix>.tsv` always; `<prefix>.json` with `--json`; `<prefix>.sorted.bam` when `--longread-fastq` triggers auto-mapping.
