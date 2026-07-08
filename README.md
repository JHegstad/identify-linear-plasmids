# identify-linear-plasmids

A multi-evidence scoring pipeline for identifying **linear plasmids** in
bacterial genome assemblies — in particular pELF1-type linear plasmids
(asymmetric hairpin/invertron ends), as described in:

- Hashimoto et al. (2019) *Front. Microbiol.* — original pELF1 description (asymmetric ends, invertron type, IDR/TATA hairpin, coverage drop)
- Boumamoud et al. (2022) *mBio* — VREfm pELF linear plasmids
- Hashimoto et al. (2023) *AAC* — Enterococcal linear plasmids

Most bacterial plasmid-finding tools assume every plasmid is circular. This
script instead combines several independent lines of evidence — assembler
metadata, terminal sequence structure, gene content, BAM coverage, assembly
graph topology, and homology search — into a single weighted score and a
HIGH/MEDIUM/LOW/NONE confidence call per contig, so genuinely linear
replicons aren't silently mis-assembled or discarded as "incomplete."

## Installation

Requires [conda](https://docs.conda.io/) (or [mamba](https://mamba.readthedocs.io/), recommended for speed):

```bash
git clone git@github.com:JHegstad/identify-linear-plasmids.git
cd identify-linear-plasmids
conda env create -f environment.yml
conda activate linear_plasmids_env
```

### Dependencies

Installed automatically via `environment.yml`:

| Package | Required for |
|---|---|
| `python=3.10`, `biopython`, `numpy`, `pandas` | Core pipeline (always required) |
| `pysam` | BAM-based evidence: coverage drop, copy number, `--bam`/auto-mapping |
| `matplotlib` | `--visualize` terminal-structure plots |
| `minimap2`, `samtools` | Auto-mapping via `--longread-fastq` / `--shortread-*` |
| `blast` | `--blast-db` homology search against a plasmid database (e.g. PLSDB) |
| `skani` | `--skani-db` ANI-based plasmid database search |

Only `biopython`/`numpy`/`pandas` are strictly required — the script degrades
gracefully when the other optional inputs/tools aren't supplied, using
whatever evidence is available.

## Usage

Sequence-only analysis (no reads, no annotation — uses assembler header
metadata, sequence structure, and size/GC alone):

```bash
python identify_linear_plasmids.py -i assembly.fasta
```

With long-read auto-mapping (recommended — enables coverage-drop and IDR
evidence; mapped automatically with minimap2, no manual BAM needed):

```bash
python identify_linear_plasmids.py -i assembly.fasta \
    --longread-fastq reads.fastq.gz --longread-preset map-ont --longread-threads 8
```

With short-read mapping (interleaved paired-end):

```bash
python identify_linear_plasmids.py -i assembly.fasta \
    --shortread-fastq reads.interleaved.fastq.gz --shortread-threads 8
```

With short-read mapping (separate R1/R2):

```bash
python identify_linear_plasmids.py -i assembly.fasta \
    --shortread-r1 R1.fastq.gz --shortread-r2 R2.fastq.gz
```

With every optional input (fullest evidence set):

```bash
python identify_linear_plasmids.py -i assembly.fasta \
    --bam reads.bam --gfa assembly.gfa \
    --annot prokka.gff --blast-db plsdb.fasta \
    --skani-db skani_db/sketches \
    --ref-gc 37.5 --chromosome-contigs chr1 \
    --json
```

### Batch mode: Hybracter output directories

If your assemblies come from [Hybracter](https://github.com/gbouras13/hybracter),
point the script at the Hybracter output root instead of a single FASTA. It
auto-discovers every sample's final assembly, per-contig circularity/typing
info, Flye/Plassembler assembly graph, and Hybracter's own QC'd long reads
(for automatic BAM mapping), and runs the full pipeline once per sample:

```bash
python identify_linear_plasmids.py --hybracter-dir hybracter_out/ \
    --annot prokka.gff --skani-db skani_db/sketches --json -o batch_report
```

This mode is mutually exclusive with `-i`/`--input`, and with
`--bam`/`--longread-fastq`/`--shortread-*` (ambiguous across multiple
samples — reads are auto-discovered per sample instead). Output is one
combined `<prefix>.tsv`/`.json` covering every sample, with a leading
`sample` column.

## Key options

Run `python identify_linear_plasmids.py --help` for the full list. The most
commonly used:

| Option | Purpose |
|---|---|
| `-i / --input` | Input assembly FASTA (mutually exclusive with `--hybracter-dir`) |
| `--hybracter-dir` | Hybracter output root — batch mode, see above |
| `-o / --output` | Output prefix (default: `linear_plasmid_report`) |
| `--bam` | Pre-mapped BAM (skip auto-mapping) |
| `--longread-fastq` | Long reads to auto-map with minimap2 (Nanopore/PacBio) |
| `--shortread-fastq` / `--shortread-r1`/`--shortread-r2` | Illumina reads to auto-map |
| `--gfa` | Assembly graph (Unicycler/Flye/Plassembler) for topology evidence |
| `--annot` | GFF3/Prokka TSV annotation for gene-based evidence |
| `--blast-db` | Plasmid database (e.g. PLSDB) for BLAST homology search |
| `--skani-db` | Pre-sketched or FASTA plasmid database for SKANI ANI search |
| `--chromosome-contigs` | Chromosome contig ID(s), for copy-number normalisation |
| `--ref-gc` | Reference chromosome GC% (auto-computed if omitted) |
| `--annotate-gfa-hairpins` | Write a Bandage-visualization GFA copy with hairpin links annotated (diagnostic only) |
| `--json` | Also write per-contig JSON detail alongside the TSV |
| `--visualize` | Generate terminal-structure PNGs for HIGH/MEDIUM/IDR-flagged contigs |

## How scoring works

Each contig accumulates points from independent evidence categories —
assembler header metadata, terminal hairpin/inverted-repeat structure,
asymmetric (hairpin + invertron) end pattern, size, GC deviation,
plasmid-associated genes (partition systems, replication genes,
toxin-antitoxin systems, IS elements), BLAST/SKANI homology to known linear
plasmids, BAM coverage drop at contig ends, assembly graph topology, and copy
number. The weighted total maps to a confidence call:

| Score | Confidence |
|---|---|
| ≥ 35 | **HIGH** |
| 20–34 | **MEDIUM** |
| 8–19 | **LOW** |
| < 8 | NONE (not reported by default; adjust with `--min-score`) |

The maximum achievable score is 100. See the `SCORING_WEIGHTS` and
`CONFIDENCE_THRESHOLDS` dicts at the top of `identify_linear_plasmids.py` for
the exact per-category weights and rationale.

## Output

- `<prefix>.tsv` — always written; one row per reported contig with score,
  confidence, and a column per evidence category.
- `<prefix>.json` — with `--json`; full per-contig evidence detail.
- `<prefix>.sorted.bam` (+ `.bai`) — written when `--longread-fastq` or
  `--shortread-*` triggers auto-mapping.
- `<prefix>_<contig>_terminal.png` — with `--visualize`; terminal-structure
  plot (coverage, soft/hard clips, IDR arms) for flagged contigs.
- `<gfa-stem>.hairpins.gfa` — with `--annotate-gfa-hairpins`.
