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
| `ncbi-amrfinderplus` | `--amrfinder` antimicrobial resistance gene screening (run `amrfinder -u` once after install to fetch its database) |

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
    --amrfinder --amrfinder-organism Enterococcus_faecium \
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

### Batch wrapper scripts

`--hybracter-dir` handles a single Hybracter *batch* run (one shared
`FINAL_OUTPUT/{complete,incomplete}/` across samples). Two other common
layouts aren't a single Hybracter batch, so they're driven by standalone
wrapper scripts instead — each just loops the core script (and, where
needed, Bakta) once per sample directory and merges the results:

| Script | For | Per-sample input |
|---|---|---|
| `run_hybracter_bakta_batch.sh` | A collection of **independent single-sample** Hybracter runs, each with its own Bakta annotation already run against that sample's Hybracter contig IDs | `SAMPLE/hybracter/`, `SAMPLE/bakta/SAMPLE.gff3` |
| `run_bakta_autocycler_batch.sh` | Annotating [Autocycler](https://github.com/rrwick/Autocycler) consensus assemblies with Bakta, run *directly* against Autocycler's own numeric contig IDs | `SAMPLE_autocycler/autocycler_out/consensus_assembly.fasta` → `SAMPLE_autocycler/bakta/SAMPLE.gff3` |
| `run_autocycler_batch.sh` | The Autocycler assemblies themselves — plain `-i` mode per sample (not `--hybracter-dir`, since there's no Hybracter output at all here) | `SAMPLE_autocycler/autocycler_out/consensus_assembly.{fasta,gfa}`, `SAMPLE_autocycler/SAMPLE.fastq.gz`, and `SAMPLE_autocycler/bakta/SAMPLE.gff3` if present |

```bash
./run_hybracter_bakta_batch.sh [SRC_DIR] [THREADS]
./run_bakta_autocycler_batch.sh [SRC_DIR] [BAKTA_DB]
./run_autocycler_batch.sh [SRC_DIR] [THREADS]
```

All three take the sample-collection root as their first argument (each
has a project-specific default baked in — edit it, or just pass a path).
Every sample is run with the fullest evidence set the wrapper can find
(GFA topology, long-read auto-mapping, AMRFinderPlus, PLSDB BLAST, skani,
`--annotate-gfa-hairpins`, `--visualize`, `--json`), and each wrapper ends
by calling `combine_batch_results.py` to merge every sample's report into
one `<SRC_DIR>/linear_plasmid_batch_summary.{tsv,json}`.

**Why annotation is a separate step for Autocycler, and why the contig IDs
matter:** `screen_genes()` matches a GFF3's seqid against each contig by
*substring*, so an annotation only gives correct gene-based evidence when
its seqids are exactly the assembly's own contig IDs — a Bakta run against
a *different* assembly of the same isolate (e.g. Hybracter's
`chromosome00001`/`plasmid00001` naming) will silently produce wrong
matches against Autocycler's numeric IDs (`1`, `2`, `3`, ...) rather than
a clean skip, since e.g. contig `"1"` is a substring of `"chromosome00001"`.
`run_bakta_autocycler_batch.sh` avoids this by running Bakta straight
against each Autocycler `consensus_assembly.fasta` with
`--keep-contig-headers`, so its output seqids are Autocycler's own IDs and
match cleanly. Run it before `run_autocycler_batch.sh` — the latter picks
up `SAMPLE_autocycler/bakta/SAMPLE.gff3` automatically via `--annot` when
present, and otherwise runs without it. (It also passes `--skip-plot` to
Bakta: the circular-plot step divides by a contig-length-derived
`step_size` and crashes on the short junk contigs a multi-assembler
consensus sometimes leaves in — the GFF3 is already fully written by the
time that step runs, so this only skips an unused cosmetic image.)

`combine_batch_results.py SRC_DIR [-o OUT_PREFIX]` can also be run
standalone to (re-)merge `SRC_DIR/*/linear_plasmid/*.tsv` reports from any
of the wrappers above into one combined summary — it discovers exactly one
non-`.amrfinder.tsv` report per `linear_plasmid/` directory and adds a
leading `sample` column (taken from the report's own filename) if one
isn't already present.

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
| `--amrfinder` | Screen for AMR genes with NCBI AMRFinderPlus (reported per contig, not scored) |
| `--amrfinder-organism` | Organism for AMRFinderPlus point-mutation screening (e.g. `Enterococcus_faecium`) |
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

**Hard gates** are applied last, after every category above has scored, and
can override everything else back down to 0:

- **Chromosomes never score as linear plasmids**, regardless of evidence — a
  contig is treated as the chromosome if it's over 500 kb (the same ceiling
  `classify_size()` already treats as implausible for a linear plasmid), or
  `"chromosome"` appears in the contig id. This applies regardless of the
  assembler's own `circular=` call: a poorly-resolved assembly that fails to
  close the chromosome into a circle emits `circular=false` on a multi-Mb
  contig, and that shouldn't earn it linear-plasmid credit — if anything
  it's assembly incompleteness, not evidence of a genuine open replicon.
- **A plasmid-sized contig the assembler already calls `circular=true` only
  scores if it carries independently-detected hairpin/telomere structure**
  (a localized terminal fold-back, or a pELF1-type asymmetric hairpin +
  invertron end). Indirect evidence alone — a skani/BLAST hit, GFA
  topology, coverage drop, gene content, size/GC/copy-number — isn't
  sufficient on its own to override an assembler's circular call, since
  those signals are also seen on ordinary circular replicons.

When a gate fires, the contig's score and breakdown are zeroed and the
reason is recorded in the `gate_reason` column (TSV) / field (JSON) — check
it before assuming a 0/NONE result means "no evidence found" rather than
"gated."

## Output

- `<prefix>.tsv` — always written; one row per reported contig with score,
  confidence, `gate_reason` (non-empty only when a hard gate zeroed the
  contig), and a column per evidence category. With every contig gated or
  below `--min-score`, this is a valid header-less near-empty file rather
  than an error.
- `<prefix>.json` — with `--json`; full per-contig evidence detail.
- `<prefix>.sorted.bam` (+ `.bai`) — written when `--longread-fastq` or
  `--shortread-*` triggers auto-mapping.
- `<prefix>_<contig>_terminal.png` — with `--visualize`; terminal-structure
  plot (coverage, soft/hard clips, IDR arms) for flagged contigs.
- `<gfa-stem>.hairpins.gfa` — with `--annotate-gfa-hairpins`.
- `<prefix>.amrfinder.tsv` — with `--amrfinder`; raw AMRFinderPlus output for
  the whole assembly (per-contig summary columns `amr_hit_count`/`amr_genes`/
  `amr_classes` are folded into the main TSV/JSON).

Each batch wrapper script (see above) additionally writes
`<SRC_DIR>/linear_plasmid_batch_summary.{tsv,json}` — every sample's report
concatenated with a leading `sample` column — and
`<SRC_DIR>/_linear_plasmid_batch_logs/<SAMPLE>.log` per-sample run logs.
