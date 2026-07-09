#!/usr/bin/env bash
#
# Batch wrapper: runs identify_linear_plasmids.py with the fullest available
# evidence set on every sample under an Autocycler output root, where each
# sample is its own SAMPLE_autocycler/ directory containing:
#   SAMPLE_autocycler/SAMPLE.fastq.gz                    (raw ONT long reads)
#   SAMPLE_autocycler/autocycler_out/consensus_assembly.fasta
#   SAMPLE_autocycler/autocycler_out/consensus_assembly.gfa
#
# Uses SAMPLE_autocycler/bakta/SAMPLE.gff3 for --annot when present (see
# run_bakta_autocycler_batch.sh — annotates consensus_assembly.fasta
# directly with --keep-contig-headers, so its seqids are Autocycler's own
# numeric contig IDs and match cleanly). Earlier revisions of this batch
# deliberately ran without --annot: the only Bakta annotations available
# at the time were run against the *Hybracter* assembly's contig IDs
# (chromosome00001, plasmid00001, ...), and since screen_genes() matches
# contig IDs by substring, feeding those in against Autocycler's numeric
# IDs would have risked spurious matches (contig "1" substring-matches
# "chromosome00001") rather than a clean skip.
#
# Usage:
#   ./run_autocycler_batch.sh [SRC_DIR] [THREADS]
#
# Outputs land in SAMPLE_autocycler/linear_plasmid/SAMPLE.{tsv,json,...}.
# Logs land in SRC_DIR/_linear_plasmid_batch_logs/.

set -uo pipefail

SRC_DIR="${1:-/home/joachim/NGS/Projects/Kres-LRE/Linear_plasmids/new-putative-linear-plasmids/longreads/ont_reads/AUTOCYCLER_OUT_020726}"
THREADS="${2:-8}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IDENTIFY="$SCRIPT_DIR/identify_linear_plasmids.py"
BLAST_DB="/space/Databases/PLSDB/plsdb.fasta"
SKANI_DB="$SCRIPT_DIR/skani_db/sketches"
AMR_ORGANISM="Enterococcus_faecium"

LOG_DIR="$SRC_DIR/_linear_plasmid_batch_logs"
mkdir -p "$LOG_DIR"

# Activate the conda env if not already active
if [ "${CONDA_DEFAULT_ENV:-}" != "linear_plasmids_env" ]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate linear_plasmids_env
fi

echo "[BATCH] Source dir : $SRC_DIR"
echo "[BATCH] Script     : $IDENTIFY"
echo "[BATCH] BLAST db   : $BLAST_DB"
echo "[BATCH] skani db   : $SKANI_DB"
echo "[BATCH] Threads    : $THREADS"
echo "[BATCH] Logs       : $LOG_DIR"
echo

FAILED=()
SKIPPED=()
DONE=()

for sample_dir in "$SRC_DIR"/*_autocycler/; do
    [ -d "$sample_dir" ] || continue
    dirname_base="$(basename "$sample_dir")"
    sample="${dirname_base%_autocycler}"

    fasta="${sample_dir}autocycler_out/consensus_assembly.fasta"
    gfa="${sample_dir}autocycler_out/consensus_assembly.gfa"
    fastq="${sample_dir}${sample}.fastq.gz"
    gff3="${sample_dir}bakta/${sample}.gff3"

    if [ ! -f "$fasta" ]; then
        echo "[SKIP] $sample: no consensus_assembly.fasta"
        SKIPPED+=("$sample")
        continue
    fi

    gfa_args=()
    [ -f "$gfa" ] && gfa_args=(--gfa "$gfa")

    annot_args=()
    if [ -f "$gff3" ]; then
        annot_args=(--annot "$gff3")
    else
        echo "[WARN] $sample: no Bakta annotation at $gff3, running without --annot"
    fi

    fastq_args=()
    if [ -f "$fastq" ]; then
        fastq_args=(--longread-fastq "$fastq" --longread-preset map-ont --longread-threads "$THREADS")
    else
        echo "[WARN] $sample: no ${sample}.fastq.gz found, running without long-read mapping"
    fi

    out_dir="${sample_dir}linear_plasmid"
    mkdir -p "$out_dir"
    out_prefix="$out_dir/$sample"

    echo "==================================================================="
    echo "[SAMPLE] $sample  ($(date '+%Y-%m-%d %H:%M:%S'))"
    echo "==================================================================="

    python "$IDENTIFY" \
        -i "$fasta" \
        "${gfa_args[@]}" \
        "${annot_args[@]}" \
        "${fastq_args[@]}" \
        --amrfinder --amrfinder-organism "$AMR_ORGANISM" --amrfinder-threads "$THREADS" \
        --blast-db "$BLAST_DB" \
        --skani-db "$SKANI_DB" \
        --annotate-gfa-hairpins \
        --visualize \
        --json \
        -o "$out_prefix" \
        > "$LOG_DIR/${sample}.log" 2>&1
    status=$?

    if [ "$status" -ne 0 ]; then
        echo "[ERROR] $sample failed (exit $status) — see $LOG_DIR/${sample}.log"
        FAILED+=("$sample")
    else
        echo "[OK] $sample → $out_prefix.tsv"
        DONE+=("$sample")
    fi
    echo
done

echo "==================================================================="
echo "[BATCH] Complete: ${#DONE[@]} ok, ${#FAILED[@]} failed, ${#SKIPPED[@]} skipped"
[ "${#FAILED[@]}" -gt 0 ]  && echo "[BATCH] FAILED : ${FAILED[*]}"
[ "${#SKIPPED[@]}" -gt 0 ] && echo "[BATCH] SKIPPED: ${SKIPPED[*]}"
echo "==================================================================="

# Combine all per-sample TSV/JSON into one batch summary
python "$SCRIPT_DIR/combine_batch_results.py" "$SRC_DIR" \
    -o "$SRC_DIR/linear_plasmid_batch_summary"

[ "${#FAILED[@]}" -eq 0 ]
