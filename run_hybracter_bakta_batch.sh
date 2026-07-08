#!/usr/bin/env bash
#
# Batch wrapper: runs identify_linear_plasmids.py with the fullest available
# evidence set on every sample under a Hybracter+Bakta collection directory,
# where each sample is its own single-sample Hybracter run (SAMPLE/hybracter/)
# with a matching whole-assembly Bakta annotation (SAMPLE/bakta/SAMPLE.gff3).
#
# This directory layout is NOT a single Hybracter batch output (one
# FINAL_OUTPUT/{complete,incomplete}/ shared across samples) — it's one
# independent Hybracter run per sample. So identify_linear_plasmids.py's own
# --hybracter-dir batch mode is invoked once per sample (each discovers
# exactly one sample), which gives us per-sample auto-GFA, auto-BAM-mapping
# from Hybracter's QC'd long reads, and auto chromosome/plasmid contig
# typing from per_contig_stats.tsv — plus a per-sample Bakta annotation and
# AMRFinderPlus/BLAST/skani evidence layered on top.
#
# Usage:
#   ./run_hybracter_bakta_batch.sh [SRC_DIR] [THREADS]
#
# Outputs land in SAMPLE_DIR/linear_plasmid/SAMPLE.{tsv,json,...}, colocated
# with the other per-sample tool outputs (bakta/, amrfinderplus/, mlst/, ...).
# Logs land in SRC_DIR/_linear_plasmid_batch_logs/.

set -uo pipefail

SRC_DIR="${1:-/home/joachim/NGS/Projects/Kres-LRE/Linear_plasmids/hybracter_bakta}"
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

for sample_dir in "$SRC_DIR"/*/; do
    sample="$(basename "$sample_dir")"
    hyb_dir="${sample_dir}hybracter"

    if [ ! -d "$hyb_dir" ]; then
        echo "[SKIP] $sample: no hybracter/ subdirectory"
        SKIPPED+=("$sample")
        continue
    fi

    gff3="${sample_dir}bakta/${sample}.gff3"
    annot_args=()
    if [ -f "$gff3" ]; then
        annot_args=(--annot "$gff3")
    else
        echo "[WARN] $sample: no Bakta annotation at $gff3, running without --annot"
    fi

    out_dir="${sample_dir}linear_plasmid"
    mkdir -p "$out_dir"
    out_prefix="$out_dir/$sample"

    echo "==================================================================="
    echo "[SAMPLE] $sample  ($(date '+%Y-%m-%d %H:%M:%S'))"
    echo "==================================================================="

    python "$IDENTIFY" \
        --hybracter-dir "$hyb_dir" \
        "${annot_args[@]}" \
        --amrfinder --amrfinder-organism "$AMR_ORGANISM" --amrfinder-threads "$THREADS" \
        --blast-db "$BLAST_DB" \
        --skani-db "$SKANI_DB" \
        --longread-threads "$THREADS" \
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
