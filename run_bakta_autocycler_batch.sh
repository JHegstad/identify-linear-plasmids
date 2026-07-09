#!/usr/bin/env bash
#
# Batch wrapper: annotates every SAMPLE_<suffix>/autocycler_out/consensus_assembly.fasta
# with Bakta (anno_env conda environment), writing to SAMPLE_<suffix>/bakta/.
#
# Works for any Autocycler-shaped consensus-assembly output, not just runs
# literally named "*_autocycler" — e.g. "*_bcapa" (a reduced Flye+Plassembler-
# only assembler ensemble run through the same Autocycler consensus step).
# SUFFIX selects which one: SRC_DIR/*_<suffix>/autocycler_out/... is globbed.
#
# Uses --keep-contig-headers so the GFF3 seqids stay exactly "1", "2", "3", ...
# matching the consensus assembly's own contig IDs — this is what lets
# identify_linear_plasmids.py's --annot matching work correctly afterwards
# (unlike the pre-existing Bakta annotations for some of these isolates,
# which were run against Hybracter's chromosome00001/plasmid00001 contig IDs
# and can't be safely reused here).
#
# --skip-plot: Bakta's circular-plot step divides by a step_size derived from
# contig length and crashes (ValueError: range() arg 3 must not be zero) on
# the very short junk contigs a multi-assembler consensus sometimes leaves in
# (a handful of bp long). The GFF3/TSV/etc. are already fully written by that
# point regardless — this only skips the cosmetic plot image, which
# identify_linear_plasmids.py doesn't use anyway.
#
# Usage:
#   ./run_bakta_autocycler_batch.sh [SRC_DIR] [BAKTA_DB] [SUFFIX]

set -uo pipefail

SRC_DIR="${1:-/home/joachim/NGS/Projects/Kres-LRE/Linear_plasmids/new-putative-linear-plasmids/longreads/ont_reads/AUTOCYCLER_OUT_020726}"
BAKTA_DB_PATH="${2:-/space/Databases/Bakta-DB/db-light}"
SUFFIX="${3:-autocycler}"

if [ "${CONDA_DEFAULT_ENV:-}" != "anno_env" ]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate anno_env
fi

LOG_DIR="$SRC_DIR/_bakta_batch_logs"
mkdir -p "$LOG_DIR"

echo "[BATCH] Source dir : $SRC_DIR"
echo "[BATCH] Bakta db   : $BAKTA_DB_PATH"
echo "[BATCH] Suffix     : *_$SUFFIX"
echo "[BATCH] Logs       : $LOG_DIR"
echo

FAILED=()
SKIPPED=()
DONE=()

for sample_dir in "$SRC_DIR"/*_"$SUFFIX"/; do
    [ -d "$sample_dir" ] || continue
    dirname_base="$(basename "$sample_dir")"
    sample="${dirname_base%_$SUFFIX}"

    fasta="${sample_dir}autocycler_out/consensus_assembly.fasta"
    if [ ! -f "$fasta" ]; then
        echo "[SKIP] $sample: no consensus_assembly.fasta"
        SKIPPED+=("$sample")
        continue
    fi

    out_dir="${sample_dir}bakta"

    echo "==================================================================="
    echo "[SAMPLE] $sample  ($(date '+%Y-%m-%d %H:%M:%S'))"
    echo "==================================================================="

    bakta --db "$BAKTA_DB_PATH" \
        --keep-contig-headers \
        --skip-plot \
        --prefix "$sample" \
        --output "$out_dir" \
        --force \
        "$fasta" \
        > "$LOG_DIR/${sample}.log" 2>&1
    status=$?

    if [ "$status" -ne 0 ]; then
        echo "[ERROR] $sample failed (exit $status) — see $LOG_DIR/${sample}.log"
        FAILED+=("$sample")
    else
        echo "[OK] $sample → $out_dir/$sample.gff3"
        DONE+=("$sample")
    fi
    echo
done

echo "==================================================================="
echo "[BATCH] Complete: ${#DONE[@]} ok, ${#FAILED[@]} failed, ${#SKIPPED[@]} skipped"
[ "${#FAILED[@]}" -gt 0 ]  && echo "[BATCH] FAILED : ${FAILED[*]}"
[ "${#SKIPPED[@]}" -gt 0 ] && echo "[BATCH] SKIPPED: ${SKIPPED[*]}"
echo "==================================================================="

[ "${#FAILED[@]}" -eq 0 ]
