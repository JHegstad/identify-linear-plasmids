#!/usr/bin/env python3
"""
identify_linear_plasmids.py
============================
Identifies linear plasmids in assembled sequences using multiple evidence
layers derived from the literature:

  Hashimoto et al. (2019) Front. Microbiol. — pELF1 original description (asymmetric ends,
                                              invertron type, IDR/TATA hairpin, coverage drop)
  Boumamoud et al. (2022) mBio            — VREfm pELF linear plasmids
  Hashimoto et al. (2023) AAC             — Enterococcal linear plasmids

Usage
-----
  python identify_linear_plasmids.py -i assembly.fasta [options]
  python identify_linear_plasmids.py -i assembly.fasta --longread-fastq reads.fastq.gz
  python identify_linear_plasmids.py -i assembly.fasta --longread-fastq reads.fastq.gz \\
      --longread-preset map-ont --longread-threads 8
  python identify_linear_plasmids.py -i assembly.fasta --bam reads.bam --gfa assembly.gfa
  python identify_linear_plasmids.py -i assembly.fasta --blast-db plsdb.fasta --annot prokka.tsv

Long-read mapping
-----------------
  --longread-fastq maps reads automatically with minimap2 and writes
  <output-prefix>.sorted.bam (+ .bai index).  Requires minimap2 and samtools
  on PATH.  The BAM is then used for IDR palindrome detection in reads,
  coverage-drop analysis, and copy-number estimation.

  Presets:
    map-ont   Nanopore (default)
    map-pb    PacBio CLR
    map-hifi  PacBio HiFi / Revio

Output: TSV report + per-contig JSON detail.

Dependencies (install with pip):
  biopython, numpy, pandas, pysam (optional, for BAM analysis)
External tools (optional, for long-read mapping):
  minimap2, samtools
"""

import argparse
import copy
import glob
import json
import math
import os
import re
import subprocess
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS  (derived from Hashimoto 2019, Boumamoud 2022, Hashimoto 2023)
# ─────────────────────────────────────────────────────────────────────────────

# Size ranges (bp) for known linear plasmid families
SIZE_RANGES = {
    "enterococcal_pELF":  (70_000, 170_000),  # Hashimoto 2023 / Boumamoud 2022: 76–164 kb
    "general":            (10_000, 500_000),  # permissive catch-all
}

# Invertron end parameters (Hashimoto 2019)
# Terminal proteins (TPs) covalently bound at 5' end — detected by:
#   lambda exonuclease resistance / exonuclease III sensitivity (wet lab only)
# In-silico proxy: absence of perfect self-complement at that end
INVERTRON_TP_GENE_KEYWORDS = [
    "terminal protein", "tp", "tpg", "tap", "telomere associated protein",
    "telomere-associated protein",
]

# Toxin-antitoxin system gene families (Hashimoto / Boumamoud)
TAS_GENES = {
    "relBE":  ["relB", "relE"],
    "vapBC":  ["vapB", "vapC"],
    "parDE":  ["parD", "parE"],
    "dinJ_yoeB": ["dinJ", "yoeB"],
    "higBA":  ["higB", "higA"],
    "mazEF":  ["mazE", "mazF"],
}

# Partition system Pfam accessions
PAR_PFAM = {
    "parA": ["PF13614", "PF01656"],
    "parB": ["PF01672"],
}

# Key gene keywords for linear plasmid identification
LINEAR_PLASMID_GENE_KEYWORDS = [
    # Partition
    "parA", "parB", "soj", "sopA", "sopB",
    # Replication — repB/Rep_2 superfamily specifically found on pELF1 (Hashimoto 2019)
    "repB", "repA", "rep", "rep_2", "rep2", "replication initiation",
    # Structural / maintenance
    "ftsk", "spoiiiie", "ftsk/spoiiiie",
    # Terminal protein (invertron end, Hashimoto 2019)
    "terminal protein", "tpg", "tap", "telomere associated",
    # TAS toxins
    "rele", "vape", "pare", "yoeb", "higb", "mazf",
    # TAS antitoxins
    "relb", "vapb", "pard", "dinj", "higa", "maze",
    # Enterococcal pELF — IS1216E flanks vanM cluster in same direction (Hashimoto 2019)
    "is1216", "is1216e", "is1216v", "is1542", "isefm1", "isefa11",
    # Conjugation / transfer machinery (ftsK confirmed on pELF1, Hashimoto 2019)
    "traN", "traD", "virB", "ftsk",
    # Resistance
    "vana", "vanb", "vanc", "vand", "vanm",
    # Aminoglycoside cluster on pELF1: aadE-sat4-aphA-3 (Hashimoto 2019)
    "aade", "sat4", "apha", "apha-3", "aphA3",
]

# Resistance genes associated with linear plasmids (Hashimoto 2019 / Boumamoud 2022 / Hashimoto 2023)
RESISTANCE_GENES = {
    "vanA":  "Glycopeptide (vancomycin, high-level)",
    "vanB":  "Glycopeptide (vancomycin, moderate)",
    "vanC":  "Glycopeptide (vancomycin, low-level)",
    "vanD":  "Glycopeptide (vancomycin)",
    "vanM":  "Glycopeptide (vancomycin, VanM-type; pELF1 Hashimoto 2019)",
    "blaZ":  "Beta-lactam (ampicillin)",
    "ermB":  "Macrolide (erythromycin; confirmed on pELF1)",
    "aac6":  "Aminoglycoside (gentamicin)",
    "optrA": "Oxazolidinone (linezolid)",
    "dfrG":  "Trimethoprim",
    # pELF1-specific resistance cluster: aadE-sat4-aphA-3 (Hashimoto 2019)
    "aadE":  "Aminoglycoside (streptomycin; pELF1 cluster)",
    "sat4":  "Streptothricin (pELF1 cluster)",
    "aphA3": "Aminoglycoside (kanamycin; pELF1 cluster)",
}

# IS elements associated with linear plasmid rearrangements
# Hashimoto 2019: IS1216V + IS1542 in Tn1546-like vanA; IS1216E flanking vanM (same direction)
# ISEfm1, ISEfa11 interrupting sat4
LINEAR_PLASMID_IS = [
    "IS1216", "IS1216E", "IS1216V", "IS1542",
    "ISEfm1", "ISEfa11",
]

# Scoring weights for each evidence category (Hashimoto 2019 / Boumamoud 2022 / Hashimoto 2023)
#
# Rescaled 2026-07 so the realistic achievable ceiling (accounting for the
# mutually-exclusive blast_hit/blast_hit_linear_db and skani_hit/skani_hit_linear_db
# tiers, see compute_score) is exactly 100. Two categories were removed outright
# during revalidation against known-positive linear plasmids (pELF1 reference,
# 51525510):
#   - circular_flag_absent: fired on every test run regardless of true topology —
#     it's true whenever assembler-header metadata is simply absent (e.g. any
#     non-Unicycler input), so it rewards missing metadata as if it were evidence.
#   - self_complement_end: tests whether the whole first/last 500bp windows are
#     mutual reverse-complements, i.e. a *symmetric* hairpin telomere. pELF1-type
#     plasmids have *asymmetric* ends (one hairpin, one invertron — see
#     asymmetric_ends/detect_asymmetric_ends) and scored only 30.8% identity
#     (chance level) on the literal pELF1 reference sequence — this metric tests
#     the wrong biological hypothesis for the plasmids this tool targets.
SCORING_WEIGHTS = {
    "hairpin_end":               8,   # Hairpin/palindromic end structure
    "asymmetric_ends":           6,   # One hairpin + one invertron end (pELF1-type, Hashimoto 2019)
    "invertron_tp_gene":         5,   # Terminal protein gene present (invertron end)
    "size_range":                4,   # Size consistent with known linear plasmids
    "gc_deviation":              4,   # GC% lower than typical chromosome
    "par_system":                4,   # ParA/ParB partition system
    "repb_rep2":                 4,   # repB/Rep_2 superfamily replication gene (Hashimoto 2019)
    "ftsK_parA_repB_combo":      5,   # ftsK + parA + repB co-occurrence (pELF1 signature)
    "tas_system":                4,   # Toxin-antitoxin system

    "is_elements":               2,   # IS elements associated with linear plasmids

    "blast_hit":                 6,   # BLAST hit to known linear plasmid db (PLSDB etc.)
    "blast_hit_linear_db":       8,   # BLAST hit explicitly to a *linear* plasmid sequence
    "skani_hit":                 6,   # SKANI hit to known linear plasmid db
    "skani_hit_linear_db":       8,   # SKANI hit explicitly to a *linear* plasmid sequence
    "plasmid_finder_no_hit":     4,   # No PlasmidFinder hit (novel rep = typical of linear, H.2019)
    "coverage_drop_ends":        4,   # Coverage drop at contig ends (hairpin inaccessibility, Hashimoto 2019)
    "assembly_graph_linear":     6,   # Assembly graph linear topology
    "enterococcal_markers":      5,   # pELF-specific markers
    "copy_number_low":           2,   # ~1 copy/cell (characteristic)
    # ── FASTA header metadata (Unicycler / Plassembler / Hybracter output) ──
    "assembler_not_circular":   13,   # Explicit circular=false in header
    "header_copy_number":        4,   # Copy number in header consistent with plasmid (~0.3–4x)
}

# Scores that are ONLY valid when the annotation is contig-specific (not shared).
# These are zeroed out when annotation falls back to assembly-wide features.
CONTIG_SPECIFIC_SCORES = {
    "par_system", "repb_rep2", "ftsK_parA_repB_combo", "tas_system",
    "is_elements", "enterococcal_markers", "invertron_tp_gene",
}

# Contigs with circular=true in header are disqualified from gene-based scoring.
# Structural evidence (TIR, hairpin) still runs in case of mis-assembly,
# but the gene module returns empty to prevent annotation cross-contamination.
CIRCULAR_DISQUALIFIES_GENE_SCORING = True

# Confidence thresholds (against the 100-point achievable ceiling, see compute_score)
CONFIDENCE_THRESHOLDS = {
    "HIGH":   35,   # ≥35 points → likely linear plasmid
    "MEDIUM": 20,   # 20–34 points → possible linear plasmid
    "LOW":     8,   # 8–19 points → weak evidence
}

# Loop sequences flanked by IDR arms in hairpin telomeres.
# Primary: TATA (pELF1, Hashimoto 2019). Variants observed or expected:
#   TATATA / ATATAT — AT-repeat extension of the TATA motif
#   ATAT            — same dinucleotide, reverse orientation (detected after
#                     right-end reversal in detect_idr_tata_ends)
#   TTAA / AATT     — AT-rich palindromic alternatives
#   TAAT / ATTA     — TA-dinucleotide palindromes
#   TATAAT          — Pribnow-box-like extension; seen in distantly related
#                     linear replicons (Borrelia telomere variants)
# An empty loop (perfect fold-back hairpin, no spacer) is also accepted.
HAIRPIN_LOOP_MOTIFS = [
    "TATA",
    "TATATA",
    "ATAT",
    "ATATAT",
    "ATATA",
    "TATAT",
    "TTAA",
    "AATT",
    "TAAT",
    "ATTA",
    "TATAAT",
    "TATAAA",
]


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 0: FASTA HEADER METADATA PARSER
# ─────────────────────────────────────────────────────────────────────────────

def parse_fasta_header(description: str) -> dict:
    """
    Extract assembler-set metadata from common FASTA header formats.

    Handles Unicycler, Plassembler, Hybracter, and Trycycler headers, e.g.:
      >chromosome00001 len=2802778 circular=true
      >plasmid00002 length=107461 plasmid_copy_number_short=2.01x plasmid_copy_number_long=1.11x
      >contig_1 length=143316 depth=1.23x circular=false
      >NODE_1_length_143316_cov_34.5

    Returns dict with:
      circular         : bool | None   (None = not stated)
      copy_number      : float | None  (mean of short/long if both present)
      depth            : float | None  (Unicycler depth= tag)
      assembler_length : int | None
    """
    meta = {
        "circular":          None,
        "copy_number":       None,
        "depth":             None,
        "assembler_length":  None,
        "raw":               description,
    }

    desc = description.lower()

    # ── Circular flag ─────────────────────────────────────────────────────────
    if "circular=true" in desc:
        meta["circular"] = True
    elif "circular=false" in desc or "linear=true" in desc:
        meta["circular"] = False
    # No flag → None (ambiguous; do not penalise)

    # ── Length ────────────────────────────────────────────────────────────────
    for pat in [r"length=(\d+)", r"len=(\d+)"]:
        m = re.search(pat, desc)
        if m:
            meta["assembler_length"] = int(m.group(1))
            break
    # SPAdes NODE format: NODE_X_length_N_cov_D
    m = re.search(r"_length_(\d+)_cov_", desc)
    if m and meta["assembler_length"] is None:
        meta["assembler_length"] = int(m.group(1))

    # ── Copy number / depth ───────────────────────────────────────────────────
    # Plassembler/Hybracter: plasmid_copy_number_short=2.01x  plasmid_copy_number_long=1.11x
    cns = []
    for pat in [r"plasmid_copy_number_short=([\d.]+)x?",
                r"plasmid_copy_number_long=([\d.]+)x?",
                r"copy_number=([\d.]+)x?"]:
        m = re.search(pat, desc)
        if m:
            cns.append(float(m.group(1)))
    if cns:
        meta["copy_number"] = round(sum(cns) / len(cns), 2)

    # Unicycler: depth=34.50x
    m = re.search(r"depth=([\d.]+)x?", desc)
    if m:
        meta["depth"] = float(m.group(1))

    return meta


def assess_header_metadata(meta: dict, chromosome_depth: float = None,
                           seq_len: int = 0) -> dict:
    """
    Interpret header metadata for linear plasmid evidence.

    Assembler topology signal:
      assembler_not_circular  (weight 13): explicit circular=false in header
      circular=true           → disqualifies gene-based scoring; not a positive signal

    circular_flag_absent is still recorded (diagnostic / TSV output) but is NOT
    scored: it is true whenever no circular= flag is present at all (e.g. any
    non-Unicycler input, including plain reference FASTAs), so it fired on every
    single evidence run examined during 2026-07 revalidation regardless of the
    contig's true topology — it measures missing metadata, not linearity.

    header_copy_number (weight 4): CN 0.3–4x consistent with plasmid
    """
    circ = meta["circular"]
    result = {
        "circular_flag":          circ,
        "assembler_not_circular": circ is False,
        # Absent flag is positive evidence only for plasmid-sized contigs
        "circular_flag_absent":   circ is None and 0 < seq_len < 500_000,
        "is_circular":            circ is True,   # used to block gene scoring
        "copy_number":            meta["copy_number"],
        "depth":                  meta["depth"],
        "header_cn_consistent":   False,
        "high_copy_circular":     False,
    }

    cn = meta["copy_number"]
    if cn is None and meta["depth"] is not None and chromosome_depth:
        cn = round(meta["depth"] / chromosome_depth, 2)
        result["copy_number"] = cn

    if cn is not None:
        result["header_cn_consistent"] = 0.3 <= cn <= 4.0
        result["high_copy_circular"]   = cn > 10.0

    return result


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 0b: HYBRACTER OUTPUT DISCOVERY
# ─────────────────────────────────────────────────────────────────────────────
#
# Hybracter (gbouras13/hybracter) is a single-assembler-per-sample pipeline
# (Flye for long reads, Plassembler/Unicycler for small plasmids) — unlike
# Autocycler's multi-assembler ensemble, there's no per-cluster voting to do;
# each sample just has one final assembly. Confirmed directory layout (2026-07,
# via gbouras13/hybracter and gbouras13/plassembler source):
#   <dir>/FINAL_OUTPUT/{complete,incomplete}/{sample}_final.fasta
#   <dir>/FINAL_OUTPUT/{complete,incomplete}/{sample}_per_contig_stats.tsv
#     (complete only: columns contig_name, contig_type [chromosome|plasmid],
#      length, gc, circular [True/False]; incomplete has no contig_type/circular
#      — an incomplete assembly carries no circularity signal at all)
#   <dir>/processing/assemblies/{sample}/assembly_graph.gfa   (Flye)
#   <dir>/processing/plassembler/{sample}/**/*.gfa            (best-effort —
#     plassembler runs Unicycler internally and the nested output dir name
#     isn't stable across versions, so this is a recursive glob)
#   <dir>/processing/qc/{sample}_filt_trim.fastq.gz           (Hybracter's own
#     QC'd long reads — reusable for auto-BAM-mapping)
#
# Note: Hybracter's complete/{sample}_final.fasta already has circular=true/
# circular=True embedded in the FASTA description for chromosome/plasmid
# records (written by plassembler's select_best_lib.py), which the existing
# parse_fasta_header/assess_header_metadata already reads correctly (it
# lowercases before matching) — so no new scoring pathway is needed for the
# circular flag. parse_hybracter_contig_stats() below is used only for
# reliable chromosome/plasmid separation (--chromosome-contigs auto-population)
# and as a diagnostic cross-check, not as a scoring input.

def discover_hybracter_samples(hybracter_dir: str) -> list:
    """
    Discover per-sample assembly outputs under a Hybracter output directory.

    Returns a list of dicts: {sample, complete, fasta, contig_stats_tsv,
    flye_gfa, plassembler_gfa, longread_fastq} — the latter four are None
    when not found.
    """
    base = Path(hybracter_dir)
    samples = []
    for status in ("complete", "incomplete"):
        final_dir = base / "FINAL_OUTPUT" / status
        if not final_dir.is_dir():
            continue
        for fasta in sorted(final_dir.glob("*_final.fasta")):
            sample = fasta.name[: -len("_final.fasta")]

            contig_stats = final_dir / f"{sample}_per_contig_stats.tsv"

            flye_gfa = base / "processing" / "assemblies" / sample / "assembly_graph.gfa"

            plassembler_dir = base / "processing" / "plassembler" / sample
            plassembler_gfas = (sorted(plassembler_dir.glob("**/*.gfa"))
                                if plassembler_dir.is_dir() else [])

            longread_fastq = base / "processing" / "qc" / f"{sample}_filt_trim.fastq.gz"

            samples.append({
                "sample":           sample,
                "complete":         status == "complete",
                "fasta":            str(fasta),
                "contig_stats_tsv": str(contig_stats) if contig_stats.exists() else None,
                "flye_gfa":         str(flye_gfa) if flye_gfa.exists() else None,
                "plassembler_gfa":  str(plassembler_gfas[0]) if plassembler_gfas else None,
                "longread_fastq":   str(longread_fastq) if longread_fastq.exists() else None,
            })
    return samples


def parse_hybracter_contig_stats(tsv_path: str) -> dict:
    """
    Parse a Hybracter {sample}_per_contig_stats.tsv.

    Returns {contig_name: {"contig_type": "chromosome"|"plasmid"|"", "length":
    int|None, "gc": float|None, "circular": bool}}. Missing/absent file
    returns {}.
    """
    if not tsv_path or not os.path.exists(tsv_path):
        return {}
    df = pd.read_csv(tsv_path, sep="\t")
    out = {}
    for _, row in df.iterrows():
        circ_raw = str(row.get("circular", "")).strip().lower()
        out[row["contig_name"]] = {
            "contig_type": row.get("contig_type", "") if pd.notna(row.get("contig_type", "")) else "",
            "length":      int(row["length"]) if pd.notna(row.get("length")) else None,
            "gc":          float(row["gc"]) if pd.notna(row.get("gc")) else None,
            "circular":    circ_raw == "true",
        }
    return out


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 2: TERMINAL HAIRPIN DETECTION (localized per-end fold-back)
# ─────────────────────────────────────────────────────────────────────────────
#
# Ported from linear-plasmid-hairpin-tools' find_hairpins.py (analyse_end),
# 2026-07. Replaces the previous detect_self_complement_ends, which compared
# the WHOLE first/last `window` bp to each other (testing whether the two ends
# are mutual reverse-complements — a *symmetric* hairpin telomere hypothesis).
# pELF1-type plasmids have *asymmetric* ends (one hairpin, one invertron —
# Hashimoto 2019) and that whole-window compare scored only 30.8% identity
# (chance level) on the literal pELF1 reference sequence: it was testing the
# wrong hypothesis. This detector instead looks for a LOCALIZED fold-back
# within a terminal window, independently at each end — it finds the
# inverted-repeat symmetry axis and checks whether the outer arm reaches the
# contig terminus, which is the actual structural signature of a hairpin
# telomere and works correctly regardless of what the opposite end looks like.

_HAIRPIN_COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def _hairpin_rev_comp(s: bytes) -> bytes:
    return s.translate(_HAIRPIN_COMP)[::-1]


def _analyse_hairpin_end(seq: bytes, k: int, window: int, end: str, edge_tol: int,
                         max_kmer_hits: int = 200):
    """
    Look for a reverse-complement inverted repeat in the terminal window at the
    given end ('5' or '3'). Locate the dominant symmetry axis and measure the
    outer arm.

    Returns None or dict:
      support  : number of supporting RC k-mer pairs at the dominant axis
      tir_span : bp spanned by the inverted repeat (both arms + loop)
      touches  : True if the outer arm reaches the contig terminus (<= edge_tol)
      gap      : bp from the terminus to the outer arm (0 = flush with the end)
      arm1/arm2: (start, end) global 0-based half-open coords of each arm
      axis     : global 0-based symmetry-axis position
    """
    L = len(seq)
    if L < k:
        return None
    if end == "3":
        off = max(0, L - window)
        w = seq[off:]
    else:  # "5"
        off = 0
        w = seq[:min(window, L)]
    n = len(w)
    if n < k:
        return None

    pos = defaultdict(list)
    for i in range(n - k + 1):
        pos[w[i:i + k]].append(i)

    # centre (axis) -> list of (i, j) supporting pairs
    centres = defaultdict(list)
    for i in range(n - k + 1):
        js = pos.get(_hairpin_rev_comp(w[i:i + k]))
        if not js or len(js) > max_kmer_hits:
            continue
        for j in js:
            if j < i:
                continue
            centres[(i + j + k - 1) // 2].append((i, j))
    if not centres:
        return None

    c0 = max(centres, key=lambda c: len(centres[c]))
    pairs = []
    for c in (c0 - 1, c0, c0 + 1):          # allow ±1 for rounding
        pairs += centres.get(c, [])
    positions = set()
    for i, j in pairs:
        positions.add(i)
        positions.add(j)
    support = len(pairs)
    minp, maxp = min(positions), max(positions)
    tir_span = (maxp - minp) + k

    lefts = [i for i, j in pairs]
    rights = [j for i, j in pairs]
    # arm1 = 5'-side arm, arm2 = 3'-side arm (global 0-based half-open)
    arm1 = (off + min(lefts), off + max(lefts) + k)
    arm2 = (off + min(rights), off + max(rights) + k)

    outer_global = off + maxp + k          # 3' edge of outermost k-mer
    inner_global = off + minp              # 5' edge of innermost k-mer
    if end == "3":
        gap = max(0, L - outer_global)
        touches = gap <= edge_tol
    else:
        gap = max(0, inner_global)
        touches = gap <= edge_tol
    return {"support": support, "tir_span": tir_span,
            "touches": touches, "gap": gap,
            "arm1": arm1, "arm2": arm2, "axis": off + c0}


def detect_terminal_hairpins(seq: str, k: int = 31, window: int = 50_000,
                             min_shared: int = 25, edge_tol: int = 100) -> dict:
    """
    Detect a terminal hairpin/TIR fold-back independently at each end of `seq`.

    Returns {"left": {...} | None, "right": {...} | None} — an end's entry is
    populated only if a fold-back with >= min_shared supporting k-mer pairs is
    found AND it reaches the contig terminus (within edge_tol bp).
    """
    seq_b = seq.upper().encode("ascii", "replace")
    out = {"left": None, "right": None}
    for end, key in (("5", "left"), ("3", "right")):
        r = _analyse_hairpin_end(seq_b, k, window, end, edge_tol)
        if r is not None and r["support"] >= min_shared and r["touches"]:
            out[key] = r
    return out


def detect_asymmetric_ends(hairpin_ends: dict,
                           gene_hits: dict,
                           cov_result: dict = None) -> dict:
    """
    pELF1 has ASYMMETRIC ends: one end is a hairpin, the other an invertron
    (terminal protein). This is a unique indicator of the pELF lineage
    (Hashimoto 2019).

    Two detection paths — either is sufficient:

    Sequence + gene:
      - Terminal hairpin fold-back at exactly one end (detect_terminal_hairpins)
        — not both, which would be a symmetric telomere type instead — AND
      - Terminal protein gene annotated

    BAM coverage (most diagnostic with Illumina reads):
      - Left end coverage drop (hairpin blocks adapter ligation) AND
      - Right end coverage NOT dropped (invertron end remains accessible)
      A symmetric drop at both ends does NOT count — that could be any palindromic
      structure, not the pELF1-specific asymmetry.
    """
    left_hairpin  = hairpin_ends.get("left") is not None
    right_hairpin = hairpin_ends.get("right") is not None
    has_hairpin = left_hairpin or right_hairpin
    has_tp_gene = bool(gene_hits.get("terminal_protein_genes"))
    seq_asymmetric = (left_hairpin != right_hairpin) and has_tp_gene

    cov = cov_result or {}
    bam_asymmetric = (cov.get("available") and
                      cov.get("left_drop") and
                      not cov.get("right_drop"))

    return {
        "asymmetric_pelf_type": seq_asymmetric or bam_asymmetric,
        "has_hairpin_end":      has_hairpin,
        "has_tp_gene":          has_tp_gene,
        "bam_asymmetric":       bam_asymmetric,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 2b: IDR / TATA HAIRPIN TELOMERE DETECTION
# ─────────────────────────────────────────────────────────────────────────────

def detect_idr_tata_ends(seq: str, search_window: int = 150,
                          min_arm: int = 6, max_arm: int = 40,
                          max_loop: int = 16, min_identity: float = 0.80) -> dict:
    """
    Search for an IDR-loop-IDR hairpin telomere structure at each sequence end.

    Structure searched in the terminal window: 5'--[arm]--[loop]--[RC(arm)]--3'

    For the left end, the window is seq[:search_window].
    For the right end, the window is seq[-search_window:] reversed (simple
    reversal, not RC) so the 3' tip is at index 0; this preserves the palindrome
    property — if the right end is RC(arm)--loop--arm, reversing gives
    arm_rev--loop_rev--RC(arm)_rev and the RC test arm_rev ≈ RC(RC(arm)_rev)
    simplifies to arm_rev ≈ complement(arm), which passes for AT-symmetric arms.
    TATA reversed is ATAT, which is included in HAIRPIN_LOOP_MOTIFS.

    A loop that is:
      - present in HAIRPIN_LOOP_MOTIFS, OR
      - ≥ 75% AT-rich, OR
      - empty (perfect fold-back hairpin)
    is considered a confirmed IDR/TATA motif.  An IDR without such a loop is
    still reported but does NOT count as "confirmed".

    Returns:
      left  : scan results for left end
      right : scan results for right end
      either_end     : True if IDR found at either end
      confirmed_idr_tata : True if IDR + valid loop found at either end
    """
    if len(seq) < 2 * min_arm:
        return {"left": {"found": False}, "right": {"found": False},
                "either_end": False, "confirmed_idr_tata": False}

    def _scan_tip(tip: str) -> dict:
        best = {"found": False, "arm_len": 0, "loop_seq": "", "identity": 0.0,
                "loop_motif_match": False, "at_rich_loop": False,
                "confirmed": False}
        tip_len = len(tip)
        for arm_len in range(min_arm, min(max_arm + 1,
                                          (tip_len - 0) // 2 + 1)):
            for loop_len in range(0, max_loop + 1):
                total = arm_len * 2 + loop_len
                if total > tip_len:
                    break
                arm1 = tip[:arm_len]
                loop = tip[arm_len:arm_len + loop_len]
                arm2 = tip[arm_len + loop_len:arm_len + loop_len + arm_len]
                rc_arm2 = str(Seq(arm2).reverse_complement())
                identity = sum(a == b for a, b in zip(arm1, rc_arm2)) / arm_len
                if identity < min_identity:
                    continue
                loop_upper = loop.upper()
                loop_motif = (not loop) or any(
                    m in loop_upper or loop_upper in m
                    for m in HAIRPIN_LOOP_MOTIFS
                )
                at_rich = (not loop) or (
                    sum(c in "AT" for c in loop_upper) / len(loop) >= 0.75
                )
                confirmed = loop_motif or at_rich
                # Keep best by: confirmed first, then identity
                if identity > best["identity"] or (
                        confirmed and not best["confirmed"] and identity >= min_identity):
                    best = {
                        "found":            True,
                        "arm_len":          arm_len,
                        "loop_seq":         loop,
                        "identity":         round(identity, 3),
                        "loop_motif_match": loop_motif,
                        "at_rich_loop":     at_rich,
                        "confirmed":        confirmed,
                    }
        return best

    left_tip  = seq[:search_window].upper()
    right_tip = seq[-search_window:].upper()[::-1]   # simple reversal; see docstring

    left_result  = _scan_tip(left_tip)
    right_result = _scan_tip(right_tip)

    return {
        "left":               left_result,
        "right":              right_result,
        "either_end":         left_result["found"] or right_result["found"],
        "confirmed_idr_tata": left_result.get("confirmed") or right_result.get("confirmed"),
    }


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 2c: COVERAGE DROP AT ENDS  (Hashimoto 2019)
# ─────────────────────────────────────────────────────────────────────────────

def detect_coverage_drop_ends(bam_file, contig_id: str,
                               contig_len: int, end_window: int = 5000,
                               tip_window: int = 10) -> dict:
    """
    Detect coverage drop at contig ends, consistent with linear plasmid terminal structure.

    Two complementary checks are run:

    Wide window (<end_window> bp, default 5 kb), threshold < 50% of body depth:
      Fires when the assembled hairpin arm is present in the contig and Illumina
      reads cannot be generated from it (adapter ligation blocked by the hairpin
      fold, Hashimoto 2019, Fig 3A). Also fires with ONT reads when mappability
      is reduced across the palindromic terminal region.

    Narrow tip window (<tip_window> bp, default 10 bp), threshold < 10% of body depth:
      Fires when the assembly was TRUNCATED at the palindrome centre — the hairpin
      arm was never assembled, so the contig starts at the hairpin tip (position 1
      = TATA / palindrome centre). In that case the very first few bases are only
      reachable by reads genuinely starting there; no reads from "before" position 1
      exist (the hairpin prevented their generation). Circular assemblies are not
      confounded: soft-clipped reads from near the genome end keep position-1 depth
      close to body depth, well above the 10% threshold.
    """
    try:
        import pysam
    except ImportError:
        return {"available": False, "message": "pysam not installed"}

    close_when_done = False
    if isinstance(bam_file, str):
        if not bam_file or not os.path.exists(bam_file):
            return {"available": False, "message": "BAM not found"}
        try:
            bam_file = pysam.AlignmentFile(bam_file, "rb")
            close_when_done = True
        except Exception as e:
            return {"available": False, "message": str(e)}

    try:
        bam = bam_file

        def region_depths(start, end):
            return [col.nsegments
                    for col in bam.pileup(contig_id, start, end,
                                          min_mapping_quality=20)]

        body_start = end_window
        body_end   = max(contig_len - end_window, end_window + 1)

        left_depths  = region_depths(0, min(end_window, contig_len))
        right_depths = region_depths(max(0, contig_len - end_window), contig_len)
        body_depths  = region_depths(body_start, body_end)

        avg_left  = mean(left_depths)  if left_depths  else 0
        avg_right = mean(right_depths) if right_depths else 0
        avg_body  = mean(body_depths)  if body_depths  else 1

        left_ratio  = avg_left  / avg_body if avg_body > 0 else 1
        right_ratio = avg_right / avg_body if avg_body > 0 else 1

        # Narrow tip check: very first / last <tip_window> bp
        tip_left_depths  = region_depths(0, min(tip_window, contig_len))
        tip_right_depths = region_depths(max(0, contig_len - tip_window), contig_len)
        avg_tip_left  = mean(tip_left_depths)  if tip_left_depths  else 0
        avg_tip_right = mean(tip_right_depths) if tip_right_depths else 0
        tip_left_ratio  = avg_tip_left  / avg_body if avg_body > 0 else 1
        tip_right_ratio = avg_tip_right / avg_body if avg_body > 0 else 1
        tip_left_drop  = tip_left_ratio  < 0.10
        tip_right_drop = tip_right_ratio < 0.10

        result = {
            "available":          True,
            "left_depth_ratio":   round(left_ratio,  3),
            "right_depth_ratio":  round(right_ratio, 3),
            "body_mean_depth":    round(avg_body, 1),
            "left_drop":          left_ratio  < 0.50,
            "right_drop":         right_ratio < 0.50,
            "left_tip_ratio":     round(tip_left_ratio,  3),
            "right_tip_ratio":    round(tip_right_ratio, 3),
            "left_tip_drop":      tip_left_drop,
            "right_tip_drop":     tip_right_drop,
            "consistent_with_linear": (
                left_ratio  < 0.50 or right_ratio  < 0.50 or
                tip_left_drop       or tip_right_drop
            ),
        }
        if close_when_done:
            bam.close()
        return result
    except Exception as e:
        if close_when_done:
            try:
                bam_file.close()
            except Exception:
                pass
        return {"available": False, "message": str(e)}


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 3: GC CONTENT ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def gc_content(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    total = len(seq) - seq.count("N")
    return round(100 * gc / total, 2) if total > 0 else 0.0


def assess_gc(seq_gc: float, ref_gc: float = None) -> dict:
    """
    Linear plasmids typically have lower GC% than the host chromosome
    (Hashimoto 2023). If a reference GC is provided, check deviation;
    otherwise use species-agnostic heuristic ranges.
    """
    result = {"gc": seq_gc, "deviation": None, "low_gc_flag": False}
    if ref_gc is not None:
        deviation = ref_gc - seq_gc
        result["deviation"] = round(deviation, 2)
        result["low_gc_flag"] = deviation >= 3.0   # ≥3 % drop is significant
    else:
        # Without a reference, flag GC outside typical bacterial range but
        # not definitively diagnostic
        result["low_gc_flag"] = seq_gc < 40.0
    return result


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 4: SIZE CLASSIFICATION
# ─────────────────────────────────────────────────────────────────────────────

def classify_size(length: int) -> dict:
    """Return size category and whether it falls in a known linear plasmid range."""
    families = []
    for family, (lo, hi) in SIZE_RANGES.items():
        if lo <= length <= hi:
            families.append(family)
    return {
        "length_bp": length,
        "in_known_range": len(families) > 0,
        "matching_families": families,
        "size_kb": round(length / 1000, 1),
    }


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 5: GENE CONTENT ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def parse_annotation(annot_file: str) -> pd.DataFrame:
    """
    Parse an annotation file in several common formats:

      1. Prokka TSV  (.tsv)  — columns: locus_tag, ftype, length_bp, gene,
                               EC_number, COG, product
                               (contig ID is the locus_tag prefix before the
                               last underscore+digits, e.g. 'AP022343_00001'
                               → contig 'AP022343')
      2. GFF3        (.gff)  — 9-column tab-separated; gene/product in
                               attributes field
      3. Simple TSV          — any tab-separated file that already has 'gene'
                               and 'product' columns; 'contig' optional

    Always returns a DataFrame with at least 'contig', 'gene', 'product' columns.
    If 'contig' cannot be determined, a blank string is used (all features are
    then matched against every contig, which is safe for small assemblies).
    """
    if not annot_file or not os.path.exists(annot_file):
        return pd.DataFrame()

    # ── peek at the first non-comment line to detect format ──────────────────
    header_line = ""
    first_data  = ""
    with open(annot_file, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line.strip():
                continue
            if not header_line:
                header_line = line
            else:
                first_data = line
                break

    cols_lower = [c.strip().lower() for c in header_line.split("\t")]

    # ── Format 1: Prokka TSV (has 'locus_tag' header) ────────────────────────
    if "locus_tag" in cols_lower:
        try:
            df = pd.read_csv(annot_file, sep="\t", comment="#", dtype=str)
            df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
            # Derive contig from locus_tag: strip trailing _NNNNN
            if "locus_tag" in df.columns:
                df["contig"] = df["locus_tag"].str.replace(
                    r"_\d+$", "", regex=True)
            else:
                df["contig"] = ""
            for col in ("gene", "product"):
                if col not in df.columns:
                    df[col] = ""
            df["gene"]    = df["gene"].fillna("").astype(str)
            df["product"] = df["product"].fillna("").astype(str)
            return df[["contig", "gene", "product"]]
        except Exception as e:
            print(f"[WARN] Prokka TSV parse error: {e}", file=sys.stderr)
            return pd.DataFrame()

    # ── Format 2: GFF3 (9 columns, no readable header) ───────────────────────
    # Detected by: first data line has 9 tab-separated fields and field[8]
    # contains '=' (GFF attributes)
    first_fields = first_data.split("\t") if first_data else []
    if len(first_fields) == 9 and "=" in first_fields[8]:
        try:
            df = pd.read_csv(annot_file, sep="\t", comment="#",
                             header=None,
                             names=["contig", "source", "feature", "start",
                                    "end", "score", "strand", "phase",
                                    "attributes"],
                             dtype=str)
            df = df[df["feature"].isin(["CDS", "gene", "rRNA", "tRNA"])].copy()
            df["gene"]    = df["attributes"].str.extract(
                r"gene=([^;]+)", expand=False).fillna("")
            df["product"] = df["attributes"].str.extract(
                r"product=([^;]+)", expand=False).fillna("")
            df["contig"]  = df["contig"].fillna("").astype(str)
            return df[["contig", "gene", "product"]]
        except Exception as e:
            print(f"[WARN] GFF3 parse error: {e}", file=sys.stderr)
            return pd.DataFrame()

    # ── Format 3: Generic TSV with header ────────────────────────────────────
    try:
        df = pd.read_csv(annot_file, sep="\t", comment="#", dtype=str)
        df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
        if "contig"  not in df.columns: df["contig"]  = ""
        if "gene"    not in df.columns: df["gene"]    = ""
        if "product" not in df.columns: df["product"] = ""
        df["gene"]    = df["gene"].fillna("").astype(str)
        df["product"] = df["product"].fillna("").astype(str)
        df["contig"]  = df["contig"].fillna("").astype(str)
        return df[["contig", "gene", "product"]]
    except Exception as e:
        print(f"[WARN] Generic TSV parse error: {e}", file=sys.stderr)
        return pd.DataFrame()


def screen_genes(annot_df: pd.DataFrame, contig_id: str) -> dict:
    """
    Screen annotation for genes associated with linear plasmids.
    Returns counts per category and named hits.

    Contig filtering: tries to match contig_id against the 'contig' column.
    If the column is missing, blank, or yields no rows (e.g. Prokka TSV where
    locus_tag prefix doesn't match FASTA header), falls back to using ALL
    features — safe for single-replicon files and gives a conservative over-
    estimate for multi-replicon files.
    """
    if annot_df.empty:
        return {}

    # Ensure required columns exist
    for col in ("contig", "gene", "product"):
        if col not in annot_df.columns:
            annot_df = annot_df.copy()
            annot_df[col] = ""

    # Try contig-specific subset first
    contig_col = annot_df["contig"].astype(str)
    n_unique_contigs = contig_col.nunique()
    subset = annot_df[contig_col.str.contains(
        re.escape(str(contig_id)), na=False, case=False)]

    # Fallback logic:
    # - Single-contig annotation (Prokka on one sequence): safe to use all
    # - Multi-contig annotation with no match: Prokka TSV prefix ≠ FASTA header.
    #   Using all genes would contaminate every contig with the whole assembly's
    #   gene set — INCORRECT. Return empty and let compute_score handle it.
    if subset.empty and not annot_df.empty:
        if n_unique_contigs <= 1:
            # Single contig in annotation — safe fallback
            subset = annot_df
        else:
            # Multi-contig: cannot assign genes without per-contig mapping
            print(f"[WARN] No annotation match for '{contig_id}' "
                  f"(annotation has {n_unique_contigs} distinct contig IDs). "
                  f"Provide a GFF3 file (--annot) for per-contig gene mapping. "
                  f"Gene scoring skipped for this contig.", file=sys.stderr)
            return {"_annotation_mismatch": True}

    hits = defaultdict(list)
    all_genes    = subset["gene"].fillna("").astype(str).str.lower().tolist()
    all_products = subset["product"].fillna("").astype(str).str.lower().tolist()
    all_text     = all_genes + all_products

    def _word_match(keyword, texts):
        """True if keyword appears as a whole word (case-insensitive)."""
        pat = r'(?<![a-z0-9])' + re.escape(keyword.lower()) + r'(?![a-z0-9])'
        return any(re.search(pat, t) for t in texts)

    # Linear plasmid gene keywords — word-boundary to prevent substring FP
    # (e.g. "rep" matching "separate", "para" matching "paracoccus")
    for kw in LINEAR_PLASMID_GENE_KEYWORDS:
        if _word_match(kw, all_text):
            hits["linear_plasmid_keywords"].append(kw)

    # Toxin-antitoxin systems
    for system, genes in TAS_GENES.items():
        found = [g for g in genes if _word_match(g, all_text)]
        if len(found) >= 1:
            hits["tas_systems"].append(system)

    # Partition systems
    for par, pfams in PAR_PFAM.items():
        if _word_match(par, all_text):
            hits["partition_genes"].append(par)

    # Resistance genes
    for gene, description in RESISTANCE_GENES.items():
        if _word_match(gene, all_text):
            hits["resistance_genes"].append(f"{gene} ({description})")

    # IS elements — substring OK: IS1216 intentionally matches IS1216E/V variants
    for is_elem in LINEAR_PLASMID_IS:
        if any(is_elem.lower() in t for t in all_text):
            hits["is_elements"].append(is_elem)

    # Enterococcal pELF markers (Hashimoto 2019 / Boumamoud 2022)
    enterococcal_m = ["is1216", "ftsK", "soj", "ftsk", "vana", "vanb"]
    hits["enterococcal_markers"] = [m for m in enterococcal_m
                                    if _word_match(m, all_text)]

    # Terminal protein genes — word-boundary to avoid FP from "C-terminal"
    hits["terminal_protein_genes"] = [
        kw for kw in INVERTRON_TP_GENE_KEYWORDS
        if _word_match(kw, all_text)
    ]

    # repB / Rep_2 superfamily (Hashimoto 2019: pELF1 replication initiation)
    rep2_keywords = ["repb", "rep_2", "rep2", "replication initiation"]
    hits["repB_rep2"] = [kw for kw in rep2_keywords
                         if _word_match(kw, all_text)]

    # ftsK + parA + repB co-occurrence signature (Hashimoto 2019 confirmed trio)
    trio = {"ftsk": _word_match("ftsk", all_text),
            "para": _word_match("parA", all_text),
            "repb": _word_match("repB", all_text)}
    hits["ftsK_parA_repB_trio"] = [k for k, v in trio.items() if v]

    # IS1216E flanking in same direction (vanM cluster structure, Hashimoto 2019)
    is_variants = [e for e in ["IS1216E", "IS1216V", "IS1216"]
                   if any(e.lower() in t for t in all_text)]
    hits["is1216_variants"] = is_variants

    return {k: list(set(v)) for k, v in hits.items()}


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 5b: AMRFINDERPLUS SCREENING (antimicrobial resistance genes)
# ─────────────────────────────────────────────────────────────────────────────
#
# Reports AMR gene content per contig using NCBI's AMRFinderPlus (`amrfinder`
# CLI + its curated reference database — bioconda: ncbi-amrfinderplus).
# Opt-in via --amrfinder; run once on the whole assembly (nucleotide mode,
# -n) and attributed back to each contig by AMRFinderPlus's own "Contig id"
# column, rather than re-run per contig.
#
# Informational only — NOT scored. This follows the same precedent as the
# keyword-based `resistance_genes` category (commit 863d4b2, "Remove
# antibiotic resistance genes from scoring"): AMR content doesn't discriminate
# linear vs circular topology, so it's reported in the TSV/JSON but does not
# contribute to SCORING_WEIGHTS/compute_score.

# Column-name variants seen across AMRFinderPlus versions — resolved
# case-insensitively by substring so minor header changes don't break parsing.
_AMR_COLUMN_ALIASES = {
    "contig":   ["contig id", "contig identifier"],
    "gene":     ["gene symbol", "element symbol"],
    "product":  ["sequence name", "element name"],
    "type":     ["element type", "type"],
    "subtype":  ["element subtype", "subtype"],
    "class":    ["class"],
    "subclass": ["subclass"],
    "identity": ["% identity to reference sequence", "% identity to reference"],
    "coverage": ["% coverage of reference sequence", "% coverage of reference"],
    "start":    ["start"],
    "stop":     ["stop"],
}


def _resolve_amr_columns(columns) -> dict:
    """Map our canonical field names to whatever AMRFinderPlus actually named
    the column in this run (case-insensitive exact match against known
    aliases). Missing fields are simply absent from the returned dict."""
    lower_cols = {c.lower(): c for c in columns}
    resolved = {}
    for field, aliases in _AMR_COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in lower_cols:
                resolved[field] = lower_cols[alias]
                break
    return resolved


def run_amrfinder(assembly_fasta: str, out_file: str, db: str = None,
                  organism: str = None, threads: int = 4) -> pd.DataFrame:
    """
    Run AMRFinderPlus in nucleotide mode on the whole assembly FASTA.
    Returns a DataFrame with canonical columns (contig, gene, product, type,
    subtype, class, subclass, identity, coverage, start, stop) — a subset may
    be missing if AMRFinderPlus's output didn't include them. Empty DataFrame
    if the tool isn't on PATH, the run fails, or there are no hits.
    """
    if not os.path.exists(assembly_fasta):
        return pd.DataFrame()

    cmd = ["amrfinder", "-n", assembly_fasta, "-o", out_file,
           "--threads", str(threads)]
    if db:
        cmd += ["-d", db]
    if organism:
        cmd += ["-O", organism]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"[WARN] amrfinder failed or not found on PATH ({e}); "
              f"skipping AMR screening.", file=sys.stderr)
        return pd.DataFrame()

    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        return pd.DataFrame()

    try:
        raw = pd.read_csv(out_file, sep="\t")
    except Exception:
        return pd.DataFrame()
    if raw.empty:
        return pd.DataFrame()

    cols = _resolve_amr_columns(raw.columns)
    if "contig" not in cols:
        print("[WARN] amrfinder output missing a recognisable Contig id "
              "column; skipping AMR screening.", file=sys.stderr)
        return pd.DataFrame()

    df = pd.DataFrame({field: raw[col] for field, col in cols.items()})
    return df


def interpret_amr_hits(amr_df: pd.DataFrame, contig_id: str) -> dict:
    """Summarise AMRFinderPlus hits for one contig."""
    if amr_df is None or amr_df.empty:
        return {"available": amr_df is not None, "hits": 0,
                "genes": [], "classes": [], "elements": []}

    subset = amr_df[amr_df["contig"].astype(str) == str(contig_id)]
    if subset.empty:
        return {"available": True, "hits": 0, "genes": [], "classes": [],
                "elements": []}

    elements = []
    for _, row in subset.iterrows():
        elements.append({
            "gene":     row.get("gene", ""),
            "product":  row.get("product", ""),
            "type":     row.get("type", ""),
            "subtype":  row.get("subtype", ""),
            "class":    row.get("class", ""),
            "subclass": row.get("subclass", ""),
            "identity": row.get("identity", None),
            "coverage": row.get("coverage", None),
            "start":    row.get("start", None),
            "stop":     row.get("stop", None),
        })

    return {
        "available": True,
        "hits":      len(subset),
        "genes":     sorted(subset["gene"].dropna().astype(str).unique().tolist())
                    if "gene" in subset else [],
        "classes":   sorted(subset["class"].dropna().astype(str).unique().tolist())
                    if "class" in subset else [],
        "elements":  elements,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 6: COPY NUMBER ESTIMATION (from BAM coverage)
# ─────────────────────────────────────────────────────────────────────────────

def estimate_copy_number(bam_file, contig_id: str,
                         chromosome_contigs: list = None) -> dict:
    """
    Estimate copy number of a contig relative to the chromosome by comparing
    mean read depths (~1 copy/cell characteristic of linear plasmids).

    bam_file: path string OR open pysam.AlignmentFile (shared handle).
    Requires pysam. Falls back gracefully if not installed.
    """
    try:
        import pysam
    except ImportError:
        return {"available": False, "message": "pysam not installed"}

    close_when_done = False
    if isinstance(bam_file, str):
        if not bam_file or not os.path.exists(bam_file):
            return {"available": False, "message": "BAM file not found"}
        try:
            bam_file = pysam.AlignmentFile(bam_file, "rb")
            close_when_done = True
        except Exception as e:
            return {"available": False, "message": str(e)}

    try:
        bam = bam_file

        def mean_depth(ctg):
            depths = [col.nsegments for col in bam.pileup(ctg, min_mapping_quality=20)]
            return np.mean(depths) if depths else 0

        plasmid_depth = mean_depth(contig_id)

        chrom_depth = None
        if chromosome_contigs:
            depths = [mean_depth(c) for c in chromosome_contigs]
            chrom_depth = np.mean(depths) if depths else None

        copy_number = None
        if chrom_depth and chrom_depth > 0:
            copy_number = round(plasmid_depth / chrom_depth, 2)

        if close_when_done:
            bam.close()
        return {
            "available": True,
            "plasmid_depth": round(plasmid_depth, 1),
            "chromosome_depth": round(chrom_depth, 1) if chrom_depth else None,
            "estimated_copy_number": copy_number,
            "consistent_with_linear": (copy_number is not None and 0.5 <= copy_number <= 3.0),
        }
    except Exception as e:
        if close_when_done:
            try:
                bam_file.close()
            except Exception:
                pass
        return {"available": False, "message": str(e)}


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 7: ASSEMBLY GRAPH TOPOLOGY (GFA parsing)
# ─────────────────────────────────────────────────────────────────────────────

# Strand-aware GFA graph classifier, ported from linear-plasmid-hairpin-tools'
# autocycler_dotplot_classify.py (2026-07). Not actually Autocycler-specific —
# it only needs standard S/L lines, so it works on Flye's assembly_graph.gfa
# and Unicycler/Plassembler's assembly.gfa too. Replaces the previous
# link-degree-counting parser, which had a real bug: a genuinely circular
# single-contig self-loop (same-strand self-link, `L 1 + 1 +`) has 2
# link-table entries with 1 partner, landing in the same "tir_like" bucket as
# a true open-ended terminal-inverted-repeat linear contig — i.e. it could
# misclassify a circular replicon as linear evidence. The component-walk
# below distinguishes a closed circular loop from a hairpin self-link
# (opposite-strand, `L 1 + 1 -`) from a genuinely open/fragmented graph by
# actually walking the strand-aware links, not just counting them.

def _gfa_parse_segments_links(gfa_file: str):
    """Parse S/L lines only. Returns (segments: {id: length}, links: set of
    (from_seg, from_strand, to_seg, to_strand))."""
    segments = {}
    links = set()
    with open(gfa_file) as fh:
        for line in fh:
            if not line:
                continue
            t = line[0]
            if t not in "SL":
                continue
            f = line.rstrip("\n").split("\t")
            if t == "S" and len(f) >= 3:
                seg = f[1]
                seq = f[2]
                length = 0
                if seq and seq != "*":
                    length = len(seq)
                else:
                    for tag in f[3:]:
                        if tag.startswith("LN:i:"):
                            try:
                                length = int(tag[5:])
                            except ValueError:
                                pass
                segments[seg] = length
            elif t == "L" and len(f) >= 5:
                links.add((f[1], f[2], f[3], f[4]))
    for (a, _, b, _) in links:
        segments.setdefault(a, 0)
        segments.setdefault(b, 0)
    return segments, links


def _gfa_build_adjacency(links):
    fnext = defaultdict(list)   # leaving u on '+'  -> list of (to, to_strand)
    rnext = defaultdict(list)   # leaving u on '-'  -> list of (to, to_strand)
    fprev = defaultdict(int)    # links arriving at u on '+'
    rprev = defaultdict(int)    # links arriving at u on '-'
    undirected = defaultdict(set)
    for (a, sa, b, sb) in links:
        (fnext if sa == "+" else rnext)[a].append((b, sb))
        if sb == "+":
            fprev[b] += 1
        else:
            rprev[b] += 1
        undirected[a].add(b)
        undirected[b].add(a)
    return fnext, rnext, fprev, rprev, undirected


def _gfa_connected_components(segments, undirected):
    visited, comps = set(), []
    for seg in segments:
        if seg in visited:
            continue
        stack, comp = [seg], []
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            comp.append(cur)
            for nb in undirected.get(cur, ()):
                if nb not in visited:
                    stack.append(nb)
        comps.append(sorted(comp))
    return sorted(comps)


def _gfa_component_is_circular_loop(component, fnext, rnext, fprev, rprev) -> bool:
    """Start at the lowest-numbered unitig on the forward strand and walk
    forward links; every unitig on the loop must have exactly one link on
    each of its four sides, and the walk must return to the start after
    covering the whole component."""
    if not component:
        return False
    first = component[0]
    num, strand = first, "+"
    visited = set()
    while num != first or not visited:
        if num in visited:
            return False
        visited.add(num)
        if (len(fnext.get(num, [])) != 1 or len(rnext.get(num, [])) != 1 or
                fprev.get(num, 0) != 1 or rprev.get(num, 0) != 1):
            return False
        nxt = fnext[num][0] if strand == "+" else rnext[num][0]
        num, strand = nxt
    return len(visited) == len(component)


def _gfa_dominant_segment_is_open(comp, segments, min_main_len: int = 5_000,
                                  max_main_len: int = SIZE_RANGES["general"][1],
                                  min_dominance_frac: float = 0.7,
                                  min_dominance_ratio: float = 5.0) -> bool:
    """
    True if one segment in a multi-segment component dominates it in length
    and every other segment is small by comparison — the signature Autocycler
    leaves for a real hairpin telomere: the main body plus a cluster of tiny
    "satellite" segments representing assembler ambiguity at the fold-back
    repeat, linked to one or both of the main segment's ends without closing
    a loop. Distinguishes that from genuine multi-contig fragmentation (e.g.
    an unresolved chromosome sitting in several comparably-large pieces),
    which won't have one segment towering over the rest.

    The same graph signature also shows up around large (multi-hundred-kb to
    Mb-scale) unresolved chromosome fragments — an isolated big piece with
    its own small tangle of assembly-ambiguity segments at its ends is not
    distinguishable from a real telomere by topology alone. max_main_len
    caps the main segment at SIZE_RANGES["general"]'s upper bound (500 kb,
    the same "plausible linear plasmid" ceiling classify_size() already
    uses) so this only fires in a size range a linear plasmid could
    plausibly occupy.
    """
    lens = sorted((segments.get(s, 0) for s in comp), reverse=True)
    if not lens:
        return False
    main_len, total_len = lens[0], sum(lens)
    second_len = lens[1] if len(lens) > 1 else 0
    return (min_main_len <= main_len <= max_main_len
            and main_len >= min_dominance_frac * total_len
            and (second_len == 0 or main_len >= min_dominance_ratio * second_len))


def _gfa_classify_components(segments, links):
    """Classify every connected component of the graph.
    Returns list of dicts: {'segs','length','topology','hairpin'}
    topology in {'circular','linear','fragmented'}."""
    fnext, rnext, fprev, rprev, undirected = _gfa_build_adjacency(links)
    out = []
    for comp in _gfa_connected_components(segments, undirected):
        cset = set(comp)
        length = sum(segments.get(s, 0) for s in comp)
        hp = any(a == b and sa != sb for (a, sa, b, sb) in links if a in cset)
        touches = any((a in cset or b in cset) for (a, _, b, _) in links)
        if _gfa_component_is_circular_loop(comp, fnext, rnext, fprev, rprev):
            topo = "circular"
        elif len(comp) == 1 and not touches:
            topo = "linear"          # isolated unitig, open ends
        elif len(comp) == 1:
            topo = "linear"          # single unitig with a hairpin link
        elif _gfa_dominant_segment_is_open(comp, segments):
            topo = "linear"          # dominant body + small satellite fragments
        else:
            topo = "fragmented"      # several comparably-sized unitigs, not a clean loop
        out.append({"segs": comp, "length": length, "topology": topo,
                    "hairpin": hp})
    return out


def _closest_contig_by_length(seg_len, length_map: dict, used: set,
                              rel_tol: float = 0.002, abs_tol: int = 50):
    """
    Find the final contig whose length is closest to seg_len, within a small
    tolerance. Exact equality is too strict for a raw pre-polish assembly
    graph (e.g. Flye's) matched against Hybracter's final polished contigs —
    polishing (medaka/polypolish/pypolca) typically shifts a contig's length
    by a handful of bases, not a meaningful fraction, so exact-length lookup
    silently drops almost every segment. Returns None if nothing is within
    tolerance, if two candidates tie for closest (ambiguous), or if the best
    candidate is already claimed by an earlier segment in this same file.
    """
    if not seg_len:
        return None
    tol = max(abs_tol, int(seg_len * rel_tol))
    candidates = sorted(
        (abs(seg_len - ln), ctg) for ctg, ln in length_map.items()
        if abs(seg_len - ln) <= tol and ctg not in used)
    if not candidates:
        return None
    if len(candidates) > 1 and candidates[0][0] == candidates[1][0]:
        return None  # tie — ambiguous
    return candidates[0][1]


def _parse_gfa_topology_single(gfa_file: str, length_map: dict = None) -> dict:
    """Parse one GFA file — see parse_gfa_topology for the multi-file wrapper."""
    if not gfa_file or not os.path.exists(gfa_file):
        return {}

    segments, links = _gfa_parse_segments_links(gfa_file)
    seg_lengths = {seg: ln for seg, ln in segments.items() if ln}
    components = _gfa_classify_components(segments, links)

    topology = {}
    for comp in components:
        for seg in comp["segs"]:
            topology[seg] = {"topology": comp["topology"], "hairpin": comp["hairpin"]}

    if length_map:
        remapped = {}
        used_ctgs = set()
        for seg, topo in topology.items():
            if seg in length_map:
                remapped[seg] = topo          # name already matches a contig ID
                used_ctgs.add(seg)
            else:
                ctg_id = _closest_contig_by_length(seg_lengths.get(seg), length_map, used_ctgs)
                if ctg_id:
                    remapped[ctg_id] = topo
                    used_ctgs.add(ctg_id)
                else:
                    remapped[seg] = topo      # left unmapped — no confident match
        topology = remapped

    return topology


def parse_gfa_topology(gfa_file, length_map: dict = None) -> dict:
    """
    Parse a GFA assembly graph (Flye/Unicycler/Plassembler output) and
    classify each segment's connected component as circular / linear /
    fragmented, walking the strand-aware graph rather than just counting
    link degree (see module comment above).

    gfa_file: a single path, or a list of paths (e.g. Hybracter's Flye
    assembly_graph.gfa *and* its Plassembler GFA — Flye's graph typically
    only covers the chromosome and any large plasmids it managed to
    resolve, while smaller/low-coverage plasmids are Plassembler-only and
    never appear in Flye's graph at all; passing both lets each contribute
    topology for the contigs it actually covers). Results are merged,
    first file wins on any contig id collision.

    Returns a dict: contig_id → {"topology": ..., "hairpin": bool}.

    length_map: optional {contig_id: seq_length} built from FASTA records.
    When GFA segment names differ from FASTA contig IDs (e.g. Plassembler's
    numeric names, or Flye's own edge_N/contig_N naming against Hybracter's
    final chromosome00001/plasmid00001 ids) segments are re-keyed by
    matching their length against the provided lengths, within a small
    tolerance (see _closest_contig_by_length).
    """
    files = [gfa_file] if isinstance(gfa_file, str) else list(gfa_file or [])
    merged: dict = {}
    for f in files:
        for ctg_id, topo in _parse_gfa_topology_single(f, length_map).items():
            merged.setdefault(ctg_id, topo)
    return merged


def is_linear_in_gfa(contig_id: str, gfa_topology: dict) -> dict:
    """Interpret GFA topology for a specific contig."""
    entry = gfa_topology.get(contig_id)
    topo = entry["topology"] if entry else "unknown"
    hairpin = entry["hairpin"] if entry else False
    is_linear = topo == "linear"
    desc = {
        "circular":   "Component forms a closed circular loop",
        "linear":     "Open-ended component" + (" with a hairpin self-link" if hairpin else ""),
        "fragmented": "Multiple unitigs that do not form a clean loop — ambiguous",
        "unknown":    "Contig not in GFA",
    }.get(topo, "unknown")
    return {
        "topology": topo,
        "hairpin": hairpin,
        "consistent_with_linear": is_linear,
        "description": desc,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 7b: GFA HAIRPIN-EDGE ANNOTATION (diagnostic, Bandage-facing)
# ─────────────────────────────────────────────────────────────────────────────
#
# Ported from linear-plasmid-hairpin-tools' add_hairpin_edges.py (2026-07).
# Detects a terminal fold-back on each segment of a raw assembler GFA (Flye's
# assembly_graph.gfa, Unicycler/Plassembler's assembly.gfa) and appends the
# same hairpin link Autocycler uses to represent one:
#   3' fold-back  ->  L <seg> + <seg> - 0M
#   5' fold-back  ->  L <seg> - <seg> + 0M
# This is an ANNOTATION only — it adds the link but does not alter segment
# sequences or split the contig at the fold apex, so it is correct for
# detection/visualization (Bandage, this script's own gfa_topology module)
# but not a length-accurate assembly-graph rebuild. It does not feed scoring.

def annotate_gfa_hairpins(gfa_path: str, out_path: str = None,
                          k: int = 31, window: int = 50_000,
                          min_shared: int = 25, edge_tol: int = 100,
                          overlap: bool = False) -> dict:
    """
    Detect terminal hairpins in `gfa_path`'s segments and write a copy with
    Autocycler-style hairpin links added.

    Returns {"out_path": str | None, "n_added": int, "links": [...]}.
    out_path is None if no hairpins were found (nothing written).
    """
    lines = []
    seqs = {}
    existing_links = set()
    with open(gfa_path) as fh:
        for line in fh:
            lines.append(line.rstrip("\n"))
            if line[:1] == "S":
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and f[2] != "*":
                    seqs[f[1]] = f[2].encode("ascii", "replace").upper()
            elif line[:1] == "L":
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5:
                    existing_links.add((f[1], f[2], f[3], f[4]))

    new_links = []   # (seg, from_strand, to_strand, overlap_str, info_dict)
    for seg, seq in seqs.items():
        for end in ("5", "3"):
            r = _analyse_hairpin_end(seq, k, window, end, edge_tol)
            if r is None or r["support"] < min_shared or not r["touches"]:
                continue
            # 3' fold -> L seg + seg - ; 5' fold -> L seg - seg +
            sa, sb = ("+", "-") if end == "3" else ("-", "+")
            if (seg, sa, seg, sb) in existing_links:
                continue
            arm1_s, arm1_e = r["arm1"]
            arm2_s, arm2_e = r["arm2"]
            arm_len = min(arm1_e - arm1_s, arm2_e - arm2_s)
            ov = f"{arm_len}M" if overlap else "0M"
            new_links.append((seg, sa, sb, ov,
                              {"segment": seg, "end": f"{end}'", "arm_len": arm_len,
                               "gap": r["gap"], "support": r["support"],
                               "seg_len": len(seq)}))
            existing_links.add((seg, sa, seg, sb))

    if not new_links:
        return {"out_path": None, "n_added": 0, "links": []}

    out = out_path or (os.path.splitext(gfa_path)[0] + ".hairpins.gfa")
    with open(out, "w") as fh:
        for ln in lines:
            fh.write(ln + "\n")
        for seg, sa, sb, ov, _info in new_links:
            fh.write(f"L\t{seg}\t{sa}\t{seg}\t{sb}\t{ov}\n")

    return {"out_path": out, "n_added": len(new_links),
            "links": [info for *_, info in new_links]}


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 8: BLAST-BASED DETECTION
# ─────────────────────────────────────────────────────────────────────────────

def run_blast(query_fasta: str, db: str, out_file: str,
              mode: str = "blastn", identity: float = 70.0,
              coverage: float = 30.0) -> pd.DataFrame:
    """
    Run BLASTn or BLASTp against a plasmid database (e.g. PLSDB).
    Returns DataFrame of significant hits above identity/coverage thresholds.
    """
    if not db or not os.path.exists(query_fasta):
        return pd.DataFrame()

    cmd = [
        mode, "-query", query_fasta, "-db", db,
        "-out", out_file, "-outfmt",
        "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle",
        "-max_target_seqs", "5", "-num_threads", "4",
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return pd.DataFrame()

    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        return pd.DataFrame()

    cols = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
            "qcovs", "evalue", "bitscore", "stitle"]
    try:
        df = pd.read_csv(out_file, sep="\t", names=cols)
        df = df[(df["pident"] >= identity) & (df["qcovs"] >= coverage)]
        return df
    except Exception:
        return pd.DataFrame()


def interpret_blast_hits(blast_df: pd.DataFrame) -> dict:
    """Summarise BLAST hits relevant to linear plasmid identification.

    Two tiers:
      linear_plasmid_hit     — any hit to a sequence with linear-plasmid keywords
      linear_plasmid_db_hit  — hit to a sequence explicitly named as a linear plasmid
                               (pELF, pBSSB, lp28, etc.) — higher scoring weight
    """
    if blast_df.empty:
        return {"hits": 0, "best_identity": 0, "linear_plasmid_hit": False,
                "linear_plasmid_db_hit": False}

    # Tier 1: broad — any topology or family keyword
    broad_keywords = [
        "linear plasmid", "linear chromosome", "linear replicon",
        "linear_plasmid", "telomere",
    ]
    # Tier 2: specific named linear plasmid families
    # Use prefix anchors only (\bpelf not \bpelf\b) so variants like
    # pELF_AA290 or pBSSB1 are matched despite the underscore/digit breaking \b.
    specific_keywords = [
        r"\bpelf",
        r"\bpbssb",
    ]

    titles_lower = blast_df["stitle"].str.lower().fillna("")

    broad_mask    = titles_lower.str.contains("|".join(broad_keywords), na=False)
    specific_mask = titles_lower.str.contains("|".join(specific_keywords),
                                               na=False, regex=True)

    broad_hits    = blast_df[broad_mask | specific_mask]
    specific_hits = blast_df[specific_mask]

    best = blast_df.iloc[0]
    return {
        "hits":                  len(blast_df),
        "best_identity":         float(blast_df["pident"].max()),
        "best_coverage":         float(blast_df["qcovs"].max()),
        "linear_plasmid_hit":    len(broad_hits) > 0,
        "linear_plasmid_db_hit": len(specific_hits) > 0,
        "top_hit":               str(best["stitle"]),
        "top_identity":          float(best["pident"]),
    }


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 8b: SKANI-BASED DETECTION
# ─────────────────────────────────────────────────────────────────────────────

def run_skani(query_fasta: str, db: str, out_file: str,
              min_ani: float = 90.0, min_af: float = 0.05,
              n_results: int = 5) -> pd.DataFrame:
    """Run skani against a plasmid database and return a normalised DataFrame.

    db may be a pre-sketched directory (uses `skani search`) or a FASTA file
    (uses `skani dist`).  For large multi-sequence FASTA databases, pre-sketch
    first: skani sketch -l db.fna -o db_dir/

    Returns columns matching interpret_blast_hits expectations:
    pident (ANI %), qcovs (align_fraction_query * 100), stitle (ref basename).
    """
    if not db or not os.path.exists(query_fasta):
        return pd.DataFrame()

    if os.path.isdir(db):
        cmd = ["skani", "search", "-d", db, "-q", query_fasta,
               "-o", out_file, "-n", str(n_results)]
    else:
        cmd = ["skani", "dist", query_fasta, db, "-o", out_file]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return pd.DataFrame()

    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        return pd.DataFrame()

    try:
        df = pd.read_csv(out_file, sep="\t")
        df["pident"] = df["ANI"]
        # Align_fraction_query is already in percent in skani ≥0.3
        df["qcovs"]  = df["Align_fraction_query"]
        # Prefer the full sequence header (Ref_name) for keyword matching;
        # fall back to the filename stem if the column is absent (older skani)
        if "Ref_name" in df.columns:
            df["stitle"] = df["Ref_name"].fillna(
                df["Ref_file"].apply(lambda p: os.path.splitext(os.path.basename(p))[0])
            )
        else:
            df["stitle"] = df["Ref_file"].apply(
                lambda p: os.path.splitext(os.path.basename(p))[0]
            )
        df = df[(df["pident"] >= min_ani) & (df["qcovs"] >= min_af * 100)]
        return df.sort_values("pident", ascending=False).reset_index(drop=True)
    except Exception:
        return pd.DataFrame()


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 9: COMPOSITE SCORING
# ─────────────────────────────────────────────────────────────────────────────

def compute_score(evidence: dict) -> dict:
    """
    Combine all evidence modules into a composite score.
    """
    score = 0
    breakdown = {}

    # 0. FASTA header metadata
    hdr = evidence.get("header", {})
    is_confirmed_circular = hdr.get("is_circular", False)

    if hdr.get("assembler_not_circular"):
        score += SCORING_WEIGHTS["assembler_not_circular"]
        breakdown["assembler_not_circular"] = SCORING_WEIGHTS["assembler_not_circular"]
    # circular_flag_absent is NOT scored (see assess_header_metadata docstring) —
    # it fires whenever header metadata is simply missing, not when linearity is
    # actually confirmed, so it carries no discriminating evidence value.

    if hdr.get("header_cn_consistent") and not hdr.get("high_copy_circular"):
        score += SCORING_WEIGHTS["header_copy_number"]
        breakdown["header_copy_number"] = SCORING_WEIGHTS["header_copy_number"]

    # 1b. Hairpin telomere at either end — either signal is sufficient:
    #   - confirmed_idr_tata: specific (arm + known loop motif, e.g. TATA)
    #   - hairpin_ends (detect_terminal_hairpins): general localized fold-back
    #     detector (ported from find_hairpins.py, 2026-07) — finds the
    #     inverted-repeat symmetry axis per end and requires the outer arm to
    #     reach the contig terminus. Correctly handles asymmetric ends (only
    #     one end needs to fold back), unlike the old whole-window compare.
    idr = evidence.get("idr_tata", {})
    hairpin_ends = evidence.get("hairpin_ends", {})
    if idr.get("confirmed_idr_tata") or hairpin_ends.get("left") or hairpin_ends.get("right"):
        score += SCORING_WEIGHTS["hairpin_end"]
        breakdown["hairpin_end"] = SCORING_WEIGHTS["hairpin_end"]

    # 2b. Asymmetric ends (pELF1-type: one hairpin + one invertron)
    asym = evidence.get("asymmetric_ends", {})
    if asym.get("asymmetric_pelf_type"):
        score += SCORING_WEIGHTS["asymmetric_ends"]
        breakdown["asymmetric_ends"] = SCORING_WEIGHTS["asymmetric_ends"]

    # Invertron terminal protein gene
    if asym.get("has_tp_gene"):
        score += SCORING_WEIGHTS["invertron_tp_gene"]
        breakdown["invertron_tp_gene"] = SCORING_WEIGHTS["invertron_tp_gene"]

    # 3. Size
    if evidence.get("size", {}).get("in_known_range"):
        score += SCORING_WEIGHTS["size_range"]
        breakdown["size_range"] = SCORING_WEIGHTS["size_range"]

    # 4. GC content
    if evidence.get("gc", {}).get("low_gc_flag"):
        score += SCORING_WEIGHTS["gc_deviation"]
        breakdown["gc_deviation"] = SCORING_WEIGHTS["gc_deviation"]

    # 5. Gene-based scoring.
    # Suppressed when annotation could not be matched to this contig
    # (Prokka-TSV fallback: all genes smeared across all contigs), OR when the
    # header already confirms circular=true (CIRCULAR_DISQUALIFIES_GENE_SCORING):
    # partition systems, TA systems, and IS elements are common chromosomal
    # background too, so on a contig the assembler itself calls circular they
    # aren't discriminating evidence for a *linear* plasmid. Structural/
    # sequence evidence (hairpin, GFA topology, blast/skani, coverage drop)
    # still runs independently and is what should override a mistaken
    # assembler circular= call (see Hashimoto 2019 hairpin-fools-assembler).
    genes = evidence.get("genes", {})
    if genes.get("_annotation_mismatch") or (CIRCULAR_DISQUALIFIES_GENE_SCORING and is_confirmed_circular):
        genes = {}   # suppress gene scoring

    if genes.get("partition_genes"):
        score += SCORING_WEIGHTS["par_system"]
        breakdown["par_system"] = SCORING_WEIGHTS["par_system"]

    # 5b. repB / Rep_2 superfamily (Hashimoto 2019)
    if genes.get("repB_rep2"):
        score += SCORING_WEIGHTS["repb_rep2"]
        breakdown["repb_rep2"] = SCORING_WEIGHTS["repb_rep2"]

    # 5c. ftsK + parA + repB trio co-occurrence (pELF1 signature, Hashimoto 2019)
    trio_found = genes.get("ftsK_parA_repB_trio", [])
    if len(trio_found) >= 3:   # all three present
        score += SCORING_WEIGHTS["ftsK_parA_repB_combo"]
        breakdown["ftsK_parA_repB_combo"] = SCORING_WEIGHTS["ftsK_parA_repB_combo"]

    # 6. TAS
    if genes.get("tas_systems"):
        score += SCORING_WEIGHTS["tas_system"]
        breakdown["tas_system"] = SCORING_WEIGHTS["tas_system"]

    # 7. IS elements
    if genes.get("is_elements"):
        score += SCORING_WEIGHTS["is_elements"]
        breakdown["is_elements"] = SCORING_WEIGHTS["is_elements"]

    # 9. BLAST — two tiers
    blast = evidence.get("blast", {})
    if blast.get("linear_plasmid_db_hit"):
        # Named linear plasmid family (pELF, pBSSB, lp28, lp36, lp54, pELF_USZ...)
        score += SCORING_WEIGHTS["blast_hit_linear_db"]
        breakdown["blast_hit_linear_db"] = SCORING_WEIGHTS["blast_hit_linear_db"]
    elif blast.get("linear_plasmid_hit"):
        # Generic linear-plasmid keyword hit
        score += SCORING_WEIGHTS["blast_hit"]
        breakdown["blast_hit"] = SCORING_WEIGHTS["blast_hit"]

    # 9b. SKANI — two tiers (independent of BLAST; same weights for comparison)
    skani = evidence.get("skani", {})
    if skani.get("linear_plasmid_db_hit"):
        score += SCORING_WEIGHTS["skani_hit_linear_db"]
        breakdown["skani_hit_linear_db"] = SCORING_WEIGHTS["skani_hit_linear_db"]
    elif skani.get("linear_plasmid_hit"):
        score += SCORING_WEIGHTS["skani_hit"]
        breakdown["skani_hit"] = SCORING_WEIGHTS["skani_hit"]

    # 9c. PlasmidFinder no-hit (novel rep = typical of linear plasmids, Hashimoto 2019)
    if evidence.get("plasmid_finder_no_hit"):
        score += SCORING_WEIGHTS["plasmid_finder_no_hit"]
        breakdown["plasmid_finder_no_hit"] = SCORING_WEIGHTS["plasmid_finder_no_hit"]

    # 10. Coverage drop at ends (hairpin inaccessibility, Hashimoto 2019)
    cov = evidence.get("coverage_drop", {})
    if cov.get("available") and cov.get("consistent_with_linear"):
        score += SCORING_WEIGHTS["coverage_drop_ends"]
        breakdown["coverage_drop_ends"] = SCORING_WEIGHTS["coverage_drop_ends"]

    # 11. Assembly graph
    gfa = evidence.get("gfa_topology", {})
    if gfa.get("consistent_with_linear"):
        score += SCORING_WEIGHTS["assembly_graph_linear"]
        breakdown["assembly_graph_linear"] = SCORING_WEIGHTS["assembly_graph_linear"]

    # 12. Organism-specific markers
    if genes.get("enterococcal_markers"):
        score += SCORING_WEIGHTS["enterococcal_markers"]
        breakdown["enterococcal_markers"] = SCORING_WEIGHTS["enterococcal_markers"]

    # 13. Copy number
    cn = evidence.get("copy_number", {})
    if cn.get("available") and cn.get("consistent_with_linear"):
        score += SCORING_WEIGHTS["copy_number_low"]
        breakdown["copy_number_low"] = SCORING_WEIGHTS["copy_number_low"]

    # ── Hard gates ──────────────────────────────────────────────────────────
    # Applied last, after every category above has had a chance to score, so
    # no future evidence category can silently bypass them.
    gate_reason = None

    # 1. Chromosomes can never be linear plasmids, regardless of evidence.
    if evidence.get("is_chromosome"):
        gate_reason = ("chromosome (>500 kb and/or 'chromosome' in contig id) "
                       "— excluded from linear-plasmid scoring regardless of "
                       "its circular= flag")
        score = 0
        breakdown = {}

    # 2. A plasmid-sized contig the assembler already calls circular=true can
    # only be scored as a linear-plasmid candidate if it carries genuine,
    # independently-detected telomere structure (hairpin fold-back or a
    # pELF1-type asymmetric hairpin+invertron end) — the actual structural
    # signature Hashimoto 2019 describes an assembler being fooled by.
    # Indirect evidence (skani/blast hit, GFA topology, coverage-drop, gene
    # content, size/GC/copy-number) is not on its own strong enough to
    # override a circular call: those signals are also seen on ordinary
    # circular replicons, so allowing them through here would make every
    # circular=true plasmid a potential false positive.
    elif is_confirmed_circular:
        has_telomere_evidence = bool(
            idr.get("confirmed_idr_tata")
            or hairpin_ends.get("left") or hairpin_ends.get("right")
            or asym.get("asymmetric_pelf_type") or asym.get("has_tp_gene"))
        if not has_telomere_evidence:
            gate_reason = ("circular=true in header with no independently-detected "
                           "hairpin/telomere evidence — indirect evidence alone "
                           "cannot override an assembler circular call")
            score = 0
            breakdown = {}

    # Confidence level
    if score >= CONFIDENCE_THRESHOLDS["HIGH"]:
        confidence = "HIGH"
    elif score >= CONFIDENCE_THRESHOLDS["MEDIUM"]:
        confidence = "MEDIUM"
    elif score >= CONFIDENCE_THRESHOLDS["LOW"]:
        confidence = "LOW"
    else:
        confidence = "NONE"

    # blast_hit/blast_hit_linear_db and skani_hit/skani_hit_linear_db are each
    # mutually exclusive (elif above) — only the higher tier of each pair can
    # ever be reached, so the lower tier's weight must not count toward the
    # achievable ceiling.
    max_score = (sum(SCORING_WEIGHTS.values())
                 - min(SCORING_WEIGHTS["blast_hit"], SCORING_WEIGHTS["blast_hit_linear_db"])
                 - min(SCORING_WEIGHTS["skani_hit"], SCORING_WEIGHTS["skani_hit_linear_db"]))
    return {
        "total_score": score,
        "max_possible": max_score,
        "percent": round(100 * score / max_score, 1),
        "confidence": confidence,
        "breakdown": breakdown,
        "gate_reason": gate_reason,
    }


# ─────────────────────────────────────────────────────────────────────────────
# LONG-READ MAPPING HELPER
# ─────────────────────────────────────────────────────────────────────────────

def _run_mapping_pipeline(mm2_cmd: list, sort_cmd: list,
                          output_bam: str) -> str:
    """
    Shared helper: pipe minimap2 stdout into samtools sort, then index.
    Returns output_bam on success, '' on failure.
    """
    try:
        mm2  = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        sort = subprocess.Popen(sort_cmd, stdin=mm2.stdout,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        mm2.stdout.close()          # allow mm2 to receive SIGPIPE if sort exits
        _, sort_err = sort.communicate()
        mm2.wait()

        if mm2.returncode != 0:
            _, mm2_err = mm2.communicate()
            print(f"[ERROR] minimap2 failed (exit {mm2.returncode}):\n"
                  f"{mm2_err.decode()}", file=sys.stderr)
            return ""
        if sort.returncode != 0:
            print(f"[ERROR] samtools sort failed (exit {sort.returncode}):\n"
                  f"{sort_err.decode()}", file=sys.stderr)
            return ""

    except Exception as e:
        print(f"[ERROR] Mapping pipeline failed: {e}", file=sys.stderr)
        return ""

    if not os.path.exists(output_bam) or os.path.getsize(output_bam) == 0:
        print(f"[ERROR] BAM file not created or empty: {output_bam}",
              file=sys.stderr)
        return ""

    try:
        subprocess.run(["samtools", "index", output_bam], check=True,
                       capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] samtools index failed: {e}", file=sys.stderr)
        return ""

    bam_size = os.path.getsize(output_bam)
    print(f"[MAP] Done → {output_bam}  ({bam_size / 1e6:.1f} MB)")
    return output_bam


def _check_mapping_deps() -> bool:
    """Return True if minimap2 and samtools are on PATH."""
    missing = []
    for tool in ("minimap2", "samtools"):
        try:
            subprocess.run([tool, "--version"], check=True, capture_output=True)
        except (FileNotFoundError, subprocess.CalledProcessError):
            missing.append(tool)
    if missing:
        print(f"[ERROR] Required tool(s) not found: {', '.join(missing)}. "
              f"Install with: conda install -c bioconda minimap2 samtools",
              file=sys.stderr)
        return False
    return True


def map_longreads(fastq: str, assembly_fasta: str, output_bam: str,
                  preset: str = "map-ont", threads: int = 4) -> str:
    """
    Map long reads to the assembly with minimap2, sort and index the result.

    Runs:
        minimap2 -ax <preset> -t <threads> assembly.fasta reads.fastq \\
            | samtools sort -@ <threads> -o output.bam
        samtools index output.bam

    Parameters
    ----------
    fastq         : path to long-read FASTQ (or FASTQ.gz)
    assembly_fasta: path to the assembly FASTA (the -i input)
    output_bam    : path for the sorted BAM output
    preset        : minimap2 preset — "map-ont" (Nanopore) or "map-pb" (PacBio)
    threads       : threads for both minimap2 and samtools

    Returns
    -------
    Path to the sorted, indexed BAM file, or empty string on failure.
    """
    if not _check_mapping_deps():
        return ""
    if not os.path.exists(fastq):
        print(f"[ERROR] Long-read FASTQ not found: {fastq}", file=sys.stderr)
        return ""
    if not os.path.exists(assembly_fasta):
        print(f"[ERROR] Assembly FASTA not found: {assembly_fasta}", file=sys.stderr)
        return ""

    print(f"[MAP] minimap2 {preset}  {fastq}  →  {output_bam}")

    mm2_cmd = [
        "minimap2", "-ax", preset,
        "-t", str(threads),
        "--secondary=no",
        assembly_fasta, fastq,
    ]
    sort_cmd = [
        "samtools", "sort", "-@", str(threads), "-o", output_bam,
    ]
    return _run_mapping_pipeline(mm2_cmd, sort_cmd, output_bam)


def _deinterleave_to_tempfiles(interleaved_path: str) -> tuple:
    """
    Split an interleaved paired-end FASTQ (plain or gzipped) into two
    temporary files containing R1 and R2 reads respectively.

    Returns (r1_tmp_path, r2_tmp_path) on success, or ('', '') on failure.
    The caller is responsible for deleting the temp files after use.
    """
    import gzip as _gz
    import tempfile as _tf

    open_fn = _gz.open if interleaved_path.endswith(".gz") else open
    try:
        r1_fh = _tf.NamedTemporaryFile(mode="w", suffix="_R1.fastq",
                                        delete=False)
        r2_fh = _tf.NamedTemporaryFile(mode="w", suffix="_R2.fastq",
                                        delete=False)
        r1_path, r2_path = r1_fh.name, r2_fh.name

        with open_fn(interleaved_path, "rt") as src:
            read_idx = 0
            buf = []
            for line in src:
                buf.append(line)
                if len(buf) == 4:          # one complete FASTQ record
                    dest = r1_fh if read_idx % 2 == 0 else r2_fh
                    dest.writelines(buf)
                    buf = []
                    read_idx += 1

        r1_fh.close()
        r2_fh.close()
        return r1_path, r2_path

    except Exception as exc:
        print(f"[ERROR] Deinterleaving failed: {exc}", file=sys.stderr)
        for p in (r1_fh.name, r2_fh.name):
            try:
                os.unlink(p)
            except OSError:
                pass
        return "", ""


def map_shortreads(assembly_fasta: str, output_bam: str,
                   r1: str = None, r2: str = None,
                   interleaved: str = None,
                   threads: int = 4) -> str:
    """
    Map Illumina short reads to the assembly with minimap2 (sr preset),
    sort and index the result.

    Accepts either:
      · Separate R1 + R2 FASTQ files  (r1= and r2=)
      · A single interleaved paired-end FASTQ  (interleaved=)
        The interleaved file is split into temp R1/R2 files in Python before
        calling minimap2 (minimap2 ≤2.30 lacks a native --interleaved flag).

    Runs:
        minimap2 -ax sr -t <threads> --secondary=no assembly.fasta R1.fq R2.fq \\
            | samtools sort -@ <threads> -o output.bam
        samtools index output.bam

    Parameters
    ----------
    assembly_fasta : path to the assembly FASTA (-i input)
    output_bam     : path for the sorted BAM output
    r1             : R1 FASTQ path (separate-files mode)
    r2             : R2 FASTQ path (separate-files mode)
    interleaved    : interleaved paired-end FASTQ path (plain or .gz)
    threads        : threads for minimap2 and samtools

    Returns
    -------
    Path to the sorted, indexed BAM file, or empty string on failure.
    """
    if not _check_mapping_deps():
        return ""
    if not os.path.exists(assembly_fasta):
        print(f"[ERROR] Assembly FASTA not found: {assembly_fasta}", file=sys.stderr)
        return ""

    tmp_r1 = tmp_r2 = None   # temp paths to clean up after mapping

    if interleaved:
        if not os.path.exists(interleaved):
            print(f"[ERROR] Interleaved FASTQ not found: {interleaved}", file=sys.stderr)
            return ""
        print(f"[MAP] Deinterleaving  {interleaved}  …", flush=True)
        tmp_r1, tmp_r2 = _deinterleave_to_tempfiles(interleaved)
        if not tmp_r1:
            return ""
        r1, r2 = tmp_r1, tmp_r2
        print(f"[MAP] minimap2 sr (paired)  {interleaved}  →  {output_bam}")

    elif r1 and r2:
        for path, label in [(r1, "R1"), (r2, "R2")]:
            if not os.path.exists(path):
                print(f"[ERROR] Short-read {label} FASTQ not found: {path}",
                      file=sys.stderr)
                return ""
        print(f"[MAP] minimap2 sr  {r1}  {r2}  →  {output_bam}")

    else:
        print("[ERROR] map_shortreads: supply either interleaved= or both r1= and r2=",
              file=sys.stderr)
        return ""

    mm2_cmd = [
        "minimap2", "-ax", "sr",
        "-t", str(threads),
        "--secondary=no",
        assembly_fasta, r1, r2,
    ]
    sort_cmd = [
        "samtools", "sort", "-@", str(threads), "-o", output_bam,
    ]
    result = _run_mapping_pipeline(mm2_cmd, sort_cmd, output_bam)

    # Clean up temp deinterleave files
    for p in (tmp_r1, tmp_r2):
        if p:
            try:
                os.unlink(p)
            except OSError:
                pass

    return result




# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZATION  (Figure-3-style terminal structure plot)
# ─────────────────────────────────────────────────────────────────────────────

def _extract_clip_profile(bam, contig_id: str, contig_len: int,
                           window: int = 5000) -> dict:
    """
    Scan terminal windows and count per-position soft- and hard-clipped reads.

    For each alignment that overlaps either terminal window:
      left panel  — reads whose 5′ alignment start is inside the left window,
                    split by whether the read has a soft clip (CIGAR starts S)
                    or hard clip (CIGAR starts H) at that start position.
      right panel — reads whose 3′ alignment end is inside the right window,
                    split by soft clip (CIGAR ends S) or hard clip (CIGAR ends H).

    A pile-up of 5′-soft-clipped reads near position 0 is the classic Illumina
    signature of a hairpin telomere: the hairpin blocks adapter ligation so all
    reads that would cross it terminate at the tip (Hashimoto 2019, Fig 3A).
    Hard clips appear when the aligner truncates alignment at the hairpin boundary.

    Returns dict with keys:
      left_soft  : {pos → count}   5′-soft-clipped reads by alignment start
      left_hard  : {pos → count}   5′-hard-clipped reads by alignment start
      right_soft : {pos → count}   3′-soft-clipped reads by alignment end
      right_hard : {pos → count}   3′-hard-clipped reads by alignment end
      left_window_end  : int (x upper bound for left plot)
      right_window_start : int (x lower bound for right plot)
    """
    from collections import defaultdict
    SOFT, HARD = 4, 5

    left_soft:  dict = defaultdict(int)
    left_hard:  dict = defaultdict(int)
    right_soft: dict = defaultdict(int)
    right_hard: dict = defaultdict(int)

    left_end   = min(window, contig_len)
    right_start = max(0, contig_len - window)

    seen: set = set()
    for region in [(0, left_end), (right_start, contig_len)]:
        try:
            for read in bam.fetch(contig_id, region[0], region[1]):
                if read.is_unmapped or not read.cigartuples:
                    continue
                uid = (read.query_name, read.reference_start, read.flag)
                if uid in seen:
                    continue
                seen.add(uid)

                cigar     = read.cigartuples
                ref_start = read.reference_start
                ref_end   = read.reference_end or ref_start

                # 5′-end clip (at the start of the alignment)
                if ref_start < left_end:
                    if cigar[0][0] == SOFT:
                        left_soft[ref_start] += 1
                    elif cigar[0][0] == HARD:
                        left_hard[ref_start] += 1

                # 3′-end clip (at the end of the alignment)
                if ref_end > right_start:
                    if cigar[-1][0] == SOFT:
                        right_soft[ref_end] += 1
                    elif cigar[-1][0] == HARD:
                        right_hard[ref_end] += 1
        except Exception:
            pass

    return {
        "left_soft":           dict(left_soft),
        "left_hard":           dict(left_hard),
        "right_soft":          dict(right_soft),
        "right_hard":          dict(right_hard),
        "left_window_end":     left_end,
        "right_window_start":  right_start,
    }


def _draw_clip_panel(ax, soft_dict: dict, hard_dict: dict,
                     x_start: int, x_end: int,
                     title: str, bin_size: int = 50) -> None:
    """
    Draw a stacked-bar histogram of soft-clip (teal) + hard-clip (orange)
    read counts in [x_start, x_end), binned into <bin_size>-bp windows.
    """
    # Build bin edges and compute counts
    edges = list(range(x_start, x_end, bin_size))
    if not edges:
        ax.set_visible(False)
        return

    soft_counts = []
    hard_counts = []
    centers     = []
    for b in edges:
        b_end = b + bin_size
        sc = sum(soft_dict.get(p, 0) for p in range(b, b_end))
        hc = sum(hard_dict.get(p, 0) for p in range(b, b_end))
        soft_counts.append(sc)
        hard_counts.append(hc)
        centers.append(b + bin_size / 2)

    w = bin_size * 0.85
    ax.bar(centers, soft_counts, width=w,
           color="#00ACC1", alpha=0.85, label="Soft clip", linewidth=0)
    ax.bar(centers, hard_counts, width=w, bottom=soft_counts,
           color="#FF7043", alpha=0.85, label="Hard clip", linewidth=0)

    ax.set_xlim(x_start, x_end)
    ax.set_title(title, fontsize=8.5, pad=3)
    ax.set_xlabel("Position (bp)", fontsize=7.5)
    ax.set_ylabel("Clipped reads", fontsize=7.5)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.7)

    # Mark position 0 / contig end with a dashed line
    boundary = x_start if x_start == 0 else x_end
    ax.axvline(boundary, color="black", linewidth=0.8, linestyle="--", alpha=0.5)


def visualize_terminal_structure(contig_id: str, seq: str, evidence: dict,
                                   bam_file: str, output_prefix: str) -> str:
    """
    Generate a multi-panel PNG for one contig:

    Panel 1 (if BAM): full-contig coverage depth with terminal-window shading
                      and coverage-drop ratio annotations.
    Panel 2 (if BAM): soft-clip and hard-clip read counts at each end.
                      Left subplot — 5′-clipped reads by alignment-start position
                      in the first 5 kb.  Pile-up near position 0 indicates a
                      hairpin telomere (reads cannot extend past the tip).
                      Right subplot — 3′-clipped reads by alignment-end position
                      in the last 5 kb; same interpretation at the right telomere.
                      Soft clips: teal (#00ACC1).  Hard clips: orange (#FF7043).
    Panel 3: contig schematic with IDR arm regions colour-coded (left end blue,
             right end deep-orange) and loop sequence labelled.

    Returns the PNG path on success, empty string on failure.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        import matplotlib.patches as mpatches
    except ImportError:
        print("[WARN] matplotlib not installed; skipping visualization. "
              "Install with: pip install matplotlib", file=sys.stderr)
        return ""

    seq_len = len(seq)
    idr  = evidence.get("idr_tata", {})
    cov  = evidence.get("coverage_drop", {})
    hp   = evidence.get("hairpin_ends", {})
    sc_r = evidence.get("score", {})

    has_bam = bool(bam_file and os.path.exists(bam_file))

    # ── Open BAM once for all data collection ────────────────────────────────
    positions: list = []
    depths: list    = []
    clip_data: dict = {}

    if has_bam:
        try:
            import pysam
            _bam = pysam.AlignmentFile(bam_file, "rb")
            # Coverage depth
            for col in _bam.pileup(contig_id, 0, seq_len, min_mapping_quality=0):
                positions.append(col.reference_pos)
                depths.append(col.nsegments)
            # Clipping profiles
            clip_data = _extract_clip_profile(_bam, contig_id, seq_len)
            _bam.close()
        except Exception as _e:
            has_bam = False
            print(f"[WARN] BAM read failed for {contig_id}: {_e}", file=sys.stderr)

    has_coverage = has_bam and bool(positions)
    has_clips    = has_bam and bool(clip_data)

    # ── Figure layout: up to 3 rows, 2 columns ───────────────────────────────
    # Coverage and structure rows span both columns; clip row uses 2 columns.
    panels = []
    if has_coverage:
        panels.append("coverage")
    if has_clips:
        panels.append("clips")
    panels.append("structure")

    n_rows   = len(panels)
    fig_h    = 3.2 * n_rows + 0.6
    fig      = plt.figure(figsize=(14, fig_h))
    gs       = GridSpec(n_rows, 2, figure=fig, hspace=0.55, wspace=0.25)

    conf  = sc_r.get("confidence", "?")
    title = (f"{contig_id}  |  {seq_len:,} bp  |  "
             f"score {sc_r.get('total_score', '?')}  [{conf}]")
    fig.suptitle(title, fontsize=10, y=1.01)

    row = 0

    # ── Panel 1: Coverage depth ───────────────────────────────────────────────
    ax_cov = None
    if "coverage" in panels:
        ax_cov = fig.add_subplot(gs[row, :])   # span both columns
        row += 1

        max_d = max(depths) if depths else 1
        ax_cov.fill_between(positions, depths, alpha=0.55, color="#4C72B0",
                            linewidth=0)
        ax_cov.set_xlabel("Position (bp)", fontsize=8)
        ax_cov.set_ylabel("Read depth", fontsize=8)
        ax_cov.set_xlim(0, seq_len)
        ax_cov.tick_params(labelsize=7)
        body_depth = cov.get("body_mean_depth", "?")
        ax_cov.set_title(f"Coverage depth  (body mean: {body_depth}×)", fontsize=9)

        end_w = 5000
        for x0, x1 in [(0, min(end_w, seq_len)),
                        (max(0, seq_len - end_w), seq_len)]:
            ax_cov.axvspan(x0, x1, alpha=0.12, color="red", linewidth=0)

        for side, ratio, drop, xpos, ha in [
            ("L", cov.get("left_depth_ratio"),  cov.get("left_drop"),
             200, "left"),
            ("R", cov.get("right_depth_ratio"), cov.get("right_drop"),
             seq_len - 200, "right"),
        ]:
            if ratio is not None:
                colour = "red" if drop else "#2E7D32"
                ax_cov.text(xpos, max_d * 0.88,
                            f"{side}: {ratio:.2f}×",
                            fontsize=8, color=colour, ha=ha)

    # ── Panel 2: Clipping profiles ────────────────────────────────────────────
    if "clips" in panels:
        ax_cl = fig.add_subplot(gs[row, 0])   # left end
        ax_cr = fig.add_subplot(gs[row, 1])   # right end
        row += 1

        lwe = clip_data.get("left_window_end",    min(5000, seq_len))
        rws = clip_data.get("right_window_start", max(0, seq_len - 5000))

        # Adapt bin size: aim for ~60 bars per panel regardless of window size
        bin_size = max(10, (lwe) // 60)

        _draw_clip_panel(
            ax_cl,
            clip_data.get("left_soft",  {}),
            clip_data.get("left_hard",  {}),
            0, lwe,
            title=f"Left end  (0 – {lwe:,} bp)\n5′-clipped read starts",
            bin_size=bin_size,
        )
        _draw_clip_panel(
            ax_cr,
            clip_data.get("right_soft", {}),
            clip_data.get("right_hard", {}),
            rws, seq_len,
            title=f"Right end  ({rws:,} – {seq_len:,} bp)\n3′-clipped read ends",
            bin_size=bin_size,
        )

    # ── Panel 3: IDR structure schematic ─────────────────────────────────────
    ax_str = fig.add_subplot(gs[row, :])   # span both columns

    ax_str.set_xlim(-seq_len * 0.02, seq_len * 1.02)
    ax_str.set_ylim(-0.8, 1.6)
    ax_str.axis("off")
    ax_str.set_title("Terminal IDR / TATA structure", fontsize=9)

    ax_str.plot([0, seq_len], [0.5, 0.5], color="black", linewidth=4,
                solid_capstyle="butt", zorder=1)

    arm_colours = {"left": "#1565C0", "right": "#E64A19"}
    arm_y  = 0.5
    arm_h  = 0.28

    for side, idr_side in [("left",  idr.get("left",  {})),
                            ("right", idr.get("right", {}))]:
        if not idr_side.get("found"):
            continue
        arm_len   = idr_side.get("arm_len", 0)
        loop_seq  = idr_side.get("loop_seq", "")
        identity  = idr_side.get("identity", 0.0)
        loop_len  = len(loop_seq)
        confirmed = idr_side.get("confirmed", False)
        col       = arm_colours[side]

        if side == "left":
            x_a1 = (0,        arm_len)
            x_lp = (arm_len,  arm_len + loop_len)
            x_a2 = (arm_len + loop_len, arm_len + loop_len + arm_len)
        else:
            x_a2 = (seq_len - arm_len,              seq_len)
            x_lp = (seq_len - arm_len - loop_len,   seq_len - arm_len)
            x_a1 = (seq_len - 2*arm_len - loop_len, x_lp[0])

        ymin_frac = (arm_y - arm_h / 2 + 0.5) / 2
        ymax_frac = (arm_y + arm_h / 2 + 0.5) / 2
        for x0, x1, alpha in [(x_a1[0], x_a1[1], 0.55),
                               (x_a2[0], x_a2[1], 0.30)]:
            ax_str.axvspan(x0, x1, ymin=ymin_frac, ymax=ymax_frac,
                           alpha=alpha, color=col, zorder=2)

        if loop_len:
            ax_str.axvspan(x_lp[0], x_lp[1], ymin=ymin_frac, ymax=ymax_frac,
                           alpha=0.75, color="gold", zorder=3)
            mid_lp = (x_lp[0] + x_lp[1]) / 2
            ax_str.text(mid_lp, arm_y + 0.58, loop_seq,
                        ha="center", va="bottom", fontsize=8,
                        fontweight="bold", color="darkgoldenrod")

        mid_arm = (x_a1[0] + x_a2[1]) / 2
        marker  = "✓" if confirmed else "~"
        ax_str.text(mid_arm, arm_y + 0.95,
                    f"IDR {side}  {arm_len} bp  {identity:.0%}  {marker}",
                    ha="center", va="bottom", fontsize=8, color=col)

    hp_left  = hp.get("left")
    hp_right = hp.get("right")
    hp_desc = (f"hairpin left: {'support=' + str(hp_left['support']) if hp_left else 'no'}   "
              f"hairpin right: {'support=' + str(hp_right['support']) if hp_right else 'no'}")
    ax_str.text(seq_len / 2, arm_y - 0.68,
                f"{hp_desc}   "
                f"IDR left: {idr.get('left', {}).get('confirmed', False)}   "
                f"IDR right: {idr.get('right', {}).get('confirmed', False)}",
                ha="center", va="center", fontsize=7.5, color="dimgray")

    legend_patches = [
        mpatches.Patch(color=arm_colours["left"],  alpha=0.6, label="IDR arm (left)"),
        mpatches.Patch(color=arm_colours["right"], alpha=0.6, label="IDR arm (right)"),
        mpatches.Patch(color="gold",               alpha=0.8, label="Loop / TATA"),
    ]
    ax_str.legend(handles=legend_patches, loc="lower right",
                  fontsize=7, framealpha=0.7)

    # ── Save ──────────────────────────────────────────────────────────────────
    safe_id  = re.sub(r"[^\w.-]", "_", contig_id)
    out_path = f"{output_prefix}_{safe_id}_terminal.png"
    try:
        plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
    except Exception as e:
        print(f"[WARN] Could not save visualization: {e}", file=sys.stderr)
        plt.close(fig)
        return ""
    plt.close(fig)
    print(f"[VIZ] {contig_id}  →  {out_path}")
    return out_path


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def analyse_contig(record, args, annot_df, gfa_topology, ref_gc=None,
                   chromosome_depth: float = None,
                   bam_handle=None, amr_df: pd.DataFrame = None) -> dict:
    """Run all modules on a single contig/sequence record.

    bam_handle: optional open pysam.AlignmentFile shared across contigs.
    When provided it is used for coverage/copy-number analysis without
    re-opening the BAM per contig; the caller is responsible for closing it.

    amr_df: optional AMRFinderPlus hits for the whole assembly (from
    run_amrfinder, run once per assembly by run_single_assembly), filtered
    here to this contig's hits. None means AMR screening wasn't requested.
    """
    seq = str(record.seq)
    cid = record.id

    evidence = {}

    # FASTA header metadata (circular flag, copy number) — Module 0
    header_meta = parse_fasta_header(record.description)
    evidence["header"] = assess_header_metadata(
        header_meta, chromosome_depth, seq_len=len(seq))

    # Chromosome identification for the hard-exclude gate below. Any contig
    # above the plausible-linear-plasmid ceiling (SIZE_RANGES["general"],
    # 500 kb — the same bound _gfa_dominant_segment_is_open() already uses)
    # is excluded regardless of its circular= flag. Originally this also
    # required circular_flag is True, mirroring the chromosome-detection
    # heuristic used elsewhere for copy-number normalisation (which
    # legitimately wants a confidently-closed reference replicon) — but that
    # made the gate *weaker* exactly when it mattered most: a poorly-
    # resolved assembly that fails to close the chromosome into a circle
    # emits circular=false on a multi-Mb contig, which used to sail straight
    # through ungated *and* collect assembler_not_circular's 13 points for
    # "not circular" on top. Confirmed in practice on a reduced-assembler
    # (Flye+Plassembler-only) consensus run: two full chromosomes scored
    # 42 HIGH and 31 MEDIUM as linear-plasmid candidates before this fix.
    evidence["is_chromosome"] = len(seq) > 500_000 or "chromosome" in cid.lower()

    # Structural (hairpin, invertron end)
    evidence["hairpin_ends"] = detect_terminal_hairpins(
        seq,
        k=getattr(args, "hairpin_k", 31),
        window=getattr(args, "hairpin_window", 50_000),
        min_shared=getattr(args, "hairpin_min_shared", 25),
        edge_tol=getattr(args, "hairpin_edge_tol", 100))
    evidence["idr_tata"]        = detect_idr_tata_ends(seq)
    evidence["size"]            = classify_size(len(seq))

    # GC
    contig_gc = gc_content(seq)
    evidence["gc"] = assess_gc(contig_gc, ref_gc)

    # Gene content
    evidence["genes"] = screen_genes(annot_df, cid) if not annot_df.empty else {}

    # AMRFinderPlus hits (informational only — not scored, see module comment)
    evidence["amr"] = interpret_amr_hits(amr_df, cid)

    # PlasmidFinder no-hit flag (set externally; default False unless --no-plasmid-finder-hit)
    evidence["plasmid_finder_no_hit"] = getattr(args, "plasmid_finder_no_hit", False)

    # BAM-based — use shared handle if provided, else open from path
    bam_src = bam_handle if bam_handle is not None else args.bam
    if bam_src:
        evidence["copy_number"]   = estimate_copy_number(
            bam_src, cid, args.chromosome_contigs)
        evidence["coverage_drop"] = detect_coverage_drop_ends(
            bam_src, cid, len(seq))
    else:
        evidence["copy_number"]   = {"available": False}
        evidence["coverage_drop"] = {"available": False}

    # Asymmetric end analysis (pELF1-type, Hashimoto 2019)
    # Runs after BAM so coverage_drop is available for the left-only drop path
    evidence["asymmetric_ends"] = detect_asymmetric_ends(
        evidence["hairpin_ends"],
        evidence.get("genes", {}),
        evidence.get("coverage_drop"))

    # GFA topology
    evidence["gfa_topology"] = is_linear_in_gfa(cid, gfa_topology)

    import tempfile as _tf

    # BLAST — write query to a temp file, clean up after use
    if args.blast_db:
        with _tf.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp_fa:
            tmp_fa.write(f">{cid}\n{seq}\n")
            tmp_fa_path = tmp_fa.name
        with _tf.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp_out:
            tmp_out_path = tmp_out.name
        try:
            blast_df = run_blast(tmp_fa_path, args.blast_db, tmp_out_path,
                                 mode="blastn", identity=args.blast_identity)
            evidence["blast"] = interpret_blast_hits(blast_df)
        finally:
            for p in (tmp_fa_path, tmp_out_path):
                try:
                    os.unlink(p)
                except OSError:
                    pass
    else:
        evidence["blast"] = {"hits": 0, "linear_plasmid_hit": False}

    # SKANI — same temp-file pattern; reuse query FA if already written
    if args.skani_db:
        with _tf.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp_fa:
            tmp_fa.write(f">{cid}\n{seq}\n")
            tmp_fa_path = tmp_fa.name
        with _tf.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp_out:
            tmp_out_path = tmp_out.name
        try:
            skani_df = run_skani(tmp_fa_path, args.skani_db, tmp_out_path,
                                 min_ani=args.skani_ani, min_af=args.skani_af)
            evidence["skani"] = interpret_blast_hits(skani_df)
        finally:
            for p in (tmp_fa_path, tmp_out_path):
                try:
                    os.unlink(p)
                except OSError:
                    pass
    else:
        evidence["skani"] = {"hits": 0, "linear_plasmid_hit": False}

    # Composite score
    evidence["score"] = compute_score(evidence)

    return {"contig": cid, "length": len(seq), "gc": contig_gc, "evidence": evidence}


def run_single_assembly(args) -> list:
    """
    Run the full pipeline on one assembly FASTA (args.input) and return the
    filtered/sorted list of per-contig result dicts, as consumed by the TSV/
    JSON writers and visualize_terminal_structure in main().

    Extracted from main() (2026-07) so --hybracter-dir batch mode can call it
    once per discovered sample.
    """
    # Load sequences
    records = list(SeqIO.parse(args.input, "fasta"))
    print(f"[INFO] Loaded {len(records)} sequences from {args.input}")

    # Filter by minimum length
    records = [r for r in records if len(r.seq) >= args.min_length]
    print(f"[INFO] {len(records)} sequences ≥ {args.min_length} bp")

    # Parse annotation
    annot_df = parse_annotation(args.annot) if args.annot else pd.DataFrame()
    if not annot_df.empty:
        print(f"[INFO] Loaded annotation: {len(annot_df)} features")

    # Parse GFA — pass length map so numeric Plassembler/Flye segment names
    # are remapped to FASTA contig IDs via sequence-length matching.
    # args.gfa may be a single path or a list (e.g. Hybracter batch mode
    # passes both Flye's and Plassembler's GFA — see parse_gfa_topology).
    length_map = {r.id: len(r.seq) for r in records}
    gfa_topology = parse_gfa_topology(args.gfa, length_map) if args.gfa else {}
    if gfa_topology:
        print(f"[INFO] Loaded GFA topology: {len(gfa_topology)} segments")

    if args.gfa and getattr(args, "annotate_gfa_hairpins", False):
        gfa_files = [args.gfa] if isinstance(args.gfa, str) else list(args.gfa)
        for gfa_file in gfa_files:
            ann = annotate_gfa_hairpins(
                gfa_file, k=getattr(args, "hairpin_k", 31),
                window=getattr(args, "hairpin_window", 50_000),
                min_shared=getattr(args, "hairpin_min_shared", 25),
                edge_tol=getattr(args, "hairpin_edge_tol", 100))
            if ann["n_added"]:
                print(f"[INFO] Annotated {ann['n_added']} hairpin link(s) → {ann['out_path']}")
            else:
                print(f"[INFO] --annotate-gfa-hairpins: no terminal hairpins detected in {gfa_file}")

    # AMRFinderPlus — run once on the whole assembly (nucleotide mode), hits
    # attributed back to each contig via its own Contig id column
    amr_df = None
    if getattr(args, "amrfinder", False):
        amr_out = f"{args.output}.amrfinder.tsv"
        amr_df = run_amrfinder(
            args.input, amr_out,
            db=getattr(args, "amrfinder_db", None),
            organism=getattr(args, "amrfinder_organism", None),
            threads=getattr(args, "amrfinder_threads", 4))
        if not amr_df.empty:
            print(f"[INFO] AMRFinderPlus: {len(amr_df)} hit(s) across "
                  f"{amr_df['contig'].nunique()} contig(s) → {amr_out}")
        else:
            print("[INFO] AMRFinderPlus: no hits (or tool unavailable)")

    # Parse header metadata for all records (needed for GC reference and
    # chromosome identification before per-contig analysis)
    all_headers = {r.id: parse_fasta_header(r.description) for r in records}

    # Compute reference GC — prefer chromosomal contig(s); fall back to median
    if args.ref_gc is None:
        chrom_gc = [gc_content(str(r.seq)) for r in records
                    if (all_headers[r.id]["circular"] is True and len(r.seq) > 500_000)
                    or "chromosome" in r.id.lower()]
        if chrom_gc:
            ref_gc = float(np.mean(chrom_gc))
            print(f"[INFO] Reference GC (chromosome mean) = {ref_gc:.2f}%")
        else:
            all_gc = [gc_content(str(r.seq)) for r in records]
            ref_gc = float(np.median(all_gc)) if all_gc else None
            print(f"[INFO] Reference GC (median all contigs) = {ref_gc:.2f}%")
    else:
        ref_gc = args.ref_gc

    # Identify chromosome contig(s) for copy-number normalisation
    chromosome_depth = None
    for r in sorted(records, key=lambda x: len(x.seq), reverse=True):
        hm = all_headers[r.id]
        is_chrom = (hm["circular"] is True and len(r.seq) > 500_000) or \
                   "chromosome" in r.id.lower()
        if is_chrom and hm["depth"] is not None:
            chromosome_depth = hm["depth"]
            print(f"[INFO] Chromosome depth from header: {chromosome_depth:.1f}x "
                  f"({r.id})")
            break

    # Auto-detect chromosome contigs for copy-number normalisation if not given
    if not args.chromosome_contigs:
        for r in records:
            hm = all_headers[r.id]
            if ((hm["circular"] is True and len(r.seq) > 500_000)
                    or "chromosome" in r.id.lower()):
                args.chromosome_contigs.append(r.id)
        if args.chromosome_contigs:
            print(f"[INFO] Auto-detected chromosome: "
                  f"{', '.join(args.chromosome_contigs)}")

    # Print header metadata summary
    print(f"\n[INFO] FASTA header metadata:")
    for r in records:
        hm = all_headers[r.id]
        circ = {True: "circular=true", False: "circular=false",
                None: "circular=?"}[hm["circular"]]
        cn = f"  CN={hm['copy_number']}x" if hm["copy_number"] else ""
        depth = f"  depth={hm['depth']}x" if hm["depth"] else ""
        print(f"  {r.id:20s}  {len(r.seq):>10,} bp  {circ}{cn}{depth}")
    print()

    # Open BAM once for the whole assembly (avoids per-contig open/close overhead)
    _bam_handle = None
    if args.bam:
        try:
            import pysam as _pysam
            _bam_handle = _pysam.AlignmentFile(args.bam, "rb")
        except Exception as _e:
            print(f"[WARN] Cannot pre-open BAM ({_e}); will open per-contig.",
                  file=sys.stderr)

    # Analyse each contig
    results = []
    for rec in records:
        print(f"  → Analysing {rec.id} ({len(rec.seq):,} bp) ...", end=" ", flush=True)
        res = analyse_contig(rec, args, annot_df, gfa_topology, ref_gc,
                             chromosome_depth=chromosome_depth,
                             bam_handle=_bam_handle, amr_df=amr_df)
        score = res["evidence"]["score"]["total_score"]
        conf  = res["evidence"]["score"]["confidence"]
        circ  = all_headers[rec.id]["circular"]
        circ_tag = " [circular=true]" if circ else (
                   " [circular=false]" if circ is False else "")
        print(f"score={score}  [{conf}]{circ_tag}")
        results.append(res)

    if _bam_handle is not None:
        _bam_handle.close()

    # Filter and sort
    results = sorted(results, key=lambda x: x["evidence"]["score"]["total_score"], reverse=True)
    results = [r for r in results if r["evidence"]["score"]["total_score"] >= args.min_score]
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Identify linear plasmids in assembled sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-i", "--input",   default=None,
                        help="Input FASTA file (assembled sequences). "
                             "Mutually exclusive with --hybracter-dir.")
    parser.add_argument("--hybracter-dir", dest="hybracter_dir", default=None,
                        help="Hybracter output directory root. Auto-discovers every "
                             "sample's FINAL_OUTPUT FASTA, per_contig_stats.tsv "
                             "(for chromosome/plasmid separation), Flye/Plassembler "
                             "GFA, and Hybracter's own QC'd long reads (for "
                             "auto-BAM-mapping), and runs the pipeline once per "
                             "sample. Mutually exclusive with -i/--input; --bam/"
                             "--longread-fastq/--shortread-* are not supported "
                             "alongside it (ambiguous across samples).")
    parser.add_argument("-o", "--output",  default="linear_plasmid_report",
                        help="Output prefix (default: linear_plasmid_report)")
    parser.add_argument("--bam",           help="Pre-mapped BAM file of reads aligned to assembly")
    parser.add_argument("--longread-fastq", dest="longread_fastq",
                        help="Long-read FASTQ (or .fastq.gz) — mapped automatically "
                             "with minimap2 to produce a BAM for IDR and coverage analysis")
    parser.add_argument("--longread-preset", dest="longread_preset",
                        default="map-ont",
                        choices=["map-ont", "map-pb", "map-hifi"],
                        help="minimap2 preset for long reads "
                             "(map-ont Nanopore, map-pb PacBio CLR, map-hifi PacBio HiFi; "
                             "default: map-ont)")
    parser.add_argument("--longread-threads", dest="longread_threads",
                        type=int, default=4,
                        help="Threads for minimap2 / samtools (default: 4)")
    parser.add_argument("--shortread-fastq", dest="shortread_fastq",
                        metavar="INTERLEAVED.fastq.gz",
                        help="Interleaved paired-end Illumina FASTQ — mapped with "
                             "minimap2 sr preset; produces <output-prefix>.short.sorted.bam")
    parser.add_argument("--shortread-r1", dest="shortread_r1",
                        metavar="R1.fastq.gz",
                        help="Illumina R1 FASTQ (paired-end, separate files)")
    parser.add_argument("--shortread-r2", dest="shortread_r2",
                        metavar="R2.fastq.gz",
                        help="Illumina R2 FASTQ (paired-end, separate files)")
    parser.add_argument("--shortread-threads", dest="shortread_threads",
                        type=int, default=4,
                        help="Threads for minimap2 / samtools short-read mapping "
                             "(default: 4)")
    parser.add_argument("--gfa",           help="Assembly graph GFA file (Unicycler/Flye/Plassembler)")
    parser.add_argument("--hairpin-k", dest="hairpin_k", type=int, default=31,
                        help="k-mer size for terminal hairpin detection (default: 31)")
    parser.add_argument("--hairpin-window", dest="hairpin_window", type=int, default=50_000,
                        help="Terminal window (bp) scanned at each contig end for "
                             "hairpin fold-backs (default: 50000)")
    parser.add_argument("--hairpin-min-shared", dest="hairpin_min_shared", type=int, default=25,
                        help="Min supporting reverse-complement k-mer pairs to call "
                             "a terminal hairpin (default: 25)")
    parser.add_argument("--hairpin-edge-tol", dest="hairpin_edge_tol", type=int, default=100,
                        help="Max gap (bp) from the terminus to still call a fold-back "
                             "TERMINAL (default: 100)")
    parser.add_argument("--annotate-gfa-hairpins", action="store_true",
                        dest="annotate_gfa_hairpins",
                        help="Write <gfa-stem>.hairpins.gfa alongside --gfa with "
                             "Autocycler-style hairpin links (L n + n -) added for "
                             "any segment with a terminal fold-back. Diagnostic "
                             "Bandage-visualization aid; does not affect scoring.")
    parser.add_argument("--annot",         help="GFF3/TSV annotation file (Prokka)")
    parser.add_argument("--amrfinder", action="store_true",
                        help="Screen for antimicrobial resistance genes with NCBI "
                             "AMRFinderPlus (run once on the whole assembly in "
                             "nucleotide mode; requires 'amrfinder' on PATH and its "
                             "database installed). Reported per contig in the TSV/"
                             "JSON output — informational only, not scored.")
    parser.add_argument("--amrfinder-db", dest="amrfinder_db", default=None,
                        help="Custom AMRFinderPlus database directory (passed as -d). "
                             "Default: whatever database amrfinder finds on its own.")
    parser.add_argument("--amrfinder-organism", dest="amrfinder_organism", default=None,
                        help="Organism for AMRFinderPlus point-mutation screening "
                             "(passed as -O), e.g. Enterococcus_faecium. "
                             "See `amrfinder -l` for the supported organism list.")
    parser.add_argument("--amrfinder-threads", dest="amrfinder_threads",
                        type=int, default=4,
                        help="Threads for amrfinder (default: 4)")
    parser.add_argument("--blast-db",      help="BLAST database path (e.g. PLSDB)")
    parser.add_argument("--blast-identity", type=float, default=70.0,
                        help="Minimum BLAST %%identity (default: 70)")
    parser.add_argument("--skani-db",
                        help="SKANI database: pre-sketched directory or FASTA file. "
                             "For large FASTA DBs, pre-sketch with: "
                             "skani sketch -l db.fna -o db_dir/")
    parser.add_argument("--skani-ani", type=float, default=90.0,
                        help="Minimum SKANI ANI %%%% (default: 90)")
    parser.add_argument("--skani-af",  type=float, default=0.05,
                        help="Minimum SKANI query alignment fraction (default: 0.05)")
    parser.add_argument("--chromosome-contigs", nargs="*", default=[],
                        help="Contig IDs of chromosome (for copy number normalisation)")
    parser.add_argument("--ref-gc", type=float, default=None,
                        help="Reference GC%% of chromosome (for GC deviation scoring)")
    parser.add_argument("--min-length", type=int, default=5_000,
                        help="Minimum contig length to analyse (default: 5000 bp)")
    parser.add_argument("--min-score", type=int, default=CONFIDENCE_THRESHOLDS["LOW"],
                        help=f"Minimum score to report (default: {CONFIDENCE_THRESHOLDS['LOW']}, "
                             f"the LOW confidence threshold)")
    parser.add_argument("--no-plasmid-finder-hit", action="store_true",
                        dest="plasmid_finder_no_hit",
                        help="Flag: PlasmidFinder found no replicon hit "
                             "(positive linear plasmid indicator, Hashimoto 2019)")
    parser.add_argument("--json", action="store_true",
                        help="Also write detailed JSON output")
    parser.add_argument("--visualize", action="store_true",
                        help="Generate Figure-3-style left-end visualization (PNG) for "
                             "contigs with IDR/TATA evidence (requires matplotlib)")
    args = parser.parse_args()

    # ── Validate -i / --hybracter-dir mutual exclusivity ──────────────────────
    if bool(args.input) == bool(args.hybracter_dir):
        parser.error("Exactly one of -i/--input or --hybracter-dir is required.")

    # ── Validate short-read arguments ─────────────────────────────────────────
    _sr_interleaved = getattr(args, "shortread_fastq", None)
    _sr_r1          = getattr(args, "shortread_r1", None)
    _sr_r2          = getattr(args, "shortread_r2", None)

    if _sr_interleaved and (_sr_r1 or _sr_r2):
        parser.error("--shortread-fastq (interleaved) and --shortread-r1/r2 "
                     "are mutually exclusive.")
    if bool(_sr_r1) != bool(_sr_r2):
        parser.error("--shortread-r1 and --shortread-r2 must be supplied together.")

    batch_mode = bool(args.hybracter_dir)
    seq_maps = {}    # sample (or None in single mode) -> {contig_id: seq}, for --visualize
    sample_bams = {}  # sample (or None in single mode) -> bam path used, for --visualize

    if batch_mode:
        if args.bam or args.longread_fastq or _sr_interleaved or _sr_r1:
            parser.error("--bam/--longread-fastq/--shortread-* are not supported "
                         "alongside --hybracter-dir (ambiguous across multiple "
                         "samples) — Hybracter's own QC'd long reads are "
                         "auto-mapped per sample instead.")

        samples = discover_hybracter_samples(args.hybracter_dir)
        if not samples:
            print(f"[ERROR] No Hybracter samples found under {args.hybracter_dir} "
                  f"(expected FINAL_OUTPUT/{{complete,incomplete}}/*_final.fasta)",
                  file=sys.stderr)
            sys.exit(1)
        print(f"[INFO] Discovered {len(samples)} Hybracter sample(s) under "
              f"{args.hybracter_dir}\n")

        results = []
        for s in samples:
            print(f"{'='*60}\n[SAMPLE] {s['sample']}  "
                  f"({'complete' if s['complete'] else 'incomplete'})\n{'='*60}")

            sample_args = copy.deepcopy(args)
            sample_args.input = s["fasta"]
            if not sample_args.gfa:
                # Pass both — Flye's graph typically only covers the
                # chromosome and any large plasmids it resolved; smaller/
                # low-coverage plasmids are often Plassembler-only and never
                # appear in Flye's graph at all. parse_gfa_topology merges
                # topology from whichever file(s) actually cover each contig.
                sample_args.gfa = [g for g in (s["flye_gfa"], s["plassembler_gfa"]) if g]

            if s["contig_stats_tsv"]:
                stats = parse_hybracter_contig_stats(s["contig_stats_tsv"])
                if stats and not sample_args.chromosome_contigs:
                    sample_args.chromosome_contigs = [
                        c for c, v in stats.items() if v["contig_type"] == "chromosome"]
                    if sample_args.chromosome_contigs:
                        print(f"[INFO] Chromosome contig(s) from per_contig_stats.tsv: "
                              f"{', '.join(sample_args.chromosome_contigs)}")

            if s["longread_fastq"]:
                bam_out = f"{args.output}_{s['sample']}.sorted.bam"
                mapped = map_longreads(
                    fastq          = s["longread_fastq"],
                    assembly_fasta = s["fasta"],
                    output_bam     = bam_out,
                    preset         = args.longread_preset,
                    threads        = args.longread_threads,
                )
                if mapped:
                    sample_args.bam = mapped
                else:
                    print(f"[WARN] Long-read mapping failed for {s['sample']}; "
                          f"continuing without BAM.", file=sys.stderr)

            sample_results = run_single_assembly(sample_args)
            for r in sample_results:
                r["sample"] = s["sample"]
            results.extend(sample_results)
            seq_maps[s["sample"]] = {rec.id: str(rec.seq)
                                     for rec in SeqIO.parse(s["fasta"], "fasta")}
            sample_bams[s["sample"]] = sample_args.bam or ""
            print()

        # Re-sort the combined multi-sample results by score
        results = sorted(results, key=lambda x: x["evidence"]["score"]["total_score"], reverse=True)

    else:
        # ── Auto-map long reads if --longread-fastq supplied ─────────────────
        if args.longread_fastq:
            if args.bam:
                print("[WARN] Both --bam and --longread-fastq supplied; "
                      "--bam takes precedence, skipping mapping.", file=sys.stderr)
            else:
                bam_out = f"{args.output}.sorted.bam"
                mapped = map_longreads(
                    fastq          = args.longread_fastq,
                    assembly_fasta = args.input,
                    output_bam     = bam_out,
                    preset         = args.longread_preset,
                    threads        = args.longread_threads,
                )
                if mapped:
                    args.bam = mapped
                else:
                    print("[WARN] Long-read mapping failed; continuing without BAM.",
                          file=sys.stderr)

        # ── Auto-map short reads if --shortread-fastq / --shortread-r1/r2 given
        if _sr_interleaved or _sr_r1:
            if args.bam:
                print("[WARN] --bam already set; skipping short-read mapping "
                      "(--bam takes precedence).", file=sys.stderr)
            else:
                bam_out = f"{args.output}.short.sorted.bam"
                mapped = map_shortreads(
                    assembly_fasta = args.input,
                    output_bam     = bam_out,
                    interleaved    = _sr_interleaved,
                    r1             = _sr_r1,
                    r2             = _sr_r2,
                    threads        = args.shortread_threads,
                )
                if mapped:
                    args.bam = mapped
                else:
                    print("[WARN] Short-read mapping failed; continuing without BAM.",
                          file=sys.stderr)

        results = run_single_assembly(args)
        seq_maps[None] = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.input, "fasta")}
        sample_bams[None] = args.bam or ""

    # Build TSV output
    rows = []
    for r in results:
        ev = r["evidence"]
        sc = ev["score"]
        genes = ev.get("genes", {})
        hdr = ev.get("header", {})
        row = {}
        if batch_mode:
            row["sample"] = r.get("sample", "")
        row.update({
            "contig":            r["contig"],
            "length_bp":         r["length"],
            "gc_pct":            r["gc"],
            "score":             sc["total_score"],
            "max_score":         sc["max_possible"],
            "pct_max":           sc["percent"],
            "confidence":        sc["confidence"],
            "gate_reason":       sc.get("gate_reason") or "",
            # Header metadata columns
            "circular_flag":     hdr.get("circular_flag", ""),
            "header_copy_number": hdr.get("copy_number", ""),
            "header_depth":      hdr.get("depth", ""),
            "hairpin_left_support":  (ev.get("hairpin_ends", {}).get("left") or {}).get("support", ""),
            "hairpin_right_support": (ev.get("hairpin_ends", {}).get("right") or {}).get("support", ""),
            "idr_confirmed":     ev.get("idr_tata", {}).get("confirmed_idr_tata", ""),
            "idr_left_arm_bp":   ev.get("idr_tata", {}).get("left", {}).get("arm_len", ""),
            "idr_left_loop":     ev.get("idr_tata", {}).get("left", {}).get("loop_seq", ""),
            "idr_right_arm_bp":  ev.get("idr_tata", {}).get("right", {}).get("arm_len", ""),
            "idr_right_loop":    ev.get("idr_tata", {}).get("right", {}).get("loop_seq", ""),
            "asymmetric_ends":   ev.get("asymmetric_ends", {}).get("asymmetric_pelf_type", ""),
            "invertron_tp_gene": ev.get("asymmetric_ends", {}).get("has_tp_gene", ""),
            "size_families":     "|".join(ev.get("size", {}).get("matching_families", [])),
            "gc_deviation":      ev.get("gc", {}).get("deviation", ""),
            "partition_genes":   "|".join(genes.get("partition_genes", [])),
            "tas_systems":       "|".join(genes.get("tas_systems", [])),
            "resistance_genes":  "|".join(genes.get("resistance_genes", [])),
            "is_elements":       "|".join(genes.get("is_elements", [])),
            "enterococcal_m":    "|".join(genes.get("enterococcal_markers", [])),
            "amr_hit_count":     ev.get("amr", {}).get("hits", ""),
            "amr_genes":         "|".join(ev.get("amr", {}).get("genes", [])),
            "amr_classes":       "|".join(ev.get("amr", {}).get("classes", [])),
            "blast_hit":         ev.get("blast", {}).get("linear_plasmid_db_hit", False)
                                 or ev.get("blast", {}).get("linear_plasmid_hit", False),
            "blast_top_hit":     ev.get("blast", {}).get("top_hit", ""),
            "blast_identity":    ev.get("blast", {}).get("best_identity", ""),
            "blast_coverage":    ev.get("blast", {}).get("best_coverage", ""),
            "skani_hit":         ev.get("skani", {}).get("linear_plasmid_db_hit", False)
                                 or ev.get("skani", {}).get("linear_plasmid_hit", False),
            "skani_top_hit":     ev.get("skani", {}).get("top_hit", ""),
            "skani_identity":    ev.get("skani", {}).get("best_identity", ""),
            "skani_coverage":    ev.get("skani", {}).get("best_coverage", ""),
            "gfa_topology":      ev.get("gfa_topology", {}).get("topology", ""),
            "gfa_linear":        ev.get("gfa_topology", {}).get("consistent_with_linear", ""),
            "copy_number":       ev.get("copy_number", {}).get("estimated_copy_number", ""),
            "score_breakdown":   str(sc["breakdown"]),
        })
        rows.append(row)

    tsv_path = f"{args.output}.tsv"
    df_out = pd.DataFrame(rows)
    df_out.to_csv(tsv_path, sep="\t", index=False)
    print(f"\n[OUTPUT] TSV report → {tsv_path}")
    print(f"[OUTPUT] {len(rows)} contigs reported (score ≥ {args.min_score})")

    # Summary to stdout. df_out can be empty with no columns at all (e.g. every
    # contig gated to 0 by compute_score's hard gates and filtered by
    # --min-score) — pd.DataFrame([]) has no "confidence" column to index.
    for conf in ["HIGH", "MEDIUM", "LOW"] if not df_out.empty else []:
        subset = df_out[df_out["confidence"] == conf]
        if len(subset) > 0:
            print(f"\n  {conf} confidence ({len(subset)} contigs):")
            for _, row in subset.iterrows():
                sample_tag = f"{row['sample']}: " if batch_mode else ""
                print(f"    {sample_tag}{row['contig']:30s}  {row['length_bp']:>10,} bp  "
                      f"score={row['score']:>3}  GC={row['gc_pct']}%  "
                      f"hairpin(L/R)={row['hairpin_left_support']}/{row['hairpin_right_support']}  "
                      f"asymmetric={row['asymmetric_ends']}")

    # JSON output
    if args.json:
        json_path = f"{args.output}.json"
        with open(json_path, "w") as fh:
            json.dump(results, fh, indent=2, default=str)
        print(f"[OUTPUT] JSON detail    → {json_path}")

    if args.visualize:
        viz_contigs = [r for r in results
                       if r["evidence"]["score"]["confidence"] in ("HIGH", "MEDIUM")
                       or r["evidence"].get("idr_tata", {}).get("either_end")]
        if viz_contigs:
            print(f"\n[VIZ] Generating terminal-structure plots for "
                  f"{len(viz_contigs)} contig(s)…")
            for r in viz_contigs:
                sample_key = r.get("sample")  # None in single mode
                seq_map = seq_maps.get(sample_key, {})
                out_prefix = f"{args.output}_{sample_key}" if sample_key else args.output
                visualize_terminal_structure(
                    contig_id     = r["contig"],
                    seq           = seq_map.get(r["contig"], ""),
                    evidence      = r["evidence"],
                    bam_file      = sample_bams.get(sample_key, ""),
                    output_prefix = out_prefix,
                )
        else:
            print("[VIZ] No contigs with IDR/TATA evidence or HIGH/MEDIUM confidence "
                  "— no plots generated.")

    print("\n[DONE]")


if __name__ == "__main__":
    main()
