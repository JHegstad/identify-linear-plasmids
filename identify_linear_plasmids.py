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

# Hairpin central motif (Boumamoud 2022 / Hashimoto 2019)
# pELF1 left end: ~5 kb inverted tandem repeats interspaced by 5'-TATA-3'
HAIRPIN_CENTRAL_MOTIF = "TATA"
HAIRPIN_IDR_WINDOW    = 5_000   # bp; inverted direct repeat flanking the TATA loop (Hashimoto 2019)
HAIRPIN_IDR_MIN_LEN   = 100     # minimum IDR arm length to count as evidence

# Terminal TATA window: if TATA is within this many bp of the contig start/end,
# treat it as "tata_at_terminus" — the assembler stopped at the palindrome centre
# (Hashimoto 2019, Fig 3A: coverage drops to zero at the hairpin tip).
TATA_TERMINUS_WINDOW = 10   # bp from contig end to still count as terminal TATA

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
SCORING_WEIGHTS = {
    "hairpin_end":              20,   # Hairpin/palindromic end structure
    "idr_tata_detected":        18,   # Inverted Direct Repeat + TATA loop (Hashimoto 2019 pELF1)
    "tata_at_terminus":         15,   # TATA at contig tip = assembler stopped at palindrome centre
    "idr_in_reads":             22,   # Full arm1—TATA—arm2 palindrome found in a raw long read
    "asymmetric_ends":          15,   # One hairpin + one invertron end (pELF1-type, Hashimoto 2019)
    "invertron_tp_gene":        12,   # Terminal protein gene present (invertron end)
    "size_range":               10,   # Size consistent with known linear plasmids
    "gc_deviation":              8,   # GC% lower than typical chromosome
    "par_system":               10,   # ParA/ParB partition system
    "repb_rep2":                 8,   # repB/Rep_2 superfamily replication gene (Hashimoto 2019)
    "ftsK_parA_repB_combo":     12,   # ftsK + parA + repB co-occurrence (pELF1 signature)
    "tas_system":                8,   # Toxin-antitoxin system

    "is_elements":               5,   # IS elements associated with linear plasmids

    "blast_hit":                15,   # BLAST hit to known linear plasmid db (PLSDB etc.)
    "blast_hit_linear_db":      20,   # BLAST hit explicitly to a *linear* plasmid sequence
    "skani_hit":                15,   # SKANI hit to known linear plasmid db
    "skani_hit_linear_db":      20,   # SKANI hit explicitly to a *linear* plasmid sequence
    "plasmid_finder_no_hit":     8,   # No PlasmidFinder hit (novel rep = typical of linear, H.2019)
    "coverage_drop_ends":       10,   # Coverage drop at contig ends (hairpin inaccessibility, Hashimoto 2019)
    "assembly_graph_linear":    15,   # Assembly graph linear topology
    "enterococcal_markers":     12,   # pELF-specific markers
    "copy_number_low":           5,   # ~1 copy/cell (characteristic)
    "self_complement_end":      18,   # Reverse complement overlap at ends
    # ── FASTA header metadata (Unicycler / Plassembler / Hybracter output) ──
    "assembler_not_circular":   30,   # Explicit circular=false in header
    "circular_flag_absent":     20,   # No circular=true flag (Unicycler only sets it when confirmed)
    "header_copy_number":        8,   # Copy number in header consistent with plasmid (~0.3–4x)
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

# Confidence thresholds
CONFIDENCE_THRESHOLDS = {
    "HIGH":   70,   # ≥70 points → likely linear plasmid
    "MEDIUM": 40,   # 40–69 points → possible linear plasmid
    "LOW":    15,   # 15–39 points → weak evidence
}


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

    Three tiers of assembler topology signal:
      assembler_not_circular  (weight 30): explicit circular=false in header
      circular_flag_absent    (weight 20): no circular=true AND size < 500 kb
                                           (Unicycler only writes circular=true when
                                            the contig circularises; absence = linear)
      circular=true           → disqualifies gene-based scoring; not a positive signal

    header_copy_number (weight 8): CN 0.3–4x consistent with plasmid
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
# MODULE 2: SELF-COMPLEMENT END ANALYSIS (hairpin / inverse duplication)
# ─────────────────────────────────────────────────────────────────────────────

def detect_self_complement_ends(seq: str, window: int = 500, min_match: int = 40) -> dict:
    """
    Check whether the terminal window of a sequence contains self-complementary
    structure (inverted duplication / hairpin telomere).

    pELF plasmids have 5'-TATA-3' at the centre of the hairpin (Boumamoud 2022 /
    Hashimoto 2019). This function checks for a shorter diagnostic window.
    """
    if len(seq) < 2 * window:
        return {"found": False, "score": 0}

    left  = seq[:window].upper()
    right = seq[-window:].upper()
    rc_right = str(Seq(right).reverse_complement())

    # Count matching positions
    matches = sum(1 for a, b in zip(left, rc_right) if a == b)
    identity = matches / window

    # Check for TATA motif very close to the ends (within 30 bp of start/end)
    # Requiring proximity avoids false positives from random TATA occurrences
    tata_left  = HAIRPIN_CENTRAL_MOTIF in left[:30]
    tata_right = HAIRPIN_CENTRAL_MOTIF in right[-30:]

    return {
        "found": identity >= 0.70,   # 70 % identity threshold for hairpin evidence
        "identity": round(identity, 3),
        "tata_motif_left":  tata_left,
        "tata_motif_right": tata_right,
        "score": round(identity * 100, 1),
    }


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 2b: INVERTED DIRECT REPEAT (IDR) + TATA LOOP DETECTION  (Hashimoto 2019)
# ─────────────────────────────────────────────────────────────────────────────

def detect_idr_tata(seq: str, window: int = 15_000, min_arm: int = 500,
                    identity_threshold: float = 0.75) -> dict:
    """
    Detect ~5 kb Inverted Direct Repeats (IDRs) interspaced by a 5'-TATA-3' motif
    at the left end of pELF1-type plasmids (Hashimoto 2019, Fig 3).

    Structure (Hashimoto 2019 Fig 3A):
        [→ IDR_arm_1 (~5 kb)] [TATA] [← IDR_arm_2 (~5 kb) = rc(IDR_arm_1)]

    The two arms are "completely identical inverted direct repeats" — arm2 is the
    reverse complement of arm1 and they point toward each other (converging arrows
    in Fig 3A).  The TATA tetranucleotide sits at the centre of the palindrome.

    Because "TATA" is a common tetranucleotide, the naive approach of taking the
    first occurrence gives a trivially short arm.  Instead this function:
      1. Examines the first (and last) <window> bp (default 15 kb to cover 2×5 kb arms).
      2. Iterates over ALL TATA positions in that region.
      3. For each candidate, computes arm identity = fraction matching between
         rc(arm1) and arm2 of the same length.
      4. Returns the candidate with the longest arm that meets the identity threshold.

    Parameters
    ----------
    window           : bp to examine at each end (default 15 000)
    min_arm          : minimum IDR arm length to consider (default 500 bp)
    identity_threshold: minimum fraction identity to call a match (default 0.75)

    Returns dict: found, arm_length, tata_position, identity, end (left/right/both).
    """
    motif     = HAIRPIN_CENTRAL_MOTIF          # "TATA"
    motif_len = len(motif)

    def _best_idr_in_region(region: str) -> dict:
        """
        Try every TATA occurrence; return the best (longest arm, ≥ identity_threshold).

        Two detection modes:
          A. Full IDR: arm1 (before TATA) vs rc(arm2) (after TATA), both ≥ min_arm bp.
             This is the ideal case when the whole palindrome is assembled.
          B. Terminal TATA: TATA is within TATA_TERMINUS_WINDOW bp of the region start.
             This happens when the assembler stopped at the palindrome centre — the
             mirror arm before the TATA was never sequenced (coverage drops to zero
             at the hairpin tip, Hashimoto 2019 Fig 3A).  Reported separately so it
             can be scored even without a full arm comparison.
        """
        best = {"found": False, "arm_length": 0, "identity": 0.0, "tata_pos": -1,
                "tata_at_terminus": False}

        # Mode B: check whether TATA is right at the contig tip
        first_tata = region.find(motif)
        if 0 <= first_tata <= TATA_TERMINUS_WINDOW:
            best["tata_at_terminus"] = True
            best["tata_pos"]         = first_tata
            # Still try Mode A from this position (arm_len = first_tata, likely 0-10 bp,
            # so will fail min_arm — that's expected; terminus flag carries the evidence)

        # Mode A: full palindrome search over all TATA positions
        start = 0
        while True:
            pos = region.find(motif, start)
            if pos == -1:
                break
            start = pos + 1

            arm_len = pos            # arm1 = region[:pos]
            if arm_len < min_arm:
                continue

            arm2_start = pos + motif_len
            arm2_end   = arm2_start + arm_len
            if arm2_end > len(region):
                continue

            arm1    = region[:arm_len]
            arm2    = region[arm2_start:arm2_end]
            rc_arm1 = str(Seq(arm1).reverse_complement())

            matches  = sum(1 for a, b in zip(rc_arm1, arm2) if a == b)
            identity = matches / arm_len

            if identity >= identity_threshold and arm_len > best["arm_length"]:
                best.update({
                    "found":      True,
                    "arm_length": arm_len,
                    "identity":   round(identity, 4),
                    "tata_pos":   pos,
                })

        return best

    n      = len(seq)
    chk    = min(window, n // 2)
    left_r  = seq[:chk].upper()
    right_r = seq[max(0, n - chk):].upper()

    left_result  = _best_idr_in_region(left_r)
    right_result = _best_idr_in_region(right_r)

    results    = {"left": left_result, "right": right_result}
    found_ends    = [e for e, r in results.items() if r["found"]]
    terminus_ends = [e for e, r in results.items() if r.get("tata_at_terminus")]

    # Best-arm summary across both ends (for scoring)
    best_arm = max((r["arm_length"] for r in results.values() if r["found"]),
                   default=0)

    return {
        "found":            len(found_ends) > 0,
        "found_ends":       found_ends,
        "best_arm_bp":      best_arm,
        # Terminal TATA: assembler stopped at palindrome centre; no arm comparison possible
        "tata_at_terminus": len(terminus_ends) > 0,
        "terminus_ends":    terminus_ends,
        "details":          results,
    }


def detect_asymmetric_ends(sc_result: dict, idr_result: dict,
                           gene_hits: dict) -> dict:
    """
    pELF1 has ASYMMETRIC ends: left = hairpin, right = invertron (terminal protein).
    This is a unique indicator of the pELF lineage (Hashimoto 2019).

    Detects asymmetry when:
      - one end shows IDR/TATA structure AND
      - terminal protein gene is annotated
    """
    has_hairpin = idr_result.get("found") or sc_result.get("found")
    has_tp_gene = bool(gene_hits.get("terminal_protein_genes"))

    # Asymmetric: one hairpin end + one invertron end (no symmetric TIR expected for pELF)
    asymmetric = has_hairpin and has_tp_gene

    return {
        "asymmetric_pelf_type": asymmetric,
        "has_hairpin_end":      has_hairpin,
        "has_tp_gene":          has_tp_gene,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 2c: COVERAGE DROP AT ENDS  (Hashimoto 2019)
# ─────────────────────────────────────────────────────────────────────────────

def detect_coverage_drop_ends(bam_file: str, contig_id: str,
                               contig_len: int, end_window: int = 5000) -> dict:
    """
    Hairpin structures at linear plasmid ends reduce Illumina short-read coverage
    because the hairpin fold prevents adapter ligation and PCR amplification
    (Hashimoto 2019, Fig 3A).

    Detects: mean depth in terminal <end_window> bp < 50% of plasmid body depth.
    """
    try:
        import pysam
    except ImportError:
        return {"available": False, "message": "pysam not installed"}

    if not bam_file or not os.path.exists(bam_file):
        return {"available": False, "message": "BAM not found"}

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")

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

        bam.close()
        return {
            "available":          True,
            "left_depth_ratio":   round(left_ratio,  3),
            "right_depth_ratio":  round(right_ratio, 3),
            "body_mean_depth":    round(avg_body, 1),
            # Drop below 50% at either end = evidence of hairpin inaccessibility
            "left_drop":  left_ratio  < 0.50,
            "right_drop": right_ratio < 0.50,
            "consistent_with_linear": (left_ratio < 0.50 or right_ratio < 0.50),
        }
    except Exception as e:
        return {"available": False, "message": str(e)}


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 2d: IDR/TATA PALINDROME DETECTION IN RAW LONG READS  (Hashimoto 2019)
# ─────────────────────────────────────────────────────────────────────────────

def detect_idr_in_reads(bam_file: str, contig_id: str, contig_len: int,
                         end_window: int = 15_000, min_arm: int = 500,
                         max_reads: int = 500) -> dict:
    """
    Scan long reads mapping near the contig ends for the full arm1—TATA—arm2
    palindrome (Hashimoto 2019, Fig 3).

    Why this is needed
    ------------------
    Short-read assembly truncates at the palindrome centre: the hairpin fold
    prevents adapter ligation so read coverage drops to zero at the tip, and the
    assembled contig starts at the TATA itself (tata_at_terminus).  A single long
    read that begins at the physical left end of the plasmid carries:

        [IDR arm1 (~5 kb)] [TATA] [IDR arm2 = rc(arm1) (~5 kb)] [plasmid body ...]

    detect_idr_tata() finds this structure correctly when the TATA is at position
    ~5000 in the read rather than at position 0.

    Strategy
    --------
    Fetch reads whose alignment START is within end_window of contig position 0
    (left end) or whose alignment END is within end_window of contig_len (right
    end).  For each qualifying read call detect_idr_tata() on the full query
    sequence.  Report the best (longest-arm) hit found.

    Parameters
    ----------
    end_window  : how far from the contig end to look for reads (default 15 000 bp)
    min_arm     : minimum IDR arm length to call a match (default 500 bp)
    max_reads   : cap on reads examined per end to bound runtime (default 500)

    Returns
    -------
    dict with keys:
      available        : bool  (False when pysam missing or BAM absent)
      found            : bool  (True if full palindrome detected in ≥1 read)
      best_arm_bp      : int   arm length of best hit (0 if not found)
      best_identity    : float identity of best hit
      best_read        : str   read name of best hit
      end              : str   "left" | "right" | None
      reads_checked    : int
    """
    try:
        import pysam
    except ImportError:
        return {"available": False, "message": "pysam not installed"}

    if not bam_file or not os.path.exists(bam_file):
        return {"available": False, "message": "BAM file not found"}

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")

        best = {
            "available":     True,
            "found":         False,
            "best_arm_bp":   0,
            "best_identity": 0.0,
            "best_read":     None,
            "end":           None,
            "reads_checked": 0,
        }

        def _check_reads_near(fetch_start: int, fetch_end: int, end_label: str):
            checked = 0
            for read in bam.fetch(contig_id, fetch_start, fetch_end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                qs = read.query_sequence
                if not qs or len(qs) < 2 * min_arm + 4:
                    continue

                # Run full palindrome detection on the raw read sequence.
                # Pass window = full read length so the palindrome is not
                # clipped if it lies away from the read ends.
                r = detect_idr_tata(qs, window=len(qs), min_arm=min_arm)

                if r["found"] and r["best_arm_bp"] > best["best_arm_bp"]:
                    # Extract identity from whichever end had the best match
                    best_det = max(
                        (r["details"][e] for e in r["found_ends"]),
                        key=lambda d: d["arm_length"],
                    )
                    best.update({
                        "found":         True,
                        "best_arm_bp":   r["best_arm_bp"],
                        "best_identity": best_det["identity"],
                        "best_read":     read.query_name,
                        "end":           end_label,
                    })

                checked += 1
                best["reads_checked"] += 1
                if checked >= max_reads:
                    break

        # Left end: reads that start within end_window of contig position 0
        _check_reads_near(0, end_window, "left")
        # Right end: reads that end within end_window of contig_len
        _check_reads_near(max(0, contig_len - end_window), contig_len, "right")

        bam.close()
        return best

    except Exception as e:
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

    # Linear plasmid gene keywords
    for kw in LINEAR_PLASMID_GENE_KEYWORDS:
        matches = [t for t in all_text if kw.lower() in t]
        if matches:
            hits["linear_plasmid_keywords"].append(kw)

    # Toxin-antitoxin systems
    for system, genes in TAS_GENES.items():
        found = [g for g in genes if any(g.lower() in t for t in all_text)]
        if len(found) >= 1:
            hits["tas_systems"].append(system)

    # Partition systems
    for par, pfams in PAR_PFAM.items():
        if any(par.lower() in t for t in all_text):
            hits["partition_genes"].append(par)

    # Resistance genes
    for gene, description in RESISTANCE_GENES.items():
        if any(gene.lower() in t for t in all_text):
            hits["resistance_genes"].append(f"{gene} ({description})")

    # IS elements
    for is_elem in LINEAR_PLASMID_IS:
        if any(is_elem.lower() in t for t in all_text):
            hits["is_elements"].append(is_elem)

    # Organism-specific markers
    def _word_match(keyword, texts):
        """True if keyword appears as a whole word in any text (case-insensitive)."""
        pat = r'(?<![a-z0-9])' + re.escape(keyword.lower()) + r'(?![a-z0-9])'
        return any(re.search(pat, t) for t in texts)

    # Enterococcal pELF markers (Hashimoto 2019 / Boumamoud 2022)
    enterococcal_m = ["is1216", "ftsK", "soj", "ftsk", "vana", "vanb"]
    hits["enterococcal_markers"] = [m for m in enterococcal_m
                                    if any(m.lower() in t for t in all_text)]

    # Terminal protein genes — word-boundary to avoid FP from e.g. "C-terminal"
    hits["terminal_protein_genes"] = [
        kw for kw in INVERTRON_TP_GENE_KEYWORDS
        if _word_match(kw, all_text)
    ]

    # repB / Rep_2 superfamily (Hashimoto 2019: pELF1 replication initiation)
    rep2_keywords = ["repb", "rep_2", "rep2", "replication initiation"]
    hits["repB_rep2"] = [kw for kw in rep2_keywords
                         if any(kw in t for t in all_text)]

    # ftsK + parA + repB co-occurrence signature (Hashimoto 2019 confirmed trio)
    trio = {"ftsk": False, "para": False, "repb": False}
    for t in all_text:
        if "ftsk" in t or "ftsK" in t: trio["ftsk"] = True
        if "para" in t or "parA" in t: trio["para"] = True
        if "repb" in t or "repB" in t: trio["repb"] = True
    hits["ftsK_parA_repB_trio"] = [k for k, v in trio.items() if v]

    # IS1216E flanking in same direction (vanM cluster structure, Hashimoto 2019)
    # Heuristic: both IS1216 variants present = flanking likely
    is_variants = [e for e in ["IS1216E", "IS1216V", "IS1216"]
                   if any(e.lower() in t for t in all_text)]
    hits["is1216_variants"] = is_variants

    return {k: list(set(v)) for k, v in hits.items()}


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 6: COPY NUMBER ESTIMATION (from BAM coverage)
# ─────────────────────────────────────────────────────────────────────────────

def estimate_copy_number(bam_file: str, contig_id: str,
                         chromosome_contigs: list = None) -> dict:
    """
    Estimate copy number of a contig relative to the chromosome by comparing
    mean read depths (~1 copy/cell characteristic of linear plasmids).

    Requires pysam. Falls back gracefully if not installed.
    """
    try:
        import pysam
    except ImportError:
        return {"available": False, "message": "pysam not installed"}

    if not bam_file or not os.path.exists(bam_file):
        return {"available": False, "message": "BAM file not found"}

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")

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

        bam.close()
        return {
            "available": True,
            "plasmid_depth": round(plasmid_depth, 1),
            "chromosome_depth": round(chrom_depth, 1) if chrom_depth else None,
            "estimated_copy_number": copy_number,
            "consistent_with_linear": (copy_number is not None and 0.5 <= copy_number <= 3.0),
        }
    except Exception as e:
        return {"available": False, "message": str(e)}


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 7: ASSEMBLY GRAPH TOPOLOGY (GFA parsing)
# ─────────────────────────────────────────────────────────────────────────────

def parse_gfa_topology(gfa_file: str, length_map: dict = None) -> dict:
    """
    Parse a GFA assembly graph (Unicycler/Plassembler output) to detect linear
    contig topology signatures (open ends, isolated contigs).

    Returns a dict: contig_id → topology_info.

    length_map: optional {contig_id: seq_length} built from FASTA records.
    When GFA segment names differ from FASTA contig IDs (e.g. Plassembler
    uses numeric names like "1", "2", ...) segments are re-keyed by matching
    their LN tag against the provided lengths. Ambiguous lengths are skipped.
    """
    if not gfa_file or not os.path.exists(gfa_file):
        return {}

    links = defaultdict(list)   # node → [(other_node, orientation)]
    segments = set()
    seg_lengths = {}             # seg_name → int length

    with open(gfa_file) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if parts[0] == "S":
                seg = parts[1]
                segments.add(seg)
                seq = parts[2] if len(parts) > 2 else "*"
                ln = len(seq) if seq != "*" else 0
                for tag in parts[3:]:
                    if tag.startswith("LN:i:"):
                        ln = int(tag[5:])
                        break
                if ln:
                    seg_lengths[seg] = ln
            elif parts[0] == "L":
                if len(parts) >= 5:
                    src, src_o, dst, dst_o = parts[1], parts[2], parts[3], parts[4]
                    links[src].append((dst, src_o, dst_o))
                    links[dst].append((src, dst_o, src_o))

    topology = {}
    for seg in segments:
        connected = links[seg]
        n_links = len(connected)
        if n_links == 0:
            topology[seg] = "isolated"
        elif n_links == 1:
            topology[seg] = "linear_end"
        elif n_links == 2:
            partners = {c[0] for c in connected}
            if len(partners) == 1:
                topology[seg] = "tir_like"
            else:
                topology[seg] = "linear_middle_or_circular"
        else:
            topology[seg] = "complex"

    if length_map:
        # Build inverse: length → contig_id, skipping any ambiguous lengths
        len_to_ctg: dict = {}
        for ctg_id, ctg_len in length_map.items():
            if ctg_len in len_to_ctg:
                len_to_ctg[ctg_len] = None  # collision — mark ambiguous
            else:
                len_to_ctg[ctg_len] = ctg_id

        remapped = {}
        for seg, topo in topology.items():
            if seg in length_map:
                remapped[seg] = topo          # name already matches a contig ID
            else:
                ln = seg_lengths.get(seg)
                ctg_id = len_to_ctg.get(ln) if ln else None
                remapped[ctg_id if ctg_id else seg] = topo
        topology = remapped

    return topology


def is_linear_in_gfa(contig_id: str, gfa_topology: dict) -> dict:
    """Interpret GFA topology for a specific contig."""
    topo = gfa_topology.get(contig_id, "unknown")
    is_linear = topo in ("linear_end", "isolated", "tir_like")
    return {
        "topology": topo,
        "consistent_with_linear": is_linear,
        "description": {
            "linear_end":     "Contig has open end(s) — strong linear indicator",
            "isolated":       "No graph links — may be linear or very distinct",
            "tir_like":       "Both ends link to same partner — TIR structure detected",
            "complex":        "Multiple connections — branching / repeat region",
            "linear_middle_or_circular": "Two distinct partners — could be circular",
            "unknown":        "Contig not in GFA",
        }.get(topo, "unknown"),
    }


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
    elif hdr.get("circular_flag_absent"):
        # No circular=true on a plasmid-sized contig — Unicycler/Plassembler would
        # have written it if the contig circularised; absence is a linear indicator
        score += SCORING_WEIGHTS["circular_flag_absent"]
        breakdown["circular_flag_absent"] = SCORING_WEIGHTS["circular_flag_absent"]

    if hdr.get("header_cn_consistent") and not hdr.get("high_copy_circular"):
        score += SCORING_WEIGHTS["header_copy_number"]
        breakdown["header_copy_number"] = SCORING_WEIGHTS["header_copy_number"]

    # 1. Hairpin / self-complement
    sc = evidence.get("self_complement", {})
    if sc.get("found"):
        score += SCORING_WEIGHTS["self_complement_end"]
        breakdown["self_complement_end"] = SCORING_WEIGHTS["self_complement_end"]
    if sc.get("tata_motif_left") or sc.get("tata_motif_right"):
        score += SCORING_WEIGHTS["hairpin_end"]
        breakdown["hairpin_end"] = SCORING_WEIGHTS["hairpin_end"]

    # 2b. IDR + TATA loop (Hashimoto 2019 pELF1 signature)
    idr = evidence.get("idr_tata", {})
    if idr.get("found"):
        score += SCORING_WEIGHTS["idr_tata_detected"]
        breakdown["idr_tata_detected"] = SCORING_WEIGHTS["idr_tata_detected"]
    # Terminal TATA: assembler stopped at palindrome centre (hairpin tip assembly artefact)
    if idr.get("tata_at_terminus"):
        score += SCORING_WEIGHTS["tata_at_terminus"]
        breakdown["tata_at_terminus"] = SCORING_WEIGHTS["tata_at_terminus"]

    # Full palindrome in a raw long read — strongest IDR evidence possible
    idr_reads = evidence.get("idr_reads", {})
    if idr_reads.get("available") and idr_reads.get("found"):
        score += SCORING_WEIGHTS["idr_in_reads"]
        breakdown["idr_in_reads"] = SCORING_WEIGHTS["idr_in_reads"]

    # 2c. Asymmetric ends (pELF1-type: one hairpin + one invertron)
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
    # Suppressed only when annotation could not be matched to this contig
    # (Prokka-TSV fallback: all genes smeared across all contigs). A circular
    # flag alone does not make gene evidence unreliable when the annotation is
    # a properly-matched GFF3 with per-contig features.
    genes = evidence.get("genes", {})
    if genes.get("_annotation_mismatch"):
        genes = {}   # suppress gene scoring — annotation not contig-specific

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

    # Confidence level
    if score >= CONFIDENCE_THRESHOLDS["HIGH"]:
        confidence = "HIGH"
    elif score >= CONFIDENCE_THRESHOLDS["MEDIUM"]:
        confidence = "MEDIUM"
    elif score >= CONFIDENCE_THRESHOLDS["LOW"]:
        confidence = "LOW"
    else:
        confidence = "NONE"

    max_score = sum(SCORING_WEIGHTS.values())
    return {
        "total_score": score,
        "max_possible": max_score,
        "percent": round(100 * score / max_score, 1),
        "confidence": confidence,
        "breakdown": breakdown,
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
# MODULE 10: LEFT-END VISUALIZATION  (Fig. 3 style, Hashimoto 2019)
# ─────────────────────────────────────────────────────────────────────────────

def visualize_left_end_structure(
    seq: str,
    contig_id: str,
    idr_result: dict,
    bam_file: str = None,
    output_prefix: str = "left_end",
    window: int = 15_000,
) -> str:
    """
    Figure-3-style visualization of the pELF1 left-end IDR/TATA hairpin
    (Hashimoto et al. 2019, Front. Microbiol., Fig. 3).

    Panel A  Coverage-depth pileup with:
             · Two converging grey arrows: IDR arm 1 (→) and arm 2 (←)
             · Bounding box at the IDR junction (as in Fig. 3A)
             · Vertical dashed line at the TATA centre
             · ~80 bp sequence snippet with TATA highlighted in red
    Panel B  Schematic stem–loop (hairpin):
             5′ ─── arm 1 ─── [TATA] ─── arm 2 ───  (as in Fig. 3B)

    Parameters
    ----------
    seq            : full contig sequence
    contig_id      : contig name (used for BAM fetch and output filename)
    idr_result     : return value of detect_idr_tata()
    bam_file       : sorted, indexed BAM — actual read depth if present,
                     otherwise a schematic depth profile is drawn
    output_prefix  : output filename prefix
    window         : bp to examine at the left end (default 15 000)

    Returns the saved PNG path, or '' on failure.
    Requires matplotlib (conda install -c conda-forge matplotlib).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.ticker import FuncFormatter
    except ImportError:
        print("[WARN] matplotlib not installed — skipping left-end visualization. "
              "Install with: conda install -c conda-forge matplotlib", file=sys.stderr)
        return ""

    # ── Left-end region ──────────────────────────────────────────────────────
    chk    = min(window, len(seq) // 2)
    region = seq[:chk].upper()
    n      = len(region)

    # ── IDR / TATA positions from detect_idr_tata() ──────────────────────────
    left_idr         = idr_result.get("details", {}).get("left", {})
    tata_pos         = left_idr.get("tata_pos", -1)
    arm_len          = left_idr.get("arm_length", 0)
    idr_identity     = left_idr.get("identity", 0.0)
    tata_at_terminus = idr_result.get("tata_at_terminus", False)

    if tata_pos < 0:
        tata_pos = region.find(HAIRPIN_CENTRAL_MOTIF)   # fallback

    # ── Coverage depth ───────────────────────────────────────────────────────
    depths    = None
    schematic = True

    if bam_file and os.path.exists(bam_file):
        try:
            import pysam
            bam = pysam.AlignmentFile(bam_file, "rb")
            arr = np.zeros(n, dtype=np.float32)
            for col in bam.pileup(contig_id, 0, n,
                                  min_mapping_quality=0, stepper="all"):
                pos = col.reference_pos
                if 0 <= pos < n:
                    arr[pos] = col.nsegments
            bam.close()
            depths    = arr
            schematic = False
        except Exception as exc:
            print(f"[WARN] BAM coverage extraction failed ({exc}); "
                  "using schematic depth.", file=sys.stderr)

    if depths is None:
        # Schematic: full depth outside IDR region, half inside IDR arms,
        # near-zero at the TATA hairpin tip (Hashimoto 2019 Fig. 3A pattern)
        depths = np.full(n, 50.0, dtype=np.float32)
        if arm_len > 0 and tata_pos > 0:
            idr_end = min(n, tata_pos + 4 + arm_len)
            depths[:idr_end] = 25.0          # half depth across IDR arms
            hw = 80                           # near-zero window at TATA tip
            depths[max(0, tata_pos - hw): min(n, tata_pos + 4 + hw)] = 3.0

    # Rolling-mean smoothing for display
    kw = max(1, n // 400)
    sm = np.convolve(depths, np.ones(kw) / kw, mode="same") if kw > 1 else depths.copy()
    max_d = float(sm.max()) or 1.0

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 9), facecolor="white")
    gs  = fig.add_gridspec(
        4, 1, height_ratios=[0.09, 0.46, 0.09, 0.36],
        hspace=0.06, top=0.92, bottom=0.06,
    )
    ax_arr  = fig.add_subplot(gs[0])
    ax_cov  = fig.add_subplot(gs[1], sharex=ax_arr)
    ax_seq  = fig.add_subplot(gs[2], sharex=ax_arr)
    ax_hair = fig.add_subplot(gs[3])

    x = np.arange(n)

    # ─────────────────────────────────────────────────────────────────────────
    # PANEL A LABEL
    # ─────────────────────────────────────────────────────────────────────────
    fig.text(0.005, 0.95, "A", fontsize=13, fontweight="bold", va="top")

    # ─────────────────────────────────────────────────────────────────────────
    # ARROW TRACK  (IDR arm 1 → …TATA… ← IDR arm 2)
    # ─────────────────────────────────────────────────────────────────────────
    ax_arr.set_xlim(0, n)
    ax_arr.set_ylim(0, 1)
    ax_arr.axis("off")
    ax_arr.set_facecolor("#EFEFEF")
    ax_arr.patch.set_visible(True)

    if arm_len > 0 and tata_pos > 0:
        arm2_end = min(n - 1, tata_pos + 4 + arm_len)

        # Arm 1: rightward → from position 0 to tata_pos
        ax_arr.annotate(
            "", xy=(tata_pos, 0.50), xytext=(0, 0.50),
            arrowprops=dict(arrowstyle="-|>", color="#444444",
                            lw=4.0, mutation_scale=22),
            annotation_clip=False, zorder=5,
        )
        # Arm 2: leftward ← from arm2_end back to tata_pos+4
        ax_arr.annotate(
            "", xy=(tata_pos + 4, 0.50), xytext=(arm2_end, 0.50),
            arrowprops=dict(arrowstyle="-|>", color="#444444",
                            lw=4.0, mutation_scale=22),
            annotation_clip=False, zorder=5,
        )

        # Arm labels
        ax_arr.text(tata_pos / 2, 0.82, "IDR arm 1",
                    ha="center", va="center", fontsize=8.5, color="#222222")
        ax_arr.text((tata_pos + 4 + arm2_end) / 2, 0.82, "IDR arm 2",
                    ha="center", va="center", fontsize=8.5, color="#222222")

        # Bounding box at the IDR junction (as in Hashimoto 2019 Fig. 3A)
        box_half = max(200, arm_len // 6)
        bx0 = max(0, tata_pos - box_half)
        bw  = box_half * 2 + 4
        ax_arr.add_patch(mpatches.FancyBboxPatch(
            (bx0, 0.05), bw, 0.90,
            boxstyle="square,pad=0", linewidth=2,
            edgecolor="black", facecolor="none", zorder=7,
        ))

    elif tata_at_terminus:
        ax_arr.text(
            n // 2, 0.50,
            "← Contig starts at TATA — assembler truncated at hairpin tip",
            ha="center", va="center", fontsize=8.5, color="#333333",
        )

    # ─────────────────────────────────────────────────────────────────────────
    # COVERAGE / PILEUP TRACK
    # ─────────────────────────────────────────────────────────────────────────
    norm = sm / max_d   # normalise to [0, 1]

    # Forward reads — coral/red (top half)
    ax_cov.fill_between(x, norm, color="#D4604A", alpha=0.82, label="Forward reads")
    # Reverse reads — steel blue (mirrored below zero)
    ax_cov.fill_between(x, -norm, color="#4A78D4", alpha=0.82, label="Reverse reads")
    ax_cov.axhline(0, color="black", lw=0.7)

    # IDR region light shading
    if arm_len > 0 and tata_pos > 0:
        ax_cov.axvspan(0, min(n, tata_pos + 4 + arm_len),
                       alpha=0.07, color="gray", zorder=0)

    # TATA vertical marker
    if tata_pos >= 0:
        ax_cov.axvline(tata_pos + 2, color="black", lw=1.5, ls="--",
                       alpha=0.85, zorder=8)
        ax_cov.text(tata_pos + 2, 0.90, " TATA",
                    transform=ax_cov.get_xaxis_transform(),
                    fontsize=8, va="bottom", color="black")

    depth_label = "Read depth  [schematic]" if schematic else f"Read depth  (max {max_d:.0f}×)"
    ax_cov.set_ylabel(depth_label, fontsize=9)
    ax_cov.set_ylim(-1.12, 1.12)
    ax_cov.set_yticks([-1, 0, 1])
    ax_cov.set_yticklabels(["rev", "0", "fwd"], fontsize=7.5)
    ax_cov.spines["top"].set_visible(False)
    ax_cov.spines["right"].set_visible(False)
    ax_cov.legend(loc="upper right", fontsize=8, framealpha=0.8)

    # ─────────────────────────────────────────────────────────────────────────
    # SEQUENCE SNIPPET TRACK
    # ─────────────────────────────────────────────────────────────────────────
    ax_seq.axis("off")

    if tata_pos >= 0:
        flank   = 36
        s_start = max(0, tata_pos - flank)
        s_end   = min(n, tata_pos + 4 + flank)

        before = region[s_start:tata_pos]
        after  = region[tata_pos + 4:s_end]

        prefix_str = f"{s_start + 1}-{before}"
        tata_str   = region[tata_pos:tata_pos + 4]    # "TATA"
        suffix_str = f"{after}-{s_end}"

        full_len   = len(prefix_str) + len(tata_str) + len(suffix_str)
        x0_tata    = len(prefix_str) / full_len
        x0_after   = (len(prefix_str) + len(tata_str)) / full_len

        # Three text objects share the same y row; monospace → equal char widths
        ax_seq.text(x0_tata, 0.55, prefix_str,
                    transform=ax_seq.transAxes,
                    ha="right", va="center",
                    fontsize=7.5, fontfamily="monospace", color="black")
        ax_seq.text(x0_tata, 0.55, tata_str,
                    transform=ax_seq.transAxes,
                    ha="left", va="center",
                    fontsize=7.5, fontfamily="monospace",
                    color="white", fontweight="bold",
                    bbox=dict(facecolor="crimson", edgecolor="none",
                              boxstyle="round,pad=0.18"))
        ax_seq.text(x0_after, 0.55, suffix_str,
                    transform=ax_seq.transAxes,
                    ha="left", va="center",
                    fontsize=7.5, fontfamily="monospace", color="black")
        ax_seq.text(x0_tata + len(tata_str) / (2 * full_len), 0.05,
                    "*1", transform=ax_seq.transAxes,
                    ha="center", va="bottom", fontsize=7, color="crimson")

    # ─────────────────────────────────────────────────────────────────────────
    # PANEL B LABEL
    # ─────────────────────────────────────────────────────────────────────────
    fig.text(0.005, 0.38, "B", fontsize=13, fontweight="bold", va="top")

    # ─────────────────────────────────────────────────────────────────────────
    # HAIRPIN SCHEMATIC  (Fig. 3B style)
    # ─────────────────────────────────────────────────────────────────────────
    ax_hair.set_xlim(-0.8, 13.0)
    ax_hair.set_ylim(-0.7, 2.6)
    ax_hair.axis("off")

    if arm_len > 0 and tata_pos > 0:
        arm_d  = 5.0    # display units per arm
        y_mid  = 1.4    # vertical centre
        lw_ds  = 2.8    # double-strand linewidth

        def _draw_dsdna(ax, x0, x1, y=y_mid, color="royalblue", lw=lw_ds):
            """Two parallel lines + evenly spaced tick marks = double-stranded DNA."""
            for dy in (-0.10, +0.10):
                ax.plot([x0, x1], [y + dy, y + dy],
                        color=color, lw=lw, solid_capstyle="round")
            for xt in np.arange(x0 + 0.45, x1, 0.60):
                ax.plot([xt, xt], [y - 0.10, y + 0.10],
                        color=color, lw=1.0, alpha=0.45)

        # Arm 1
        _draw_dsdna(ax_hair, 0.0, arm_d)

        # Hairpin loop (semicircle) — the TATA is at the tip
        loop_r  = 0.55
        loop_cx = arm_d + loop_r
        theta   = np.linspace(np.pi / 2, -np.pi / 2, 200)
        ax_hair.plot(loop_cx + loop_r * np.cos(theta),
                     y_mid + loop_r * np.sin(theta),
                     color="royalblue", lw=lw_ds)

        # TATA label inside loop (red box, matching Fig. 3B *1 annotation)
        ax_hair.text(loop_cx + loop_r * 0.55, y_mid, "*1\nTATA",
                     ha="center", va="center", fontsize=8,
                     fontfamily="monospace", color="crimson",
                     fontweight="bold", linespacing=1.3)

        # Arm 2
        arm2_x0 = arm_d + 2 * loop_r
        _draw_dsdna(ax_hair, arm2_x0, arm2_x0 + arm_d)

        # 5' label at the left end
        ax_hair.text(-0.20, y_mid, "5′",
                     ha="right", va="center",
                     fontsize=10, fontstyle="italic", color="black")

        # Arm-length brackets + annotations
        for xa, xb, lbl in [
            (0.0,    arm_d,
             f"IDR arm 1  (~{arm_len // 1000} kb, {idr_identity:.0%} identity)"),
            (arm2_x0, arm2_x0 + arm_d,
             "IDR arm 2  (rc of arm 1)"),
        ]:
            y_br = y_mid - 0.60
            ax_hair.annotate("", xy=(xb, y_br), xytext=(xa, y_br),
                             arrowprops=dict(arrowstyle="<->",
                                             color="#888888", lw=1.0))
            ax_hair.text((xa + xb) / 2, y_br - 0.16, lbl,
                         ha="center", va="top", fontsize=7.5, color="#555555")

        # Caption footnote
        ax_hair.text(
            (arm2_x0 + arm_d) / 2, -0.50,
            "*1  5′-TATA-3′ is the sequence of the hairpin loop "
            "at the left end of pELF1 (Hashimoto et al. 2019, Fig. 3B)",
            ha="center", va="top", fontsize=7.5, color="#666666",
            fontstyle="italic",
        )

    else:
        msg = ("TATA detected at terminus — full IDR arm not assembled."
               if tata_at_terminus else
               "No IDR/TATA structure detected in this window.")
        ax_hair.text(6.0, 1.4, msg, ha="center", va="center",
                     fontsize=10, color="#888888")

    # ─────────────────────────────────────────────────────────────────────────
    # X-AXIS FORMATTING (coverage panel only)
    # ─────────────────────────────────────────────────────────────────────────
    def _fmt_bp(v, _):
        return f"{v / 1000:.0f} kb" if v >= 1000 else f"{int(v)} bp"

    ax_cov.xaxis.set_major_formatter(FuncFormatter(_fmt_bp))
    plt.setp(ax_arr.get_xticklabels(), visible=False)
    plt.setp(ax_seq.get_xticklabels(), visible=False)

    # ─────────────────────────────────────────────────────────────────────────
    # TITLE AND FIGURE CAPTION
    # ─────────────────────────────────────────────────────────────────────────
    depth_src = "schematic depth" if schematic else "BAM read depth"
    arm_info  = (f"{arm_len:,} bp arm, TATA at pos {tata_pos + 1}"
                 if arm_len > 0 and tata_pos >= 0
                 else ("TATA at terminus" if tata_at_terminus else "no IDR detected"))
    fig.suptitle(
        f"Left-end IDR/TATA hairpin structure — {contig_id}  "
        f"({depth_src}; {arm_info})",
        fontsize=10.5, y=0.975,
    )
    fig.text(
        0.50, 0.005,
        "FIGURE 3 | Visualization of the left-end sequence of pELF1.  "
        "(A) Coverage of short reads mapped to the left-end sequence.  "
        "Dark arrows represent completely identical inverted direct repeats (IDR).  "
        "These IDRs confront each other at the region enclosed in a square.  "
        "(B) The structure of the left end of pELF1 was predicted using the hairpin model.",
        ha="center", va="bottom", fontsize=7.5, color="#666666", fontstyle="italic",
    )

    # ─────────────────────────────────────────────────────────────────────────
    # SAVE
    # ─────────────────────────────────────────────────────────────────────────
    safe_id  = re.sub(r"[^\w.\-]", "_", contig_id)
    out_path = f"{output_prefix}_left_end_{safe_id}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[OUTPUT] Left-end visualization → {out_path}")
    return out_path


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def analyse_contig(record, args, annot_df, gfa_topology, ref_gc=None,
                   chromosome_depth: float = None) -> dict:
    """Run all modules on a single contig/sequence record."""
    seq = str(record.seq)
    cid = record.id

    evidence = {}

    # FASTA header metadata (circular flag, copy number) — Module 0
    header_meta = parse_fasta_header(record.description)
    evidence["header"] = assess_header_metadata(
        header_meta, chromosome_depth, seq_len=len(seq))

    # Structural (Hashimoto 2019: hairpin, IDR/TATA, invertron end)
    evidence["self_complement"] = detect_self_complement_ends(seq)
    evidence["idr_tata"]        = detect_idr_tata(seq)
    evidence["size"]            = classify_size(len(seq))

    # GC
    contig_gc = gc_content(seq)
    evidence["gc"] = assess_gc(contig_gc, ref_gc)

    # Gene content
    evidence["genes"] = screen_genes(annot_df, cid) if not annot_df.empty else {}

    # Asymmetric end analysis (pELF1-type, Hashimoto 2019)
    evidence["asymmetric_ends"] = detect_asymmetric_ends(
        evidence["self_complement"], evidence["idr_tata"],
        evidence.get("genes", {}))

    # PlasmidFinder no-hit flag (set externally; default False unless --no-plasmid-finder-hit)
    evidence["plasmid_finder_no_hit"] = getattr(args, "plasmid_finder_no_hit", False)

    # BAM-based
    if args.bam:
        evidence["copy_number"]   = estimate_copy_number(
            args.bam, cid, args.chromosome_contigs)
        evidence["coverage_drop"] = detect_coverage_drop_ends(
            args.bam, cid, len(seq))
        evidence["idr_reads"]     = detect_idr_in_reads(
            args.bam, cid, len(seq))
    else:
        evidence["copy_number"]   = {"available": False}
        evidence["coverage_drop"] = {"available": False}
        evidence["idr_reads"]     = {"available": False}

    # GFA topology
    evidence["gfa_topology"] = is_linear_in_gfa(cid, gfa_topology)

    # BLAST
    if args.blast_db:
        tmp_fa  = f"/tmp/{cid}_query.fa"
        tmp_out = f"/tmp/{cid}_blast.tsv"
        with open(tmp_fa, "w") as fh:
            fh.write(f">{cid}\n{seq}\n")
        blast_df = run_blast(tmp_fa, args.blast_db, tmp_out,
                             mode="blastn", identity=args.blast_identity)
        evidence["blast"] = interpret_blast_hits(blast_df)
    else:
        evidence["blast"] = {"hits": 0, "linear_plasmid_hit": False}

    # SKANI
    if args.skani_db:
        tmp_fa   = f"/tmp/{cid}_query.fa"
        tmp_out  = f"/tmp/{cid}_skani.tsv"
        with open(tmp_fa, "w") as fh:
            fh.write(f">{cid}\n{seq}\n")
        skani_df = run_skani(tmp_fa, args.skani_db, tmp_out,
                             min_ani=args.skani_ani, min_af=args.skani_af)
        evidence["skani"] = interpret_blast_hits(skani_df)
    else:
        evidence["skani"] = {"hits": 0, "linear_plasmid_hit": False}

    # Composite score
    evidence["score"] = compute_score(evidence)

    return {"contig": cid, "length": len(seq), "gc": contig_gc, "evidence": evidence}


def main():
    parser = argparse.ArgumentParser(
        description="Identify linear plasmids in assembled sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-i", "--input",   required=True,
                        help="Input FASTA file (assembled sequences)")
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
    parser.add_argument("--gfa",           help="Assembly graph GFA file (Unicycler)")
    parser.add_argument("--annot",         help="GFF3/TSV annotation file (Prokka)")
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
    parser.add_argument("--min-score", type=int, default=15,
                        help="Minimum score to report (default: 15)")
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

    # ── Validate short-read arguments ─────────────────────────────────────────
    _sr_interleaved = getattr(args, "shortread_fastq", None)
    _sr_r1          = getattr(args, "shortread_r1", None)
    _sr_r2          = getattr(args, "shortread_r2", None)

    if _sr_interleaved and (_sr_r1 or _sr_r2):
        parser.error("--shortread-fastq (interleaved) and --shortread-r1/r2 "
                     "are mutually exclusive.")
    if bool(_sr_r1) != bool(_sr_r2):
        parser.error("--shortread-r1 and --shortread-r2 must be supplied together.")

    # ── Auto-map long reads if --longread-fastq supplied ─────────────────────
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

    # ── Auto-map short reads if --shortread-fastq / --shortread-r1/r2 given ──
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

    # Parse GFA — pass length map so numeric Plassembler segment names are
    # remapped to FASTA contig IDs via sequence-length matching.
    length_map = {r.id: len(r.seq) for r in records}
    gfa_topology = parse_gfa_topology(args.gfa, length_map) if args.gfa else {}
    if gfa_topology:
        print(f"[INFO] Loaded GFA topology: {len(gfa_topology)} segments")

    # Compute reference GC if not provided (median of all contigs)
    if args.ref_gc is None:
        all_gc = [gc_content(str(r.seq)) for r in records]
        ref_gc = float(np.median(all_gc)) if all_gc else None
        print(f"[INFO] Reference GC (median) = {ref_gc:.2f}%")
    else:
        ref_gc = args.ref_gc

    # Parse header metadata for all records; identify chromosome depth for
    # copy-number normalisation (longest contig marked circular=true, or
    # contig whose ID contains "chromosome")
    all_headers = {r.id: parse_fasta_header(r.description) for r in records}
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

    # Analyse each contig
    results = []
    for rec in records:
        print(f"  → Analysing {rec.id} ({len(rec.seq):,} bp) ...", end=" ", flush=True)
        res = analyse_contig(rec, args, annot_df, gfa_topology, ref_gc,
                             chromosome_depth=chromosome_depth)
        score = res["evidence"]["score"]["total_score"]
        conf  = res["evidence"]["score"]["confidence"]
        circ  = all_headers[rec.id]["circular"]
        circ_tag = " [circular=true]" if circ else (
                   " [circular=false]" if circ is False else "")
        print(f"score={score}  [{conf}]{circ_tag}")
        results.append(res)

    # Filter and sort
    results = sorted(results, key=lambda x: x["evidence"]["score"]["total_score"], reverse=True)
    results = [r for r in results if r["evidence"]["score"]["total_score"] >= args.min_score]

    # Build TSV output
    rows = []
    for r in results:
        ev = r["evidence"]
        sc = ev["score"]
        genes = ev.get("genes", {})
        hdr = ev.get("header", {})
        rows.append({
            "contig":            r["contig"],
            "length_bp":         r["length"],
            "gc_pct":            r["gc"],
            "score":             sc["total_score"],
            "max_score":         sc["max_possible"],
            "pct_max":           sc["percent"],
            "confidence":        sc["confidence"],
            # Header metadata columns
            "circular_flag":     hdr.get("circular_flag", ""),
            "header_copy_number": hdr.get("copy_number", ""),
            "header_depth":      hdr.get("depth", ""),
            "hairpin_identity":  ev.get("self_complement", {}).get("identity", ""),
            "tata_motif":        (ev.get("self_complement", {}).get("tata_motif_left") or
                                  ev.get("self_complement", {}).get("tata_motif_right")),
            "idr_tata_found":    ev.get("idr_tata", {}).get("found", ""),
            "tata_at_terminus":  ev.get("idr_tata", {}).get("tata_at_terminus", ""),
            "idr_arm_bp":        ev.get("idr_tata", {}).get("best_arm_bp", ""),
            "idr_in_reads":      ev.get("idr_reads", {}).get("found", ""),
            "idr_read_arm_bp":   ev.get("idr_reads", {}).get("best_arm_bp", ""),
            "idr_read_identity": ev.get("idr_reads", {}).get("best_identity", ""),
            "idr_read_name":     ev.get("idr_reads", {}).get("best_read", ""),
            "asymmetric_ends":   ev.get("asymmetric_ends", {}).get("asymmetric_pelf_type", ""),
            "invertron_tp_gene": ev.get("asymmetric_ends", {}).get("has_tp_gene", ""),
            "size_families":     "|".join(ev.get("size", {}).get("matching_families", [])),
            "gc_deviation":      ev.get("gc", {}).get("deviation", ""),
            "partition_genes":   "|".join(genes.get("partition_genes", [])),
            "tas_systems":       "|".join(genes.get("tas_systems", [])),
            "resistance_genes":  "|".join(genes.get("resistance_genes", [])),
            "is_elements":       "|".join(genes.get("is_elements", [])),
            "enterococcal_m":    "|".join(genes.get("enterococcal_markers", [])),
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

    tsv_path = f"{args.output}.tsv"
    df_out = pd.DataFrame(rows)
    df_out.to_csv(tsv_path, sep="\t", index=False)
    print(f"\n[OUTPUT] TSV report → {tsv_path}")
    print(f"[OUTPUT] {len(rows)} contigs reported (score ≥ {args.min_score})")

    # Summary to stdout
    for conf in ["HIGH", "MEDIUM", "LOW"]:
        subset = df_out[df_out["confidence"] == conf]
        if len(subset) > 0:
            print(f"\n  {conf} confidence ({len(subset)} contigs):")
            for _, row in subset.iterrows():
                print(f"    {row['contig']:30s}  {row['length_bp']:>10,} bp  "
                      f"score={row['score']:>3}  GC={row['gc_pct']}%  "
                      f"hairpin={row['tata_motif']}  "
                      f"idr_tata={row['idr_tata_found']}  "
                      f"tata_terminus={row['tata_at_terminus']}  "
                      f"asymmetric={row['asymmetric_ends']}")

    # JSON output
    if args.json:
        json_path = f"{args.output}.json"
        with open(json_path, "w") as fh:
            json.dump(results, fh, indent=2, default=str)
        print(f"[OUTPUT] JSON detail    → {json_path}")

    # Figure-3-style left-end visualization for IDR/TATA-positive contigs
    if args.visualize:
        records_seqs = {r.id: str(r.seq) for r in records}
        viz_count = 0
        for res in results:
            idr = res["evidence"].get("idr_tata", {})
            if idr.get("found") or idr.get("tata_at_terminus"):
                visualize_left_end_structure(
                    seq           = records_seqs[res["contig"]],
                    contig_id     = res["contig"],
                    idr_result    = idr,
                    bam_file      = args.bam,
                    output_prefix = args.output,
                )
                viz_count += 1
        if viz_count == 0:
            print("[INFO] --visualize: no contigs with IDR/TATA evidence to plot.")

    print("\n[DONE]")


if __name__ == "__main__":
    main()
