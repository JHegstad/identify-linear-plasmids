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
    """
    if len(seq) < 2 * window:
        return {"found": False, "score": 0}

    left  = seq[:window].upper()
    right = seq[-window:].upper()
    rc_right = str(Seq(right).reverse_complement())

    matches = sum(1 for a, b in zip(left, rc_right) if a == b)
    identity = matches / window

    return {
        "found": identity >= 0.70,   # 70 % identity threshold for hairpin evidence
        "identity": round(identity, 3),
        "score": round(identity * 100, 1),
    }


def detect_asymmetric_ends(sc_result: dict,
                           gene_hits: dict,
                           cov_result: dict = None) -> dict:
    """
    pELF1 has ASYMMETRIC ends: left = hairpin, right = invertron (terminal protein).
    This is a unique indicator of the pELF lineage (Hashimoto 2019).

    Two detection paths — either is sufficient:

    Sequence + gene:
      - Self-complementary hairpin at one end (detect_self_complement_ends) AND
      - Terminal protein gene annotated

    BAM coverage (most diagnostic with Illumina reads):
      - Left end coverage drop (hairpin blocks adapter ligation) AND
      - Right end coverage NOT dropped (invertron end remains accessible)
      A symmetric drop at both ends does NOT count — that could be any palindromic
      structure, not the pELF1-specific asymmetry.
    """
    has_hairpin = sc_result.get("found")
    has_tp_gene = bool(gene_hits.get("terminal_protein_genes"))
    seq_asymmetric = has_hairpin and has_tp_gene

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
# MODULE 2c: COVERAGE DROP AT ENDS  (Hashimoto 2019)
# ─────────────────────────────────────────────────────────────────────────────

def detect_coverage_drop_ends(bam_file: str, contig_id: str,
                               contig_len: int, end_window: int = 5000) -> dict:
    """
    Detect coverage drop at contig ends, consistent with linear plasmid terminal structure.

    With Illumina reads: hairpin ends prevent adapter ligation and PCR amplification,
    causing a sharp drop to near-zero at the terminal region (Hashimoto 2019, Fig 3A).
    This is the most specific signal.

    With long reads (ONT/PacBio): reads can span hairpin loops, so a drop is not
    guaranteed by the same mechanism. A drop may still occur due to assembly truncation
    at the terminus or reduced mappability of the terminal palindromic sequence, but
    is less specific than with Illumina.

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

    # Structural (hairpin, invertron end)
    evidence["self_complement"] = detect_self_complement_ends(seq)
    evidence["size"]            = classify_size(len(seq))

    # GC
    contig_gc = gc_content(seq)
    evidence["gc"] = assess_gc(contig_gc, ref_gc)

    # Gene content
    evidence["genes"] = screen_genes(annot_df, cid) if not annot_df.empty else {}

    # PlasmidFinder no-hit flag (set externally; default False unless --no-plasmid-finder-hit)
    evidence["plasmid_finder_no_hit"] = getattr(args, "plasmid_finder_no_hit", False)

    # BAM-based
    if args.bam:
        evidence["copy_number"]   = estimate_copy_number(
            args.bam, cid, args.chromosome_contigs)
        evidence["coverage_drop"] = detect_coverage_drop_ends(
            args.bam, cid, len(seq))
    else:
        evidence["copy_number"]   = {"available": False}
        evidence["coverage_drop"] = {"available": False}

    # Asymmetric end analysis (pELF1-type, Hashimoto 2019)
    # Runs after BAM so coverage_drop is available for the left-only drop path
    evidence["asymmetric_ends"] = detect_asymmetric_ends(
        evidence["self_complement"],
        evidence.get("genes", {}),
        evidence.get("coverage_drop"))

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
                      f"hairpin_id={row['hairpin_identity']}  "
                      f"asymmetric={row['asymmetric_ends']}")

    # JSON output
    if args.json:
        json_path = f"{args.output}.json"
        with open(json_path, "w") as fh:
            json.dump(results, fh, indent=2, default=str)
        print(f"[OUTPUT] JSON detail    → {json_path}")

    if args.visualize:
        print("[INFO] --visualize: IDR/TATA visualization has been removed (TATA recognition disabled).")

    print("\n[DONE]")


if __name__ == "__main__":
    main()
