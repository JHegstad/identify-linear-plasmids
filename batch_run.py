#!/usr/bin/env python3
"""
Batch run of identify_linear_plasmids.py over the Andreas filesender collection.
Uses Illumina R1/R2 where available; falls back to ONT for 51565239.
"""
import os
import subprocess
import sys

BASE   = "/home/joachim/NGS/Projects/Kres-LRE/Linear_plasmids/filesender_Andreas"
FA_DIR = f"{BASE}/fasta/final_fasta"
R1_DIR = f"{BASE}/illumina/illumina_r1"
R2_DIR = f"{BASE}/illumina/illumina_r2"
ONT_DIR= f"{BASE}/Longread/ont_reads"
OUT_BASE = f"{BASE}/results"
SKANI_DB = "skani_db/sketches"
THREADS  = 8

# ── Sample manifest ──────────────────────────────────────────────────────────
# (sample_id, fasta, r1, r2, ont)
# r1/r2=None → use ont instead; ont always provided for reference
SAMPLES = [
    ("51525510",   "51525510_final.fasta",
     "51525510_lib649800_merged_1.fastq.gz",
     "51525510_lib649800_merged_2.fastq.gz",
     "51525510.fastq.gz"),

    ("51535247",   "51535247_final.fasta",
     "51535247_S2_L001_R1_001.fastq.gz",
     "51535247_S2_L001_R2_001.fastq.gz",
     "51535247.fastq.gz"),

    ("51558664",   "51558664_final.fasta",
     "51558664_S17_L001_R1_001.fastq.gz",
     "51558664_S17_L001_R2_001.fastq.gz",
     "51558664.fastq.gz"),

    ("51560653",   "51560653_final.fasta",
     "51560653_S4_L001_R1_001.fastq.gz",
     "51560653_S4_L001_R2_001.fastq.gz",
     "51560653.fastq.gz"),

    ("51561019",   "51561019_final.fasta",
     "51561019_S13_L001_R1_001.fastq.gz",
     "51561019_S13_L001_R2_001.fastq.gz",
     "51561019.fastq.gz"),

    # No Illumina — run with ONT
    ("51565239",   "51565239_final.fasta",
     None, None,
     "51565239.fastq.gz"),

    # Two Illumina pairs; prefer _S11_L001
    ("51575324",   "51575324_final.fasta",
     "51575324_S11_L001_R1_001.fastq.gz",
     "51575324_S11_L001_R2_001.fastq.gz",
     "51575324_filt.fastq.gz"),

    ("51575456",   "51575456_final.fasta",
     "51575456_S6_L001_R1_001.fastq.gz",
     "51575456_S6_L001_R2_001.fastq.gz",
     "51575456_filt.fastq.gz"),

    ("51579387",   "51579387_final.fasta",
     "51579387-924M00186_SS1012_S21_L001_R1_001.fastq.gz",
     "51579387-924M00186_SS1012_S21_L001_R2_001.fastq.gz",
     "51579387_filt.fastq.gz"),

    ("51581625",   "51581625_final.fasta",
     "51581625-924K00395-SS00007068_S16_L001_R1_001.fastq.gz",
     "51581625-924K00395-SS00007068_S16_L001_R2_001.fastq.gz",
     "51581625_filt.fastq.gz"),

    ("51586111",   "51586111_final.fasta",
     "51586111_S3_L001_R1_001.fastq.gz",
     "51586111_S3_L001_R2_001.fastq.gz",
     "51586111.fastq.gz"),

    ("51587501",   "51587501_final.fasta",
     "51587501_S5_L001_R1_001.fastq.gz",
     "51587501_S5_L001_R2_001.fastq.gz",
     "51587501.fastq.gz"),

    ("51600266",   "51600266_final.fasta",
     "51600266-925604525_S3_L001_R1_001.fastq.gz",
     "51600266-925604525_S3_L001_R2_001.fastq.gz",
     "51600266.fastq.gz"),

    ("KresLRE103", "KresLRE103_final.fasta",
     "KresLRE-103_lib649801_10136_3_1.fastq.gz",
     "KresLRE-103_lib649801_10136_3_2.fastq.gz",
     "ET_KresLRE-103_28_F_long.fastq.gz"),

    ("KresLRE105", "KresLRE105_final.fasta",
     "KresLRE-105_lib649796_10136_3_1.fastq.gz",
     "KresLRE-105_lib649796_10136_3_2.fastq.gz",
     "ET_KresLRE105_30_R1_long.fastq.gz"),

    ("KresLRE-119","KresLRE-119_final.fasta",
     "KresLRE-119_S5_L001_R1_001.fastq.gz",
     "KresLRE-119_S5_L001_R2_001.fastq.gz",
     "ET_KresLRE119_32_R1_long.fastq.gz"),

    ("KresLRE-157","KresLRE-157_final.fasta",
     "KresLRE-157_S22_L001_R1_001.fastq.gz",
     "KresLRE-157_S22_L001_R2_001.fastq.gz",
     "ET_KresLRE-157_38_F_long.fastq.gz"),

    ("KresLRE-28", "KresLRE-28_final.fasta",
     "KresLRE-28_S33_R1_001.fastq.gz",
     "KresLRE-28_S33_R2_001.fastq.gz",
     "ET_KresLRE28_6_R1_long.fastq.gz"),

    ("KresLRE-35", "KresLRE-35_final.fasta",
     "KresLRE-35_S1_R1_001.fastq.gz",
     "KresLRE-35_S1_R2_001.fastq.gz",
     "ET_KresLRE-35_8_F_long.fastq.gz"),

    ("KresVRE0042","KresVRE0042_final.fasta",
     "KresVRE0042_S1_R1_001.fastq.gz",
     "KresVRE0042_S1_R2_001.fastq.gz",
     "ET_KresVRE0042_14_F_long.fastq.gz"),
]


def run_sample(sample_id, fasta_name, r1_name, r2_name, ont_name):
    fasta = f"{FA_DIR}/{fasta_name}"
    out_dir = f"{OUT_BASE}/{sample_id}"
    os.makedirs(out_dir, exist_ok=True)
    out_prefix = f"{out_dir}/report"

    cmd = [
        "python", "identify_linear_plasmids.py",
        "-i", fasta,
        "--chromosome-contigs", "chromosome00001",
        "--skani-db", SKANI_DB,
        "--json",
        "-o", out_prefix,
    ]

    if r1_name and r2_name:
        r1 = f"{R1_DIR}/{r1_name}"
        r2 = f"{R2_DIR}/{r2_name}"
        cmd += ["--shortread-r1", r1, "--shortread-r2", r2,
                "--shortread-threads", str(THREADS)]
    else:
        ont = f"{ONT_DIR}/{ont_name}"
        cmd += ["--longread-fastq", ont,
                "--longread-preset", "map-ont",
                "--longread-threads", str(THREADS)]

    print(f"\n{'='*60}")
    print(f"  {sample_id}  ({'Illumina' if r1_name else 'ONT'})")
    print(f"{'='*60}")
    log_path = f"{out_dir}/run.log"
    with open(log_path, "w") as log:
        result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
        log.write(result.stdout)
        print(result.stdout)
    return result.returncode


if __name__ == "__main__":
    os.makedirs(OUT_BASE, exist_ok=True)
    failed = []
    for sid, fa, r1, r2, ont in SAMPLES:
        rc = run_sample(sid, fa, r1, r2, ont)
        if rc != 0:
            failed.append(sid)

    print(f"\n{'='*60}")
    print(f"Batch complete. {len(SAMPLES) - len(failed)}/{len(SAMPLES)} succeeded.")
    if failed:
        print(f"Failed: {failed}")
        sys.exit(1)
