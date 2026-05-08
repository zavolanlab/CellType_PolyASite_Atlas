#!/usr/bin/env python3
"""
Dynamic umi_tools dedup pipeline
=================================
- Uses umi_tools dedup --method=unique --per-gene
- Chromosome pre-filtered with samtools view into a temp indexed BAM
- Each job gets a private TMPDIR to avoid parallel umi_tools temp file conflicts
- Small chromosomes: one job directly
- Large chromosomes: split into CB+UB sub-chunks, dedup each, keep separate
- Output: per-chunk coordinate-sorted indexed BAMs ready for umi_tools count
- No final merge — count runs on each BAM separately

USAGE:
    python run_umi_dedup_pipeline.py

ONLY CHANGE THE PARAMETERS BELOW.
"""

import subprocess
import math
import os
import sys
import json

# ══════════════════════════════════════════════════════════════════════════════
# PARAMETERS — change these
# ══════════════════════════════════════════════════════════════════════════════
BAM               = "/scicore/home/zavolan/GROUP/CellType_PolyASite_Atlas/scqpas_internal/data/Liver_TLHJuly13S6_bam/10X_SingleCellLiverLandscape_TLHJuly13S6.bam"
OUTDIR            = "/scicore/home/zavolan/rados0001/Projects/CellType_PolyASite_Atlas/wf_runs/10x_liver_sample_test/sanity_test/umi_dedup_splitted_samtools_vFC"
CONDA_ENV_SAMTOOLS = "umi_tools_env"   # env with samtools
CONDA_ENV_UMITOOLS = "umi_tools_env"   # env with umi_tools (can be same)
MAX_READS_PER_JOB = 20_000_000        # reads threshold — above this, split into sub-chunks
MIN_READS         = 100_000           # skip chromosomes with fewer reads
THREADS           = 16
MEM_GB            = 6                 # RAM per dedup job
TIME_DEDUP        = "00:30:00"
TIME_SPLIT        = "01:30:00"
QOS               = "1day"
PARTITION         = "scicore"
# ══════════════════════════════════════════════════════════════════════════════

def run(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running: {cmd}\n{result.stderr}")
        sys.exit(1)
    return result.stdout.strip()

def submit(script_path, dependency=None):
    cmd = "sbatch --parsable"
    if dependency:
        cmd += f" --dependency=afterok:{dependency}"
    cmd += f" {script_path}"
    job_id = run(cmd)
    print(f"  Submitted {os.path.basename(script_path)} → job {job_id}")
    return job_id

def write_script(path, content):
    with open(path, "w") as f:
        f.write(content)
    os.chmod(path, 0o755)

# ── step 1: read idxstats ─────────────────────────────────────────────────────
print("Reading BAM index stats...")
idxstats_out = run(f"conda run -n {CONDA_ENV_SAMTOOLS} samtools idxstats {BAM}")

chromosomes = {}
for line in idxstats_out.splitlines():
    parts = line.split("\t")
    if len(parts) < 4:
        continue
    chrom, length, mapped, unmapped = parts[0], parts[1], int(parts[2]), int(parts[3])
    if mapped >= MIN_READS:
        chromosomes[chrom] = mapped

print(f"\nFound {len(chromosomes)} chromosomes with >= {MIN_READS:,} reads:\n")
print(f"  {'Chromosome':<12} {'Reads':>15}   {'Strategy'}")
print(f"  {'-'*12} {'-'*15}   {'-'*30}")
for chrom, reads in sorted(chromosomes.items(), key=lambda x: -x[1]):
    n_chunks = max(1, math.ceil(reads / MAX_READS_PER_JOB))
    strategy = "direct (1 job)" if n_chunks == 1 else f"split into {n_chunks} sub-chunks"
    print(f"  {chrom:<12} {reads:>15,}   {strategy}")

total_dedup_jobs = sum(max(1, math.ceil(r / MAX_READS_PER_JOB)) for r in chromosomes.values())
print(f"\nTotal dedup jobs: {total_dedup_jobs}")
print(f"Output directory: {OUTDIR}\n")

confirm = input("Proceed? (y/n): ").strip().lower()
if confirm != "y":
    print("Aborted.")
    sys.exit(0)

# ── create output directories ─────────────────────────────────────────────────
for d in [OUTDIR, f"{OUTDIR}/scripts", f"{OUTDIR}/logs",
          f"{OUTDIR}/chunks", f"{OUTDIR}/deduped", f"{OUTDIR}/tmp"]:
    os.makedirs(d, exist_ok=True)

with open(f"{OUTDIR}/chromosome_plan.json", "w") as f:
    plan = {c: {"reads": r, "n_chunks": max(1, math.ceil(r / MAX_READS_PER_JOB))}
            for c, r in chromosomes.items()}
    json.dump(plan, f, indent=2)

# ── bash variable expressions — expanded once at Python time ──────────────────
# These become literal bash strings in the generated scripts
SAMTOOLS_CMD = f"$(conda run -n {CONDA_ENV_SAMTOOLS} which samtools)"
UMITOOLS_CMD = f"$(conda run -n {CONDA_ENV_UMITOOLS} which umi_tools)"

# ── generate scripts per chromosome ──────────────────────────────────────────
all_dedup_job_ids = []

for chrom, reads in chromosomes.items():
    n_chunks = max(1, math.ceil(reads / MAX_READS_PER_JOB))

    if n_chunks == 1:
        script_path = f"{OUTDIR}/scripts/dedup_{chrom}.sh"
        script = f"""#!/bin/bash
#SBATCH --job-name=dedup_{chrom}
#SBATCH --cpus-per-task={THREADS}
#SBATCH --mem={MEM_GB}G
#SBATCH --time={TIME_DEDUP}
#SBATCH --output={OUTDIR}/logs/dedup_{chrom}_%j.out
#SBATCH --error={OUTDIR}/logs/dedup_{chrom}_%j.err
#SBATCH --qos={QOS}
#SBATCH --partition={PARTITION}

set -e  # exit on any error

source ~/.bashrc
conda activate {CONDA_ENV_UMITOOLS}
SAMTOOLS={SAMTOOLS_CMD}
UMITOOLS={UMITOOLS_CMD}

# private temp dir per job — prevents parallel umi_tools jobs from
# interfering with each other via shared /tmp files
export TMPDIR={OUTDIR}/tmp/$SLURM_JOB_ID
mkdir -p $TMPDIR

INPUT={OUTDIR}/deduped/{chrom}_input.bam

echo "[$(date)] Extracting {chrom} reads to temp BAM..."
$SAMTOOLS view -@ {THREADS} -b {BAM} {chrom} -o $INPUT
$SAMTOOLS index -@ {THREADS} $INPUT

echo "[$(date)] Deduplicating {chrom}..."
$UMITOOLS dedup \\
    --method=unique \\
    --extract-umi-method=tag \\
    --multimapping-detection-method=NH \\
    --cell-tag=CB \\
    --umi-tag=UB \\
    --per-gene \\
    --gene-tag=XT \\
    --assigned-status-tag=XS \\
    --per-cell \\
    -I $INPUT \\
    -S {OUTDIR}/deduped/{chrom}_dedup.bam \\
    -L {OUTDIR}/logs/dedup_{chrom}.log

echo "[$(date)] Sorting and indexing output..."
$SAMTOOLS sort -@ {THREADS} -m 1G \\
    {OUTDIR}/deduped/{chrom}_dedup.bam \\
    -o {OUTDIR}/deduped/{chrom}_dedup.sorted.bam
$SAMTOOLS index -@ {THREADS} {OUTDIR}/deduped/{chrom}_dedup.sorted.bam

echo "[$(date)] Cleaning up temp files..."
rm $INPUT ${{INPUT}}.bai {OUTDIR}/deduped/{chrom}_dedup.bam
rm -rf $TMPDIR

echo "Reads after dedup: $($SAMTOOLS view -c {OUTDIR}/deduped/{chrom}_dedup.sorted.bam)"
echo "[$(date)] Done: {chrom}"
"""
        write_script(script_path, script)
        job_id = submit(script_path)
        all_dedup_job_ids.append(job_id)

    else:
        # ── large chromosome: split by CB+UB first ────────────────────────────

        split_script_path = f"{OUTDIR}/scripts/split_{chrom}.sh"
        split_script = f"""#!/bin/bash
#SBATCH --job-name=split_{chrom}
#SBATCH --cpus-per-task={THREADS}
#SBATCH --mem=32G
#SBATCH --time={TIME_SPLIT}
#SBATCH --output={OUTDIR}/logs/split_{chrom}_%j.out
#SBATCH --error={OUTDIR}/logs/split_{chrom}_%j.err
#SBATCH --qos={QOS}
#SBATCH --partition={PARTITION}

set -e  # exit on any error

source ~/.bashrc
conda activate {CONDA_ENV_SAMTOOLS}
SAMTOOLS={SAMTOOLS_CMD}

export TMPDIR={OUTDIR}/tmp/$SLURM_JOB_ID
mkdir -p $TMPDIR

echo "[$(date)] Extracting and sorting {chrom} by CB+UB..."
$SAMTOOLS view -@ {THREADS} -b {BAM} {chrom} | \\
    $SAMTOOLS sort -@ {THREADS} -m 1G -t CB -t UB - \\
    -o {OUTDIR}/chunks/{chrom}_sorted_tags.bam

TOTAL=$($SAMTOOLS view -@ {THREADS} -c {OUTDIR}/chunks/{chrom}_sorted_tags.bam)
CHUNK_SIZE=$(( TOTAL / {n_chunks} ))
echo "[$(date)] {chrom}: $TOTAL reads → {n_chunks} chunks of ~$CHUNK_SIZE each"

$SAMTOOLS view -H {OUTDIR}/chunks/{chrom}_sorted_tags.bam > {OUTDIR}/chunks/{chrom}_header.sam

$SAMTOOLS view -@ {THREADS} {OUTDIR}/chunks/{chrom}_sorted_tags.bam | awk \\
    -v chunk_size=$CHUNK_SIZE \\
    -v outdir="{OUTDIR}/chunks" \\
    -v chrom="{chrom}" \\
    -v header="{OUTDIR}/chunks/{chrom}_header.sam" \\
'BEGIN {{
    chunk=0; count=0; prev_key=""
    outfile=outdir "/" chrom "_chunk_" sprintf("%04d", chunk) ".sam"
    while ((getline line < header) > 0) print line > outfile
    close(header)
}}
{{
    cb=""; ub=""
    for(i=12;i<=NF;i++) {{
        if($i ~ /^CB:Z:/) cb=$i
        if($i ~ /^UB:Z:/) ub=$i
    }}
    key=cb"\\t"ub
    if(key != prev_key && count >= chunk_size && prev_key != "") {{
        close(outfile); chunk++; count=0
        outfile=outdir "/" chrom "_chunk_" sprintf("%04d", chunk) ".sam"
        while ((getline line < header) > 0) print line > outfile
        close(header)
    }}
    print > outfile
    prev_key=key; count++
}}
END {{ close(outfile); print "Created " chunk+1 " chunks for {chrom}" }}'

# convert SAM chunks to coordinate-sorted BAM
for sam in {OUTDIR}/chunks/{chrom}_chunk_*.sam; do
    base="${{sam%.sam}}"
    $SAMTOOLS view -@ {THREADS} -b $sam | \\
        $SAMTOOLS sort -@ {THREADS} -m 1G - -o $base.bam
    $SAMTOOLS index -@ {THREADS} $base.bam
    rm $sam
done

rm {OUTDIR}/chunks/{chrom}_header.sam
rm {OUTDIR}/chunks/{chrom}_sorted_tags.bam
rm -rf $TMPDIR

echo "[$(date)] Done splitting {chrom}"
echo "Chunks created: $(ls {OUTDIR}/chunks/{chrom}_chunk_*.bam | wc -l)"
"""
        write_script(split_script_path, split_script)
        split_id = submit(split_script_path)

        # step B: dedup array on each sub-chunk
        dedup_array_path = f"{OUTDIR}/scripts/dedup_{chrom}_array.sh"
        dedup_array = f"""#!/bin/bash
#SBATCH --job-name=dedup_{chrom}
#SBATCH --array=0-{n_chunks - 1}
#SBATCH --cpus-per-task={THREADS}
#SBATCH --mem={MEM_GB}G
#SBATCH --time={TIME_DEDUP}
#SBATCH --output={OUTDIR}/logs/dedup_{chrom}_%A_%a.out
#SBATCH --error={OUTDIR}/logs/dedup_{chrom}_%A_%a.err
#SBATCH --qos={QOS}
#SBATCH --partition={PARTITION}

set -e  # exit on any error

source ~/.bashrc
conda activate {CONDA_ENV_UMITOOLS}
SAMTOOLS={SAMTOOLS_CMD}
UMITOOLS={UMITOOLS_CMD}

export TMPDIR={OUTDIR}/tmp/$SLURM_JOB_ID
mkdir -p $TMPDIR

CHUNK=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
INPUT={OUTDIR}/chunks/{chrom}_chunk_${{CHUNK}}.bam
OUTPUT={OUTDIR}/deduped/{chrom}_chunk_${{CHUNK}}_dedup.bam
FINAL={OUTDIR}/deduped/{chrom}_chunk_${{CHUNK}}_dedup.sorted.bam

echo "[$(date)] Deduplicating {chrom} chunk $CHUNK"
echo "Input: $($SAMTOOLS view -c $INPUT) reads"

$UMITOOLS dedup \\
    --method=unique \\
    --extract-umi-method=tag \\
    --multimapping-detection-method=NH \\
    --cell-tag=CB \\
    --umi-tag=UB \\
    --per-gene \\
    --gene-tag=XT \\
    --assigned-status-tag=XS \\
    --per-cell \\
    -I $INPUT \\
    -S $OUTPUT \\
    -L {OUTDIR}/logs/dedup_{chrom}_chunk_$CHUNK.log

$SAMTOOLS sort -@ {THREADS} -m 1G $OUTPUT -o $FINAL
$SAMTOOLS index -@ {THREADS} $FINAL
rm $OUTPUT
rm -rf $TMPDIR

echo "Output: $($SAMTOOLS view -c $FINAL) reads"
echo "[$(date)] Done: {chrom} chunk $CHUNK"
"""
        write_script(dedup_array_path, dedup_array)
        dedup_id = submit(dedup_array_path, dependency=split_id)
        all_dedup_job_ids.append(dedup_id)

# ── print summary ─────────────────────────────────────────────────────────────
all_deps = ":".join(all_dedup_job_ids)
print(f"\n{'='*60}")
print(f"All jobs submitted.")
print(f"Monitor with: squeue -u $USER")
print(f"All dedup job IDs: {all_deps}")
print(f"\nOnce all complete, run umi_tools count on:")
print(f"  {OUTDIR}/deduped/*_dedup.sorted.bam")
print(f"\nCheck with:")
print(f"  ls {OUTDIR}/deduped/*.sorted.bam | wc -l")
print(f"  jobstats")