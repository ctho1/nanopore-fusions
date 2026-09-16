#!/bin/bash
# Driver: concatenate all *.fastq.gz chunks per sample subfolder, then submit
# two independent SLURM jobs per sample (JAFFAL and FusionSeeker).
#
# Usage:  ./submit_fusions.sh            # all samples under ./fastq
#         ./submit_fusions.sh S1 S2      # only the named samples
#         DRY_RUN=1 ./submit_fusions.sh  # concatenate, print sbatch commands, submit nothing
#
# Layout:
#   ./fastq/<sample>/*.fastq.gz            input chunks
#   ./fastq_concat/<sample>.fastq.gz       concatenated, untrimmed, unsplit
#   ./results_jaffal/<sample>/             JAFFAL only: bpipe work dir, tmp/, results
#   ./results_fusionseeker/<sample>/       FusionSeeker only: alignment/, tmp/, fusionseeker_out/
#   ./log/<sample>.jaffal.<jobid>.{out,err}.txt, ./log/<sample>.fusionseeker.<jobid>.{out,err}.txt
#
# NOTE: run from a directory that is NOT under a symlink (apptainer bind mounts
# silently fail on symlinked paths). If unsure: cd $(readlink -f .)

set -Eeuo pipefail
shopt -s nullglob

FASTQ_ROOT=./fastq
CONCAT_DIR=./fastq_concat
JAFFAL_RESULT_DIR=./results_jaffal
FS_RESULT_DIR=./results_fusionseeker
LOG_DIR=./log

JAFFAL_JOB=./run_jaffal_sample.sh
FS_JOB=./run_fusionseeker_sample.sh
DRY_RUN="${DRY_RUN:-0}"

if [[ "$(pwd)" != "$(readlink -f "$(pwd)")" ]]; then
    echo "ERROR: current directory is under a symlink; cd \$(readlink -f .) first." >&2
    exit 1
fi
[[ -f "$JAFFAL_JOB" ]] || { echo "ERROR: $JAFFAL_JOB not found" >&2; exit 1; }
[[ -f "$FS_JOB" ]]     || { echo "ERROR: $FS_JOB not found" >&2; exit 1; }
command -v sbatch >/dev/null || { echo "ERROR: sbatch not found" >&2; exit 1; }

mkdir -p "$CONCAT_DIR" "$JAFFAL_RESULT_DIR" "$FS_RESULT_DIR" "$LOG_DIR"

# Sample selection: arguments, or all subfolders under FASTQ_ROOT.
if (( $# > 0 )); then
    samples=("$@")
else
    samples=()
    for d in "$FASTQ_ROOT"/*/; do samples+=("$(basename "$d")"); done
fi
(( ${#samples[@]} > 0 )) || { echo "ERROR: no sample subdirectories under $FASTQ_ROOT" >&2; exit 1; }

submit() {
    # submit <sample> <tool> <job script> <fastq_abs> <out_dir>
    local sample="$1" tool="$2" script="$3" fastq_abs="$4" out_dir="$5"
    local cmd=(sbatch
        --job-name="${tool}_${sample}"
        --output="$LOG_DIR/${sample}.${tool}.%j.out.txt"
        --error="$LOG_DIR/${sample}.${tool}.%j.err.txt"
        "$script" "$sample" "$fastq_abs" "$out_dir")
    if [[ "$DRY_RUN" == "1" ]]; then
        echo "  DRY_RUN: ${cmd[*]}"
    else
        local out
        out=$("${cmd[@]}")
        echo "  $tool: $out"
    fi
}

for sample in "${samples[@]}"; do
    sample_dir="$FASTQ_ROOT/$sample"
    [[ -d "$sample_dir" ]] || { echo "WARNING: $sample_dir not found; skipping" >&2; continue; }
    chunks=("$sample_dir"/*.fastq.gz)
    (( ${#chunks[@]} > 0 )) || { echo "WARNING: no *.fastq.gz in $sample_dir; skipping" >&2; continue; }

    echo "[$(date)] ===== $sample (${#chunks[@]} chunk(s)) ====="

    ## Concatenate ####################################################################
    # cat of gzip members is a valid multi-member gzip. A single chunk is used
    # directly (no copy).
    if (( ${#chunks[@]} == 1 )); then
        fastq_abs="$(readlink -f "${chunks[0]}")"
        echo "  single chunk, using directly: $fastq_abs"
    else
        concat="$CONCAT_DIR/${sample}.fastq.gz"
        if [[ -s "$concat" ]]; then
            echo "  concatenated fastq exists, reusing: $concat"
        else
            echo "  concatenating -> $concat"
            cat "${chunks[@]}" > "${concat}.tmp"
            gzip -t "${concat}.tmp"
            mv "${concat}.tmp" "$concat"
        fi
        fastq_abs="$(readlink -f "$concat")"
    fi

    ## Submit ########################################################################
    jaffal_out="$(readlink -f "$JAFFAL_RESULT_DIR")/$sample"
    fs_out="$(readlink -f "$FS_RESULT_DIR")/$sample"
    mkdir -p "$jaffal_out" "$fs_out"

    submit "$sample" jaffal       "$JAFFAL_JOB" "$fastq_abs" "$jaffal_out"
    submit "$sample" fusionseeker "$FS_JOB"     "$fastq_abs" "$fs_out"
done

echo "[$(date)] Done. Monitor with: squeue -u \$USER"
