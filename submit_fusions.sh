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
# The driver switches to the physical repository directory before creating paths,
# which avoids Apptainer bind-mount problems caused by symlinked working paths.

set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$SCRIPT_DIR"

FASTQ_ROOT=./fastq
CONCAT_DIR=./fastq_concat
JAFFAL_RESULT_DIR=./results_jaffal
FS_RESULT_DIR=./results_fusionseeker
LOG_DIR=./log

JAFFAL_JOB=./run_jaffal_sample.sh
FS_JOB=./run_fusionseeker_sample.sh
DRY_RUN="${DRY_RUN:-0}"

[[ "$DRY_RUN" == "0" || "$DRY_RUN" == "1" ]] || { echo "ERROR: DRY_RUN must be 0 or 1" >&2; exit 2; }
[[ -f "$JAFFAL_JOB" ]] || { echo "ERROR: $JAFFAL_JOB not found" >&2; exit 1; }
[[ -f "$FS_JOB" ]]     || { echo "ERROR: $FS_JOB not found" >&2; exit 1; }
if [[ "$DRY_RUN" == "0" ]]; then
    command -v sbatch >/dev/null || { echo "ERROR: sbatch not found" >&2; exit 1; }
fi
if stat -c '%s' "$JAFFAL_JOB" >/dev/null 2>&1; then
    STAT_FLAVOR=gnu
elif stat -f '%z' "$JAFFAL_JOB" >/dev/null 2>&1; then
    STAT_FLAVOR=bsd
else
    echo "ERROR: the installed stat command cannot report file size and modification time" >&2
    exit 1
fi

mkdir -p "$CONCAT_DIR" "$JAFFAL_RESULT_DIR" "$FS_RESULT_DIR" "$LOG_DIR"

# Sample selection: arguments, or all subfolders under FASTQ_ROOT.
if (( $# > 0 )); then
    samples=("$@")
else
    samples=()
    for d in "$FASTQ_ROOT"/*/; do samples+=("$(basename "$d")"); done
fi
(( ${#samples[@]} > 0 )) || { echo "ERROR: no sample subdirectories under $FASTQ_ROOT" >&2; exit 1; }

for sample in "${samples[@]}"; do
    [[ "$sample" =~ ^[[:alnum:]][[:alnum:]_.-]*$ ]] || {
        echo "ERROR: invalid sample name '$sample' (allowed: letters, numbers, dot, underscore, hyphen)" >&2
        exit 2
    }
done

concat_tmp=""
manifest_tmp=""
cleanup() {
    [[ -z "$concat_tmp" ]] || rm -f -- "$concat_tmp"
    [[ -z "$manifest_tmp" ]] || rm -f -- "$manifest_tmp"
}
trap cleanup EXIT

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
        manifest="$CONCAT_DIR/${sample}.chunks.tsv"
        concat_tmp="${concat}.tmp"
        manifest_tmp="${manifest}.tmp"

        : > "$manifest_tmp"
        for chunk in "${chunks[@]}"; do
            if [[ "$STAT_FLAVOR" == "gnu" ]]; then
                chunk_size="$(stat -c '%s' "$chunk")"
                chunk_mtime="$(stat -c '%Y' "$chunk")"
            else
                chunk_size="$(stat -f '%z' "$chunk")"
                chunk_mtime="$(stat -f '%m' "$chunk")"
            fi
            printf '%s\t%s\t%s\n' \
                "$(readlink -f "$chunk")" \
                "$chunk_size" \
                "$chunk_mtime" >> "$manifest_tmp"
        done

        if [[ -s "$concat" && -f "$manifest" ]] && cmp -s "$manifest_tmp" "$manifest"; then
            echo "  concatenated fastq is up to date, reusing: $concat"
            rm -f -- "$manifest_tmp"
            manifest_tmp=""
        else
            echo "  concatenating -> $concat"
            cat "${chunks[@]}" > "$concat_tmp"
            gzip -t "$concat_tmp"
            mv -- "$concat_tmp" "$concat"
            concat_tmp=""
            mv -- "$manifest_tmp" "$manifest"
            manifest_tmp=""
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
