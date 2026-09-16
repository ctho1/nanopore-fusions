#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal,requeue
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --job-name=jaffal
#SBATCH --mail-type=ALL
#SBATCH --output=./log/%x_%j.out.txt
#SBATCH --error=./log/%x_%j.err.txt

# JAFFAL (long-read JAFFA pipeline via Apptainer) for ONE sample.
# Usage (normally via submit_fusions.sh):
#   sbatch run_jaffal_sample.sh <sample> <fastq_abs> <out_dir_abs>
#
# Everything this job writes stays inside <out_dir_abs>:
#   <out_dir_abs>/                 bpipe working dir (.bpipe/, intermediates, jaffa_results.csv)
#   <out_dir_abs>/tmp/             TMPDIR for apptainer/bpipe/minimap2 scratch files
# The FusionSeeker job for the same sample uses a different directory tree
# (./results_fusionseeker/<sample>/), so the two never touch the same files.
# The only shared input is the read-only concatenated fastq.

set -Eeuo pipefail
trap 'status=$?; echo "[$(date)] ERROR at line ${LINENO}: ${BASH_COMMAND} (exit ${status})" >&2; exit ${status}' ERR

(( $# == 3 )) || { echo "Usage: $0 <sample> <fastq_abs> <out_dir_abs>" >&2; exit 2; }
SAMPLE="$1"
[[ "$SAMPLE" =~ ^[[:alnum:]][[:alnum:]_.-]*$ ]] || { echo "ERROR: invalid sample name: $SAMPLE" >&2; exit 2; }
FASTQ="$(readlink -f -- "$2")"
mkdir -p -- "$3"
OUT_DIR="$(readlink -f -- "$3")"
[[ "$OUT_DIR" != "/" ]] || { echo "ERROR: refusing to use / as output directory" >&2; exit 2; }
TMP_DIR="$OUT_DIR/tmp"
DONE_MARKER="$OUT_DIR/.jaffal.complete"
INPUT_STATE="$OUT_DIR/.jaffal-input-state"

REF_DIR=/scratch/tmp/thomachr/software/JAFFA
JAFFA_SIF="$REF_DIR/jaffa_latest.sif"
ANNOTATION=gencode49          # files under REF_DIR are hg38_gencode49.* (lowercase)
THREADS="${SLURM_CPUS_PER_TASK:-16}"
[[ "$THREADS" =~ ^[0-9]+$ ]] && (( THREADS >= 1 )) || {
    echo "ERROR: SLURM_CPUS_PER_TASK must be a positive integer (got: $THREADS)" >&2
    exit 1
}

command -v module >/dev/null || { echo "ERROR: environment modules are not available" >&2; exit 1; }
module purge
module load palma/2024a Apptainer

command -v apptainer >/dev/null || { echo "ERROR: apptainer not found" >&2; exit 1; }
[[ -f "$JAFFA_SIF" ]] || { echo "ERROR: JAFFA image not found: $JAFFA_SIF" >&2; exit 1; }
[[ -d "$REF_DIR" ]]   || { echo "ERROR: JAFFA reference dir not found: $REF_DIR" >&2; exit 1; }
[[ -r "$FASTQ" ]]     || { echo "ERROR: fastq not readable: $FASTQ" >&2; exit 1; }
mkdir -p "$OUT_DIR" "$TMP_DIR"

FASTQ_STAT="$(stat -c '%s:%Y' "$FASTQ")"
IMAGE_STAT="$(stat -c '%s:%Y' "$JAFFA_SIF")"
INPUT_SIGNATURE="fastq=$FASTQ:$FASTQ_STAT|image=$JAFFA_SIF:$IMAGE_STAT|genome=hg38|annotation=$ANNOTATION"

# Scratch files (host side and inside the container) go to the sample's own tmp/.
export TMPDIR="$TMP_DIR"
export APPTAINER_TMPDIR="$TMP_DIR"
export APPTAINERENV_TMPDIR=/tmp

if [[ -s "$OUT_DIR/jaffa_results.csv" \
        && -f "$DONE_MARKER" \
        && -r "$INPUT_STATE" \
        && "$(<"$INPUT_STATE")" == "$INPUT_SIGNATURE" ]]; then
    echo "[$(date)] $SAMPLE: jaffa_results.csv already exists in $OUT_DIR; nothing to do."
    exit 0
fi

echo "[$(date)] $SAMPLE: JAFFAL start"
echo "  fastq:   $FASTQ"
echo "  out_dir: $OUT_DIR"
echo "  tmp:     $TMP_DIR"
echo "  threads: $THREADS"
apptainer --version

# Run inside OUT_DIR so bpipe's .bpipe/ dir, intermediates and results stay
# there. The sample tmp/ is bound to /tmp inside the container, so container
# scratch also lands in this sample's tree. The fastq path is absolute and
# symlink-resolved (avoids the bind-mount-on-symlink pitfall).
cd "$OUT_DIR"
apptainer run \
    -B "$REF_DIR:/ref:ro" \
    -B "$TMP_DIR:/tmp" \
    -B "$FASTQ:$FASTQ:ro" \
    "$JAFFA_SIF" \
    -p readLayout=single \
    -p genome=hg38 -p annotation="$ANNOTATION" \
    -n "$THREADS" \
    /JAFFA/JAFFAL.groovy \
    "$FASTQ"

[[ -s "$OUT_DIR/jaffa_results.csv" ]] || { echo "ERROR: JAFFAL finished but jaffa_results.csv is missing/empty" >&2; exit 1; }
printf '%s\n' "$INPUT_SIGNATURE" > "${INPUT_STATE}.tmp"
mv -- "${INPUT_STATE}.tmp" "$INPUT_STATE"
touch "$DONE_MARKER"
rm -rf -- "$TMP_DIR"
echo "[$(date)] $SAMPLE: JAFFAL finished -> $OUT_DIR/jaffa_results.csv"
