#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=normal,requeue
#SBATCH --time=6:00:00
#SBATCH --mem=80G
#SBATCH --job-name=fusionseeker
#SBATCH --mail-type=ALL
#SBATCH --output=./log/%x_%j.out.txt
#SBATCH --error=./log/%x_%j.err.txt

# Dorado aligner (minimap2 -x splice) + FusionSeeker for ONE sample.
# Usage (normally via submit_fusions.sh):
#   sbatch run_fusionseeker_sample.sh <sample> <fastq_abs> <out_dir_abs>
#
# Everything this job writes stays inside <out_dir_abs>:
#   alignment/<sample>_sorted.bam(.bai)   Dorado/minimap2 alignment
#   alignment/<sample>.dorado.log         Dorado stderr
#   alignment/read_counts.tsv             total / mapped reads (from the BAM)
#   tmp/                                  samtools sort temp, TMPDIR (removed on success)
#   fusionseeker_out/                     FusionSeeker output (confident_genefusion.txt)
# The JAFFAL job for the same sample uses a different directory tree
# (./results_jaffal/<sample>/), so the two never touch the same files.
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

export PATH="/scratch/tmp/thomachr/software:$PATH"
export PATH="/scratch/tmp/thomachr/software/FusionSeeker:$PATH"
export PATH="/scratch/tmp/thomachr/software/bsalign:$PATH"

DORADO=/scratch/tmp/thomachr/software/dorado-1.1.1-linux-x64/bin/dorado
REFERENCE=/scratch/tmp/thomachr/references/hg38/hg38.fa

TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-36}"
[[ "$TOTAL_THREADS" =~ ^[0-9]+$ ]] && (( TOTAL_THREADS >= 2 )) || {
    echo "ERROR: at least 2 CPUs are required (got: $TOTAL_THREADS)" >&2
    exit 1
}
if (( TOTAL_THREADS >= 5 )); then
    SORT_THREADS=4
else
    SORT_THREADS=1
fi
DORADO_THREADS=$(( TOTAL_THREADS - SORT_THREADS ))
FUSION_THREADS="$TOTAL_THREADS"

ALIGN_DIR="$OUT_DIR/alignment"
TMP_DIR="$OUT_DIR/tmp"
FUSION_OUT="$OUT_DIR/fusionseeker_out"
SORTED_BAM="$ALIGN_DIR/${SAMPLE}_sorted.bam"
BAM_STAGE="${SORTED_BAM}.${SLURM_JOB_ID:-$$}.tmp"
DORADO_LOG="$ALIGN_DIR/${SAMPLE}.dorado.log"
READ_COUNT_FILE="$ALIGN_DIR/read_counts.tsv"
FUSION_STAGE="$OUT_DIR/.fusionseeker_out.${SLURM_JOB_ID:-$$}.tmp"
ALIGN_STATE="$ALIGN_DIR/.input-state"
FUSION_STATE="$OUT_DIR/.fusionseeker-input-state"
FUSION_DONE="$OUT_DIR/.fusionseeker-complete"

cleanup() {
    rm -f -- "$BAM_STAGE" "${BAM_STAGE}.bai" "${READ_COUNT_FILE}.tmp" \
        "${ALIGN_STATE}.tmp" "${FUSION_STATE}.tmp"
    [[ ! -d "$FUSION_STAGE" ]] || rm -rf -- "$FUSION_STAGE"
}
trap cleanup EXIT

command -v module >/dev/null || { echo "ERROR: environment modules are not available" >&2; exit 1; }
module purge
module load palma/2022a GCC/11.3.0 SAMtools/1.16.1

[[ -r "$REFERENCE" ]]  || { echo "ERROR: reference not readable: $REFERENCE" >&2; exit 1; }
[[ -x "$DORADO" ]]     || { echo "ERROR: Dorado not executable: $DORADO" >&2; exit 1; }
[[ -r "$FASTQ" ]]      || { echo "ERROR: fastq not readable: $FASTQ" >&2; exit 1; }
command -v samtools     >/dev/null || { echo "ERROR: samtools not found" >&2; exit 1; }
command -v fusionseeker >/dev/null || { echo "ERROR: fusionseeker not found" >&2; exit 1; }
mkdir -p "$ALIGN_DIR" "$TMP_DIR"

FASTQ_STAT="$(stat -c '%s:%Y' "$FASTQ")"
REFERENCE_STAT="$(stat -c '%s:%Y' "$REFERENCE")"
ALIGN_SIGNATURE="fastq=$FASTQ:$FASTQ_STAT|reference=$REFERENCE:$REFERENCE_STAT|preset=splice"

# Any tool that honours TMPDIR (samtools, python, fusionseeker) scratches here.
export TMPDIR="$TMP_DIR"

echo "[$(date)] $SAMPLE: start"
echo "  fastq:   $FASTQ"
echo "  out_dir: $OUT_DIR"
echo "  tmp:     $TMP_DIR"
"$DORADO" --version
samtools --version | head -n 1

## 1. Alignment: dorado aligner | samtools sort (no unsorted BAM on disk) ###############
if [[ -s "$SORTED_BAM" ]] \
        && samtools quickcheck -q "$SORTED_BAM" \
        && [[ -r "$ALIGN_STATE" ]] \
        && [[ "$(<"$ALIGN_STATE")" == "$ALIGN_SIGNATURE" ]]; then
    if [[ ! -s "${SORTED_BAM}.bai" ]]; then
        echo "[$(date)] $SAMPLE: BAM index is missing; creating it"
        samtools index -@ "$SORT_THREADS" "$SORTED_BAM"
    fi
    echo "[$(date)] $SAMPLE: valid sorted BAM exists, skipping alignment"
else
    echo "[$(date)] $SAMPLE: Dorado aligner ($DORADO_THREADS threads) | samtools sort ($SORT_THREADS threads)"
    rm -f -- "$BAM_STAGE"
    if ! "$DORADO" aligner \
            --mm2-opts "-x splice" \
            --threads "$DORADO_THREADS" \
            "$REFERENCE" \
            "$FASTQ" \
            2> "$DORADO_LOG" \
        | samtools sort \
            -@ "$SORT_THREADS" \
            -m 2G \
            -T "$TMP_DIR/${SAMPLE}.sort" \
            -o "$BAM_STAGE" \
            -; then
        echo "ERROR: Dorado/samtools sort failed for $SAMPLE" >&2
        tail -n 100 "$DORADO_LOG" >&2 || true
        exit 1
    fi
    samtools quickcheck -v "$BAM_STAGE"
    samtools index -@ "$SORT_THREADS" "$BAM_STAGE" "${BAM_STAGE}.bai"
    mv -- "$BAM_STAGE" "$SORTED_BAM"
    mv -- "${BAM_STAGE}.bai" "${SORTED_BAM}.bai"
    printf '%s\n' "$ALIGN_SIGNATURE" > "${ALIGN_STATE}.tmp"
    mv -- "${ALIGN_STATE}.tmp" "$ALIGN_STATE"
fi

## 2. Read counts from the BAM ##########################################################
# -F 0x900: primary records incl. unmapped = input reads; -F 0x904: mapped reads.
total_reads=$(samtools view -c -@ "$SORT_THREADS" -F 0x900 "$SORTED_BAM")
mapped_reads=$(samtools view -c -@ "$SORT_THREADS" -F 0x904 "$SORTED_BAM")
printf 'sample\ttotal_reads\tmapped_reads\n%s\t%s\t%s\n' "$SAMPLE" "$total_reads" "$mapped_reads" > "${READ_COUNT_FILE}.tmp"
mv -- "${READ_COUNT_FILE}.tmp" "$READ_COUNT_FILE"
echo "[$(date)] $SAMPLE: total reads = $total_reads, mapped = $mapped_reads"

## 3. FusionSeeker ######################################################################
BAM_STAT="$(stat -c '%s:%Y' "$SORTED_BAM")"
FUSION_SIGNATURE="$ALIGN_SIGNATURE|bam=$SORTED_BAM:$BAM_STAT|datatype=nanopore"
if [[ -s "$FUSION_OUT/confident_genefusion.txt" \
        && -f "$FUSION_DONE" \
        && -r "$FUSION_STATE" \
        && "$(<"$FUSION_STATE")" == "$FUSION_SIGNATURE" ]]; then
    echo "[$(date)] $SAMPLE: FusionSeeker results exist, skipping"
else
    echo "[$(date)] $SAMPLE: FusionSeeker ($FUSION_THREADS threads)"
    rm -rf -- "$FUSION_STAGE"
    fusionseeker \
        --bam "$SORTED_BAM" \
        --datatype nanopore \
        --ref "$REFERENCE" \
        --thread "$FUSION_THREADS" \
        -o "$FUSION_STAGE"
    [[ -s "$FUSION_STAGE/confident_genefusion.txt" ]] || {
        echo "ERROR: FusionSeeker finished but confident_genefusion.txt is missing/empty" >&2
        exit 1
    }
    rm -rf -- "$FUSION_OUT"
    mv -- "$FUSION_STAGE" "$FUSION_OUT"
    printf '%s\n' "$FUSION_SIGNATURE" > "${FUSION_STATE}.tmp"
    mv -- "${FUSION_STATE}.tmp" "$FUSION_STATE"
    touch "$FUSION_DONE"
fi

rm -rf -- "$TMP_DIR"
echo "[$(date)] $SAMPLE: finished -> $FUSION_OUT/confident_genefusion.txt"
