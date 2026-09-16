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
FASTQ="$(readlink -f "$2")"
OUT_DIR="$(readlink -f "$3")"

export PATH="/scratch/tmp/thomachr/software:$PATH"
export PATH="/scratch/tmp/thomachr/software/FusionSeeker:$PATH"
export PATH="/scratch/tmp/thomachr/software/bsalign:$PATH"

DORADO=/scratch/tmp/thomachr/software/dorado-1.1.1-linux-x64/bin/dorado
REFERENCE=/scratch/tmp/thomachr/references/hg38/hg38.fa

TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-36}"
SORT_THREADS=4
DORADO_THREADS=$(( TOTAL_THREADS - SORT_THREADS ))
FUSION_THREADS="$TOTAL_THREADS"

ALIGN_DIR="$OUT_DIR/alignment"
TMP_DIR="$OUT_DIR/tmp"
FUSION_OUT="$OUT_DIR/fusionseeker_out"
SORTED_BAM="$ALIGN_DIR/${SAMPLE}_sorted.bam"
DORADO_LOG="$ALIGN_DIR/${SAMPLE}.dorado.log"
READ_COUNT_FILE="$ALIGN_DIR/read_counts.tsv"

ml purge
ml palma/2022a GCC/11.3.0 SAMtools/1.16.1

[[ -r "$REFERENCE" ]]  || { echo "ERROR: reference not readable: $REFERENCE" >&2; exit 1; }
[[ -x "$DORADO" ]]     || { echo "ERROR: Dorado not executable: $DORADO" >&2; exit 1; }
[[ -r "$FASTQ" ]]      || { echo "ERROR: fastq not readable: $FASTQ" >&2; exit 1; }
command -v samtools     >/dev/null || { echo "ERROR: samtools not found" >&2; exit 1; }
command -v fusionseeker >/dev/null || { echo "ERROR: fusionseeker not found" >&2; exit 1; }
mkdir -p "$ALIGN_DIR" "$TMP_DIR"

# Any tool that honours TMPDIR (samtools, python, fusionseeker) scratches here.
export TMPDIR="$TMP_DIR"

echo "[$(date)] $SAMPLE: start"
echo "  fastq:   $FASTQ"
echo "  out_dir: $OUT_DIR"
echo "  tmp:     $TMP_DIR"
"$DORADO" --version
samtools --version | head -n 1

## 1. Alignment: dorado aligner | samtools sort (no unsorted BAM on disk) ###############
if [[ -s "$SORTED_BAM" && -s "${SORTED_BAM}.bai" ]] && samtools quickcheck -q "$SORTED_BAM"; then
    echo "[$(date)] $SAMPLE: sorted BAM exists, skipping alignment"
else
    echo "[$(date)] $SAMPLE: Dorado aligner ($DORADO_THREADS threads) | samtools sort ($SORT_THREADS threads)"
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
            -o "${SORTED_BAM}.tmp" \
            -; then
        echo "ERROR: Dorado/samtools sort failed for $SAMPLE" >&2
        tail -n 100 "$DORADO_LOG" >&2 || true
        rm -f "${SORTED_BAM}.tmp"
        exit 1
    fi
    mv "${SORTED_BAM}.tmp" "$SORTED_BAM"
    samtools index -@ "$SORT_THREADS" "$SORTED_BAM"
    samtools quickcheck -v "$SORTED_BAM"
fi

## 2. Read counts from the BAM ##########################################################
# -F 0x900: primary records incl. unmapped = input reads; -F 0x904: mapped reads.
total_reads=$(samtools view -c -@ "$SORT_THREADS" -F 0x900 "$SORTED_BAM")
mapped_reads=$(samtools view -c -@ "$SORT_THREADS" -F 0x904 "$SORTED_BAM")
printf 'sample\ttotal_reads\tmapped_reads\n%s\t%s\t%s\n' "$SAMPLE" "$total_reads" "$mapped_reads" > "$READ_COUNT_FILE"
echo "[$(date)] $SAMPLE: total reads = $total_reads, mapped = $mapped_reads"

## 3. FusionSeeker ######################################################################
if [[ -s "$FUSION_OUT/confident_genefusion.txt" ]]; then
    echo "[$(date)] $SAMPLE: FusionSeeker results exist, skipping"
else
    echo "[$(date)] $SAMPLE: FusionSeeker ($FUSION_THREADS threads)"
    rm -rf "$FUSION_OUT"
    fusionseeker \
        --bam "$SORTED_BAM" \
        --datatype nanopore \
        --ref "$REFERENCE" \
        --thread "$FUSION_THREADS" \
        -o "$FUSION_OUT"
fi

rm -rf "$TMP_DIR"
echo "[$(date)] $SAMPLE: finished -> $FUSION_OUT/confident_genefusion.txt"
