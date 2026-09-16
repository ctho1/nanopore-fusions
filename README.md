# nanopore-fusions

Gene fusion detection from Nanopore cDNA reads on the PALMA-II cluster (UKM Münster).
Two independent callers per sample, each running as its own SLURM job:

- **FusionSeeker**: Dorado aligner (minimap2 `-x splice`, hg38) piped into `samtools sort`, then FusionSeeker.
- **JAFFAL**: JAFFA long-read pipeline via Apptainer (own alignment inside the container, hg38 / GENCODE 49).

No read trimming or splitting: both tools start from the concatenated raw FASTQ per sample.

## Layout

```
fastq/<sample>/*.fastq.gz          input chunks, one subfolder per sample
fastq_concat/<sample>.fastq.gz     concatenated input (created by submit_fusions.sh)
results_jaffal/<sample>/           JAFFAL only: bpipe work dir, tmp/, jaffa_results.csv
results_fusionseeker/<sample>/     FusionSeeker only: alignment/, tmp/, fusionseeker_out/
log/                               SLURM stdout/stderr per sample and tool
```

The two tools never write to the same directory tree; the only shared file is the read-only concatenated FASTQ.

## Usage

```bash
cd $(readlink -f .)          # apptainer bind mounts fail under symlinked paths
./submit_fusions.sh          # all samples under fastq/
./submit_fusions.sh S1 S2    # selected samples
DRY_RUN=1 ./submit_fusions.sh   # concatenate only, print sbatch commands
```

`submit_fusions.sh` concatenates the chunks of each sample and submits
`run_jaffal_sample.sh` (16 CPUs, 32 GB) and `run_fusionseeker_sample.sh` (36 CPUs, 80 GB),
both with `--time=6:00:00` and `--partition=normal,requeue`.

Both job scripts are idempotent: finished stages are skipped, so a rerun of
`submit_fusions.sh` after a time-out only restarts what is missing.

## Prerequisites (paths configured at the top of each job script)

- Dorado 1.1.1 and hg38 reference FASTA
- FusionSeeker and bsalign in `PATH`
- JAFFA Apptainer image plus hg38/gencode49 reference (`prepare_jaffa_hg38_reference.sh`)
- Modules: `palma/2022a GCC/11.3.0 SAMtools/1.16.1` (FusionSeeker job), `palma/2024a Apptainer` (JAFFAL job)

## Outputs

- `results_fusionseeker/<sample>/fusionseeker_out/confident_genefusion.txt`
- `results_fusionseeker/<sample>/alignment/read_counts.tsv`
- `results_jaffal/<sample>/jaffa_results.csv`
