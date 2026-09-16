# 🧬 Nanopore Fusion Detection

Detect gene fusions from Oxford Nanopore cDNA reads on the PALMA-II cluster
(UKM Münster). Each sample is analysed by two independent callers in separate
SLURM jobs:

- 🔎 **FusionSeeker** — Dorado/minimap2 splice alignment against hg38, followed
  by FusionSeeker.
- 🧪 **JAFFAL** — the JAFFA long-read pipeline in an Apptainer container, using
  hg38 and GENCODE 49.

Both workflows start from the concatenated, untrimmed FASTQ. They use separate
output trees and only share the read-only input file.

## 🔄 Workflow

```text
fastq/<sample>/*.fastq.gz
            │
            ▼
fastq_concat/<sample>.fastq.gz
            │
       ┌────┴───────────┐
       ▼                ▼
  JAFFAL job      FusionSeeker job
       │                │
       ▼                ▼
jaffa_results.csv  confident_genefusion.txt
```

## 📁 Repository layout

```text
fastq/<sample>/*.fastq.gz          Input chunks, one directory per sample
fastq_concat/<sample>.fastq.gz     Concatenated input created by the driver
results_jaffal/<sample>/           JAFFAL work directory and results
results_fusionseeker/<sample>/     Alignment and FusionSeeker results
log/                               SLURM stdout/stderr per sample and tool
```

Generated FASTQs, logs, and result files are ignored by Git; the directory
structure is retained through `.gitkeep` files.

## ✅ Prerequisites

The paths near the top of the job scripts must match the cluster installation:

- Dorado 1.1.1 and an hg38 reference FASTA
- FusionSeeker and `bsalign` available through `PATH`
- A JAFFA Apptainer image and the prepared hg38/GENCODE 49 reference
- PALMA modules:
  - `palma/2022a GCC/11.3.0 SAMtools/1.16.1` for FusionSeeker
  - `palma/2024a Apptainer` for JAFFAL
- SLURM commands such as `sbatch` and `squeue`

## 🖥️ SLURM resources

| Caller | CPUs | Memory | Time limit | Partition |
| --- | ---: | ---: | ---: | --- |
| JAFFAL | 16 | 32 GB | 6 hours | `normal,requeue` |
| FusionSeeker | 36 | 80 GB | 6 hours | `normal,requeue` |

## 🚀 Quick start

1. Place each sample's gzipped FASTQ chunks in its own directory:

   ```text
   fastq/
   ├── S1/
   │   ├── chunk_01.fastq.gz
   │   └── chunk_02.fastq.gz
   └── S2/
       └── reads.fastq.gz
   ```

2. Submit all samples or select specific sample IDs:

   ```bash
   ./submit_fusions.sh
   ./submit_fusions.sh S1 S2
   ```

3. Inspect the queue and logs:

   ```bash
   squeue -u "$USER"
   tail -f log/S1.jaffal.*.out.txt
   ```

Sample IDs may contain letters, numbers, dots, underscores, and hyphens.

### 🧭 Preview without submitting

The dry run prepares concatenated input files and prints the `sbatch` commands,
but does not require SLURM and does not submit jobs:

```bash
DRY_RUN=1 ./submit_fusions.sh
```

## ♻️ Safe reruns

The workflow is restart-friendly:

- Concatenated FASTQs are rebuilt when the chunk list, file size, or modification
  time changes.
- BAM and caller results are reused only when their recorded FASTQ/reference
  signature still matches; changed inputs trigger the affected stages again.
- A missing BAM index is recreated without realigning when the BAM itself is
  still current.
- FusionSeeker publishes its output only after a successful run and a non-empty
  result file, then writes a completion marker.
- JAFFAL writes a completion marker only after `jaffa_results.csv` has been
  verified.

After a timeout or failure, run `./submit_fusions.sh` again for the affected
sample. Incomplete stages are resumed or rebuilt.

## 📊 Main outputs

| Caller | Result |
| --- | --- |
| FusionSeeker | `results_fusionseeker/<sample>/fusionseeker_out/confident_genefusion.txt` |
| FusionSeeker QC | `results_fusionseeker/<sample>/alignment/read_counts.tsv` |
| JAFFAL | `results_jaffal/<sample>/jaffa_results.csv` |

## ⚠️ Notes

- Run the scripts from a real path rather than a symlinked checkout; Apptainer
  bind mounts can otherwise fail. The driver automatically switches to the
  physical repository directory.
- Do not launch the two per-sample job scripts manually unless you provide the
  required absolute FASTQ and output paths. The driver handles this normally.
- Cluster software and reference paths are site-specific and should be reviewed
  before the first run.
