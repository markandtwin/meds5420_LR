# RNA-seq Analysis Pipeline

This repository contains a pipeline for processing long-read RNA-seq data, performing quality control (QC), alignment, and gene quantification. The pipeline is designed to run on a High-Performance Computing (HPC) cluster using Slurm job scheduler, leveraging tools like NanoPlot, minimap2, and Rsubread::featureCounts. The output is suitable for downstream differential gene expression analysis (e.g., with DESeq2).

## Overview

The pipeline processes FASTQ files from long-read RNA-seq (e.g., Oxford Nanopore) through the following steps:
1. **FASTQ Quality Control**: Generate QC reports for FASTQ files using NanoPlot.
2. **Alignment**: Align FASTQ reads to a reference genome using minimap2.
3. **BAM Quality Control**: Assess alignment quality with NanoPlot on BAM files.
4. **Gene Quantification**: Count reads per gene using Rsubread::featureCounts.

The pipeline was developed for samples like `WT_D0_1.chr21` and `WT_D7_2.chr21`, as used in a differential expression analysis comparing wild-type conditions at day 0 and day 7.

## Prerequisites

- **HPC Cluster**: Access to a Slurm-based HPC cluster (e.g., with modules for bioinformatics tools).
- **Software**:
  - NanoPlot (`NanoPlot`)
  - minimap2 (`minimap2`)
  - samtools (`samtools`)
  - R with Rsubread (`Rsubread`)
- **Input Files**:
  - FASTQ files (e.g., `/home/FCAM/meds5420/Zhang_LR/fastq/WT_D0_1.chr21.fastq`)
  - Reference genome (e.g., `hg38_chr21.fa`)
  - Annotation file (e.g., `hg38_chr21.gtf`)
- **Storage**: Ensure sufficient disk space for FASTQ, BAM, and output files (e.g., `/home/FCAM/meds5420/Zhang_LR/`).

## Directory Structure

```
├── scripts/
│   ├── 01_nanoplot_fastq_qc.sh
│   ├── 02_minimap2.sh
│   ├── 03_nanoplot_bam_qc.sh
│   ├── 04_featureCounts.sh
├── data/
│   ├── fastq/
│   ├── reference/
│   ├── output/
└── README.md
```

- `scripts/`: Contains Slurm job scripts for each pipeline step.
- `data/`: Stores input FASTQ files, reference genome, and output directories.
- `output/`: Contains NanoPlot reports, BAM files, and featureCounts results.

## Usage

Run the pipeline steps sequentially on the HPC cluster using Slurm. Ensure all dependencies are loaded and input files are in place.

### Step 1: FASTQ Quality Control
Run NanoPlot to generate QC reports for FASTQ files.

```bash
sbatch 01_nanoplot_fastq_qc.sh
```

**Output**: HTML and text reports in `/home/FCAM/meds5420/Zhang_LR/NanoPlot_fastq_QC/WT_D0_1.chr21/NanoPlot-report.html`.

### Step 2: Alignment with minimap2
Align FASTQ reads to the reference genome (e.g., hg38 chromosome 21) using minimap2.

```bash
sbatch 02_minimap2.sh
```

**Output**: BAM files in `/home/FCAM/meds5420/Zhang_LR/minimap2_bam/WT_D0_1.chr21.bam`.

### Step 3: BAM Quality Control
Run NanoPlot on BAM files to assess alignment quality.

```bash
sbatch 03_nanoplot_bam_qc.sh
```

**Output**: QC reports in `/home/FCAM/meds5420/Zhang_LR/NanoPlot_bam_QC/WT_D0_1.chr21/`.

### Step 4: Gene Quantification with featureCounts
Quantify gene expression by counting reads per gene using Rsubread::featureCounts.

```bash
sbatch 04_featureCounts.sh
```

**Output**: Count matrix in `/home/FCAM/meds5420/Zhang_LR/featureCounts/hg38_chr21_quant_name`.

## Notes

- **File Paths**: Update paths in the scripts to match your HPC directory structure (e.g., `/home/FCAM/meds5420/Zhang_LR/`).
- **BAM File Validation**: If `featureCounts` fails with errors like `sequence length in the BAM record is out of the expected region: 2427, 1497`, validate BAM files:
  ```bash
  samtools quickcheck /home/FCAM/meds5420/Zhang_LR/minimap2_bam/*.bam
  ```
  Check NanoPlot FASTQ reports for read quality issues (e.g., short reads) and realign if needed.

- **Viewing NanoPlot Reports**: HTML reports (e.g., `NanoPlot-report.html`) cannot be opened directly on a headless HPC. Transfer to a local machine:
  ```bash
  scp username@hpc.server:/home/FCAM/meds5420/Zhang_LR/NanoPlot_fastq_QC/WT_D0_1.chr21/NanoPlot-report.html .
  ```
  Or use a text-based browser:
  ```bash
  lynx /home/FCAM/meds5420/Zhang_LR/NanoPlot_fastq_QC/WT_D0_1.chr21/NanoPlot-report.html
  ```

- **Downstream Analysis**: Use the count matrix for differential expression analysis with DESeq2. See the accompanying R Markdown file (`DESeq2_Analysis.Rmd`) for an example.

- **Alternative Tools**: For more accurate long-read quantification, consider [IsoQuant](https://github.com/ablab/IsoQuant) or [Bambu](https://github.com/GoekeLab/bambu) instead of `featureCounts`.

## Troubleshooting

- **Job Failures**: Check Slurm logs (`slurm-*.out`) for errors. Ensure modules are loaded (e.g., `module load R minimap2 samtools`).
- **Missing Files**: Verify FASTQ, reference, and GTF files exist. Update script paths if necessary.
- **featureCounts Error**: If you encounter errors in `04_featureCounts.sh`, ensure the reference genome and GTF match the BAM files. Re-run `minimap2` with a consistent reference if needed.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue for suggestions, bug fixes, or additional features.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
