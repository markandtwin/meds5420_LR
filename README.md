# RNA-seq Analysis Pipeline

This repository contains a pipeline for processing long-read RNA-seq data, performing quality control (QC), alignment, and gene quantification. The pipeline is designed to run on a High-Performance Computing (HPC) cluster using Slurm job scheduler, leveraging tools like NanoPlot, minimap2, and Rsubread::featureCounts. The output is suitable for downstream differential gene expression analysis (e.g., with DESeq2).

## Overview

The pipeline processes FASTQ files from long-read RNA-seq (e.g., Oxford Nanopore) through the following steps:
1. **FASTQ Quality Control**: Generate QC reports for FASTQ files using NanoPlot.
2. **Alignment**: Align FASTQ reads to a reference genome using minimap2.
3. **BAM Quality Control**: Assess alignment quality with NanoPlot on BAM files.
4. **Gene Quantification**: Count reads per gene using featureCounts in subread.

The pipeline was developed for samples from two conditions like `WT_D0_*` and `WT_D7_*`, as used in a differential expression analysis comparing wild-type conditions at day 0 and day 7.

## Prerequisites

- **HPC Cluster**: Access to a Slurm-based HPC cluster (e.g., with modules for bioinformatics tools).
- **Software**:
  - NanoPlot (`NanoPlot`)
  - minimap2 (`minimap2`)
  - samtools (`samtools`)
  - featureCounts within subread (`subread`)
- **Input Files**:
  - FASTQ files (e.g., `/home/FCAM/meds5420/Zhang_LR/fastq/WT_D0_1.chr21.fastq`)
  - Reference genome (e.g., `hg38_chr21.fa`)
  - Annotation file (e.g., `hg38_chr21.gtf`)


## Directory Structure
Prepare a folder in your user diretory for today's long-read RNA-seq data processing, with the substucture like this:
```
├── scripts/
│   ├── 01_nanoplot_fastq_qc.sh
│   ├── 02_minimap2.sh
│   ├── 03_nanoplot_bam_qc.sh
│   ├── 04_featureCounts.sh
└── eofiles/
```
For the scripts for each step, you can either copy it from the instruction below or from `/home/FCAM/meds5420/Zhang_LR/scripts/`.

Data that you need access from `/home/FCAM/meds5420/Zhang_LR/` directory (you don't have to make your own copy):
```
├── genome_chr21/
│   ├── hg38_chr21.gtf
│   ├── hg38_chr21.bed
│   ├── hg38_chr21.fa
└── short_read_bam/
│   ├── iNeu-D7-1.chr21.bam.bai
│   ├── iNeu-D7-1.chr21.bam
│   ├── H9-D1-1.chr21.bam.bai
│   ├── H9-D1-1.chr21.bam
└── fastq/
    ├── WT_D0_1.chr21.fastq
    ├── WT_D0_2.chr21.fastq
    ├── WT_D0_3.chr21.fastq
    ├── WT_D7_1.chr21.fastq
    ├── WT_D7_2.chr21.fastq
    └── WT_D7_3.chr21.fastq
```


- `scripts/`: Contains Slurm job scripts for each pipeline step.
- `eofile/`: Contains job error and out files.
- `genome_chr21/`: Contains genome sequence, gene annotation in bed and bed files.
- `fastq/`: Stores input FASTQ files, reference genome, and output directories.
- `short_read_bam`: Short-read RNA-seq data from same samples as reference.

## Usage

Run the pipeline steps sequentially on the HPC cluster using Slurm. Ensure all dependencies are loaded and input files are in place.

### Step 1: FASTQ Quality Control
Before you start, make sure you navigate to the `scripts` directory.
Check the `01_nanoplot_fastq_qc.sh` file to make sure the files are available for all the commands. 
```bash
#!/bin/bash
#BATCH --job-name=NanoPlot
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 3   
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=END
#SBATCH --mem=2G
#SBATCH --mail-user=yourname@uchc.edu
#SBATCH --output=../eofiles/%x.%j.out  # Standard output log
#SBATCH --error=../eofiles/%x.%j.err   # Standard error log

# Print job information

set -e  # Exit immediately if any command fails (important for debugging)

# Record start time
start_time=$(date +%s)  # Get timestamp in seconds


date
echo "Hostname: $(hostname)"

# Load required modules
module load NanoPlot

# Directory containing the FASTQ files
FASTQ_DIR="/home/FCAM/meds5420/Zhang_LR/fastq"

# Directory to store the output of NanoPlot
OUTPUT_DIR="../NanoPlot_fastq_QC/"


# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all the .fastq files in the directory
for sample in "$FASTQ_DIR"/*.fastq; do

	# Extract the sample name from the file name
	SAMPLE_NAME=$(basename "$sample" ".fastq")
    
	# Run NanoPlot
	NanoPlot --fastq "$sample" --loglength -t 3 \
    			-o "$OUTPUT_DIR/$SAMPLE_NAME" --plots dot

done

echo "NanoPlot quantification completed for all samples."

# Record end time
end_time=$(date +%s)  # Get timestamp in seconds

# Calculate runtime duration
runtime=$((end_time - start_time))
echo "Pipeline completed successfully."
echo "Total runtime: $((runtime / 60)) minutes and $((runtime % 60)) seconds."

```


Run NanoPlot to generate QC reports for FASTQ files.

```bash
sbatch 01_nanoplot_fastq_qc.sh
```

**Output**: HTML and text reports in `/home/FCAM/meds5420/YourUsrName/NanoPlot_fastq_QC/WT_D0_1.chr21/NanoPlot-report.html`.

### Step 2: Alignment with minimap2
Check the `02_minimap2.sh` file to make sure the files are available for all the commands. For this part, we submit it `general` partition to get more resource.
```bash
#BATCH --job-name=minimap2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10   
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=8G
#SBATCH --mail-user=yourname@uchc.edu
#SBATCH --output=../eofiles/%x.%j.out  # Standard output log
#SBATCH --error=../eofiles/%x.%j.err   # Standard error log

# Print job information

set -e  # Exit immediately if any command fails (important for debugging)

# Record start time
start_time=$(date +%s)  # Get timestamp in seconds


date
echo "Hostname: $(hostname)"

# Load required modules
module load minimap2
module load samtools

# Directory containing the FASTQ files
FASTQ_DIR="/home/FCAM/meds5420/Zhang_LR/fastq"

# Directory to store the output of minimap2
OUTPUT_DIR="../minimap2_bam/"

# Directory to store the reference file
GENOME_DIR="/home/FCAM/meds5420/Zhang_LR/genome_chr21/"



# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all the .bam files in the directory
for sample in "$FASTQ_DIR"/*.fastq; do

	# Extract the sample name from the file name
    SAMPLE_NAME=$(basename "$sample" ".fastq")
    
    # Run minimap2
    minimap2 -ax splice -t 10 -B 3 -O 3,20 --junc-bed  "$GENOME_DIR"/*.bed \
		"$GENOME_DIR"/*.fasta  $sample  > "$OUTPUT_DIR"/$SAMPLE_NAME.sam
		
	# Run samtools
	samtools sort "$OUTPUT_DIR"/$SAMPLE_NAME.sam > "$OUTPUT_DIR"/$SAMPLE_NAME.bam
	samtools index "$OUTPUT_DIR"/$SAMPLE_NAME.bam
	rm "$OUTPUT_DIR"/$SAMPLE_NAME.sam
	
done

echo "minimap2 alignment completed for all samples."

# Record end time
end_time=$(date +%s)  # Get timestamp in seconds

# Calculate runtime duration
runtime=$((end_time - start_time))
echo "Pipeline completed successfully."
echo "Total runtime: $((runtime / 60)) minutes and $((runtime % 60)) seconds."

```


Align FASTQ reads to the reference genome (e.g., hg38 chromosome 21) using minimap2.

```bash
sbatch 02_minimap2.sh
```

**Output**: BAM files in `/home/FCAM/meds5420/YourUsrName/minimap2_bam/WT_D0_1.chr21.bam`.

### Step 3: BAM Quality Control
Check the `03_nanoplot_bam_qc.sh` file to make sure the files are available for all the commands. 
```bash
#!/bin/bash
#BATCH --job-name=NanoPlot
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 3   
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=END
#SBATCH --mem=2G
#SBATCH --mail-user=yourname@uchc.edu
#SBATCH --output=../eofiles/%x.%j.out  # Standard output log
#SBATCH --error=../eofiles/%x.%j.err   # Standard error log

# Print job information

set -e  # Exit immediately if any command fails (important for debugging)

# Record start time
start_time=$(date +%s)  # Get timestamp in seconds


date
echo "Hostname: $(hostname)"

# Load required modules
module load NanoPlot

# Directory containing the bam files
bam_DIR="../minimap2_bam"

# Directory to store the output of NanoPlot
OUTPUT_DIR="../NanoPlot_bam_QC/"


# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all the .bam files in the directory
for sample in "$bam_DIR"/*.bam; do

	# Extract the sample name from the file name
	SAMPLE_NAME=$(basename "$sample" ".bam")
	    
	# Run NanoPlot
	NanoPlot --bam "$sample" --loglength -t 3 -o "$OUTPUT_DIR/$SAMPLE_NAME" --plots dot

done

echo "NanoPlot QC completed for all samples."

# Record end time
end_time=$(date +%s)  # Get timestamp in seconds

# Calculate runtime duration
runtime=$((end_time - start_time))
echo "Pipeline completed successfully."
echo "Total runtime: $((runtime / 60)) minutes and $((runtime % 60)) seconds."

```

Run NanoPlot on BAM files to assess alignment quality.

```bash
sbatch 03_nanoplot_bam_qc.sh
```

**Output**: QC reports in `/home/FCAM/meds5420/YourUsrName/NanoPlot_bam_QC/WT_D0_1.chr21/`.

Mount or copy (using `scp`) the `/home/FCAM/meds5420/YourUsrName/` directory to your local machine. Load the short-read and long-read data side by side to appreciate the difference. Navigate to genes `PAXBP1` `ITSN1` `TRAPPC10` and see what you can find from the data.


### Step 4: Gene Quantification with featureCounts
Check the `04_featureCounts.sh` file to make sure the files are available for all the commands. 
```bash
#!/bin/bash
#BATCH --job-name=NanoPlot
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1   
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=END
#SBATCH --mem=2G
#SBATCH --mail-user=yourname@uchc.edu
#SBATCH --output=../eofiles/%x.%j.out  # Standard output log
#SBATCH --error=../eofiles/%x.%j.err   # Standard error log

# Print job information

set -e  # Exit immediately if any command fails (important for debugging)

# Record start time
start_time=$(date +%s)  # Get timestamp in seconds


date
echo "Hostname: $(hostname)"

# Load required modules
module load subread

# Directory containing the bam files
bam_DIR="../minimap2_bam"

# Directory to store the reference file
GENOME_DIR="/home/FCAM/meds5420/Zhang_LR/genome_chr21/"

# Directory to store the output of NanoPlot
OUTPUT_DIR="../featureCounts/"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"


# Run featureCounts from subread module, gene_id could be changed to gene_name
featureCounts 	-a $GENOME_DIR/hg38_chr21.gtf -L  -g gene_name -T 1 \
    			-o $OUTPUT_DIR/hg38_chr21_quant_name  $bam_DIR/*.bam

echo "featureCounts completed for all samples."

# Record end time
end_time=$(date +%s)  # Get timestamp in seconds

# Calculate runtime duration
runtime=$((end_time - start_time))
echo "Pipeline completed successfully."
echo "Total runtime: $((runtime / 60)) minutes and $((runtime % 60)) seconds."

```

Quantify gene expression by counting reads per gene using Rsubread::featureCounts.

```bash
sbatch 04_featureCounts.sh
```

**Output**: Count matrix at `/home/FCAM/meds5420/YourUsrName/featureCounts/hg38_chr21_quant_name`.

## Notes

- **File Paths**: Update paths in the scripts to match your HPC directory structure (e.g., `/home/FCAM/meds5420/YourUsrName/`). Check all the path with `YourUsrName` `yourname` and `YourAccount`.

- **Viewing NanoPlot Reports**: HTML reports (e.g., `NanoPlot-report.html`) cannot be opened directly on a headless HPC. Transfer to a local machine:
  ```bash
  scp -r YourAccount@transfer.cam.uchc.edu:/home/FCAM/meds5420/YourUsrName/ .
  ```

- **Downstream Analysis**: Use the count matrix for differential expression analysis with DESeq2. See the accompanying R Markdown file (`Deseq2_featureCounts_LR.Rmd`) for an example.

- **Alternative Tools**: For more accurate long-read quantification, consider [IsoQuant](https://github.com/ablab/IsoQuant) or [Bambu](https://github.com/GoekeLab/bambu) instead of `featureCounts`.

## Troubleshooting

- **Job Failures**: Check Slurm logs (`eofile`) for errors. Ensure modules are loaded (e.g., `module load  minimap2 `).
- **Missing Files**: Verify FASTQ, reference, and GTF files exist. Update script paths if necessary.
- **featureCounts Error**: If you encounter errors in `04_featureCounts.sh` but the job was completed, ignore the error as the data was extracted for chromosome 21 only. Otherwise, ensure the reference genome and GTF match the BAM files. Re-run `minimap2` with a consistent reference if needed.


## License

This project is licensed under the GLP-3 License. See the [LICENSE](LICENSE) file for details.
