---
title: \sf MEDS5420 - Long Read (LR) RNA-seq data format, QC, and processing.
header-includes:
- \usepackage{color}
- \usepackage{float}
- \DeclareUnicodeCharacter{2212}{-}
date: "Apr 23, 2025" 
output:
  bookdown::html_document2:
    toc: true
fontsize: 14pt
geometry: margin=1in
---


# Background
* LR vs short read 
* LR DNA sequencing
* LR RNA sequencing

```{r  out.width = "100%", echo=F, fig.align = "center", fig.cap="pipeline flow chart"}
#library(knitr)
knitr::include_graphics("./XXXXX.png") 
```


# Practice: QC and processing of LR RNA-seq data:
* Check the .fastq files
* Using `NanoPlot` to check the quality of LR raw data.
* Using `minimap2` to align the data to genome sequence.
* Using `NanoPlot` to check the quality of aligned LR data.
* Using `featureCounts` to quantify the reads aligned to each gene.



# Check the data

## navigate to the directory with the fastq data
```{r engine='bash', eval=F, echo=TRUE}
cd /home/FCAM/meds5420/Your_folder/
ls -lth ../Zhang_LR/fastq/
```

## check the data

```{r engine='bash', eval=F, echo=TRUE}
head ../Zhang_LR/fastq/WT_D0_1.chr21.fastq
```

## `blat` the second sequence in UCSC genome browser


# Check the Quality of LR RNA-seq raw data

## Take a quick look at NanoPlot-the tool we're gonna use for this purpose

You can find the introduction of NanoPlot here: \

https://github.com/wdecoster/NanoPlot \

We simply load the module on HPC and check the help information to get a quick idead how it works.

```{r engine='bash', eval=F, echo=TRUE}
module load NanoPlot
NanoPlot -h
```
