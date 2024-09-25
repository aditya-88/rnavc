# RNA-VC #

[![DOI](https://zenodo.org/badge/582600349.svg)](https://zenodo.org/doi/10.5281/zenodo.13837677)

## Variant calling pipeline for RNA-Seq data ##

### Introduction ###

The script is based on the GATK best practices workflow for RNA-Seq data.
https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
Script developed by Aditya Singh
https://github.com/aditya-88
dr.singhaditya@hotmail.com

The script takes in the following arguments:
1. Path to the reference genome
2. Path to the input BAM file aligned using STAR two-pass mode
3. Path to the output directory
4. Known sites for GATK BaseRecalibrator

The script performs the following steps:
1. Mark duplicates
2. Split'N'Trim and reassign mapping qualities
3. Base quality score recalibration
4. Variant calling
5. Variant filtering

The script requires the following tools to be installed:
1. GATK
2. Picard tools (for MarkDuplicates, is part of GATK)
3. STAR (for alignment, not required for this script)

