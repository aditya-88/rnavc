#!/bin/bash

# The script is based on the GATK best practices workflow for RNA-Seq data.
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
# Script developed by Aditya Singh
# https://github.com/aditya-88
# dr.singhaditya@hotmail.com

# The script takes in the following arguments:
# 1. Path to the reference genome
# 2. Path to the input BAM file aligned using STAR two-pass mode
# 3. Path to the output directory
# 4. Known sites for GATK BaseRecalibrator

# The script performs the following steps:
# 1. Mark duplicates
# 2. Split'N'Trim and reassign mapping qualities
# 3. Base quality score recalibration
# 4. Variant calling
# 5. Variant filtering

# The script requires the following tools to be installed:
# 1. GATK
# 2. Picard tools (for MarkDuplicates, is part of GATK)
# 3. STAR (for alignment, not required for this script)

# Arguments and options
REF=$1
BAM=$2
KNOWN=$3 #optional. Defaults to 1000G GRCh38 known sites from GATK resource bundle: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
OUT=$4 #optional
gatk=$5 #optional

# Get number of processors
cpus=$(nproc)
# Set number of threads to all but one processor
cpus=$(($cpus - 1))

# Get OS type
OSTYPE=$(uname)

# Check if Mac or Linux
if [[ "$OSTYPE" == "Darwin"* ]]; then
    # Mac OSX
    # Get total memory
    mem=$(sysctl -n hw.memsize)
    # Convert to GB
    mem=$(($mem / 1024 / 1024 / 1024))
else
    # Linux
    # Get total memory
    mem=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    # Convert to GB
    mem=$(($mem / 1024 / 1024))
fi
# Set memory to 85% of total memory
mem=$(($mem * 85 / 100))

# Get complete address of the current directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Check if no arguments are provided, if yes print usage
if [ $# -eq 0 ]; then
    echo "Usage: bash run.sh <reference genome> <input BAM file> <known sites> <output directory> <GATK executable>"
    exit 1
fi

# Check if the reference genome is provided, if not exit
if [ ! $REF ]; then
    echo "Please provide the path to the reference genome"
    exit 1
fi

# Check if the input BAM file is provided, if not exit
if [ ! $BAM ]; then
    echo "Please provide the path to the input BAM file"
    exit 1
fi

# Get sample name
sample=$(basename $BAM | cut -d. -f1)

# Check if the output directory is provided, if not get sample directory
if [ ! $OUT ]; then
    OUT=$(dirname $BAM)/$sample
fi

# Check if the output directory exists, if not create it
if [ ! -d $OUT ]; then
    mkdir -p $OUT
fi
# Check if the known sites are provided, if not use the default
if [ ! $KNOWN ]; then
    KNOWN=$DIR/db/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
fi
# Check if the GATK executable is provided, if not use the default
if [ ! $gatk ]; then
    gatk=gatk
fi

# Saving parameters and system information
echo $(date) > $OUT/$sample.log
echo "::System information::" >> $OUT/$sample.log
echo "Number of processors: "$cpus >> $OUT/$sample.log
echo "Total memory: "$mem"G" >> $OUT/$sample.log
echo "::Parameters::" >> $OUT/$sample.log
echo "Reference genome: "$REF >> $OUT/$sample.log
echo "Input BAM file: "$BAM >> $OUT/$sample.log
echo "Known sites: "$KNOWN >> $OUT/$sample.log
echo "Output directory: "$OUT >> $OUT/$sample.log
echo "GATK executable: "$gatk >> $OUT/$sample.log
echo "---" >> $OUT/$sample.log

# Mark duplicates
echo "Marking duplicates" >> $OUT/$sample.log
$gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" MarkDuplicates -I $BAM -O $OUT/$sample.markdup.bam -M $OUT/$sample.markdup.metrics.txt 1>> $OUT/$sample.log 2>> $OUT/$sample.err

# Split'N'Trim and reassign mapping qualities
echo "Split'N'Trim and reassign mapping qualities" >> $OUT/$sample.log
$gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" SplitNCigarReads -R $REF -I $OUT/$sample.markdup.bam -O $OUT/$sample.split.bam 1>> $OUT/$sample.log 2>> $OUT/$sample.err

# Base quality score recalibration
echo "Base quality score recalibration" >> $OUT/$sample.log
$gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" BaseRecalibrator -R $REF -I $OUT/$sample.split.bam --known-sites $KNOWN -O $OUT/$sample.recal_data.table 1>> $OUT/$sample.log 2>> $OUT/$sample.err
$gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" ApplyBQSR -R $REF -I $OUT/$sample.split.bam --bqsr-recal-file $OUT/$sample.recal_data.table -O $OUT/$sample.recal.bam 1>> $OUT/$sample.log 2>> $OUT/$sample.err

# Variant calling
echo "Variant calling" >> $OUT/$sample.log
$gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" HaplotypeCaller -R $REF -I $OUT/$sample.recal.bam -O $OUT/$sample.raw.vcf 1>> $OUT/$sample.log 2>> $OUT/$sample.err

# Variant filtering
echo "Variant filtering" >> $OUT/$sample.log
$gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" VariantFiltration -R $REF -V $OUT/$sample.raw.vcf -O $OUT/$sample.filtered.vcf --filter-expression "QD < 2.0" --filter-name "QD2" --filter-expression "FS > 60.0" --filter-name "FS60" --filter-expression "MQ < 40.0" --filter-name "MQ40" --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 1>> $OUT/$sample.log 2>> $OUT/$sample.err


# Index the output VCF file
echo "Indexing the output VCF file" >> $OUT/$sample.log
$gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" IndexFeatureFile -I $OUT/$sample.filtered.vcf 1>> $OUT/$sample.log 2>> $OUT/$sample.err

echo "---" >> $OUT/$sample.log
echo "Done" >> $OUT/$sample.log