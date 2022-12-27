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

# Exit on error from any part of the script
set -eE
# Exit the entire script if any command in a pipe fails
set -o pipefail
# Exit the entire script if ctrl-c is pressed
trap "exit" INT

# Print welcome message
echo "Welcome to the RNA-Seq variant calling pipeline!"

# Arguments and options
REF=$1
BAM=$2
KNOWN=$3
OUT=$4 #optional
cpus=$5 #optional
mem=$6 #optional
gatk=$7 #optional

# Set number of threads to all available cpus if not provided
if [ ! $cpus ]; then
    cpus=$(nproc)
fi

# Check of memory provided or else get from the system
if [ ! $mem ]; then
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
    # Set memory to 90% of total memory
    mem=$(($mem * 90 / 100))
fi

# Check if no arguments are provided, if yes print usage
if [ $# -eq 0 ]; then
    echo "Usage: bash rnavc <reference genome> <input BAM file> <known sites> <output directory> <threads> <memory in GB> <GATK 4.x executable>"
    exit 1
fi
## Initial checks

echo "Checking arguments"
# Check if the reference genome is provided, if not exit
if [ ! $REF ]; then
    echo "X Please provide the path to the reference genome"
    exit 1
fi

# Check if the input BAM file is provided, if not exit
if [ ! $BAM ]; then
    echo "X Please provide the path to the input BAM file"
    exit 1
fi
# Check if the known sites are provided, if not then exit
if [ ! $KNOWN ]; then
    echo "X Please provide the path to the known sites"
    exit 1
fi

echo "> Required arguments provided"

echo "Checking files"
# Check if provided files exist
if [ ! -f $REF ]; then
    echo "X Reference genome not found"
    exit 1
fi
if [ ! -f $BAM ]; then
    echo "X Input BAM file not found"
    exit 1
fi
if [ ! -f $KNOWN ]; then
    echo "X Known sites not found"
    exit 1
fi
echo "> Provided files exist"

# Get sample name
sample=$(basename $BAM | cut -d. -f1)

# Check if the output directory is provided, if not get sample directory
if [ ! $OUT ]; then
    OUT=$(dirname $BAM)/$sample
else
    OUT=$OUT/$sample
fi

# Check if the output directory exists, if not create it
if [ ! -d $OUT ]; then
    mkdir -p $OUT
fi
# Check if the GATK executable is provided, if not use the default
if [ ! $gatk ]; then
    gatk=gatk
fi

# Delete the logs and error files if they exist
if [ -f $OUT/$sample.log ]; then
    rm $OUT/$sample.log
fi
if [ -f $OUT/$sample.err ]; then
    rm $OUT/$sample.err
fi

echo "Started: "$(date) > $OUT/$sample.log
# Delete all empty files in the output directory
echo "Housekeeping: Deleting empty files from the output folder" >>$OUT/$sample.log
find $OUT -type f -empty -delete

# Checking if the pipeline was run for the sample
if [ -f $OUT/$sample.filtered.vcf.idx ]; then
    echo -e "The pipeline was already run for the sample:\n$sample\nExiting."
    exit
fi

echo -e "Running the pipeline\nCheck the log and error files for details and progress:\n$OUT/$sample.log\n$OUT/$sample.err"

# Saving parameters and system information
echo "Allocated threads: "$cpus >> $OUT/$sample.log
echo "Allocated memory: "$mem"G" >> $OUT/$sample.log
echo "Parameters:" >> $OUT/$sample.log
echo "Reference genome: "$REF >> $OUT/$sample.log
echo "Input BAM file: "$BAM >> $OUT/$sample.log
echo "Known sites: "$KNOWN >> $OUT/$sample.log
echo "Output directory: "$OUT >> $OUT/$sample.log
echo "GATK executable: "$gatk >> $OUT/$sample.log
echo "---" >> $OUT/$sample.log

# Check if the BAM file has read groups, if not add them
echo "Checking read groups" >> $OUT/$sample.log
if [ $(samtools view -H $BAM | grep -c "@RG") -eq 0 ]; then
    echo ">Read groups not found, adding." >> $OUT/$sample.log
    # Check if the part was run before, if yes skip
    if [ -f $OUT/$sample.rg.bam ]; then
        echo ">AddOrReplaceReadGroups already run" >> $OUT/$sample.log
    else
        $gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" AddOrReplaceReadGroups -I $BAM -O $OUT/$sample.rg.bam -RGID $sample -RGLB $sample -RGPL UNKNOWN -RGPU $sample -RGSM $sample 1>> $OUT/$sample.log 2>> $OUT/$sample.err
    fi
    # Set the BAM file to the one with read groups
    BAM=$OUT/$sample.rg.bam
fi

# Mark duplicates
echo "Marking duplicates" >> $OUT/$sample.log
# Check if the part was run before, if yes skip
if [ -f $OUT/$sample.markdup.bam ]; then
    echo ">MarkDuplicates already run" >> $OUT/$sample.log
else
    $gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" MarkDuplicates -I $BAM -O $OUT/$sample.markdup.bam -M $OUT/$sample.markdup.metrics.txt 1>> $OUT/$sample.log 2>> $OUT/$sample.err
fi

# Split'N'Trim and reassign mapping qualities
echo "Split'N'Trim and reassign mapping qualities" >> $OUT/$sample.log
# Check if the part was run before, if yes skip
if [ -f $OUT/$sample.split.bam ]; then
    echo ">SplitNCigarReads already run" >> $OUT/$sample.log
else
    $gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" SplitNCigarReads -R $REF -I $OUT/$sample.markdup.bam -O $OUT/$sample.split.bam 1>> $OUT/$sample.log 2>> $OUT/$sample.err
fi

# Base quality score recalibration
echo "Base quality score recalibration" >> $OUT/$sample.log
# Check if the VCF file is indexed and if not index it
if [ ! -f $KNOWN.idx ]; then
    echo ">Indexing known sites" >> $OUT/$sample.log
    $gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" IndexFeatureFile -F $KNOWN 1>> $OUT/$sample.log 2>> $OUT/$sample.err
fi
# Check if the part was run before, if yes skip
if [ -f $OUT/$sample.recal.bam ]; then
    echo ">BaseRecalibrator and ApplyBQSR already run" >> $OUT/$sample.log
else
    $gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" BaseRecalibrator -R $REF -I $OUT/$sample.split.bam --known-sites $KNOWN -O $OUT/$sample.recal_data.table 1>> $OUT/$sample.log 2>> $OUT/$sample.err
    $gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" ApplyBQSR -R $REF -I $OUT/$sample.split.bam --bqsr-recal-file $OUT/$sample.recal_data.table -O $OUT/$sample.recal.bam 1>> $OUT/$sample.log 2>> $OUT/$sample.err
fi

# Variant calling
echo "Variant calling" >> $OUT/$sample.log
# Check if the part was run before, if yes skip
if [ -f $OUT/$sample.raw.vcf ]; then
    echo ">HaplotypeCaller already run" >> $OUT/$sample.log
else
    $gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" HaplotypeCaller -R $REF -I $OUT/$sample.recal.bam -O $OUT/$sample.raw.vcf 1>> $OUT/$sample.log 2>> $OUT/$sample.err
fi

# Variant filtering
echo "Variant filtering" >> $OUT/$sample.log
# Check if the part was run before, if yes skip
if [ -f $OUT/$sample.filtered.vcf ]; then
    echo ">VariantFiltration already run" >> $OUT/$sample.log
else
    $gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" VariantFiltration -R $REF -V $OUT/$sample.raw.vcf -O $OUT/$sample.filtered.vcf --filter-expression "QD < 2.0" --filter-name "QD2" --filter-expression "FS > 60.0" --filter-name "FS60" --filter-expression "MQ < 40.0" --filter-name "MQ40" --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 1>> $OUT/$sample.log 2>> $OUT/$sample.err
fi


# Index the output VCF file
echo "Indexing the output VCF file" >> $OUT/$sample.log
# Check if the part was run before, if yes skip
if [ -f $OUT/$sample.filtered.vcf.idx ]; then
    echo ">IndexFeatureFile already run" >> $OUT/$sample.log
else
    $gatk --java-options "-Xmx"$mem"G -XX:ParallelGCThreads="$cpus"" IndexFeatureFile -I $OUT/$sample.filtered.vcf 1>> $OUT/$sample.log 2>> $OUT/$sample.err
fi

# Clean up
echo "Cleaning up" >> $OUT/$sample.log
rm $OUT/$sample.rg.bam
rm $OUT/$sample.markdup.bam
rm $OUT/$sample.markdup.metrics.txt
rm $OUT/$sample.split.bam
rm $OUT/$sample.recal_data.table

echo "Completed: "$(date) > $OUT/$sample.log
echo "---" >> $OUT/$sample.log
echo "For any questions, please contact: " >> $OUT/$sample.log
echo "Aditya Singh" >> $OUT/$sample.log
echo "dr.singhaditya@hotmail.com" >> $OUT/$sample.log
echo "CMM, Karolinska Institutet, Stockholm, Sweden" >> $OUT/$sample.log
