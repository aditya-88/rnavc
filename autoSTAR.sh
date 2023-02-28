#/bin/bash

# AUTHOR: Aditya Singh

# STAR Aligner for NGS data
# version 0.1
# This program aligns NGS data with reference genome and annotation file
# It does the following steps:
# 1. Takes in an input directory with fastq files and an optional list of sample names to process
# 2. Finds all FASTQ/FASTQ.GZ files in the directory and selects the ones in the list, if provided
# 3. Automatically detects if the file is paired end or single end
# 4. Also detects if the sample was run on multiple lanes
# 5. Aligns the files with STAR
# 6. Creates a folder for each sample and puts the output files in it

# Usage: autostar <inputDir refDir annon.gtf outFolder threads>

# Check if the user has provided the correct number of arguments
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "autostar <inputDir refDir annon.gtf outFolder threads list>"
    exit
fi

# Check if the input directory exists
if [ ! -d "$1" ]; then
    echo "Input directory does not exist"
    exit
fi

# Check if the user provided the output directory, if not, create one using the input directory name
if [ -z "$4" ]; then
    echo "No output directory provided. Creating one using the input directory name"
    outDir="$1"_autoSTAR
    mkdir -p "$outDir"
else
    outDir="$4"
fi

# Check if the output directory exists, if not, create one
if [ ! -d "$outDir" ]; then
    echo "Output directory does not exist. Creating one"
    mkdir -p "$outDir"
fi

# Check if user provided threads, if not, detect and use all threads
if [ -z "$5" ]; then
    echo "No threads provided. Detecting and using all threads"
    threads=$(nproc)
else
    threads="$5"
fi

# Check if the reference directory exists
if [ ! -d "$2" ]; then
    echo "Reference directory does not exist"
    exit
fi

# Check if the annotation file exists
if [ ! -f "$3" ]; then
    echo "Annotation file does not exist"
    exit
fi

# Check if STAR is installed
star="star"
if ! command -v "$star" &> /dev/null
then
    echo "STAR could not be found. Please install STAR and make sure it is in your PATH"
    exit
fi

# Detect all FASTQ/FASTQ.GZ files in the input directory
echo "Detecting FASTQ/FASTQ.GZ files in the input directory"
files=$(find "$1" -type f -name "*.fastq" -o -name "*.fastq.gz")
# Check if the user provided a list of samples to process, if so, select only those files
if [ ! -z "$6" ]; then
    echo "Selecting files from the list"
    # Check if the list file exists
    if [ ! -f "$6" ]; then
        echo "List file does not exist"
        exit
    fi
    # Get the list of samples to process
    samples=$(cat "$6")
    # Select only the files that match the samples
    files=$(echo "$files" | grep -E "$samples")
fi
# Print the count of files detected divided by 2
echo "Detected $(echo "$files" | wc -l | awk '{print $1/2}') samples"

# Loop through all the files
for file in $files
do
    # Get the file name
    name=$(basename "$file" .fastq)
    name=$(basename "$name" .fastq.gz)
    name=$(basename "$name" .fq)
    name=$(basename "$name" .fq.gz)
    # Check if the file is paired end or single end
    if [[ $name == *R1* ]]; then
        # Check if the second file exists
        if [ ! -f "${file/R1/R2}" ]; then
            echo "Second file for $name does not exist. Skipping..."
            continue
        fi
        # Check if the sample has been processed already
        outFile="$outDir"/"$name"_autoSTAR_
        if ls "$outFile"* 1> /dev/null 2>&1; then
            echo "File already processed. Skipping..."
            continue
        fi
        echo "Processing $name"
        # Align the files
        # Check if compressed or not and align accordingly
        if [[ $file == *.gz ]]; then
            "$star" --runMode alignReads --genomeDir "$2" --runThreadN "$threads" --readFilesCommand zcat --readFilesIn "$file" "${file/R1/R2}" --sjdbGTFfile "$3" --outFileNamePrefix "$outFile" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
        else
            "$star" --runMode alignReads --genomeDir "$2" --runThreadN "$threads" --readFilesIn "$file" "${file/R1/R2}" --sjdbGTFfile "$3" --outFileNamePrefix "$outFile" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
        fi
    # Check if the sample has "R2" in the name, if so, skip it.
    elif [[ $name == *R2* ]]; then
        continue
    else
        # Check if the sample has been processed already
        outFile="$outDir"/"$name"_autoSTAR_
        if ls "$outFile"* 1> /dev/null 2>&1; then
            echo "File already processed. Skipping..."
            continue
        fi
        echo "Processing $name"
        # Align the files
        # Check if compressed or not and align accordingly
        if [[ $file == *.gz ]]; then
            "$star" --runMode alignReads --genomeDir "$2" --runThreadN "$threads" --readFilesCommand zcat --readFilesIn "$file" --sjdbGTFfile "$3" --outFileNamePrefix "$outFile" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
        else
            "$star" --runMode alignReads --genomeDir "$2" --runThreadN "$threads" --readFilesIn "$file" --sjdbGTFfile "$3" --outFileNamePrefix "$outFile" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
        fi
    fi
    # Delete temporary STAR files and directories
    echo "Deleting temporary files and directories"
    rm -rf "$outFile"STARtmp
    rm -rf "$outFile"STARgenome
done

# Merge all gene counts files into a single counts file
# Get all the gene count files in the output directory
geneCountFiles=$(find "$outDir" -type f -name "*ReadsPerGene.out.tab")
# Create a header
echo -e "GeneID\t$(basename "$geneCountFiles" | cut -d'_' -f1 | paste -sd '\t')" > "$outDir"/geneCounts.txt
# Get the gene names
geneNames=$(cut -f1 "$geneCountFiles" | sort | uniq)
# Loop through all the gene names
for gene in $geneNames
do
    # Get the counts for each gene
    counts=$(grep -w "$gene" "$geneCountFiles" | cut -f2 | paste -sd '\t')
    # Print the gene name and counts
    echo -e "$gene\t$counts" >> "$outDir"/geneCounts.txt
done


echo "All Done"
echo "Thanks for using autostar!"

