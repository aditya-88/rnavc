#!/bin/bash

# This scipt takes the user input folder and merges all found STAR gene counts files into a single file

# Usage: ./merge_counts.sh <input_folder>

# Check if the user provided an input folder
if [ -z "$1" ]
then
    echo "Please provide an input folder"
    exit 1
fi

# Check if the input folder exists
if [ ! -d "$1" ]
then
    echo "Input folder does not exist"
    exit 1
fi

# Variables
input_folder=$1

# Create a new folder to store the merged files in the input folder and set it as output_folder
mkdir -p $input_folder/merged
output_folder=$input_folder/merged

# Merge all STAR ReadsPerGene.out.tab files in the input folder into a single counts file removing the first four lines from each counts file
for file in $input_folder/*ReadsPerGene.out.tab
do
    # Extract the sample name from the file path
    sample=$(echo $file | rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1)

    # Create a new header file for each sample
    echo -e "Geneid\t$sample" > $output_folder/$sample.header

    # Extract the counts from the file
    cut -f 1,2 $file > $output_folder/$sample.counts

    # Merge the header and counts files
    cat $output_folder/$sample.header $output_folder/$sample.counts > $output_folder/$sample.merged
done

# Create a list of all merged files
files=$(ls $output_folder/*.merged)

# Merge all merged files into a single file
paste -d '\t' $files > $output_folder/merged_counts.tab

# Remove the header files
rm $output_folder/*.header

# Remove the counts files
rm $output_folder/*.counts

# Remove the merged files
rm $output_folder/*.merged

# Remove the first column (geneid) from the merged_counts.tab file
cut -f 2- $output_folder/merged_counts.tab > $output_folder/merged_counts_cut.tab

# Remove the first line (geneid) from the merged_counts.tab file
tail -n +2 $output_folder/merged_counts.tab > $output_folder/merged_counts_tail.tab

# Remove the merged_counts.tab file
rm $output_folder/merged_counts.tab

# Move the merged_counts_cut.tab file to the merged_counts.tab file
mv $output_folder/merged_counts_cut.tab $output_folder/merged_counts.tab

# Move the merged_counts_tail.tab file to the merged_counts.tab file
mv $output_folder/merged_counts_tail.tab $output_folder/merged_counts.tab




