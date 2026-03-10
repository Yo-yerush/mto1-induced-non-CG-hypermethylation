#!/bin/bash

n_chr_2_filter=Chr3
output_dir=/home/yoyerush/yo/filter_bam_results_by_chromosome/rnaseq

# Check if chromosome name was provided as argument
if [ $# -eq 1 ]; then
    manual_chr_name="$1"
    echo "Using manually specified chromosome name: $manual_chr_name"
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"
cd "$output_dir"

# Function to find the correct chromosome name in BAM file
find_chr3_name() {
    local bam_file=$1
    echo "Available chromosomes in BAM file:"
    local chr_names=($(samtools view -H "$bam_file" | grep "^@SQ" | cut -f2 | sed 's/SN://'))
    
    # Display all available chromosomes with numbers for easy reference
    echo "Found ${#chr_names[@]} chromosomes:"
    for i in "${!chr_names[@]}"; do
        echo "$((i+1)). ${chr_names[i]}"
    done
    echo ""
    
    # Try to find Chr3 variants (including more possibilities)
    for chr_name in "${chr_names[@]}"; do
        case "$chr_name" in
            "Chr3"|"chr3"|"3"|"NC_003074.8"|"chromosome3"|"Chromosome3"|"CHROMOSOME_III"|"III")
                echo "Found potential Chr3 match: $chr_name"
                echo "$chr_name"
                return 0
                ;;
        esac
    done
    
    echo "No obvious Chr3 match found."
    echo "Please examine the list above and identify which one corresponds to chromosome 3."
    echo "Common patterns to look for:"
    echo "  - Something with '3' in it"
    echo "  - NCBI RefSeq starting with 'NC_003074'"
    echo "  - Roman numeral 'III'"
    echo "  - Any variant of 'Chr3', 'chr3', 'chromosome3'"
    return 1
}

# Find the correct chromosome name from the first BAM file
first_sample=14
first_bam=/home/yaelh/rnaseq_2023_results/met_pipeline_for_TE_expression/results_300524/sorted_bam/met"$first_sample"_sorted.bam

if [ ! -f "$first_bam" ]; then
    echo "Error: Cannot find first BAM file: $first_bam"
    exit 1
fi

# Use manual chromosome name if provided, otherwise auto-detect
if [ ! -z "$manual_chr_name" ]; then
    actual_chr_name="$manual_chr_name"
    echo "Using manually specified chromosome: $actual_chr_name"
else
    echo "Detecting chromosome naming convention..."
    actual_chr_name=$(find_chr3_name "$first_bam")

    if [ $? -ne 0 ] || [ -z "$actual_chr_name" ]; then
        echo ""
        echo "Error: Could not automatically detect Chr3 name."
        echo ""
        echo "Please run the script again with the correct chromosome name as an argument:"
        echo "Usage: $0 [chromosome_name]"
        echo "Example: $0 'NC_003074.8'"
        echo "Example: $0 '3'"
        echo "Example: $0 'chr3'"
        echo ""
        echo "Or modify the script to include your chromosome naming pattern."
        exit 1
    fi
fi

echo "Using chromosome name: $actual_chr_name"
echo ""

# Loop through sample names
for sample_name in 14 15 16 20 22 23; do
    echo "Processing sample: $sample_name"
    
    # Input BAM file
    input_bam=/home/yaelh/rnaseq_2023_results/met_pipeline_for_TE_expression/results_300524/sorted_bam/met"$sample_name"_sorted.bam
    output_bam=met"$sample_name"_"$n_chr_2_filter"_sorted.bam

    # Check if input BAM file exists
    if [ ! -f "$input_bam" ]; then
        echo "Warning: Input BAM file not found: $input_bam"
        continue
    fi

    echo "Input BAM: $input_bam"
    echo "Output BAM: $output_bam"
    echo "Filtering for chromosome: $actual_chr_name"

    # Filter BAM file to keep only chromosome 3 (using detected name)
    echo "Filtering BAM file..."
    if samtools view -b "$input_bam" "$actual_chr_name" > temp_"$output_bam"; then
        temp_size=$(stat -f%z temp_"$output_bam" 2>/dev/null || stat -c%s temp_"$output_bam" 2>/dev/null)
        if [ "$temp_size" -gt 1000 ]; then
            echo "Successfully filtered $temp_size bytes of data"
        else
            echo "Warning: Filtered file is very small ($temp_size bytes). Check if chromosome name is correct."
        fi
    else
        echo "Error: Failed to filter BAM file"
        continue
    fi

    # Rename chromosome to Chr3 format if it's not already Chr3
    if [ "$actual_chr_name" != "Chr3" ]; then
        echo "Renaming chromosome from $actual_chr_name to Chr3..."
        samtools view -h temp_"$output_bam" | sed 's/'"$actual_chr_name"'/Chr3/g' | samtools view -b > "$output_bam"
        rm temp_"$output_bam"
    else
        echo "Chromosome name is already Chr3, no renaming needed."
        mv temp_"$output_bam" "$output_bam"
    fi

    # Index the filtered BAM file
    samtools index "$output_bam"

    echo "Filtered BAM file created: $output_bam"
done

echo ""
echo "All samples processed."