#!/bin/bash

n_chr_2_filter=Chr3
# ncbi_chr_name=NC_003071.7 # chr2
ncbi_chr_name=NC_003074.8 # chr3
output_dir=/home/yoyerush/yo/filter_bam_results_by_chromosome/
cd $output_dir

# Loop through sample names
for sample_name in 9 10 11 18 19; do
    echo "Processing sample: $sample_name"
    
    # Input BAM file
    input_bam=/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S"$sample_name"/S"$sample_name"_sorted.bam
    output_bam=S"$sample_name"_"$n_chr_2_filter"_sorted.bam

    # Filter BAM file to keep only chromosome n
    samtools view -b "$input_bam" "$ncbi_chr_name" > "$output_bam"

    # Change chromosome name from NCBI format to Chr3 format
    samtools view -h "$output_bam" | sed 's/'"$ncbi_chr_name"'/Chr3/g' | samtools view -b > temp_"$output_bam"
    mv temp_"$output_bam" "$output_bam"

    # Index the filtered BAM file
    samtools index "$output_bam"

    echo "Filtered BAM file created: $output_bam"
done

echo ""
echo "All samples processed."