#!/bin/bash
#-------------------------------------------------------------
# epiTEome tool
# https://github.com/jdaron/epiTEome
# https://doi.org/10.1186/s13059-017-1232-0
#-------------------------------------------------------------
#
# Usage: ./epiTEome_te_insertion_sites_script.sh [--dont_indx]
#   --dont_indx          Skip genome indexing step
#   --dont_concatenate   Skip concatenation of unmapped reads
#   --run_test           Use test data
#   --custom_test_list   TE list as 'te_ids_list' argument
#
#-------------------------------------------------------------
#
# # first get unmapped reads using:
# ./run_bismark_yo.sh -s ./mto1_samples_bismark.txt -g /home/yoyerush/yo/methylome_pipeline/Bismark/res_unmapped_040825/TAIR10_chr_all.fas.gz -o bismark_results -n 32 -m 16G --um
#
#-------------------------------------------------------------
#
# # there is error writing: << Unknown command 'filter' >>
# # using midofied epiTEome.pl script - use samtools=1.16 instead of the 'ngsutils' package
# # also can modify the 'bamutils' script. example in 'install_3piTEome_env.sh' script
#
#-------------------------------------------------------------

# epiTEome scripts dir - also the output files path
main_dir=/home/yoyerush/yo/methylome_pipeline/transposition_epiTEome
samples_name=("wt_1" "mto1_1" "mto1_2" "mto1_3" "wt_2")
te_ids_list=rndm_retrotransposons_list.txt
test_te_list=test_list.txt
bismark_results=/home/yoyerush/yo/methylome_pipeline/Bismark/res_unmapped_040825/bismark_results
genome_file=TAIR10_chr_all.fna
gff_file=tair10TEs.gff3
n_threads=16
batch_size=200
# MAX_JOBS=5 # can use number of samples if want

# genome indexed file with - '.epiTEome.masked.fasta' suffix
# genome_indx_file=$(basename "$genome_file" | sed 's/\.[^.]*$/.epiTEome.masked.fasta/')

# results directory suffix
res_suffix=$(basename "$te_ids_list" _list.txt)

cd $main_dir

#-------------------------------------------------------------

# '--dont_concatenate' arguments
dont_concatenate=0
for arg in "$@"; do
    if [[ "$arg" == "--dont_concatenate" ]]; then
        dont_concatenate=1
        break
    fi
done

# '--dont_indx' arguments
dont_indx=0
for arg in "$@"; do
    if [[ "$arg" == "--dont_indx" ]]; then
        dont_indx=1
        break
    fi
done

# '--run_test' arguments
run_test=0
for arg in "$@"; do
    if [[ "$arg" == "--run_test" ]]; then
        run_test=1
        break
    fi
done

# # '--custom_test_list' arguments
# test_te_list=../test_data/teid.lst
# for arg in "$@"; do
#     if [[ "$arg" == "--custom_test_list" ]]; then
#         test_te_list=../te_lists/"$te_ids_list"
#         break
#     fi
# done

#-------------------------------------------------------------


if [[ $run_test -eq 0 ]]; then
    
    mkdir -p "$main_dir"/results_"$res_suffix"
    
    
    # index the genome
    if [[ $dont_indx -eq 0 ]]; then
        mkdir -p "$main_dir"/results_"$res_suffix"/genome_indx_"$res_suffix"
        cd "$main_dir"/results_"$res_suffix"/genome_indx_"$res_suffix"
        cp "$main_dir"/genome_indx/"$genome_file" ./
        cp "$main_dir"/genome_indx/"$gff_file" ./
        perl ../../epiteome_scripts/idxEpiTEome.pl -l 150 -gff $gff_file -t ../../te_lists/"$te_ids_list" -fasta $genome_file
        cd ../
    fi
    
    # PER_JOB_THREADS=$(( n_threads / MAX_JOBS ))
    # # PER_JOB_THREADS=$((n_threads / ${#samples_name[@]})) # n_threads divide by number of samples
    # export main_dir res_suffix bismark_results gff_file te_ids_list dont_concatenate PER_JOB_THREADS
    #
    # process_one() {
    #     sample="$1"
    #     workdir="$main_dir/results_${res_suffix}/${sample}"
    #     mkdir -p "$workdir"
    #     cd "$workdir"
    #     export TMPDIR="$main_dir/tmp/${sample}"
    #     export SAMTOOLS_TMPDIR="$TMPDIR"
    #     mkdir -p "$TMPDIR"
    #
    #     if [[ $dont_concatenate -eq 0 ]]; then
    #         zcat "$bismark_results/${sample}/${sample}_unmapped_reads_1.fq.gz" \
    #         "$bismark_results/${sample}/${sample}_unmapped_reads_2.fq.gz" \
    #         > "${sample}_unmapped_reads.fq"
    #     fi
    #
    #     perl ../../epiteome_scripts/epiTEome_Yo_edit.pl -gff ../../genome_indx/"$gff_file" -t ../../te_lists/"$te_ids_list" -ref ../genome_indx_"$res_suffix"/TAIR10_chr_all.epiTEome.masked.fasta -un "${sample}_unmapped_reads.fq" -p "$PER_JOB_THREADS"
    #
    #     echo -e "\n***   Completed processing ${sample}    ***\n"
    # }
    # export -f process_one
    
    # printf "%s\n" "${samples_name[@]}" | parallel -j "$MAX_JOBS" --halt soon,fail=1 process_one {}
    
    # Loop over samples
    for sample in "${samples_name[@]}"; do
        
        mkdir -p "$main_dir"/results_"$res_suffix"/"$sample"
        cd "$main_dir"/results_"$res_suffix"/"$sample"
        
        # concatenate paired unmapped reads
        if [[ $dont_concatenate -eq 0 ]]; then
            zcat $bismark_results/"$sample"/"$sample"_unmapped_reads_1.fq.gz $bismark_results/"$sample"/"$sample"_unmapped_reads_2.fq.gz > "$sample"_unmapped_reads.fq
        fi
        
        # run epiTEome
        perl ../../epiteome_scripts/epiTEome_Yo_edit.pl -gff ../genome_indx_"$res_suffix"/"$gff_file" -t ../../te_lists/"$te_ids_list" -ref ../genome_indx_"$res_suffix"/TAIR10_chr_all.epiTEome.masked.fasta -un "$sample"_unmapped_reads.fq -p $n_threads -b $batch_size
        
        cd $main_dir
        
        echo ""
        echo "***   Completed processing $sample    ***"
        echo ""
    done
    
else
    res_suffix_test=$(basename "$test_te_list" _list.txt)
    test_data_res=test_data_res_"$res_suffix_test"
    mkdir "$main_dir"/"$test_data_res"
    mkdir "$main_dir"/"$test_data_res"/genome_indx_test
    cp "$main_dir"/test_data/unmapped.fastq.bz2 "$main_dir"/"$test_data_res"
    cp "$main_dir"/test_data/Chr2.fasta "$main_dir"/"$test_data_res"/genome_indx_test
    
    cd "$main_dir"/"$test_data_res"
    
    perl ../epiteome_scripts/idxEpiTEome.pl -l 150 -gff ../genome_indx/tair10TEs.gff3 -t ../te_lists/"$test_te_list" -fasta genome_indx_test/Chr2.fasta
    
    cd "$main_dir"/"$test_data_res"
    
    perl ../epiteome_scripts/epiTEome_Yo_edit.pl -gff ../genome_indx/tair10TEs.gff3 -t ../te_lists/"$test_te_list" -ref genome_indx_test/Chr2.epiTEome.masked.fasta -un unmapped.fastq.bz2 -p 1
    
    cd $main_dir
fi


##########################################
### R script to find overlap TEGs in mto1 vs wt (insert to met1 github page in different file)

# library(dplyr)
# 
# cg_file <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/TEG_CG_genom_annotations.csv"
# chg_file <- gsub("/CG/TEG_CG_", "/CHG/TEG_CHG_", cg_file)
# chh_file <- gsub("/CG/TEG_CG_", "/CHH/TEG_CHH_", cg_file)
# cg_data <- read.csv(cg_file)
# chg_data <- read.csv(chg_file)
# chh_data <- read.csv(chh_file)
# methylation_data <- rbind(cg_data, chg_data, chh_data) %>%
# distinct(gene_id, .keep_all = T) %>%
# select(gene_id, Derives_from)
# 
# rnaseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all_genes_results_mto1_vs_wt.csv") %>%
# filter(padj < 0.05) %>%
# select(gene_id)
# 
# overlap_tegs <- merge(methylation_data, rnaseq, by = "gene_id") %>%
# filter(Derives_from != "") %>%
# distinct(Derives_from)
# 
# write.table(overlap_tegs, "overlap_TEGs_list.txt", row.names = F, quote = F, col.names = F)