##########################################
# Create conda environment
conda create -n bedgraph_env -y
conda activate bedgraph_env
conda install -y bioconda::ucsc-bedgraphtobigwig

mkdir /home/yoyerush/yo/filter_bam_results_by_chromosome/methylaion_bigwig
cd /home/yoyerush/yo/filter_bam_results_by_chromosome/methylaion_bigwig

../cx_2_bigwig.sh --at_rename /home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S9/methylation_extractor/S9_R1_bismark_bt2_pe.CX_report.txt mto1_1 Chr3

../cx_2_bigwig.sh --at_rename /home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S10/methylation_extractor/S10_R1_bismark_bt2_pe.CX_report.txt mto1_2 Chr3

../cx_2_bigwig.sh --at_rename /home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S11/methylation_extractor/S11_R1_bismark_bt2_pe.CX_report.txt mto1_3 Chr3

../cx_2_bigwig.sh --at_rename /home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S18/methylation_extractor/S18_R1_bismark_bt2_pe.CX_report.txt wt_1 Chr3

../cx_2_bigwig.sh --at_rename /home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S19/methylation_extractor/S19_R1_bismark_bt2_pe.CX_report.txt wt_2 Chr3
