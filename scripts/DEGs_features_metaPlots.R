library(DMRcaller)
library(dplyr)
library(ggplot2)
library(GenomicFeatures)
library(plyranges)
library(parallel)

dir.create("/home/yoyerush/yo/methylome_pipeline/DEGs_metaPlots/", showWarnings = F)
dir.create("/home/yoyerush/yo/methylome_pipeline/DEGs_metaPlots/upregulated", showWarnings = F)
dir.create("/home/yoyerush/yo/methylome_pipeline/DEGs_metaPlots/downregulated", showWarnings = F)

source("/home/yoyerush/yo/methylome_pipeline/Methylome.At_test_051024/scripts/Gene_features_metaPlot_fun.R")

n.cores = 40

### mto1 and wt 'CX' file path
var_table = suppressWarnings(read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto1.txt", header = F, sep = "\t"))
vars_vector = unique(var_table[,1])
var1_path = var_table[grep(vars_vector[1], var_table[,1]), 2]
var2_path = var_table[grep(vars_vector[2], var_table[,1]), 2]

### load the data for replicates
source("/home/yoyerush/yo/methylome_pipeline/Methylome.At/scripts/load_replicates.R")
source("/home/yoyerush/yo/methylome_pipeline/Methylome.At/scripts/trimm_and_rename_seq.R")
load_vars = mclapply(list(var1_path,var2_path), function(x) load_replicates(x,n.cores,T), mc.cores = 2)
meth_var1 = trimm_and_rename(load_vars[[1]])
meth_var2 = trimm_and_rename(load_vars[[2]])

### annotation file
annotation.gr = read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/annotation_files/Methylome.At_annotations.csv.gz") %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  trimm_and_rename()

### DEGs list
RNA_up = read.csv("/home/yoyerush/yo/rnaseq_23_results_methionine/mto1_vs_wt/mto1_vs_wt.FC.up.DE.csv")[,1]
RNA_down = read.csv("/home/yoyerush/yo/rnaseq_23_results_methionine/mto1_vs_wt/mto1_vs_wt.FC.down.DE.csv")[,1]
non_DE = read.csv("/home/yoyerush/yo/rnaseq_23_results_methionine/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv")
non_DE = non_DE[non_DE$padj > 0.1, 1] %>% unique()

rna_list = list(upregulated = RNA_up,
                downregulated = RNA_down,
                non_DE = non_DE)


### main metaPlot function
for (rna in names(rna_list)) {
  
  dir.create(paste0("/home/yoyerush/yo/methylome_pipeline/DEGs_metaPlots/",rna), showWarnings = F)
  setwd(paste0("/home/yoyerush/yo/methylome_pipeline/DEGs_metaPlots/",rna))
  
  Genes_features_metaPlot(methylationPool_var1 = meth_var1,
                          methylationPool_var2 = meth_var2,
                          var1="wt",
                          var2="mto1",
                          annotations_file = annotation.gr,
                          n.random = rna_list[[rna]],
                          minReadsC = 6,
                          binSize = 10,
                          n.cores = n.cores)
  
  setwd("/home/yoyerush/yo/methylome_pipeline/DEGs_metaPlots/")
}










### test
if (F) {
  methylationPool_var1 = meth_var1
  methylationPool_var2 = meth_var2
  var1="wt"
  var2="mto1"
  annotations_file = annotation.gr
  n.random = RNA_up
  minReadsC = 6
  n.cores = 40
}


