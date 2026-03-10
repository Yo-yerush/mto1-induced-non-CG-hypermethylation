library(DMRcaller)
library(parallel)
library(ggplot2)
library(dplyr)

setwd("/home/yoyerush/yo/methylome_pipeline/total_meth_values_mto1/")
scripts_dir = "/home/yoyerush/yo/methylome_pipeline/Methylome.At/Methylome.At_scripts/"

# path for CX_report files
sample_df = read.table("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto1.txt")
wt_path = sample_df[sample_df$V1 == "wt", 2]
mto1_path = sample_df[sample_df$V1 == "mto1", 2]

# load CX files
source(paste0(scripts_dir,"load_replicates.R"))
load_vars = mclapply(list(wt_path,mto1_path), function(x) {load_replicates(x, n.cores = 6)}, mc.cores = 2)

# trimm ChrC and ChrM
trimm_Chr <- function(gr_obj) {
  remove_seqnames = c("NC_000932.1","NC_037304.1")
  gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
  seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
  return(sort(gr_obj))
}

wt_replicates = trimm_Chr(load_vars[[1]]$methylationDataReplicates)
mto1_replicates = trimm_Chr(load_vars[[2]]$methylationDataReplicates)

source(paste0(scripts_dir,"total_meth_levels.R"))
total_meth_levels(wt_replicates, mto1_replicates, "wt", "mto1")

## plot total methylation levels (5-mC%)
#source(paste0(scripts_dir,"total_meth_levels.R"))
#total_meth_levels(load_wt,load_mto1, "wt", "mto1", "gray30", "maroon3")