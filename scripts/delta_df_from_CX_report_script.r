###########################
#                         #
#  $ conda activate Renv  #
#                         #
###########################

output_dir <- "/home/yoyerush/yo/methylome_pipeline/mto1_delta_df_4_paper/"

wt_path <- c(
    "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S18/methylation_extractor/S18_R1_bismark_bt2_pe.CX_report.txt",
    "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S19/methylation_extractor/S19_R1_bismark_bt2_pe.CX_report.txt"
)
mutant_path <- c(
    "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S9/methylation_extractor/S9_R1_bismark_bt2_pe.CX_report.txt",
    "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S10/methylation_extractor/S10_R1_bismark_bt2_pe.CX_report.txt",
    "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S11/methylation_extractor/S11_R1_bismark_bt2_pe.CX_report.txt"
)

n.cores <- 6

############# ############# ############# #############

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(DMRcaller)
library(parallel)

source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/load_replicates.R")

############# ############# ############# #############

# load replicates data
n.cores.load <- ifelse(n.cores > 1, 2, 1)
load_vars <- mclapply(list(wt_path, mutant_path), function(x) load_replicates(x, n.cores, T), mc.cores = n.cores.load)

meth_wt <- trimm_and_rename(load_vars[[1]]) %>%
as.data.frame() %>%
mutate(Proportion = readsM/readsN)

meth_mutant <- trimm_and_rename(load_vars[[2]]) %>%
as.data.frame() %>%
mutate(Proportion = readsM/readsN)

############# ############# ############# #############

meth_delta <- meth_wt %>%
mutate(delta = meth_mutant$Proportion - meth_wt$Proportion) %>%
    select(seqnames, start, end, delta, context) %>%
    filter(!is.nan(delta))


CG_delta <- meth_delta[meth_delta$context == "CG", ]
CHG_delta <- meth_delta[meth_delta$context == "CHG", ]
CHH_delta <- meth_delta[meth_delta$context == "CHH", ]

write.csv(rbind(CG_delta, CHG_delta, CHH_delta),
    paste0(output_dir, "mto1_delta_df.csv.gz"),
    row.names = FALSE,
    quote = FALSE
)
#