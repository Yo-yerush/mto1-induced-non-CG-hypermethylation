#########################################################################################################
#   Calculate delta methylation in methylome '.wig' files
#
# Description:
#   This example script using input methylome from 'Stroud et al.' published data.
#
#   Stroud, H., Greenberg, M. V., Feng, S., Bernatavichute, Y. V., & Jacobsen, S. E. (2013).
#   Comprehensive analysis of silencing mutants reveals complex regulation of the Arabidopsis methylome.
#   Cell, 152(1), 352-364.

#   SRA: SRP014726
#   BioProject: PRJNA172021
#   GEO: GSE39901
#   https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP014726&o=biosample_s%3Aa%3Bacc_s%3Aa
#########################################################################################################

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)

file_path <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/GSE39901_RAW/"
output_dir <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/mutants_figs/"

source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/edit_TE_file.R")

############# ############# ############# #############

import_wt <- function(cntx) {
    cat(paste0("import 'wt' files in ", cntx, " context..."))
    wt_2 <- import(paste0(file_path, "GSM980986_WT_rep2_", cntx, ".wig.gz"), format = "wig") %>% as.data.frame()
    cat("...")
    wt_3 <- import(paste0(file_path, "GSM980987_WT_rep3_", cntx, ".wig.gz"), format = "wig") %>% as.data.frame()
    cat("...")
    wt_merged <- merge(wt_2, wt_3, by = c("seqnames", "start", "end"), suffixes = c("_2", "_3"))
    cat("...")
    wt_merged$score_wt <- (abs(wt_merged$score_2) + abs(wt_merged$score_3)) / 2
    wt_merged$seqnames <- gsub("chr", "Chr", wt_merged$seqnames)
    cat(" done\n")
    return(select(wt_merged, seqnames, start, end, score_wt))
}

wt_cg <- import_wt("CG")
wt_chg <- import_wt("CHG")
wt_chh <- import_wt("CHH")

############# ############# ############# #############

mut_list <- c("ddm1" = "GSM981009", "met1" = "GSM981031", "cmt2" = "GSM981002", "cmt3" = "GSM981003")


for (mutant_f in names(mut_list)) {
    cntx_res <- function(mutant_ff, mut_n, context_f) {
        comparison_name <- paste0(mutant_ff, "_vs_wt")

        cat(paste0("import '", mutant_f, "' file in ", context_f, " context..."))
        mut_file <- import(paste0(file_path, mut_n, "_", mutant_ff, "_", context_f, ".wig.gz"), format = "wig") %>% as.data.frame()
        mut_file$score <- abs(mut_file$score)
        mut_file$seqnames <- gsub("chr", "Chr", mut_file$seqnames)

        if (context_f == "CG") {
            wt_cntx <- wt_cg
        } else if (context_f == "CHG") {
            wt_cntx <- wt_chg
        } else if (context_f == "CHH") {
            wt_cntx <- wt_chh
        }

        cat(" delta...")
        delta_merged <- merge(mut_file, wt_cntx, by = c("seqnames", "start", "end"), suffixes = c("_mut", "_wt")) %>%
            mutate(delta = (score - score_wt)) %>%
            mutate(context = context_f) %>%
            select(seqnames, start, end, delta, context)

        cat(" done\n")
        return(delta_merged)
    }

    CG_delta <- cntx_res(mutant_f, mut_list[mutant_f], "CG")
    CHG_delta <- cntx_res(mutant_f, mut_list[mutant_f], "CHG")
    CHH_delta <- cntx_res(mutant_f, mut_list[mutant_f], "CHH")

    write.csv(rbind(CG_delta, CHG_delta, CHH_delta), 
              paste0(output_dir, mutant_f, "_delta_df.csv.gz"), 
              row.names = FALSE, 
              quote = FALSE)
}