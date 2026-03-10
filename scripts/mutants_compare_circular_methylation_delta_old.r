library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(circlize)

file_path <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/GSE39901_RAW/"
output_dir <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/mutants_figs/"

source("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/scripts/trimm_and_rename_seq.R")
source("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/scripts/edit_TE_file.R")

############# ############# ############# #############

ann.file <- read.csv("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/annotation_files/Methylome.At_annotations.csv.gz") %>%
     makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
     trimm_and_rename()

genes_type <- ann.file[which(ann.file$type == "gene")]

TE_4_dens <- edit_TE_file(read.csv("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/annotation_files/TAIR10_Transposable_Elements.txt", sep = "\t"))

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

mut_list <- c("cmt3" = "GSM981003", "ago4" = "GSM980991", "suvh8" = "GSM981062", "GSM981038" = "nrpb2")

############# start the plot #############
svg(paste0(output_dir, "DMRs_Density_mutants_data.svg"), width = 3.25 * length(mut_list), height = 3.25, family = "serif")

par(mfrow = c(1, length(mut_list)), mar = c(rep(0.5, 4)))

for (mutant_f in names(mut_list)) {
     cntx_res <- function(mutant_f, mut_n, context_f) {
          comparison_name <- paste0(mutant_f, "_vs_wt")

          cat(paste0("import '", mutant_f, "' file in ", context_f, " context..."))
          mut_file <- import(paste0(file_path, mut_n, "_", mutant_f, "_", context_f, ".wig.gz"), format = "wig") %>% as.data.frame()
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
               select(-score, -score_wt)
          cat(" done\n")
          return(delta_merged)
     }

     CG_delta <- cntx_res(mutant_f, mut_list[mutant_f], "CG")
     CHG_delta <- cntx_res(mutant_f, mut_list[mutant_f], "CHG")
     CHH_delta <- cntx_res(mutant_f, mut_list[mutant_f], "CHH")

     #####################################
     cat(paste0("plotting '", mutant_f, "'..."))
     circos.par(start.degree = 90)
     circos.genomicInitialize(as.data.frame(ann.file)[, 1:3], sector.names = paste0("Chr ", 1:5), axis.labels.cex = 0.325, labels.cex = 1.35)

     circos.genomicDensity(
          list(
               CG_delta[CG_delta$delta > 0, 1:3],
               CG_delta[CG_delta$delta < 0, 1:3]
          ),
          bg.col = "#fafcff", bg.border = NA, count_by = "number",
          col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
     )
     cat(".")
     circos.genomicDensity(
          list(
               CHG_delta[CHG_delta$delta > 0, 1:3],
               CHG_delta[CHG_delta$delta < 0, 1:3]
          ),
          bg.col = "#fafcff", bg.border = NA, count_by = "number",
          col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
     )
     cat(".")
     circos.genomicDensity(
          list(
               CHH_delta[CHH_delta$delta > 0, 1:3],
               CHH_delta[CHH_delta$delta < 0, 1:3]
          ),
          bg.col = "#fafcff", bg.border = NA, count_by = "number",
          col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
     )
     cat(".")
     circos.genomicDensity(
          list(
               as.data.frame(genes_type)[1:3],
               as.data.frame(TE_4_dens)[1:3]
          ),
          bg.col = "#fafcff", bg.border = NA, count_by = "number",
          col = c("gray80", "#fcba0320"), border = T, track.height = 0.165, track.margin = c(0, 0)
     )

     circos.clear()
     cat(" done\n\n\n")
}

dev.off()
