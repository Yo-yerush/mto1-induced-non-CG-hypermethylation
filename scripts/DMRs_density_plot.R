library(dplyr)
library(circlize)
library(GenomicRanges)

### load functions
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/edit_TE_file.R")

### load annotation file
ann.file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_annotations.csv.gz") %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  trimm_and_rename()

### load TEs annotation file
TE_4_dens <- edit_TE_file(read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10_Transposable_Elements.txt", sep = "\t"))

#####################################

chr_amount <- length(seqnames(ann.file)@values)
# genes_type <- ann.file[which(ann.file$type == "gene")]

cntx_file <- function(context) {
  ############# read DMRs file
  dmrs_file <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/DMRs_", context, "_mto1_vs_wt.csv"))
  dmrs_file <- dmrs_file[, c("seqnames", "start", "end", "log2FC")]
  return(dmrs_file)
}

CG_file <- cntx_file("CG")
CHG_file <- cntx_file("CHG")
CHH_file <- cntx_file("CHH")

#####################################
############# the plot #############
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_Density_mto1_vs_wt.svg"), width = 3.25, height = 3.25, family = "serif")

circos.par(start.degree = 90)
circos.genomicInitialize(as.data.frame(ann.file)[, 1:3], sector.names = paste0("Chr ", seq(chr_amount)), axis.labels.cex = 0.325, labels.cex = 1.35)

circos.genomicDensity(
  list(
    CG_file[CG_file$log2FC > 0, 1:3],
    CG_file[CG_file$log2FC < 0, 1:3]
  ),
  bg.col = "#fafcff", bg.border = NA, count_by = "number",
  col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
)

circos.genomicDensity(
  list(
    CHG_file[CHG_file$log2FC > 0, 1:3],
    CHG_file[CHG_file$log2FC < 0, 1:3]
  ),
  bg.col = "#fafcff", bg.border = NA, count_by = "number",
  col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
)

circos.genomicDensity(
  list(
    CHH_file[CHH_file$log2FC > 0, 1:3],
    CHH_file[CHH_file$log2FC < 0, 1:3]
  ),
  bg.col = "#fafcff", bg.border = NA, count_by = "number",
  col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
)

circos.genomicDensity(
  list(
    # as.data.frame(genes_type)[1:3],
    as.data.frame(TE_4_dens)[1:3]
  ),
  bg.col = "#fafcff", bg.border = NA, count_by = "number",
  col = c("#fcba0320"), border = T, track.height = 0.165, track.margin = c(0, 0)
) # "gray80", 

circos.clear()
dev.off()
