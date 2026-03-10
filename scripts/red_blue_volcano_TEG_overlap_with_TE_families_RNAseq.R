library(dplyr)
library(GenomicRanges)
library(ggplot2)

###   ###    ###   ###   ###
##   ##    ##   ##   ##   ##
###   ###    ###   ###   ###

with_manual_Derives_from = F

###   ###    ###   ###   ###
##   ##    ##   ##   ##   ##
###   ###    ###   ###   ###


TEG_annotations <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_annotations.csv.gz") %>%
  filter(type == "transposable_element_gene") %>%
  dplyr::select(-type, -gene_model_type)

gene_2_TE_ids <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_description_file.csv.gz") %>%
  filter(!is.na(Derives_from)) %>%
  select(gene_id, Derives_from)

TEG_2_TE = merge.data.frame(TEG_annotations, gene_2_TE_ids, by = "gene_id", all.x = T) %>%
  dplyr::relocate(gene_id, .before = Derives_from)

TEG_n_TE_ids = TEG_2_TE %>% filter(!is.na(Derives_from)) %>% select(gene_id, Derives_from)

if (with_manual_Derives_from) {
  unknown_TEG = TEG_2_TE %>%
    filter(is.na(Derives_from)) %>%
    dplyr::select(-Derives_from) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)
}
#################################################

############# find 'Derives_from' to unknown_TEG

### upload and edit TE file
TE_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt",
  sep = "\t"
) %>%
  mutate(seqnames = NA) %>% # Add a new column with NA values
  dplyr::select(seqnames, Transposon_min_Start, Transposon_max_End, orientation_is_5prime, everything()) %>%
  dplyr::rename(Derives_from = Transposon_Name)# %>% filter(grepl("Gypsy|Copia|LINE", Transposon_Super_Family))

for (i in 1:5) {
  TE_file$seqnames[grep(paste0("AT", i, "TE"), TE_file$Derives_from)] = paste0("Chr", i)
}
TE_file$orientation_is_5prime = gsub("true","+",TE_file$orientation_is_5prime)
TE_file$orientation_is_5prime = gsub("false","-",TE_file$orientation_is_5prime)
names(TE_file)[1:4] = c("seqnames", "start", "end", "strand")

if (with_manual_Derives_from) {
  TE_gr = makeGRangesFromDataFrame(TE_file, keep.extra.columns = T)

  ### overlap 'TE_file' and 'ann_file' GRanges objects by position
  m <- findOverlaps(unknown_TEG, TE_gr)
  TEG <- unknown_TEG[queryHits(m)]
  mcols(TEG) <- c(mcols(TEG), mcols(TE_gr[subjectHits(m)]))

  unknown_TEG_df <- as.data.frame(TEG) %>%
    # add column that contain all the values from each row
    mutate(tmp = paste(gene_id, Transposon_Super_Family, sep = "_")) %>%
    distinct(tmp, .keep_all = T) %>%
    dplyr::select(-tmp) %>%
    .[, -(1:5)] %>%
    # filter(grepl("Gypsy|Copia|LINE", Transposon_Super_Family)) %>%
    distinct(gene_id, .keep_all = T) %>%
    dplyr::select(-Transposon_Family)
}
#################################################
### merge with RNAseq results
RNA_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
  # filter(gene_model_type == "transposable_element_gene") %>%
  dplyr::rename(gene_id = "locus_tag") %>%
  .[, 1:4] %>%
  filter(!is.na(padj))

## bind and merge for final deTEG_n_TE_df
TEG_n_TE_ids_with_SP = merge(TEG_n_TE_ids, TE_file[, c("Derives_from", "Transposon_Super_Family")], by = "Derives_from") %>% dplyr::relocate(gene_id, .before = Derives_from)

if (with_manual_Derives_from) {
  deTEG_n_TE_df = rbind(TEG_n_TE_ids_with_SP, unknown_TEG_df) %>%
    merge(., RNA_file, by = "gene_id")
} else {
  deTEG_n_TE_df = TEG_n_TE_ids_with_SP %>%
    merge(., RNA_file, by = "gene_id")
}


#
#TEG_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_annotations.csv.gz") %>%
#  filter(type == "transposable_element_gene") %>%
#  merge.data.frame(., RNA_file, by = "gene_id") %>%
#  dplyr::select(-type, -gene_model_type) %>%
#  makeGRangesFromDataFrame(., keep.extra.columns = T)

## overlap 'TE_file' and 'ann_file' GRanges objects by position
#m = findOverlaps(TEG_file, TE_gr)
#TEG = TEG_file[queryHits(m)]
#mcols(TEG) = c(mcols(TEG), mcols(TE_gr[subjectHits(m)]))
#
#TEG_df = as.data.frame(TEG) %>%
#  # add column that contain all the values from each row
#  mutate(tmp = do.call(paste, c(., sep = "_"))) %>%
#  distinct(tmp, .keep_all = T) %>%
#  dplyr::select(-tmp) %>%
#  .[,-(1:5)] # %>% filter(padj < 0.05)


######################################
grouped_families <- function(superFamily=NULL, Family=NULL, SP_df=F) { # 'SP_df' for results data frame of super-family
  ###
  # for super-family results
  if (SP_df) {
    x = deTEG_n_TE_df %>% filter(Transposon_Super_Family == superFamily) %>% dplyr::select(-pValue) %>% arrange(padj)
  }
  
  # for family results
  if (!is.null(Family)) {
    x = deTEG_n_TE_df %>% filter(Transposon_Family == Family) %>% dplyr::select(-pValue) %>% arrange(padj)
  }
  ###
  
  ###
  # for 'annotate and significant' data frame
  if (!is.null(superFamily) & !SP_df) {
    xx = deTEG_n_TE_df %>% filter(Transposon_Super_Family == superFamily) # %>% distinct(gene_id, .keep_all = T)
    
    final_df = data.frame(families = unique(xx$Transposon_Family), annotate = NA, significant = NA)
    
    for (i.fam in 1:length(final_df$families)) {
      final_df$annotate[i.fam] = xx %>% filter(Transposon_Family == final_df$families[i.fam]) %>% nrow()
      final_df$significant[i.fam] = xx %>% filter(Transposon_Family == final_df$families[i.fam]) %>% filter(padj < 0.05) %>% nrow()
    }
    
    x = final_df %>% arrange(-annotate)
  }
  ###
  return(x)
}

######################################

#Copia = grouped_families(superFamily = "LTR/Copia")
#Gypsy = grouped_families(superFamily = "LTR/Gypsy")
#LINE = grouped_families(superFamily = "LINE/L1")

######################################

# df for volcano plot
retro_TE = rbind(
  grouped_families(superFamily = "LTR/Copia", SP_df = T),
  grouped_families(superFamily = "LTR/Gypsy", SP_df = T),
  grouped_families(superFamily = "LINE/L1", SP_df = T)
)
row.names(retro_TE) = 1:nrow(retro_TE)

# volcano plot 
mydf <- retro_TE

# Factor 'geneCat' with levels in the desired order
mydf$geneCat <- with(mydf, ifelse(padj < 0.05 & log2FoldChange > 1, "Upregulated",
                                  ifelse(padj < 0.05 & log2FoldChange < -1, "Downregulated", "nonDE")))
mydf$geneCat <- factor(mydf$geneCat, levels = c("Upregulated", "Downregulated", "nonDE"))

svg("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/red_blue_retro-TEGs_RNAseq_volcano.svg",
    width = 3.8, height = 2, family = "serif")

ggplot(mydf, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = geneCat), alpha = 0.75, size = 1) +
    xlab("log2(Fold-Change)") +
    ylab("-log10(padj)") +
    theme_bw() +
    scale_colour_manual(
      #name = "Gene Expression",
      name = "",
      values = c("nonDE" = "gray60", "Upregulated" = "#a84848", "Downregulated" = "#5d60ba")
    ) +
    guides(color = guide_legend(override.aes = list(size = 2.5))) +

  # Add a vertical line starting from y=5 at x=1
  geom_segment(aes(x = rep(1, nrow(mydf)), y = -log10(0.05), xend = 1, yend = Inf), 
               col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed") +
  # Add a horizontal line ending at x=-1 from the left
  geom_segment(aes(x = rep(-1, nrow(mydf)), y = -log10(0.05), xend = -1, yend = Inf), 
               col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed") +
  # Add a horizontal line starting from x=1 to the right
  geom_segment(aes(x = rep(1, nrow(mydf)), y = -log10(0.05), xend = Inf, yend = -log10(0.05)), 
               col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed") +
  geom_segment(aes(x = rep(-Inf, nrow(mydf)), y = -log10(0.05), xend = -1, yend = -log10(0.05)), 
               col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed")

dev.off()

####################################
