library(dplyr)
library(GenomicRanges)
library(ggplot2)

context = "all"
#context = "nonCG"
#context = "CG"
#context = "CHG"
#context = "CHH"

#################################################

TE_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt",
  sep = "\t"
) %>%
  mutate(seqnames = NA) %>% # Add a new column with NA values
  dplyr::select(seqnames, Transposon_min_Start, Transposon_max_End, orientation_is_5prime, everything())
# dplyr::rename(gene_id = Transposon_Name)

for (i in 1:5) {
  TE_file$seqnames[grep(paste0("AT", i, "TE"), TE_file$Transposon_Name)] <- paste0("Chr", i)
}
TE_file$orientation_is_5prime <- gsub("true", "+", TE_file$orientation_is_5prime)
TE_file$orientation_is_5prime <- gsub("false", "-", TE_file$orientation_is_5prime)

names(TE_file)[1:4] <- c("seqnames", "start", "end", "strand")
TE_gr <- makeGRangesFromDataFrame(TE_file, keep.extra.columns = T)

#################################################

  if (context == "all" | context == "CG") {
     x_pos = 0.6
  } else if (context == "nonCG" | context == "CHG") {
     x_pos = 0.4
  } else if (context == "CHH") {
     x_pos = 0.3
  }

  #################################################

  if (context == "all") {
    DMR_file_0 <- rbind(
      read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/TEG_CG_genom_annotations.csv"),
      read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv"),
      read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/TEG_CHH_genom_annotations.csv")
    )
  } else if (context == "nonCG") {
    DMR_file_0 <- rbind(
      read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv"),
      read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/TEG_CHH_genom_annotations.csv")
    )
  } else {
    DMR_file_0 <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/TEG_", context, "_genom_annotations.csv"))
  }

  DMR_file <- DMR_file_0 %>%
    dplyr::rename(log2FoldChange = log2FC) %>%
    .[, c("gene_id", "log2FoldChange", "pValue")]

  TEG_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/RA_costume_annotations_files/Methylome.At_annotations.csv.gz") %>%
    filter(type == "transposable_element_gene") %>%
    merge.data.frame(., DMR_file, by = "gene_id") %>%
    dplyr::select(-type, -gene_model_type) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)

  ## overlap 'TE_file' and 'ann_file' GRanges objects by position
  m <- findOverlaps(TEG_file, TE_gr)
  TEG <- TEG_file[queryHits(m)]
  mcols(TEG) <- c(mcols(TEG), mcols(TE_gr[subjectHits(m)]))

  TEG_df <- as.data.frame(TEG) %>%
    # add column that contain all the values from each row
    mutate(tmp = do.call(paste, c(., sep = "_"))) %>%
    distinct(tmp, .keep_all = T) %>%
    dplyr::select(-tmp) %>%
    .[, -(1:5)] # %>% filter(padj < 0.05)


  ######################################
  grouped_families <- function(superFamily = NULL, Family = NULL, SP_df = F) { # 'SP_df' for results data frame of super-family
    ###
    # for super-family results
    if (SP_df) {
      x <- TEG_df %>%
        filter(Transposon_Super_Family == superFamily) %>%
        arrange(pValue)
    }

    # for family results
    if (!is.null(Family)) {
      x <- TEG_df %>%
        filter(Transposon_Family == Family) %>%
        arrange(pValue)
    }
    ###

    ###
    # for 'annotate and significant' data frame
    if (!is.null(superFamily) & !SP_df) {
      xx <- TEG_df %>% filter(Transposon_Super_Family == superFamily) # %>% distinct(gene_id, .keep_all = T)

      final_df <- data.frame(families = unique(xx$Transposon_Family), annotate = NA, significant = NA)

      for (i.fam in 1:length(final_df$families)) {
        final_df$annotate[i.fam] <- xx %>%
          filter(Transposon_Family == final_df$families[i.fam]) %>%
          nrow()
        final_df$significant[i.fam] <- xx %>%
          filter(Transposon_Family == final_df$families[i.fam]) %>%
          filter(pValue < 0.05) %>%
          nrow()
      }

      x <- final_df %>% arrange(-annotate)
    }
    ###
    return(x)
  }

  ######################################

  # df for volcano plot
  retro_TE <- rbind(
    grouped_families(superFamily = "LTR/Copia", SP_df = T),
    grouped_families(superFamily = "LTR/Gypsy", SP_df = T),
    grouped_families(superFamily = "LINE/L1", SP_df = T)
  )
  row.names(retro_TE) <- 1:nrow(retro_TE)

  # volcano plot
  mydf <- retro_TE

  ## Factor 'geneCat' with levels in the desired order
  #mydf$geneCat <- with(mydf, ifelse(pValue < 0.05 & log2FoldChange > 1, "Upregulated",
  #  ifelse(pValue < 0.05 & log2FoldChange < -1, "Downregulated", "nonDE")
  #))
  #mydf$geneCat <- factor(mydf$geneCat, levels = c("Upregulated", "Downregulated", "nonDE"))

  vplot <- ggplot(mydf, aes_string(x = "log2FoldChange", y = "-log10(pValue)", color = "Transposon_Super_Family")) +
    geom_point(alpha = 0.4, size = 0.75) +
    # change axis titles
    xlab("log2(Fold-Change)") +
    ylab("-log10(padj)") +
    # theme_classic() + #
    theme_bw() +
    theme(legend.position = "none") +
    scale_colour_manual(
      name = "",
      values = c("LTR/Gypsy" = "#842dcc", "LTR/Copia" = "#159e35", "LINE/L1" = "#c79924"),
      limits = c("LTR/Gypsy", "LTR/Copia", "LINE/L1")
    ) +
    scale_y_continuous(breaks = c(0, 50, 150, 250)) +
    # guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 0.65))) +

    # Add a vertical line starting from y=5 at x=1
    geom_segment(aes(x = rep(x_pos, nrow(mydf)), y = -log10(0.05), xend = x_pos, yend = Inf),
      col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed"
    ) +
    # Add a horizontal line ending at x=-1 from the left
    geom_segment(aes(x = rep(-x_pos, nrow(mydf)), y = -log10(0.05), xend = -x_pos, yend = Inf),
      col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed"
    ) +
    # Add a horizontal line starting from x=1 to the right
    geom_segment(aes(x = rep(x_pos, nrow(mydf)), y = -log10(0.05), xend = Inf, yend = -log10(0.05)),
      col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed"
    ) +
    geom_segment(aes(x = rep(-Inf, nrow(mydf)), y = -log10(0.05), xend = -x_pos, yend = -log10(0.05)),
      col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed"
    )

  tiff(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_on_retro-TEGs_volcano_",context,".tif"), width = 1.75, height = 1.15, units = "in", res = 600, family = "serif"
  )
  (vplot)
  dev.off()

  #######################################

  # vplot_with_legend <- vplot + theme(legend.position = "right")
  # legend_plot <- cowplot::get_legend(vplot_with_legend)
  # tiff("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/retro-TEGs_volcano_legend.tif",
  #   width = 1.5, height = 1, units = "in", res = 600, family = "serif"
  # )
  # ggdraw(legend_plot)
  # dev.off()
