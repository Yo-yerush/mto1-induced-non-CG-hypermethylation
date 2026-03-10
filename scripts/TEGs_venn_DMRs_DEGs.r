library(dplyr)
library(VennDiagram)

DMR_file = rbind(
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/TEG_CG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/TEG_CHH_genom_annotations.csv")
) %>% 
distinct(gene_id)


RNA_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(gene_model_type == "transposable_element_gene") %>%
    dplyr::rename(gene_id = "locus_tag") %>%
    filter(padj < 0.05) %>%
    distinct(gene_id)


gene_sets <- list(
    DMRs = DMR_file$gene_id,
    DEGs = RNA_file$gene_id
)

#fill = c("#440154ff", '#21908dff', '#fde725ff')
venn_colors = c("yellow", "purple", "green","red", "blue")
category.position = c(0, 0)
resolution = 300
cex = 0.75

venn.diagram(
    x = gene_sets,
    category.names = c("\nunique-TEGs                    \noverlap DMRs                   ",
                       "\n\n\n\n\n\n\n\nDE-TEGs"),
    filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_deTEGs_VennDiagram.png"),
    disable.logging = T,
    output = T,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = resolution,
    compression = "lzw",
    lwd = 1,
    fill = venn_colors[1:length(gene_sets)],
    alpha = rep(0.45, length(gene_sets)),
    col = rep("white", length(gene_sets)),
    cex = cex,
    fontfamily = "serif",
    cat.cex = cex,
    cat.default.pos = "outer",
    cat.pos = category.position,
    cat.fontface = 2,
    cat.fontfamily = "serif"
    #    cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
    #    col=venn_colors,
    #    rotation = 1
)

overlap_size <- length(intersect(DMR_file$gene_id, RNA_file$gene_id))
cat("Number of overlapping genes:", overlap_size, "\n")
