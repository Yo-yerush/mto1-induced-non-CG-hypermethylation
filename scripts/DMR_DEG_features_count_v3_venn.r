library(dplyr)
library(VennDiagram)

####################################
# if want FALSE, go to 'v2' script #
#                                  #
# unique_tairs_in_DMRs <- TRUE     #
#                                  #
####################################

##################### DMRs IDs
read_DMRs <- function(f) {
    x <- NULL
    for (c in c("CG", "CHG", "CHH")) {
        x <- rbind(x, read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", c, "/", f, "_", c, "_genom_annotations.csv")))
    }
    return(select(x, gene_id)) # all DMRs
}

Promoters <- read_DMRs("Promoters")
Genes <- read_DMRs("Genes")
CDS <- read_DMRs("CDS")
Introns <- read_DMRs("Introns")
fiveUTRs <- read_DMRs("fiveUTRs")
threeUTRs <- read_DMRs("threeUTRs")

##################### DEGs IDs
DEGs <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all_genes_results_mto1_vs_wt.csv") %>%
    filter(padj < 0.05) %>%
    distinct(gene_id)

##################### genes sets
gene_sets_1 <- list(
    Promoter = Promoters$gene_id,
    Genes = Genes$gene_id,
    DEGs = DEGs$gene_id
)

gene_sets_CDS <- list(
    CDS = CDS$gene_id,
    DEGs = DEGs$gene_id
)

gene_sets_Introns <- list(
    Introns = Introns$gene_id,
    DEGs = DEGs$gene_id
)

gene_sets_fiveUTRs <- list(
    fiveUTRs = fiveUTRs$gene_id,
    DEGs = DEGs$gene_id
)

gene_sets_threeUTRs <- list(
    threeUTRs = threeUTRs$gene_id,
    DEGs = DEGs$gene_id
)

##################### venn plots
venn_colors <- c("#928e92", "#928e92", "#d69641")
resolution <- 900
cex_gb_p <- 0.65
cex <- 0.5

### Genes n promoters
venn.diagram(
    x = gene_sets_1,
    category.names = c("          Promoters", "Genes-body            ", "DEGs"),
    filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_DEGs_Gb_Pr_unique_VennDiagram.tiff"),
    disable.logging = T,
    output = T,
    imagetype = "tiff",
    height = 1.65,
    width = 1.65,
    units = "in",
    resolution = resolution,
    lwd = 1,
    fill = venn_colors[1:length(gene_sets_1)],
    alpha = rep(0.45, length(gene_sets_1)),
    col = rep("white", length(gene_sets_1)),
    cex = cex_gb_p,
    fontfamily = "serif",
    cat.cex = cex_gb_p,
    cat.default.pos = "outer",
    cat.fontface = 2,
    cat.fontfamily = "serif"
)

### CDS
ann_name <- "CDS"
venn.diagram(
    x = gene_sets_CDS,
    category.names = c(ann_name, "DEGs"),
    filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_DEGs_", ann_name, "_unique_VennDiagram.svg"),
    disable.logging = T,
    output = T,
    imagetype = "svg",
    height = 0.8,
    width = 0.8,
    resolution = resolution,
    lwd = 1,
    fill = venn_colors[2:3],
    alpha = rep(0.45, 2),
    col = rep("white", 2),
    cex = cex,
    fontfamily = "serif",
    cat.cex = cex,
    cat.default.pos = "outer",
    cat.pos = c(180, 180),
    cat.fontface = 2,
    cat.fontfamily = "serif",
    ext.pos = 0, ,
    ext.percent = rep(0.1, 3),
    ext.dist = -0.45,
    ext.length = 0.85
)

### Introns
ann_name <- "Introns"
venn.diagram(
    x = gene_sets_Introns,
    category.names = c(ann_name, "DEGs"),
    filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_DEGs_", ann_name, "_unique_VennDiagram.svg"),
    disable.logging = T,
    output = T,
    imagetype = "svg",
    height = 0.8,
    width = 0.8,
    resolution = resolution,
    lwd = 1,
    fill = venn_colors[2:3],
    alpha = rep(0.45, 2),
    col = rep("white", 2),
    cex = cex,
    fontfamily = "serif",
    cat.cex = cex,
    cat.default.pos = "outer",
    cat.pos = c(180, 180),
    cat.fontface = 2, ,
    cat.fontfamily = "serif",
    ext.pos = 0,
    ext.dist = -0.45,
    ext.length = 0.85
)

### fiveUTRs
ann_name <- "fiveUTRs"
venn.diagram(
    x = gene_sets_fiveUTRs,
    category.names = c("5'UTRs          ", "DEGs"),
    filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_DEGs_", ann_name, "_unique_VennDiagram.svg"),
    disable.logging = T,
    output = T,
    imagetype = "svg",
    height = 0.8,
    width = 0.8,
    resolution = resolution,
    lwd = 1,
    fill = venn_colors[2:3],
    alpha = rep(0.45, 2),
    col = rep("white", 2),
    cex = cex,
    fontfamily = "serif",
    cat.cex = cex,
    cat.default.pos = "outer",
    cat.pos = c(180, 180),
    cat.fontface = 2,
    cat.fontfamily = "serif",
    ext.pos = 0,
    ext.dist = c(-0.3, 0.15),
    ext.length = 0.85,
    ext.line.lwd = 0.75,
    ext.line.lty = 1
)

### threeUTRs
ann_name <- "threeUTRs"
venn.diagram(
    x = gene_sets_threeUTRs,
    category.names = c("3'UTRs        ", "DEGs"),
    filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_DEGs_", ann_name, "_unique_VennDiagram.svg"),
    disable.logging = T,
    output = T,
    imagetype = "svg",
    height = 0.8,
    width = 0.8,
    resolution = resolution,
    lwd = 1,
    fill = venn_colors[2:3],
    alpha = rep(0.45, 2),
    col = rep("white", 2),
    cex = cex,
    fontfamily = "serif",
    cat.cex = cex,
    cat.default.pos = "outer",
    cat.pos = c(180, 180),
    cat.fontface = 2,
    cat.fontfamily = "serif",
    ext.pos = 0,
    ext.dist = c(-0.3, 0.15),
    ext.length = 0.85,
    ext.line.lwd = 0.75,
    ext.line.lty = 1
)




############################################################
############################################################
############################################################
############################################################
############################################################


# venn_Gb <- function(sets, ann_name) {
#    if (ann_name == "fiveUTRs") {
#        ann_name_2 <- "5'UTR"
#    } else if (ann_name == "threeUTRs") {
#        ann_name_2 <- "3'UTR"
#    } else {
#        ann_name_2 <- ann_name
#    }
#    venn.diagram(
#        x = sets,
#        category.names = c(ann_name_2, "DEGs"),
#        filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_DEGs_", ann_name, "_unique_VennDiagram.png"),
#        disable.logging = T,
#        output = T,
#        imagetype = "png",
#        height = 300,
#        width = 300,
#        resolution = resolution,
#        lwd = 1,
#        fill = venn_colors[2:3],
#        alpha = rep(0.45, length(sets)),
#        col = rep("white", length(sets)),
#        cex = cex,
#        fontfamily = "serif",
#        cat.cex = cex,
#        cat.default.pos = "outer",
#        cat.pos = c(180, 180),
#        cat.fontface = 2,
#        cat.fontfamily = "serif"
#    )
# }
#
# venn_Gb(gene_sets_CDS, "CDS")
# venn_Gb(gene_sets_Introns, "Introns")
# venn_Gb(gene_sets_fiveUTRs, "fiveUTRs")
# venn_Gb(gene_sets_threeUTRs, "threeUTRs")

