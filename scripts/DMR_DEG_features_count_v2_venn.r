library(dplyr)
library(VennDiagram)

unique_tairs_in_DMRs <- TRUE

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

##################### DEGs IDs
DEGs <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(padj < 0.05) %>%
    rename(gene_id = locus_tag) %>%
    distinct(gene_id)


gene_sets <- list(
    Promoter = Promoters$gene_id,
    Genes = Genes$gene_id,
    DEGs = DEGs$gene_id
)

#venn_colors <- c("#a05b9c", "#a05b9c", "#71c071")
venn_colors <- c("#928e92", "#928e92", "#d69641")
resolution <- 300
cex <- 0.65

if (unique_tairs_in_DMRs) {
    venn.diagram(
        x = gene_sets,
        category.names = c("          Promoters", "Genes-body            ", "DEGs"),
        filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_DEGs_Gb_Pr_unique_VennDiagram.png"),
        disable.logging = T,
        output = T,
        imagetype = "png",
        height = 480,
        width = 480,
        resolution = resolution,
        lwd = 1,
        fill = venn_colors[1:length(gene_sets)],
        alpha = rep(0.45, length(gene_sets)),
        col = rep("white", length(gene_sets)),
        cex = cex,
        fontfamily = "serif",
        cat.cex = cex,
        cat.default.pos = "outer",
        cat.fontface = 2,
        cat.fontfamily = "serif"
    )
} else {
    # Compute total counts (including duplicates)
    n_promoters <- nrow(Promoters)
    n_genes <- nrow(Genes)
    n_degs <- nrow(DEGs)

    # Create frequency tables for each set based on gene_id
    promoters_counts <- table(Promoters$gene_id)
    genes_counts <- table(Genes$gene_id)
    degs_counts <- table(DEGs$gene_id)

    # Intersection between Promoters and DEGs
    common_promoters_degs <- intersect(names(promoters_counts), names(degs_counts))
    n_intersect_promoters_degs <- sum(sapply(common_promoters_degs, function(g) {
        min(promoters_counts[g], degs_counts[g])
    }))

    # Intersection between Promoters and Genes
    common_promoters_genes <- intersect(names(promoters_counts), names(genes_counts))
    n_intersect_promoters_genes <- sum(sapply(common_promoters_genes, function(g) {
        min(promoters_counts[g], genes_counts[g])
    }))

    # Intersection between Genes and DEGs
    common_genes_degs <- intersect(names(genes_counts), names(degs_counts))
    n_intersect_genes_degs <- sum(sapply(common_genes_degs, function(g) {
        min(genes_counts[g], degs_counts[g])
    }))

    # Triple intersection among Promoters, Genes, and DEGs
    common_all <- Reduce(intersect, list(names(promoters_counts), names(genes_counts), names(degs_counts)))
    n_intersect_all <- sum(sapply(common_all, function(g) {
        min(promoters_counts[g], genes_counts[g], degs_counts[g])
    }))

    png("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_DEGs_Gb_Pr_non.unique_VennDiagram.png", width = 480, height = 480, res = 300)
    venn.plot <- draw.triple.venn(
        area1 = n_promoters,
        area2 = n_genes,
        area3 = n_degs,
        n12 = n_intersect_promoters_genes, # Overlap between Promoters and Genes
        n23 = n_intersect_genes_degs, # Overlap between Genes and DEGs
        n13 = n_intersect_promoters_degs, # Overlap between Promoters and DEGs
        n123 = n_intersect_all, # Triple intersection
        category = c("          Promoters", "Genes-body            ", "DEGs"),
        imagetype = "png",
        height = 480,
        width = 480,
        resolution = resolution,
        lwd = 1,
        fill = venn_colors[1:length(gene_sets)],
        alpha = rep(0.45, length(gene_sets)),
        col = rep("white", length(gene_sets)),
        cex = cex,
        fontfamily = "serif",
        cat.cex = cex,
        cat.default.pos = "outer",
        # cat.pos = category.position,
        cat.fontface = 2,
        cat.fontfamily = "serif"
    )
    dev.off()
}
