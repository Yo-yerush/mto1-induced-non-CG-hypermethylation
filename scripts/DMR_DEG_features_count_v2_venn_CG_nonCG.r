library(dplyr)
library(VennDiagram)

for (context in c("CG", c("CHG", "CHH"))) {

    cntx_name = ifelse(context == "CG", "CG", "non-CG")

    ##################### DMRs IDs
    read_DMRs <- function(context, f) {
        x <- NULL
        for (c in context) {
            x <- rbind(x, read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", c, "/", f, "_", c, "_genom_annotations.csv")))
        }
        return(select(x, gene_id)) # all DMRs
    }

    Promoters <- read_DMRs(context, "Promoters")
    Genes <- read_DMRs(context, "Genes")

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

    venn.diagram(
        x = gene_sets,
        category.names = c("          Promoters", "Genes-body            ", "DEGs"),
        filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_DEGs_", cntx_name, "_Gb_Pr_unique_VennDiagram.png"),
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
}
