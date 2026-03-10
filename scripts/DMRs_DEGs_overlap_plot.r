library(dplyr)

DEGs_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(padj < 0.05) %>%
    rename(gene_id = locus_tag)

DEGs_up <- DEGs_file %>%
    filter(log2FoldChange > 0) %>%
    distinct(gene_id)
DEGs_down <- DEGs_file %>%
    filter(log2FoldChange < 0) %>%
    distinct(gene_id)


#########################################
load_csv <- function(cntx, ann, direction) {
    x <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", cntx, "/", ann, "_", cntx, "_genom_annotations.csv"))

    if (direction == "up") {
        y <- DEGs_up
    }
    if (direction == "down") {
        y <- DEGs_down
    }

    x_hyper <- x %>%
        filter(regionType == "gain") %>%
        select(gene_id) %>%
        merge(y)
    x_hypo <- x %>%
        filter(regionType == "loss") %>%
        select(gene_id) %>%
        merge(y)

    return(list(
        n_hyper = nrow(x_hyper),
        n_hypo = nrow(x_hypo),
        n_unique = rbind(x_hyper, x_hypo) %>% distinct(gene_id) %>% nrow()
    ))
}

# CG_genes = load_csv("CG", "Genes", "up")
# CG_promo = load_csv("CG", "Promoters", "up")
# CHG_genes = load_csv("CHG/Genes_CHG_genom_annotations.csv")
# CHG_promo = load_csv("CHG/Promoters_CHG_genom_annotations.csv")
# CHH_genes = load_csv("CHH/Genes_CHH_genom_annotations.csv")
# CHH_promo = load_csv("CHH/Promoters_CHH_genom_annotations.csv")

direction_col_up <- c("Hypo" = "#d96c6c", "Hyper" = "#d96c6c")
direction_col_down <- c("Hypo" = "#6c96d9", "Hyper" = "#6c96d9")

{
    svg("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_DEGs_overlap_plot.svg", width = 5.75, height = 4.25, family = "serif")

    par(mar = c(0.1, 3, 1, 0.1))

    # A 4x5 matrix: fill row by row from 1 to 20
    layout_mat <- matrix(1:20, nrow = 4, ncol = 5, byrow = TRUE)
    layout(
        mat = layout_mat,
        widths = c(0.3, 0.4, 1, 1, 1),
        heights = c(0.25, 1, 1, 0.75)
    )

    # first row (context lables)
    plot.new()
    text(0, 0, " ", cex = 1.5, font = 2)
    plot.new()
    text(0, 0, " ", cex = 1.5, font = 2)
    plot.new()
    text(-0.05, 0.35, "CG", cex = 2.5, font = 2, pos = 4)
    plot.new()
    text(-0.05, 0.35, "CHG", cex = 2.5, font = 2, pos = 4)
    plot.new()
    text(-0.05, 0.35, "CHH", cex = 2.5, font = 2, pos = 4)

    # second and third rows (plots)
    # enpty
    plot.new()
    text(0.5, 0.5, " ", cex = 1.5, font = 2, srt = 90)
    # Add a new plot with text "Upregulated"
    plot.new()
    text(0.5, 0.5, "Upregulated", cex = 1.5, font = 4, srt = 90)
    direction_colors = direction_col_up

    for (l.direction in c("up", "down")) {
        if (l.direction == "down") {
            direction_colors = direction_col_down

            plot.new()
            text(0.5, 0.5, " ", cex = 1.5, font = 2, srt = 90)
            # Add a new plot with text "Downregulated"
            plot.new()
            text(0.5, 0.5, "Downregulated", cex = 1.5, font = 4, srt = 90)
        }

        for (l.cntx in c("CG", "CHG", "CHH")) {
            l_promoters <- load_csv(l.cntx, "Promoters", l.direction)
            l_genes <- load_csv(l.cntx, "Genes", l.direction)

            plot_data <- data.frame(
                Promoters = c(l_promoters$n_hypo, l_promoters$n_hyper),
                Genes = c(l_genes$n_hypo, l_genes$n_hyper)
            )
            plot_data <- as.matrix(plot_data)
            row.names(plot_data) <- c("loss", "gain")

            barplot(
                height = plot_data,
                beside = FALSE,
                col = direction_colors,
                border = "black",
                space = 0.4,
                ylab = "",
                # main = paste0(context, " Context"),
                names.arg = c("", ""),
                # cex.names = 1.15,
                las = 2,
                cex.axis = 1.5,
                density = c(25, NA),
                angle = c(45, NA)
            )
        }

        # # Add a new plot with text "Up-/Down- regulated"
        # if (l.direction == "up") {
        #     plot.new()
        #     text(0.5, 0.5, "Upregulated", cex = 1.5, font = 2, srt = 90, pos = 2)
        # } else if (l.direction == "down") {
        #     plot.new()
        #     text(0.5, 0.5, "Downregulated", cex = 1.5, font = 2, srt = 90, pos = 2)
        # }
    }

    # forth row (labels)
    plot.new()
    text(0.5, 0.5, " ", cex = 1.5, font = 2)
    plot.new()
    text(0.5, 0.5, " ", cex = 1.5, font = 2)
    i <- 1
    while (i <= 3) {
        plot.new()
        text(0.2, 0.75, "Promoters", cex = 1.5, font = 2, srt = 45, pos = 1)
        text(0.7, 0.75, "Genes", cex = 1.5, font = 2, srt = 45, pos = 1)

        i <- i + 1
    }


    # y axis shared title
    par(fig = c(0, 1, 0, 1), new = TRUE, xpd = NA)
    plot.new()
    text(
        x = 0.01,
        y = 0.575,
        labels = "DMRs Count",
        cex = 2.5,
        font = 2,
        srt = 90
    )

    dev.off()
}

# legend
{
    svg("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_DEGs_overlap_legend.svg", width = 3, height = 3, family = "serif")
    plot.new()
    par(mar = c(0, 0, 0, 0))
    legend(
        "topleft",
        legend = rep("", 4),
        fill = c("#d96c6c", "#6c96d9", "#d96c6d", "#6c96d9"),
        density = c(15, 15, NA, NA),
        angle = c(45, 45, 0, 0),
        ncol = 2,
        bty = "n",
        x.intersp = 0,
        cex = 2
    )

    text(0.35, 0.825, "Upregulated", cex = 0.75, font = 1, srt = 0, pos = 4)
    text(0.35, 0.605, "Downregulated", cex = 0.75, font = 1, srt = 0, pos = 4)

    text(0.5, 1, "Hypo / Hyper", cex = 0.75, font = 1, srt = 0, pos = 2)
    # text(0.45, 0.4, "Downregulated", cex = 1.5, font = 1, srt = 0, pos = 4)

    dev.off()
}

# legend
#{
#    svg("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_DEGs_overlap_legend.svg", width = 5, height = 5, family = "serif")
#    plot.new()
#    par(mar = c(0, 0, 0, 0))
#    legend(
#        "left",
#        legend = rep("", 4),
#        fill = c("#d96c6c", "#6c96d9", "#d96c6d", "#6c96d9"),
#        density = c(15, 15, NA, NA),
#        angle = c(45, 45, 0, 0),
#        ncol = 2,
#        bty = "",
#        x.intersp = -0.25,
#        cex = 4,
#        box.lwd = 5
#    )
#
#    text(0.45, 0.5725, "Upregulated", cex = 1.5, font = 1, srt = 0, pos = 4)
#    text(0.45, 0.4, "Downregulated", cex = 1.5, font = 1, srt = 0, pos = 4)
#
#    text(0.05, 0.7, "Hypo / Hyper", cex = 1.5, font = 1, srt = 0, pos = 4)
#    # text(0.45, 0.4, "Downregulated", cex = 1.5, font = 1, srt = 0, pos = 4)
#
#    dev.off()
#}




# par(mgp = c(3, 1, 0))
#
## correlation legend
# text(
#    x = 0 + 0.4,
#    y = 600 * 0.9,
#    labels = "DMRs",
#    cex = 0.8,
#    pos = 4,
#    col = "black",
#    font = 2
# )
#
# legend(
#    x = 0 + 0.5,
#    y = 600 * 1,
#    legend = c("Hyper", "Hypo"),
#    fill = rev(direction_colors),
#    border = "black",
#    cex = 0.8,
#    bty = "n",
#    # title = expression(bold(" Correlation")),
#    title = "",
#    x.intersp = 0.5,
#    xjust = 0,
# )
#
# par(xpd = FALSE)