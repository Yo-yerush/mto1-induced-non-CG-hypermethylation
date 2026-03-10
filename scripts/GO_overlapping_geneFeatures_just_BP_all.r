library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(GenomicRanges)
library(topGO)

########################
GO_pValue = 0.01
########################

RNAseq_0 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
    dplyr::select(gene_id, DEG_log2FC, padj, pValue)

########################

for (beckground_list in c("DEGs", "DMRs", "all_genes")) {
    
    for (direction_i in c("upregulated", "downregulated", "both_sides")) {

        # output directory
        output_path <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/"
        dir.create(paste0(output_path, "GO_overlap"), showWarnings = F)

        if (beckground_list == "DEGs") {
            output_dir <- "GO_overlap/DEGs_beckground"
        } else if (beckground_list == "DMRs") {
            output_dir <- "GO_overlap/DMRs_beckground"
        } else {
            output_dir <- "GO_overlap"
        }

        dir.create(paste0(output_path, output_dir), showWarnings = F)

        # gene list and color
        if (direction_i == "upregulated") {
            RNAseq <- RNAseq_0 %>% filter(DEG_log2FC > 0)
            sig_color <- "#cf534c"
        } else if (direction_i == "downregulated") {
            RNAseq <- RNAseq_0 %>% filter(DEG_log2FC < 0)
            sig_color <- "#6397eb"
        } else if (direction_i == "both_sides") {
            RNAseq <- RNAseq_0
            sig_color <- "#F2A672"
        }


        overlapped_context <- function(context) {
            gene_feature <- rbind(
                read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Promoters_", context, "_genom_annotations.csv")),
                read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Genes_", context, "_genom_annotations.csv"))
            ) %>%
                dplyr::select(gene_id)

            return(gene_feature)
        }

        ################

        all_DMRs_bind <- rbind(
            overlapped_context("CG"),
            overlapped_context("CHG"),
            overlapped_context("CHH")
        ) %>%
            distinct(gene_id)

        overlapped_sig_genes <- all_DMRs_bind %>%
            merge.data.frame(RNAseq, ., by = "gene_id") %>%
            filter(padj < 0.05) %>%
            mutate(for_geneList = 1) %>%
            dplyr::select(gene_id, for_geneList)


        # differents beckground gene list
        if (beckground_list == "DEGs") {
            overlapped <- merge(overlapped_sig_genes, filter(RNAseq_0, padj < 0.05), by = "gene_id", all.y = T)
        } else if (beckground_list == "DMRs") {
            overlapped <- merge(overlapped_sig_genes, all_DMRs_bind, by = "gene_id", all.y = T)
        } else {
            overlapped <- merge(overlapped_sig_genes, RNAseq_0, by = "gene_id", all.y = T)
        }

        overlapped$for_geneList[is.na(overlapped$for_geneList)] <- 0

        geneList <- overlapped$for_geneList
        names(geneList) <- overlapped$gene_id


        res_list <- list()
        for (GO_type_loop in c("BP")) { # , "MF", "CC")) {
            myGOdata <- new("topGOdata",
                ontology = GO_type_loop,
                allGenes = geneList,
                geneSelectionFun = function(x) (x == 1),
                # description = "Test",
                annot = annFUN.org,
                # nodeSize = 5,
                mapping = "org.At.tair.db"
            )

            sg <- sigGenes(myGOdata)
            str(sg)
            numSigGenes(myGOdata)

            resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")


            allRes <- GenTable(myGOdata,
                Fisher = resultFisher,
                orderBy = "Fisher", ranksOf = "Fisher", topNodes = length(resultFisher@score)
            )
            allRes$Fisher <- as.numeric(allRes$Fisher)
            allRes$Term <- gsub(",", ";", allRes$Term)
            allRes$type <- GO_type_loop

            res_list[[GO_type_loop]] <- allRes[allRes$Fisher <= GO_pValue, ]
        }

        ### plot
        # just BP !!!!!!!!!
        res_bind <- rbind(res_list[["BP"]]) # , res_list[["CC"]], res_list[["MF"]])
        res_bind$type <- "Biological Process"

        res_bind <- res_bind[!grepl("cellular_component|biological_process|molecular_function|macromolecule biosynthetic process", res_bind$Term), ]

        write.csv(res_bind, paste0(output_path, output_dir, "/GO_all_", direction_i, "_DEGs_overlap.csv"), row.names = F)

        bubble_plot <- res_bind %>%
            ggplot(aes(Significant, reorder(Term, Significant), size = Annotated, color = Fisher)) +
            scale_color_gradient("p.value", low = sig_color, high = "black") + # theme_classic() +
            labs( # title = paste0(type_name," - ",gain_loss,"regulated transcripts"),
                x = "Significant", y = ""
            ) +
            theme_bw() +
            theme(
                # plot.title=element_text(hjust=0.5),
                # legend.key.size = unit(0.25, "cm"),
                # legend.title = element_text(size = 9.5),
                # legend.position = "right",
                # text = element_text(family = "serif")
                legend.position = "none"
            ) +
            geom_point() +
            facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
            guides(color = guide_colorbar(order = 1, barheight = 4))

        # print(paste(treatment,context,annotation, sep = "_"))

        svg(paste0(output_path, output_dir, "/GO_all_", direction_i, "_DEGs_overlap.svg"), width = 4.25, height = 3, family = "serif")
        print(bubble_plot)
        dev.off()

        svg(paste0(output_path, output_dir, "/GO_all_", direction_i, "_DEGs_overlap_Annotated_legend.svg"), width = 1, height = 1.5, family = "serif")
        legend_plot <- ggplot(res_bind, aes(x = 1, y = 1, size = Annotated)) +
            geom_point() +
            scale_size(
                # range = c(1, 6),
                name = " ",
                breaks = c(min(res_bind$Annotated), 50, 100, max(res_bind$Annotated))
            ) +
            theme_void() +
            theme(legend.position = "right")

        legend_only <- cowplot::get_legend(legend_plot)
        grid.newpage()
        grid.draw(legend_only)
        dev.off()
    }
}
