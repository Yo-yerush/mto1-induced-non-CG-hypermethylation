  library(ggplot2)
  
  genome_ann = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/"
  output_path = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/"

for (context in c("CG", "CHG", "CHH")) {
    DMRsReplicates_TE_file.0 = paste0(genome_ann, context, "/Transposable_Elements_", context, "_genom_annotations.csv")
    DMRsReplicates_Genes_file.0 = paste0(genome_ann, context, "/Genes_", context, "_genom_annotations.csv")

    # if (file.exists(DMRsReplicates_TE_file.0)) {

    DMRsReplicates_TE_file = read.csv(DMRsReplicates_TE_file.0)
    DMRsReplicates_Genes_file = read.csv(DMRsReplicates_Genes_file.0)
    DMRsReplicates_Genes_file = DMRsReplicates_Genes_file[DMRsReplicates_Genes_file$gene_model_type == "protein_coding", ]

    ######################### TE vs Genes plot
    GvsT_df = data.frame(
        x = c("Genes", "TEs"),
        y = c(nrow(DMRsReplicates_Genes_file), nrow(DMRsReplicates_TE_file))
    )
    GvsT_plot = ggplot(GvsT_df, aes(x = x, y = y, fill = x)) +
        geom_bar(stat = "identity", position = "stack") +
        geom_col(color = "black") +
        xlab("") +
        ylab("Number of DMRs") +
        ggtitle(context) +
        theme_classic() +
        theme(
            axis.text.x = element_text(face = "bold", size = 11),
            legend.position = "none"
        ) +
        scale_fill_manual(values = c("Genes" = "#6d6b6b", "TEs" = "#ab846b"))


Height_v = 1.75 # 2.1
Width_v = 1.6

    svg(file = paste0(output_path, context, "_TE.vs.ProteinCodingGenes.svg"), width = Width_v, height = Height_v, family = "serif")
    plot(GvsT_plot)
    dev.off()
    # }
}
