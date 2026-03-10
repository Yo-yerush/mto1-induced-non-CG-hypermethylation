library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(GenomicRanges)
library(topGO)

TE_overlap_beckground = F
against_sig_genes = T # only if 'TE_overlap_beckground = F'

########################

RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(gene_model_type == "protein_coding") %>%
    dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
    dplyr::select(gene_id, DEG_log2FC, padj, pValue)

########################

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
TE <- makeGRangesFromDataFrame(TE_file, keep.extra.columns = T)

########################

if (TE_overlap_beckground) {
    # filter to detected (in RNAseq analysis) protein-coding genes
    methylome_annotations <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_200225/annotation_files/Methylome.At_annotations.csv.gz") %>%
        filter(type == "gene") %>%
        filter(gene_model_type == "protein_coding") %>%
        merge.data.frame(., distinct(RNAseq, gene_id), by = "gene_id") %>%
        relocate(gene_id, .before = type) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = T)

        methylome_annotations = c(
            methylome_annotations,
            promoters(methylome_annotations)
        )

    m <- findOverlaps(methylome_annotations, TE)
    genes_overlap_TEs <- methylome_annotations[queryHits(m)]
    genes_overlap_TEs$overlapping_TE <- paste(TE$Transposon_Super_Family[subjectHits(m)], TE$Transposon_Family[subjectHits(m)], sep = ", ")

    genes_overlap_TEs <- as.data.frame(genes_overlap_TEs)
}

########################

overlapped_TE_context <- function(context) {
    #TE <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Transposable_Elements_", context, "_genom_annotations.csv")) %>%
    #    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    overlapped <- data.frame()

            gene_feature <- rbind(
                read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Promoters_", context, "_genom_annotations.csv")),
                read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Genes_", context, "_genom_annotations.csv"))
            ) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    #gene_feature <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/#yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/#Promoters_", context, "_genom_annotations.csv")) %>%
    #    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    m <- findOverlaps(gene_feature, TE)
    m_df <- gene_feature[queryHits(m)]
    m_df$overlapping_TE <- paste(TE$Transposon_Super_Family[subjectHits(m)], TE$Transposon_Family[subjectHits(m)], sep = ", ")
    # m_df$Transposon_Super_Family <- TE$Transposon_Super_Family[subjectHits(m)]
    # m_df$Transposon_Family <- TE$Transposon_Family[subjectHits(m)]

    overlapped <- rbind(overlapped, as.data.frame(m_df)) %>%
        dplyr::select(gene_id)

    return(overlapped)
}

################

overlapped_TE_sig_genes <- rbind(
    overlapped_TE_context("CG"),
    overlapped_TE_context("CHG"),
    overlapped_TE_context("CHH")
) %>%
    distinct(gene_id) %>%
    merge.data.frame(RNAseq, ., by = "gene_id") %>%
    filter(padj < 0.05) %>%
    mutate(for_geneList = 1) %>%
    dplyr::select(gene_id, for_geneList)

################

if (TE_overlap_beckground) {
    # beckgroud is overlap protein coding genes with TEs
   overlapped_TE <- merge(overlapped_TE_sig_genes, genes_overlap_TEs, by = "gene_id", all.y = T)
   is_te_beckground = "TE"

} else if (!TE_overlap_beckground) {
    # beckgroud is all protein coding genes
   overlapped_TE <- merge(overlapped_TE_sig_genes, RNAseq, by = "gene_id", all.y = T)
   is_te_beckground = "all_pcGenes"

} else if (against_sig_genes) {
    # beckgroud is just the DEGs
    overlapped_TE <- merge(overlapped_TE_sig_genes, filter(RNAseq, padj < 0.05), by = "gene_id", all.y = T)
    is_te_beckground = "DEGs"
}

overlapped_TE$for_geneList[is.na(overlapped_TE$for_geneList)] <- 0

geneList <- overlapped_TE$for_geneList
names(geneList) <- overlapped_TE$gene_id

################

res_list <- list()
for (GO_type_loop in c("BP")) { #, "MF", "CC")) {
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

    res_list[[GO_type_loop]] <- allRes[allRes$Fisher <= 0.01, ]
}

### plot
# just BP !!!!!!!!!
res_bind <- rbind(res_list[["BP"]]) # , res_list[["CC"]], res_list[["MF"]])
res_bind$type <- "Biological Process"

res_bind <- res_bind[!grepl("cellular_component|biological_process|molecular_function|macromolecule biosynthetic process", res_bind$Term), ]

write.csv(res_bind, paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_TE_overlap/GO_genes_promoters_DEGs_overlap_with_TEs_",is_te_beckground,"_beckground.csv"), row.names = F)

# gain_col = "#cf534c"
# loss_col = "#6397eb"

bubble_plot <- res_bind %>%
    ggplot(aes(Significant, reorder(Term, Significant), size = Annotated, color = Fisher)) +
    scale_color_gradient("p.value", low = "#F2A672", high = "black") + # theme_classic() +
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

svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_TE_overlap/GO_genes_promoters_DEGs_overlap_with_TEs_",is_te_beckground,"_beckground.svg"), width = 4.25, height = 3, family = "serif")
print(bubble_plot)
dev.off()

svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_TE_overlap/GO_genes_promoters_DEGs_overlap_with_TEs_Annotated_legend_",is_te_beckground,"_beckground.svg"), width = 1, height = 1.5, family = "serif")
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
