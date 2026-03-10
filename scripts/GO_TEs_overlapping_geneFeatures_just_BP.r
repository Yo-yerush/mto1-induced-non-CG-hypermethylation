library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(GenomicRanges)
library(topGO)

RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
    dplyr::select(gene_id, DEG_log2FC, padj, pValue)

overlapped_TE_context <- function(context) {
    TE <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Transposable_Elements_", context, "_genom_annotations.csv")) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    overlapped <- data.frame()
    for (i.feature in c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs")) {
        gene_feature <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/", i.feature, "_", context, "_genom_annotations.csv")) %>%
            makeGRangesFromDataFrame(keep.extra.columns = TRUE)

        m <- findOverlaps(gene_feature, TE)
        m_df <- gene_feature[queryHits(m)]
        m_df$overlapping_TE <- paste(TE$Transposon_Super_Family[subjectHits(m)], TE$Transposon_Family[subjectHits(m)], sep = ", ")
        # m_df$Transposon_Super_Family <- TE$Transposon_Super_Family[subjectHits(m)]
        # m_df$Transposon_Family <- TE$Transposon_Family[subjectHits(m)]

        overlapped <- rbind(overlapped, as.data.frame(m_df))
    }

    return(overlapped)
}

################
remove_dup_DMR <- function(y) {
    y <- as.character(unique(unlist(strsplit(y, ","))))
    paste(y, collapse = ",")
}
################

################
overlapped_TE <- rbind(
    overlapped_TE_context("CG"),
    overlapped_TE_context("CHG"),
    overlapped_TE_context("CHH")
) %>%
    dplyr::select(
        -seqnames, -start, -end, -width, -strand, -pValue, -direction, -regionType,
        -sumReadsM1, -sumReadsN1, -proportion1, -sumReadsM2, -sumReadsN2, -proportion2, -cytosinesCount
    ) %>%
    dplyr::rename(DMR_log2FC = log2FC) %>%
    mutate(
        DMR_log2FC = round(DMR_log2FC, 3),
        CG_DMRs = ifelse(grepl("CG", context), DMR_log2FC, NA),
        CHG_DMRs = ifelse(grepl("CHG", context), DMR_log2FC, NA),
        CHH_DMRs = ifelse(grepl("CHH", context), DMR_log2FC, NA)
    ) %>%
    merge.data.frame(RNAseq, ., by = "gene_id", all.y = TRUE) %>%
    mutate(tmp = paste(gene_id, type, sep = "_")) %>%
    group_by(tmp) %>%
    summarise(
        across(contains("CG_DMRs") | contains("CHG_DMRs") | contains("CHH_DMRs"), ~ remove_dup_DMR(paste(., collapse = ","))), # apply remove_dup_DMR function for DMR columns
        across(!contains("CG_DMRs") | contains("CHG_DMRs") | contains("CHH_DMRs"), dplyr::first) # for other columns
    ) %>%
    as.data.frame() %>%
    mutate(across(contains("_DMRs"), ~ gsub("NA", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",,", ",", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub("^,", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",$", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",", ", ", .))) %>%
    dplyr::relocate(CG_DMRs, CHG_DMRs, CHH_DMRs, .before = context) %>%
    dplyr::relocate(type, overlapping_TE, .after = gene_id) %>%
    dplyr::relocate(Computational_description, .before = Short_description) %>%
    dplyr::select(-context, -DMR_log2FC, -tmp, -(Short_description:last_col())) %>%
    # filter(pValue < 0.05) %>%
    arrange(pValue) %>%
    arrange(type) %>%
    mutate(across(contains("padj") | contains("pValue"), ~ gsub(" NA", NA, .)),
        pValue = as.numeric(formatC(.$pValue, format = "e", digits = 3)),
        padj = as.numeric(formatC(.$padj, format = "e", digits = 3)),
        DEG_log2FC = round(DEG_log2FC, 3)
    ) %>%
    ### add to this script:
    dplyr::select(gene_id, pValue) %>%
    distinct(gene_id, .keep_all = T)


  overlapped_TE$pValue[is.na(overlapped_TE$pValue)] <- 0.999
  overlapped_TE$pValue[overlapped_TE$pValue == 0] <- 1e-300

  geneList <- ifelse(overlapped_TE$pValue < 0.05, 1, 0)
  names(geneList) <- overlapped_TE$gene_id


res_list <- list()
for (GO_type_loop in c("BP", "MF", "CC")) {
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
    allRes$Fisher = as.numeric(allRes$Fisher)
    allRes$Term = gsub(",", ";", allRes$Term)
    allRes$type = GO_type_loop

    res_list[[GO_type_loop]] <- allRes[allRes$Fisher <= 0.01, ]
}

    ### plot
    # just BP !!!!!!!!!
    res_bind = rbind(res_list[["BP"]]) # , res_list[["CC"]], res_list[["MF"]])
    res_bind$type = "Biological Process"

    res_bind = res_bind[!grepl("cellular_component|biological_process|molecular_function|macromolecule biosynthetic process", res_bind$Term), ]

    #gain_col = "#cf534c"
    #loss_col = "#6397eb"

    bubble_plot = res_bind %>%
        ggplot(aes(Significant, reorder(Term, Significant), size = Annotated, color = Fisher)) +
        scale_color_gradient("p.value", low = "#F2A672", high = "black") + # theme_classic() +
        labs( # title = paste0(type_name," - ",gain_loss,"regulated transcripts"),
            x = "Significant", y = ""
        ) +
        theme_bw() +
        theme(
            # plot.title=element_text(hjust=0.5),
            #legend.key.size = unit(0.25, "cm"),
            #legend.title = element_text(size = 9.5),
            #legend.position = "right",
            #text = element_text(family = "serif")
            legend.position = "none"
        ) +
        geom_point() +
        facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
        guides(color = guide_colorbar(order = 1, barheight = 4))

    # print(paste(treatment,context,annotation, sep = "_"))

svg("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_DEGs_overlap_with_TEs.svg", width = 4.25, height = 4, family = "serif")
print(bubble_plot)
dev.off()

svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_DEGs_overlap_with_TEs_Annotated_legend.svg"), width = 1, height = 1.5, family = "serif")
legend_plot <- ggplot(res_bind, aes(x = 1, y = 1, size = Annotated)) +
    geom_point() +
        scale_size(
        #range = c(1, 6),
        name = " ",
        breaks = c(min(res_bind$Annotated), 50, 100, max(res_bind$Annotated))
    ) +
    theme_void() +
    theme(legend.position = "right")

legend_only <- cowplot::get_legend(legend_plot)
grid.newpage()
grid.draw(legend_only)
dev.off()




########### get genes from term
GO_term <- "GO:0006355"
genes_in_term <- genesInTerm(myGOdata, GO_term)
tair_gene_list <- data.frame(gene_id = genes_in_term[[GO_term]])

RNAseq_2 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag)

a = merge.data.frame(tair_gene_list, RNAseq_2, by = "gene_id") %>% arrange(pValue)
b = a %>% filter(pValue < 0.05)
