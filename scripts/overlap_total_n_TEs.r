library(dplyr)
#library(ggplot2)
#library(cowplot)
#library(grid)
library(GenomicRanges)

RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
    dplyr::select(gene_id, DEG_log2FC, padj) %>%
    filter(padj < 0.05)

################################################
################################################
# total overlapping

overlap_general = function(context.f.1, ann.f.1) {

    feature_file_fun <- function(context.f, ann.f) {
        feature_file <- data.frame()

        for (context in context.f) {
            for (ann in ann.f) {
                # DMR results file
                ann_DMRs <- read.csv(paste0(
                    "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/",
                    context, "/", ann, "_", context, "_genom_annotations.csv"
                )) %>%
                    dplyr::select(
                        gene_id, log2FC, context, type
                    ) %>%
                    dplyr::rename(DMR_log2FC = log2FC)

                feature_file <- rbind(feature_file, ann_DMRs)
            }
        }
        return(feature_file)
    }

    df <- data.frame()
    for (xx in context.f.1) {
        for (yy in ann.f.1) {
            df <- rbind(df, feature_file_fun(xx, yy))
        }
    }

    df = df %>%
        merge(RNAseq, ., by = "gene_id") %>%
        distinct(gene_id, .keep_all = T)
    
    return(df)
    
    ## merged_list <- list(
    #all_res <- feature_file_fun(c("CG", "CHG", "CHH"), c("Genes", "Promoters")) %>%
    #    merge(RNAseq, ., by = "gene_id") %>%
    #    distinct(gene_id, .keep_all = T) # ,
    #CG_genes <- feature_file_fun("CG", "Genes") %>%
    #    merge(RNAseq, ., by = "gene_id") %>%
    #    distinct(gene_id, .keep_all = T) # ,
    #non_CG_promoters <- feature_file_fun(c("CHG", "CHH"), "Promoters") %>%
    #    merge(RNAseq, ., by = "gene_id") %>%
    #    distinct(gene_id, .keep_all = T)
}

#)

################################################
################################################
# TEs overlapping

#upstream <- readxl::read_excel("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_near_DEGs/TEs_near_DEGs_201124.xlsx",
#    sheet = "Up.Stream"
#) %>%
#    dplyr::select(gene_id)
#
#downstream <- readxl::read_excel("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_near_DEGs/TEs_near_DEGs_201124.xlsx",
#    sheet = "Down.Stream"
#) %>%
#    dplyr::select(gene_id)


overlapped_TE_sig_genes <- function(x, y) {

    overlapped_TE_context <- function(context, ann) {
        TE <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Transposable_Elements_", context, "_genom_annotations.csv")) %>%
            makeGRangesFromDataFrame(keep.extra.columns = TRUE)

        overlapped <- data.frame()

        gene_feature <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/", ann, "_", context, "_genom_annotations.csv")) %>%
            makeGRangesFromDataFrame(keep.extra.columns = TRUE)

        m <- findOverlaps(gene_feature, TE)
        m_df <- gene_feature[queryHits(m)]
        m_df$overlapping_TE <- paste(TE$Transposon_Super_Family[subjectHits(m)], TE$Transposon_Family[subjectHits(m)], sep = ", ")
        # m_df$Transposon_Super_Family <- TE$Transposon_Super_Family[subjectHits(m)]
        # m_df$Transposon_Family <- TE$Transposon_Family[subjectHits(m)]

        overlapped <- rbind(overlapped, as.data.frame(m_df)) %>%
            dplyr::select(gene_id)

        return(overlapped)
    }

    df <- data.frame()
    for (xx in x) {
        for (yy in y) {
            df <- rbind(df, overlapped_TE_context(xx, yy))
        }
    }
    df <- df %>%
        distinct(gene_id) %>%
        merge.data.frame(RNAseq, ., by = "gene_id") %>%
        filter(padj < 0.05) %>%
        distinct(gene_id)

    return(df)
}


################################################
################################################


a = overlap_general(c("CG", "CHG", "CHH"), c("Genes", "Promoters"))
b = overlapped_TE_sig_genes(c("CG", "CHG", "CHH"), c("Genes", "Promoters"))


data.frame(
    Overlap_all = nrow(a),
    Overlap_TEs = nrow(b)
)

df_combos <- data.frame()
combos <- list(
    list(name = "CG-Genes", contexts = c("CG"), ann = c("Genes")),
    list(name = "CG-Promoters", contexts = c("CG"), ann = c("Promoters")),
    list(name = "CG-Genes-Promoters", contexts = c("CG"), ann = c("Genes", "Promoters")),
    list(name = "CHG-Genes", contexts = c("CHG"), ann = c("Genes")),
    list(name = "CHG-Promoters", contexts = c("CHG"), ann = c("Promoters")),
    list(name = "CHG-Genes-Promoters", contexts = c("CHG"), ann = c("Genes", "Promoters")),
    list(name = "CHH-Genes", contexts = c("CHH"), ann = c("Genes")),
    list(name = "CHH-Promoters", contexts = c("CHH"), ann = c("Promoters")),
    list(name = "CHH-Genes-Promoters", contexts = c("CHH"), ann = c("Genes", "Promoters")),
    list(name = "CG-CHG-Genes", contexts = c("CG", "CHG"), ann = c("Genes")),
    list(name = "CG-CHG-Promoters", contexts = c("CG", "CHG"), ann = c("Promoters")),
    list(name = "CG-CHG-Genes-Promoters", contexts = c("CG", "CHG"), ann = c("Genes", "Promoters")),
    list(name = "CG-CHH-Genes", contexts = c("CG", "CHH"), ann = c("Genes")),
    list(name = "CG-CHH-Promoters", contexts = c("CG", "CHH"), ann = c("Promoters")),
    list(name = "CG-CHH-Genes-Promoters", contexts = c("CG", "CHH"), ann = c("Genes", "Promoters")),
    list(name = "CHG-CHH-Genes", contexts = c("CHG", "CHH"), ann = c("Genes")),
    list(name = "CHG-CHH-Promoters", contexts = c("CHG", "CHH"), ann = c("Promoters")),
    list(name = "CHG-CHH-Genes-Promoters", contexts = c("CHG", "CHH"), ann = c("Genes", "Promoters")),
    list(name = "CG-CHG-CHH-Genes", contexts = c("CG", "CHG", "CHH"), ann = c("Genes")),
    list(name = "CG-CHG-CHH-Promoters", contexts = c("CG", "CHG", "CHH"), ann = c("Promoters")),
    list(name = "CG-CHG-CHH-Genes-Promoters", contexts = c("CG", "CHG", "CHH"), ann = c("Genes", "Promoters"))
)

for (combo in combos) {
    g <- overlap_general(combo$contexts, combo$ann)
    t <- overlapped_TE_sig_genes(combo$contexts, combo$ann)
    df_combos <- rbind(df_combos, data.frame(
        combination = combo$name,
        Overlap_TEs = nrow(t),
        Overlap_all = nrow(g)
    ))
}

df_combos$ratio = round(df_combos$Overlap_TEs / df_combos$Overlap_all, 2)
df_combos = arrange(df_combos, desc(ratio))

openxlsx::write.xlsx(
    df_combos,
    "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEs_n_DEGs_overlap_withis_total_all_combos.xlsx"
)
