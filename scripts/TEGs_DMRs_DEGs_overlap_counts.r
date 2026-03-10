library(dplyr)
library(VennDiagram)

DMR_file.0 <- rbind(
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/TEG_CG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/TEG_CHH_genom_annotations.csv")
) %>%
    select(gene_id, regionType, context)


RNA_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(gene_model_type == "transposable_element_gene") %>%
    dplyr::rename(gene_id = "locus_tag") %>%
    filter(padj < 0.05) %>%
    distinct(gene_id, .keep_all = T) %>%
    select(gene_id, log2FoldChange, padj)


##########
# DMRs and DE-TEGs overlap
merge(DMR_file.0, RNA_file, by = "gene_id") %>%
    group_by(gene_id) %>%
    summarise(
        context = paste(unique(context), collapse = ", "),
        regionType = paste(unique(regionType), collapse = "; "),
        log2FoldChange, padj
    ) %>%
    as.data.frame() %>%
    distinct(gene_id, .keep_all = T) %>%
    View()
