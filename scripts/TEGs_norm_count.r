library(dplyr)

TEG_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_annotations.csv.gz") %>%
    filter(type == "transposable_element_gene") %>%
    distinct(gene_id)

#######################

DEGs <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag) %>%
    select(gene_id, log2FoldChange, padj) %>%
    filter(padj < 0.05) %>%
    # filter(pValue < 0.05) %>%
    distinct(gene_id, .keep_all = T)

DEGs$log2FoldChange = round(DEGs$log2FoldChange, 2)
DEGs$padj <- formatC(DEGs$padj, format = "e", digits = 2, drop0trailing = TRUE) %>% as.numeric

#######################

RNAseq = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/norm.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = X) %>%
    mutate(across(where(is.numeric), round, digits = 1)) # round

wt_columns = grep("met20|met21|met22", names(RNAseq))
mto1_columns = grep("met14|met15|met16", names(RNAseq))

names(RNAseq)[wt_columns] <- paste0("wt_", 1:3)
names(RNAseq)[mto1_columns] <- paste0("mto1_", 1:3)

RNAseq$wt_sum = rowSums(RNAseq[wt_columns], na.rm = TRUE)
RNAseq$mto1_sum = rowSums(RNAseq[mto1_columns], na.rm = TRUE)

#######################

merged_0 = merge.data.frame(TEG_file, RNAseq, by = "gene_id")
merged = merge.data.frame(merged_0, DEGs, by = "gene_id") %>%
    arrange(padj)

#######################

openxlsx::write.xlsx(
    merged,
    "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEGs_norm_count.xlsx"
)
