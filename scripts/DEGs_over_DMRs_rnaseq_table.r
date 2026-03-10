library(dplyr)

########
# RNAseq - DEGs
RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag, log2FC = log2FoldChange, Symbol = gene) %>%
    dplyr::select(gene_id, log2FC, padj, Symbol, short_description) %>%
    filter(padj < 0.05)
RNAseq$log2FC <- round(RNAseq$log2FC, 2)
RNAseq$padj <- formatC(RNAseq$padj, format = "e", digits = 2, drop0trailing = TRUE) %>% as.numeric()
RNAseq$Symbol[is.na(RNAseq$Symbol)] <- ""

#########
# DMRs
ann_dir_path <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/"

genes_DMR <- rbind(
    read.csv(paste0(ann_dir_path, "CG/Promoters_CG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path, "CG/Genes_CG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path, "CHG/Promoters_CHG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path, "CHG/Genes_CHG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path, "CHH/Promoters_CHH_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path, "CHH/Genes_CHH_genom_annotations.csv"))
) %>%
    #filter(gene_model_type == "protein_coding") %>%
    distinct(gene_id)

merged_rna_dmr = merge.data.frame(RNAseq, genes_DMR, by = "gene_id")

write.csv(merged_rna_dmr, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DEGs_overlapped_table.csv", row.names = F)



############################################
n.total = nrow(merged_rna_dmr)
n.down = filter(merged_rna_dmr, log2FC < 0) %>% nrow()
(n.down / n.total) * 100
