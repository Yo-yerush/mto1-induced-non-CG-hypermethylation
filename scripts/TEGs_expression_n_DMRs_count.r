library(dplyr)

DMR_file <- rbind(
    #read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/TEG_CG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/TEG_CHH_genom_annotations.csv")
) %>%
    mutate(tmp = paste(seqnames, start, end, gene_id, sep = "_")) 
    
    DMRs_count = DMR_file %>% distinct(tmp, .keep_all = T) %>% select(gene_id)

    DMRs_unique_TEGs_count = DMR_file %>% distinct(gene_id)


RNA_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(gene_model_type == "transposable_element_gene") %>%
    dplyr::rename(gene_id = "locus_tag") %>%
    .[, 1:4] %>%
    filter(padj < 0.05)

RNA_DMR = merge(DMRs_unique_TEGs_count, RNA_file, by = "gene_id")

downregulated = RNA_DMR %>% filter(log2FoldChange < 0) %>% nrow()

# overlap TEGs out of all DE-TEGs
(nrow(RNA_DMR) / nrow(RNA_file)) * 100

# downregulated out of overlap TEGs
(downregulated / nrow(RNA_DMR)) * 100
