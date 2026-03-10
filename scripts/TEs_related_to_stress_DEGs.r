library(dplyr)
library(tidyr)

stress_paper <- readxl::read_excel("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/Nearby transposable elements impact plant stress gene regulatory networks/12864_2021_8215_MOESM2_ESM.xlsx", sheet = 1) %>%
    as.data.frame() %>%
    select(Stress, "TE-gene location", "Differential expression", TEF, "Genes near TEF and differentially expressed") %>%
    rename(gene_id = "Genes near TEF and differentially expressed") %>%
    separate_rows(gene_id, sep = ",")



RNAseq_2 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag)



merged = merge.data.frame(stress_paper, RNAseq_2, by = "gene_id") %>% 
    #merge.data.frame(., overlapped_TE_sig_genes, by = "gene_id") %>%
    filter(padj < 0.05) %>%
    arrange(padj)
