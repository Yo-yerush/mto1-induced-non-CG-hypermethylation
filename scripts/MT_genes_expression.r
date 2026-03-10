library(dplyr)

RNA_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = "locus_tag",
    Symbol = "gene") 

RNA_file$log2FoldChange = round(RNA_file$log2FoldChange, 2)
RNA_file$padj <- formatC(RNA_file$padj, format = "e", digits = 2, drop0trailing = TRUE) %>% as.numeric()
RNA_file$pValue <- formatC(RNA_file$pValue, format = "e", digits = 2, drop0trailing = TRUE) %>% as.numeric()



RNA_file1 <- RNA_file %>%
    filter(grepl("DMT|DRM|MET1$|CMT|DNMT", Symbol))

RNA_file2 <- RNA_file %>%
    filter(grepl("AT4G08990", gene_id))

rbind(RNA_file1, RNA_file2) %>%
arrange(padj) %>%
select(gene_id, log2FoldChange, padj, pValue, Symbol, short_description) %>%
write.csv(., "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/MT_genes.csv")
