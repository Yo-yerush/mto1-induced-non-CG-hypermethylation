library(dplyr)

TEG_DMR = rbind(
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/TEG_CG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/TEG_CHH_genom_annotations.csv")
) %>%
    select(seqnames, start, end, width, strand, gene_id, context, regionType, type) # pValue, log2FC,



# unique TEGs count
distinct(TEG_DMR, gene_id) %>% nrow()

# unique region count (DMRs in all contexts together)
TEG_DMR %>%
    mutate(tmp = paste(seqnames, start, end, sep = "_")) %>%
    distinct(tmp) %>%
    nrow()


############
# table for supp (grouped tabel)
TEG_DMR %>%
    mutate(tmp = paste(seqnames, start, end, sep = "_")) %>%
    group_by(tmp) %>%
    summarise(
        seqnames, start, end, width, strand, gene_id,
        context = paste(unique(context), collapse = ", "),
        regionType, type
    ) %>%
    as.data.frame() %>%
    distinct(tmp, .keep_all = T) %>%
    select(-tmp) %>%
    write.csv(., "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEG_DMRs_tabel.csv", row.names = F)

############





