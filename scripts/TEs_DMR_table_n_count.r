library(dplyr)

ann_dir_path = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/"

TEs_DMR = rbind(
    read.csv(paste0(ann_dir_path, "CG/Transposable_Elements_CG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path, "CHG/Transposable_Elements_CHG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path, "CHH/Transposable_Elements_CHH_genom_annotations.csv"))
) %>% 
    select(seqnames, start, end, width, strand, context, regionType, gene_id, Transposon_Family, Transposon_Super_Family) # pValue, log2FC,


# DMRs count (each context)
nrow(TEs_DMR)

# unique TEGs count
distinct(TEs_DMR, gene_id) %>% nrow()

# unique region count (DMRs in all contexts together)
TEs_DMR %>%
    mutate(tmp = paste(seqnames, start, end, sep = "_")) %>%
    distinct(tmp) %>%
    nrow()


############
# table for supp (grouped tabel)
final_df = TEs_DMR %>%
    mutate(tmp = paste(seqnames, start, end, sep = "_")) %>%
    group_by(tmp) %>%
    summarise(
        seqnames, start, end, width, strand,
        context = paste(unique(context), collapse = ", "),
        regionType,
        gene_id = paste(unique(gene_id), collapse = "; "),
        Transposon_Family = paste(unique(Transposon_Family), collapse = "; "),
        Transposon_Super_Family = paste(unique(Transposon_Super_Family), collapse = "; ")
    ) %>%
    as.data.frame() %>%
    rename(TE_id = gene_id) %>%
    distinct(tmp, .keep_all = T) %>%
    select(-tmp)

write.csv(final_df, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEs_DMRs_tabel.csv", row.names = F)

############






