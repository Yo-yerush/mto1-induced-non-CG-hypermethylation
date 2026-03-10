library(dplyr)

ann_dir_path = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/"

genes_DMR = rbind(
    read.csv(paste0(ann_dir_path,"CG/Promoters_CG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CG/CDS_CG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CG/Introns_CG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CG/fiveUTRs_CG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CG/threeUTRs_CG_genom_annotations.csv")),

    read.csv(paste0(ann_dir_path,"CHG/Promoters_CHG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CHG/CDS_CHG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CHG/Introns_CHG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CHG/fiveUTRs_CHG_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CHG/threeUTRs_CHG_genom_annotations.csv")),

    read.csv(paste0(ann_dir_path,"CHH/Promoters_CHH_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CHH/CDS_CHH_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CHH/Introns_CHH_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CHH/fiveUTRs_CHH_genom_annotations.csv")),
    read.csv(paste0(ann_dir_path,"CHH/threeUTRs_CHH_genom_annotations.csv"))
) %>% 
    select(seqnames, start, end, width, strand, context, regionType, type, gene_id, Symbol, Short_description) # pValue, log2FC,


############
# table for supp (grouped tabel)
final_df = genes_DMR %>%
    mutate(tmp = paste(seqnames, start, end, sep = "_")) %>%
    group_by(tmp) %>%
    summarise(
        seqnames, start, end, width, strand,
        context = paste(unique(context), collapse = "; "),
        regionType,
        type = paste(unique(type), collapse = "; "),
        gene_id = paste(unique(gene_id), collapse = "; "),
        Symbol = paste(unique(Symbol), collapse = "; "),
        Short_description = paste(unique(Short_description), collapse = "; ")
    ) %>%
    as.data.frame() %>%
    dplyr::rename(direction = regionType) %>%
    distinct(tmp, .keep_all = T) %>%
    select(-tmp)

final_df$Symbol = gsub("^; ", "", final_df$Symbol)
final_df$Symbol = gsub("; $", "", final_df$Symbol)
final_df$Short_description = gsub("^; ", "", final_df$Short_description)
final_df$Short_description = gsub("; $", "", final_df$Short_description)

write.csv(final_df, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/gene_features_DMRs_tabel.csv", row.names = F)