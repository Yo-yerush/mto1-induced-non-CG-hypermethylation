library(dplyr)

############################################

DMR_file <- rbind(
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/TEG_CG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/TEG_CHH_genom_annotations.csv")
) %>%
    select(seqnames, start, end, context, regionType, gene_id)

#DMR_file <- DMR_file.0 %>%
#    distinct(gene_id, .keep_all = T)

############################################

RNA_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(gene_model_type == "transposable_element_gene") %>%
    dplyr::rename(gene_id = "locus_tag") %>%
    filter(padj < 0.05) %>%
    select(gene_id, log2FoldChange, padj) %>%
    distinct(gene_id, .keep_all = T)
    
RNA_file$log2FoldChange = round(RNA_file$log2FoldChange, 2)
RNA_file$padj <- formatC(RNA_file$padj, format = "e", digits = 2, drop0trailing = TRUE) %>% as.numeric()

normCounts_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/norm.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = X) %>%
    mutate(across(where(is.numeric), round, digits = 1)) # round

############################################

# use TEGs with 'Derives_from' column
gene_2_TE_ids <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_description_file.csv.gz") %>%
    filter(!is.na(Derives_from)) %>%
    distinct(gene_id, Derives_from)

# filter to retro-TEs
TE_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt",
    sep = "\t"
) %>%
    mutate(seqnames = NA) %>% # Add a new column with NA values
    dplyr::rename(Derives_from = Transposon_Name) %>%
    dplyr::select(Derives_from, Transposon_Family, Transposon_Super_Family) %>%
    filter(grepl("Gypsy|Copia|LINE", Transposon_Super_Family)) %>%
    merge(., gene_2_TE_ids, by = "Derives_from")

############################################

retro_DMRs <- merge(DMR_file, TE_file, by = "gene_id")
retro_DEGs <- merge(RNA_file, TE_file, by = "gene_id")

############################################
### tables for supp (grouped tabel)

# DMRs over TEGs table
retro_DMRs %>%
    mutate(tmp = paste(seqnames, start, end, sep = "_")) %>%
    group_by(tmp) %>%
    summarise(
        seqnames, start, end,
        context = paste(unique(context), collapse = ", "),
        regionType, gene_id, Derives_from, Transposon_Family, Transposon_Super_Family
    ) %>%
    as.data.frame() %>%
    distinct(tmp, .keep_all = T) %>%
    select(-tmp) %>%
    write.csv(., "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEG_DMRs_tabel.csv", row.names = F)


# DEGs over TEGs (include normCounts)
wt_columns = grep("met20|met21|met22", names(normCounts_file))
mto1_columns = grep("met14|met15|met16", names(normCounts_file))

names(normCounts_file)[wt_columns] <- paste0("wt_", 1:3)
names(normCounts_file)[mto1_columns] <- paste0("mto1_", 1:3)

normCounts_file$mto1_sum = rowSums(normCounts_file[mto1_columns], na.rm = TRUE)
normCounts_file$wt_sum = rowSums(normCounts_file[wt_columns], na.rm = TRUE)

merged = merge.data.frame(normCounts_file, retro_DEGs, by = "gene_id") %>%
    arrange(padj) %>%
        write.csv(., "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEG_expression_tabel.csv", row.names = F)

############################################

message(
    "\n\n",
    "hyper-DMRs: ",
    (nrow(filter(retro_DMRs, regionType == "gain")) / nrow(retro_DMRs)) * 100,
    "\n",
    "downregulated DEGs: ",
    (nrow(filter(retro_DEGs, log2FoldChange < 0)) / nrow(retro_DEGs)) * 100,
    "\n\n",
    "overlap DMRs (unique TEGs): ",
    nrow(retro_DMRs),
    "\n",
    "overlap DMRs (not unique TEGs): ",
    DMR_file.0 %>%
        merge(., TE_file, by = "gene_id") %>%
        mutate(tmp = paste(seqnames, start, end, sep = "_")) %>%
        distinct(tmp) %>%
        nrow(),
    "\n\n",
    "DE-TEGs: ",
    nrow(retro_DEGs),
    "\n",
    "overlap DE-TEGs: ",
        length(intersect(retro_DMRs$gene_id, retro_DEGs$gene_id)),
    "\n"
)


yo = retro_DEGs
yo2 = retro_DMRs[,c("gene_id", "regionType")] %>% distinct(gene_id, .keep_all = T)

yo3 = merge(yo2, yo, by = "gene_id")
write.csv(yo3, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEG_overlapped_ids_table.csv", row.names = F)
