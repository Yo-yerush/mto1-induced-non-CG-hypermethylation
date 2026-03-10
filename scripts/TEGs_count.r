library(dplyr)
library(GenomicRanges)

##############################################################
############ TEGs expression
RNA_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    # filter(gene_model_type == "transposable_element_gene") %>%
    dplyr::rename(gene_id = "locus_tag") %>%
    .[, 1:4]

TEG_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_annotations.csv.gz") %>%
    filter(type == "transposable_element_gene") %>%
    merge.data.frame(., RNA_file, by = "gene_id") %>%
    select(-type, -gene_model_type) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)


#TE_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt",
#    sep = "\t"
#) %>%
#    mutate(seqnames = NA) %>% # Add a new column with NA values
#    dplyr::select(seqnames, Transposon_min_Start, Transposon_max_End, orientation_is_5prime, everything())
## dplyr::rename(gene_id = Transposon_Name)
#
#for (i in 1:5) {
#    TE_file$seqnames[grep(paste0("AT", i, "TE"), TE_file$Transposon_Name)] <- paste0("Chr", i)
#}
#TE_file$orientation_is_5prime <- gsub("true", "+", TE_file$orientation_is_5prime)
#TE_file$orientation_is_5prime <- gsub("false", "-", TE_file$orientation_is_5prime)
#
#names(TE_file)[1:4] <- c("seqnames", "start", "end", "strand")
#TE_gr <- makeGRangesFromDataFrame(TE_file, keep.extra.columns = T)

#################################################

RNA_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    # filter(gene_model_type == "transposable_element_gene") %>%
    dplyr::rename(gene_id = "locus_tag") %>%
    .[, 1:4]

TEG_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_annotations.csv.gz") %>%
    filter(type == "transposable_element_gene") %>%
    merge.data.frame(., RNA_file, by = "gene_id") %>%
    select(-type, -gene_model_type) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)

# overlap 'TE_file' and 'ann_file' GRanges objects by position
#m <- findOverlaps(TEG_file, TE_gr)
#TEG <- TEG_file[queryHits(m)]
#mcols(TEG) <- c(mcols(TEG), mcols(TE_gr[subjectHits(m)]))

TEG_df <- as.data.frame(TEG_file) %>%
    # add column that contain all the values from each row
    mutate(tmp = do.call(paste, c(., sep = "_"))) %>%
    distinct(tmp, .keep_all = T) %>%
    select(-tmp) %>%
    .[, -(1:5)] # %>% filter(padj < 0.05)


###
# total significants
total_sig_TEGs = TEG_df %>% filter(padj < 0.05)

# downregulated
downregulated_TEGs = total_sig_TEGs %>% filter(log2FoldChange < 0)

paste0(round((nrow(downregulated_TEGs) / nrow(total_sig_TEGs)) * 100, 2), "% downregulated TEGs")
##############################################################

##############################################################
############ overlap with DMRs
DMRs_file = rbind(
    read.csv("C:/Users/yonatany/OneDrive - Migal/Desktop/mto1_vs_wt/genome_annotation/CG/TEG_CG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/OneDrive - Migal/Desktop/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/OneDrive - Migal/Desktop/mto1_vs_wt/genome_annotation/CHH/TEG_CHH_genom_annotations.csv")
)

######### overlap DMRs with TEGs
overlap_DMRs_with_TEGs = DMRs_file %>%
    mutate(tmp = paste(seqnames, start, end, gene_id, sep = "_")) %>%
    distinct(tmp, .keep_all = T)
nrow(overlap_DMRs_with_TEGs)

gain_count <- overlap_DMRs_with_TEGs %>%
    filter(regionType == "gain") %>%
    nrow()
loss_count <- overlap_DMRs_with_TEGs %>%
    filter(regionType == "loss") %>%
    nrow()

(gain_count / nrow(overlap_DMRs_with_TEGs)) * 100
(loss_count / nrow(overlap_DMRs_with_TEGs)) * 100
#########

######### overlap DMRs with TEGs - unique TEGs count
overlap_DMRs_with_unique_TEGs = overlap_DMRs_with_TEGs %>%
    distinct(gene_id)
nrow(overlap_DMRs_with_unique_TEGs)

######### overlap DMRs with TEGs - TEGs expression
unique_rna_dmr_TEGs = DMRs_file %>%
    distinct(gene_id) %>%
    merge.data.frame(., total_sig_TEGs, by = "gene_id") %>%
nrow(unique_rna_dmr_TEGs)

upregulated_count <- unique_rna_dmr_TEGs %>%
    filter(log2FoldChange > 0) %>%
    nrow()
downregulated_count <- unique_rna_dmr_TEGs %>%
    filter(log2FoldChange < 0) %>%
    nrow()

upregulated_count
downregulated_count
(upregulated_count / nrow(unique_rna_dmr_TEGs)) * 100
(downregulated_count / nrow(unique_rna_dmr_TEGs)) * 100
#########