library(dplyr)
library(openxlsx)

##################### DEGs IDs
DEGs_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(padj < 0.05) %>%
    rename(gene_id = locus_tag)
    
up_DEGs <- DEGs_file %>%
    filter(log2FoldChange > 0) %>%
    select(gene_id)

down_DEGs <- DEGs_file %>%
    filter(log2FoldChange < 0) %>%
    select(gene_id)
#####################

features_df <- data.frame(
    Feature = c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs"),
    CG = NA,
    CHG = NA,
    CHH = NA,
    Total_DMRs = NA,
    unique_DMRs = NA,
    overlapping_upregulated_DEGs = NA,
    overlapping_downregulated_DEGs = NA
)

for (i.feature in c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs")) {
    i.pos <- grep(i.feature, features_df$Feature)

    ### by context
    for (cntx in c("CG", "CHG", "CHH")) {
        features_df[i.pos, cntx] = read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", cntx, "/", i.feature, "_", cntx, "_genom_annotations.csv")) %>%
            # distinct(., gene_id) %>%  # count by DMRs. distinct for count by genes
            nrow()
    }

    ### total DMRs
    feature_df <- rbind(
        read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/", i.feature, "_CG_genom_annotations.csv")),
        read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/", i.feature, "_CHG_genom_annotations.csv")),
        read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/", i.feature, "_CHH_genom_annotations.csv"))
    ) %>%
        select(gene_id)

    features_df$Total_DMRs[i.pos] <- feature_df %>% nrow()
    features_df$unique_DMRs[i.pos] <- feature_df %>% distinct(., gene_id) %>% nrow()
    features_df$overlapping_upregulated_DEGs[i.pos] <- feature_df %>% distinct(., gene_id) %>% merge.data.frame(., up_DEGs) %>% nrow()
    features_df$overlapping_downregulated_DEGs[i.pos] <- feature_df %>% distinct(., gene_id) %>% merge.data.frame(., down_DEGs) %>% nrow()
}

features_df$Feature <- gsub("fiveUTRs", "5'UTRs", features_df$Feature)
features_df$Feature <- gsub("threeUTRs", "3'UTRs", features_df$Feature)

xl_headers <- names(features_df)
################
# save and edit EXCEL
wb <- createWorkbook()
# Define styles
cell_n_font_style <- createStyle(fontName = "Times New Roman", borderColour = "black")
last_row_style <- createStyle(fontName = "Times New Roman", border = "Bottom", borderStyle = "thick", borderColour = "black")
header_style <- createStyle(fontName = "Times New Roman", textDecoration = "bold", border = "Bottom", borderStyle = "thick", borderColour = "black")

sheet_name <- "Gene_features_count"
df <- features_df

addWorksheet(wb, sheet_name)
writeData(wb, sheet_name, df)

addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
addStyle(wb, sheet_name, style = last_row_style, rows = nrow(df) + 1, cols = 1:ncol(df), gridExpand = TRUE)
addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)

# Remove gridlines
showGridLines(wb, sheet_name, showGridLines = FALSE)

saveWorkbook(wb, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMR_DEG_features_count.xlsx", overwrite = T)
