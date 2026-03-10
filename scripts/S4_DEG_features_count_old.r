library(dplyr)
library(openxlsx)

read_DMRs = function(c,f) {
    x = read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", c, "/", f, "_", c, "_genom_annotations.csv"))
    return(x)
}

##################### DEGs IDs
DEGs_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(padj < 0.05) %>%
    rename(gene_id = locus_tag)
    
up_DEGs <- DEGs_file %>% filter(log2FoldChange > 0) %>% select(gene_id)
down_DEGs <- DEGs_file %>% filter(log2FoldChange < 0) %>% select(gene_id)
#####################

features_df_up <- data.frame(
    Feature = c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs"),
    CG = NA,
    CHG = NA,
    CHH = NA,
    Total_DEGs = NA
)
features_df_down = features_df_up

for (i.feature in features_df_up$Feature) {
    i.pos <- grep(i.feature, features_df_up$Feature)

    ### by context
    for (cntx in c("CG", "CHG", "CHH")) {
        DMR_ann_res = read_DMRs(cntx, i.feature) %>%
            select(gene_id, regionType)
        
        gain_tairs = DMR_ann_res %>% filter(regionType == "gain") %>% distinct(gene_id)
        loss_tairs = DMR_ann_res %>% filter(regionType == "loss") %>% distinct(gene_id)
        
        gain_up = merge(gain_tairs, up_DEGs) %>% nrow()
        loss_up = merge(loss_tairs, up_DEGs) %>% nrow()
        gain_down = merge(gain_tairs, down_DEGs) %>% nrow()
        loss_down = merge(loss_tairs, down_DEGs) %>% nrow()
        
        features_df_up[i.pos, cntx] = paste0(gain_up, "/", loss_up)
        features_df_down[i.pos, cntx] = paste0(gain_down, "/", loss_down)
    }
    
    DMRs_all_cntx = rbind(read_DMRs("CG", i.feature), read_DMRs("CHG", i.feature), read_DMRs("CHH", i.feature)) %>% distinct(gene_id)

    features_df_up[i.pos, "Total_DEGs"] <- merge(DMRs_all_cntx, up_DEGs) %>% nrow()
    features_df_down[i.pos, "Total_DEGs"] <- merge(DMRs_all_cntx, down_DEGs) %>% nrow()
}

for (direction in c("up","down")) {

    if(direction == "up") {
        loop_df = features_df_up
    } else {
        loop_df = features_df_down
    }

   loop_df$Feature <- gsub("fiveUTRs", "5'UTRs", loop_df$Feature)
   loop_df$Feature <- gsub("threeUTRs", "3'UTRs", loop_df$Feature)

   xl_headers <- names(loop_df)
   ################
   # save and edit EXCEL
   wb <- createWorkbook()
   # Define styles
   cell_n_font_style <- createStyle(fontName = "Times New Roman", borderColour = "black")
   last_row_style <- createStyle(fontName = "Times New Roman", border = "Bottom", borderStyle = "thick", borderColour = "black")
   header_style <- createStyle(fontName = "Times New Roman", textDecoration = "bold", border = "Bottom", borderStyle = "thick", borderColour = "black")

   sheet_name <- "DEG_features_count"
   df <- loop_df

   addWorksheet(wb, sheet_name)
   writeData(wb, sheet_name, df)

   addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
   addStyle(wb, sheet_name, style = last_row_style, rows = nrow(df) + 1, cols = 1:ncol(df), gridExpand = TRUE)
   addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)

   # Remove gridlines
   showGridLines(wb, sheet_name, showGridLines = FALSE)

   saveWorkbook(wb, paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/", direction, "_DEG_features_count.xlsx"), overwrite = T)

}
