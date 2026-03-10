library(dplyr)
library(openxlsx)

supplementary_data <- read.xlsx("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/supplementary/supplementary_tables_090125.xlsx", sheet = "S6 - DEGs, DMRs, and TEs")
names(supplementary_data) = supplementary_data[1, ]
supplementary_data = supplementary_data[-1, ]

supplementary_data$DEG_log2FC <- as.numeric(supplementary_data$DEG_log2FC)
supplementary_data$padj <- as.numeric(supplementary_data$padj)
supplementary_data$pValue <- as.numeric(supplementary_data$pValue)

df_CG = supplementary_data %>% filter(!is.na(CG_DMRs))
df_CHG = supplementary_data %>% filter(!is.na(CHG_DMRs))
df_CHH = supplementary_data %>% filter(!is.na(CHH_DMRs))

df_CG_gain = df_CG[-grep("-", df_CG$CG_DMRs), ]
df_CG_loss = df_CG[grep("-", df_CG$CG_DMRs), ]
df_CHG_gain = df_CHG[-grep("-", df_CHG$CHG_DMRs), ]
df_CHG_loss = df_CHG[grep("-", df_CHG$CHG_DMRs), ]
df_CHH_gain = df_CHH[-grep("-", df_CHH$CHH_DMRs), ]
df_CHH_loss = df_CHH[grep("-", df_CHH$CHH_DMRs), ]

df_up = supplementary_data %>% filter(DEG_log2FC > 0)
df_down = supplementary_data %>% filter(DEG_log2FC < 0)

#####################################

features_df_up <- data.frame(
    Feature = c("promoter", "CDS", "intron", "5'UTR", "3'UTR"),
    CG_gain = NA,
    CG_loss = NA,
    CHG_gain = NA,
    CHG_loss = NA,
    CHH_gain = NA,
    CHH_loss = NA,
    Total_DEGs = NA
)
features_df_down = features_df_up

#####################################

for (i.feature in features_df_up$Feature) {
    i.pos <- grep(i.feature, features_df_up$Feature)

    for(context in c("CG", "CHG", "CHH")) {
        df_gain_loop <- get(paste0("df_", context, "_gain"))
        df_loss_loop <- get(paste0("df_", context, "_loss"))
        
        features_df_up[[paste0(context, "_gain")]][i.pos]  <- df_gain_loop %>% filter(type == i.feature, DEG_log2FC > 0) %>% nrow()
        features_df_up[[paste0(context, "_loss")]][i.pos]  <- df_loss_loop %>% filter(type == i.feature, DEG_log2FC > 0) %>% nrow()
        
        features_df_down[[paste0(context, "_gain")]][i.pos]  <- df_gain_loop %>% filter(type == i.feature, DEG_log2FC < 0) %>% nrow()
        features_df_down[[paste0(context, "_loss")]][i.pos]  <- df_loss_loop %>% filter(type == i.feature, DEG_log2FC < 0) %>% nrow()
    }

    features_df_up$Total_DEGs[i.pos] = filter(df_up, type == i.feature) %>% nrow()
    features_df_down$Total_DEGs[i.pos] = filter(df_down, type == i.feature) %>% nrow()
}


for (direction in c("up","down")) {

    if(direction == "up") {
        loop_df = features_df_up
    } else {
        loop_df = features_df_down
    }

   xl_headers <- names(loop_df)
   ################
   # save and edit EXCEL
   wb <- createWorkbook()
   # Define styles
   cell_n_font_style <- createStyle(fontName = "Times New Roman", borderColour = "black")
   last_row_style <- createStyle(fontName = "Times New Roman", border = "Bottom", borderStyle = "thick", borderColour = "black")
   header_style <- createStyle(fontName = "Times New Roman", textDecoration = "bold", border = "Bottom", borderStyle = "thick", borderColour = "black")

   sheet_name <- "TEs-DEG_features_count"
   df <- loop_df

   addWorksheet(wb, sheet_name)
   writeData(wb, sheet_name, df)

   addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
   addStyle(wb, sheet_name, style = last_row_style, rows = nrow(df) + 1, cols = 1:ncol(df), gridExpand = TRUE)
   addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)

   # Remove gridlines
   showGridLines(wb, sheet_name, showGridLines = FALSE)

   saveWorkbook(wb, paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/", direction, "_TEs_DEG_features_count.xlsx"), overwrite = T)

}
