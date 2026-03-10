library(dplyr)
#library(tidyr)
library(stringr)
library(purrr)
library(writexl)
library(openxlsx)

treatment = "mto1"
all_res = as.data.frame(readxl::read_xlsx("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_genes.xlsx",
                                          sheet = treatment))

load_corr_res <- function(cntx) {
  x = rbind(
    read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/mto1_paper_DMRs_RNA_corr/mto1/Genes/mto1_DMRs_DEGs_corr_genes_",cntx,".csv")),
    read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/mto1_paper_DMRs_RNA_corr/mto1/Promoters/mto1_DMRs_DEGs_corr_promoters_",cntx,".csv"))
  ) %>%
    select(gene_id, cor, cor_pval) %>%
    mutate(cor = round(cor, 3)) %>%
    mutate(cor_pval = as.numeric(sprintf("%.3e", cor_pval))) %>%
    merge.data.frame(., all_res, by = "gene_id")
}

remove_dup_DMR <- function(y) {
  y <- as.character(unique(unlist(strsplit(y, ","))))
  paste(y, collapse = ",")
}

################################################################################
# old group method
{
#offspring_fun <- function(go_id, xx = as.list(GO.db::GOBPOFFSPRING)) { # 'GOBPCHILDREN' for child terms
#  
#  child_terms_0 = as.character(xx[[go_id]])
#  child_terms = child_terms_0
#  
#  for (i in 1:length(child_terms_0)) {
#    child_terms = c(child_terms, as.character(xx[[child_terms[i]]]))
#  } 
#  
#  return(child_terms[!is.na(child_terms)] %>% unique()) # %>% paste(collapse = "|"))
#}
#
#grep_position <- function(x) {
#  vec = NULL
#  for (terms_l in x) {
#    vec = c(vec, grepl(terms_l, all_res$GO.biological.process))
#  }
#  return(unique(vec))
#}
## GO term offsprings
#child_terms_epigenetic = offspring_fun("GO:0040029") # epigenetic regulation of_gene expression
#child_terms_chromatin_org = offspring_fun("GO:0006325") # chromatin organization
#child_terms_chromatin_rem = offspring_fun("GO:0006338") # chromatin remodeling
}
################################################################################

cor_res = rbind(load_corr_res("CG"),
                load_corr_res("CHG"),
                load_corr_res("CHH")) %>%
  # Group by gene_id and summarize (add DMRs values for the gene, make in one row for each gene)
  group_by(gene_id) %>%
  summarise(
    across(contains("_Genes") | contains("_Promoters"), ~remove_dup_DMR(paste(., collapse = ","))),  # apply remove_dup_DMR function for DMR columns
    across(!contains("_Genes") & !contains("_Promoters"), first)  # for other columns
  ) %>%
  as.data.frame() %>%
  mutate(across(contains("_Genes") | contains("_Promoters"), ~ str_replace_all(.x, "NA", ""))) %>%
  mutate(across(contains("_Genes") | contains("_Promoters"), ~ str_replace_all(.x, ",", ", "))) %>% # change delimiter
  relocate(contains("_Genes") | contains("_Promoters"), .before = Symbol) %>%
  arrange(RNA_pvalue)
# old group method
{
# %>%
  # Adding a new column 'group' to cor_res by filtering with the appropriate conditions
  #mutate(group = case_when(
  #  grepl("^2\\.1\\.1\\.37", EC) ~ "DNA_methyltransferase",
  #  grepl("dna \\(cytosine-5\\)-methyltransferase", tolower(Protein.names)) ~ "DNA_methyltransferase",
  #  grepl("dna \\(cytosine-5\\)-methyltransferase", tolower(Short_description)) ~ "DNA_methyltransferase",
  #  grepl("dna \\(cytosine-5\\)-methyltransferase", tolower(Computational_description)) ~ "DNA_methyltransferase",
  #  grepl("dna \\(cytosine-5\\)-methyltransferase", tolower(GO.biological.process)) ~ "DNA_methyltransferase",
  #  
  #  grepl("SET domain", gsub("SET-domain", "SET domain", Short_description)) ~ "Histone_Lysine_MTs",
  #  grepl("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(Symbol)) ~ "Histone_Lysine_MTs",
  #  grepl("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(old_symbols)) ~ "Histone_Lysine_MTs",
  #  grepl("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(Protein.names)) ~ "Histone_Lysine_MTs",
  #  grepl("class v-like sam-binding methyltransferase", tolower(Protein.families)) ~ "Histone_Lysine_MTs",
  #  
  #  grepl("chromatin remodeling|chromatin remodeler|chromatin structure", tolower(GO.biological.process)) ~ "chromatin_organization_related",
  #  grepl("chromatin remodeling|chromatin remodeler|chromatin structure", tolower(Gene_description)) ~ "chromatin_organization_related",
  #  grepl("chromatin remodeling|chromatin remodeler|chromatin structure", tolower(Computational_description)) ~ "chromatin_organization_related",
  #  grepl("chromatin remodeling|chromatin remodeler|chromatin structure", tolower(Function)) ~ "chromatin_organization_related",
  #  grepl("chromatin remodeling|chromatin remodeler|chromatin structure", tolower(Short_description)) ~ "chromatin_organization_related",
  #  
  #  grepl("chromosome", tolower(GO.cellular.component)) ~ "chromatin_organization_related",
  #  grepl("chromatin", tolower(GO.cellular.component)) ~ "chromatin_organization_related",
  #  grepl("chromatin", tolower(Function)) ~ "chromatin_organization_related",
  #  grepl("nucleosome", tolower(GO.cellular.component)) ~ "chromatin_organization_related",
  #  grepl("jumonji", tolower(Short_description)) ~ "chromatin_organization_related",
  #  grepl("helicase domain", tolower(Short_description)) ~ "chromatin_organization_related",
  #  grepl("histone h3 acetylation", tolower(GO.biological.process)) ~ "chromatin_organization_related",
  #  grep_position(child_terms_chromatin_org) ~ "chromatin_organization_related",
  #  grep_position(child_terms_chromatin_rem) ~ "chromatin_organization_related",
  #  
  #  grep_position(child_terms_epigenetic) ~ "epigenetic_reg_gene_expression",
  #  
  #  TRUE ~ NA_character_  # Default to NA for rows that don't match any condition
  #))
}
################################################################################





categories <- list(
  Metabolism = c("metabol","synthesis"),
  EpigeneticRegulation = c("chromatin", "histone", "methylation", "acetylation", "epigenetic"),
  #SignalTransduction = c("kinase", "phosphorylation", "signaling", "signal"),
  ProteinModification = c("ubiquitin", "ligase", "proteasome", "degradation"),
  Transport = c("transporter", "channel", "pump"),
  StressResponse = c("stress", "defense", "resistance")
)

assign_categories <- function(row) {
  function_desc <- tolower(paste(row['GO.molecular.function'], row['Function'], row['Computational_description'], sep=" "))
  assigned_categories <- c()
  for (category_name in names(categories)) {
    keywords <- categories[[category_name]]
    if (any(sapply(keywords, function(k) grepl(k, function_desc)))) {
      assigned_categories <- c(assigned_categories, category_name)
    }
  }
  if (length(assigned_categories) == 0) {
    return("Other")
  } else {
    return(paste(assigned_categories, collapse=", "))
  }
}

cor_res$Category <- apply(cor_res, 1, assign_categories)
cor_res = relocate(cor_res, Category, .before = gene_id)

# make 'Metabolism' category if its 'Other' and there is related metabolic pathways from DBs
cor_res$Category[(!is.na(cor_res$AraCyc.Db) | !is.na(cor_res$KEGG_pathway)) & cor_res$Category == "Other"] = "Metabolism"

cor_res = arrange(cor_res, Category)

# move 'Other' categoty to the bottom
other_tmp_df = filter(cor_res, Category == "Other")
cor_res = filter(cor_res, Category != "Other") %>% rbind(., other_tmp_df)

################################################################################
clean_ASCII <- function(x) {
  x = gsub("\002", " ", x)
  x = gsub("\036", " ", x)
  return(x)
}

################

xl_headers = names(cor_res)
numeric_cols = grep("^cor|RNA_|CG_|CHG_|CHH_", xl_headers)
RNA_cols = grep("RNA_", xl_headers)
DMRs_cols = grep("CG_|CHG_|CHH_", xl_headers)
p_cols = grep("cor_p|RNA_p", xl_headers)
lfc_cols = grep("^cor$|RNA_log2FC", xl_headers)
other_cols = (grep("CHH_Promoters", xl_headers)+2) : length(xl_headers)

################

# save and edit EXCEL
wb <- createWorkbook()
# Define styles
style_up <- createStyle(bgFill = "#f59d98")
style_down <- createStyle(bgFill = "#c3ccf7")
style_p <- createStyle(bgFill = "#f7deb0")
style_other <- createStyle(bgFill = "#daf7d7")
cell_n_font_style <- createStyle(border = "TopBottomLeftRight", borderColour = "black")
header_style <- createStyle(textDecoration = "bold", border = "TopBottomLeftRight", borderStyle = "double")

DMR_style_up <- createStyle(fgFill = "#f59d96", border = "TopBottomLeftRight", borderColour = "black")
DMR_style_down <- createStyle(fgFill = "#c3ccf7", border = "TopBottomLeftRight", borderColour = "black")
DMR_style_shared <- createStyle(fgFill = "#daf7d7", border = "TopBottomLeftRight", borderColour = "black")

# for (sheet_name in names(xl_list)) {
sheet_name = "mto1"
#df <- xl_list[[sheet_name]]
df = data.frame(lapply(cor_res, clean_ASCII))

#### for make this columns as numeric
#duplicate_DMRs_values_rows = df[,DMRs_cols] %>% # find duplicate DMRs values (they cant be numeric...)
#  mutate(across(everything(), ~grepl(",", .))) %>%
#  reduce(`|`) %>%
#  which()

#non_duplicate_DMRs_values_rows = (1:nrow(df))[-duplicate_DMRs_values_rows]
#df[-duplicate_DMRs_values_rows,DMRs_cols] = sapply(df[-duplicate_DMRs_values_rows, DMRs_cols], as.numeric)

df[,RNA_cols] = sapply(df[,RNA_cols], as.numeric)
df$cor = as.numeric(df$cor)
df$cor_pval = as.numeric(df$cor_pval)

#### 

addWorksheet(wb, sheet_name)
#  setColWidths(wb, sheet_name, cols = other_cols, widths = 10)
#  setColWidths(wb, sheet_name, cols = lfc_cols, widths = 4)
#  setColWidths(wb, sheet_name, cols = p_cols, widths = 6)


writeData(wb, sheet_name, df)

addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)

# colors for DMRs values (minus as blue as plus as red)
for (i.c in DMRs_cols) {
  for (i.r in 1:nrow(df)) {
    
    # for cells with both plus and minus direction
    if ( length(grep("^-.*, [0-9]", df[i.r,i.c])) != 0 | length(grep("^[0-9].*, -[0-9]", df[i.r,i.c])) != 0 ) {
      addStyle(wb, sheet_name, style = DMR_style_shared, 
               rows = i.r + 1, # first row in excel is the headers
               cols = i.c,
               gridExpand = FALSE)
      
      # for cells with minus direction
    } else if ( length(grep("-", df[i.r,i.c])) != 0 ) {
      addStyle(wb, sheet_name, style = DMR_style_down, 
               rows = i.r + 1, # first row in excel is the headers
               cols = i.c,
               gridExpand = FALSE)
      
      # for cells with plus direction 
    } else if ( length(grep("^[0-9]", df[i.r,i.c])) != 0 ) {
      addStyle(wb, sheet_name, style = DMR_style_up, 
               rows = i.r + 1, # first row in excel is the headers
               cols = i.c,
               gridExpand = FALSE)
    }
  }     
}
#      addStyle(wb, sheet_name, style = header_style, rows = duplicate_DMRs_values_rows , cols = DMRs_cols, gridExpand = TRUE)
#      addStyle(wb, sheet_name, style = header_style, rows = non_duplicate_DMRs_values_rows , cols = DMRs_cols, gridExpand = TRUE)


#for (col in p_cols) {
  #conditionalFormatting(wb, sheet_name, cols = p_cols[1], rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
  #conditionalFormatting(wb, sheet_name, cols = p_cols[2], rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
  #conditionalFormatting(wb, sheet_name, cols = p_cols[3], rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
#}

for (col in lfc_cols) {
  conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = ">0", style = style_up)
  conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = "<0", style = style_down)
}

for (col in other_cols[other_cols %% 2 == 0]) {
  conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
  conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
}

saveWorkbook(wb, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/merge_results_with_cor_groups_mto1.xlsx",
             overwrite = T)

# csv file for upload to web stuff
write.csv(df, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/CSV_merge_results_with_cor_groups_mto1.csv",
          row.names = F)
