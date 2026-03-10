library(dplyr)
library(writexl)
library(openxlsx)

treatment = "mto1"
all_res = as.data.frame(readxl::read_xlsx("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_genes.xlsx", sheet = treatment))
all_res = all_res[,-ncol(all_res)] # %>% relocate("transcript_id", .after = "locus_tag")

################
offspring_fun <- function(go_id, xx = as.list(GO.db::GOBPOFFSPRING)) { # 'GOBPCHILDREN' for child terms
  
  child_terms_0 = as.character(xx[[go_id]])
  child_terms = child_terms_0
  
  for (i in 1:length(child_terms_0)) {
    child_terms = c(child_terms, as.character(xx[[child_terms[i]]]))
  } 
  
  return(child_terms[!is.na(child_terms)] %>% unique()) # %>% paste(collapse = "|"))
}

child_terms_response = offspring_fun("GO:0006950") # response to stress
child_terms_stimulus = offspring_fun("GO:0050896") # response to stimulus
child_terms_biotic = offspring_fun("GO:0009607") # response to biotic stimulus
child_terms_abiotic = offspring_fun("GO:0009628") # response to abiotic stimulus
child_terms_bio_p = offspring_fun("GO:0009058") # biosynthetic process
child_terms_aa = offspring_fun("GO:0008652") # amino acid biosynthetic process
child_terms_methionine = offspring_fun("GO:0009086") # methionine biosynthetic process
child_terms_chromatin_org = offspring_fun("GO:0006325") # chromatin organization
child_terms_chromatin_rem = offspring_fun("GO:0006338") # chromatin remodeling
################

grep_position <- function(x) {
  vec = NULL
  for (terms_l in x) {
    vec = c(vec, grep(terms_l, all_res$Gene.Ontology..biological.process.))
  }
  return(unique(vec))
}

################
xl_list = list(
  response_to_stress = all_res[grep_position(child_terms_response),], # %>% filter(RNA_pvalue < 0.05) # 'response to stress' child terms
  response_to_stimulus = all_res[grep_position(child_terms_stimulus),],
  response_to_biotic_stimulus = all_res[grep_position(child_terms_biotic),],
  response_to_abiotic_stimulus = all_res[grep_position(child_terms_abiotic),],
  biosynthetic_process = all_res[grep_position(child_terms_bio_p),],
  amino_acid_biosynthetic_process = all_res[grep_position(child_terms_aa),],
  methionine_biosynthetic_process = all_res[grep_position(child_terms_methionine),],
  #chromatin_organization = all_res[grep_position(child_terms_chromatin_org),]
  chromatin_remodeling = all_res[grep_position(child_terms_chromatin_rem),] # %>% filter(RNA_pvalue < 0.05) # 'defense response' related terms
)
################

################ eddit xl_list data frames
xl_list <- lapply(xl_list, function(x) {
  x = x %>% filter(RNA_pvalue < 0.05) %>%
    distinct(locus_tag, .keep_all = T) %>%
    arrange(RNA_pvalue) %>%
    select(-gene_model_type) %>% 
    relocate(c("gene", "short_description", "Function..CC."), .after = "RNA_pvalue")# %>% 
    #relocate("Function..CC.", .after = "short_description")
  
  x$Function..CC. = gsub("FUNCTION: ","", x$Function..CC.)
  names(x) = gsub("Gene.Ontology\\.\\.","GO\\.", names(x))
  names(x) = gsub("\\.\\.CC\\.","", names(x))
  
  return(x)
})
################

################
clean_ASCII <- function(x) {
  x = gsub("\002", " ", x)
  x = gsub("\036", " ", x)
  
#  x = gsub("[[:punct:]]", " ", x)
  #x = iconv(x, from = 'UTF-8', to = 'ASCII')
  return(x)
}
################
xl_headers = names(xl_list[[1]])
numeric_cols = grep("RNA_|CG_|CHG_|CHH_", xl_headers)
p_cols = grep("RNA_p", xl_headers)
lfc_cols = grep("RNA_log2FC|CG_|CHG_|CHH_", xl_headers)
other_cols = grep("Curator_summary", xl_headers) : length(xl_headers)
################
# save and edit EXCEL
wb <- createWorkbook()
# Define styles
style_up <- createStyle(fontName = "Times New Roman", bgFill = "#f59d98")
style_down <- createStyle(fontName = "Times New Roman", bgFill = "#c3ccf7")
style_p <- createStyle(fontName = "Times New Roman", bgFill = "#f7deb0")
style_other <- createStyle(fontName = "Times New Roman", bgFill = "#daf7d7")
cell_n_font_style <- createStyle(fontName = "Times New Roman", border = "TopBottomLeftRight", borderColour = "black")
header_style <- createStyle(fontName = "Times New Roman", textDecoration = "bold", border = "TopBottomLeftRight", borderStyle = "double")

for (sheet_name in names(xl_list)) {
  addWorksheet(wb, sheet_name)
#  setColWidths(wb, sheet_name, cols = other_cols, widths = 10)
#  setColWidths(wb, sheet_name, cols = lfc_cols, widths = 4)
#  setColWidths(wb, sheet_name, cols = p_cols, widths = 6)
  
  df <- xl_list[[sheet_name]]
  df <- data.frame(lapply(df, clean_ASCII))
  df[,numeric_cols] = sapply(df[, numeric_cols], as.numeric)
  
  writeData(wb, sheet_name, df)
  
  addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
  addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  
  conditionalFormatting(wb, sheet_name, cols = p_cols[1], rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
  conditionalFormatting(wb, sheet_name, cols = p_cols[2], rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
  
  for (col in lfc_cols) {
    conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = ">0", style = style_up)
    conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = "<0", style = style_down)
  }
  
  for (col in other_cols[other_cols %% 2 == 1]) {
    conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
    conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
  }
}
saveWorkbook(wb, paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/", treatment, "_GO_groups.xlsx"), overwrite = T)
################