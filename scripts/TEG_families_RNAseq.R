library(dplyr)

ann_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_description_file.csv.gz") %>%
  filter(!is.na(Derives_from)) %>%
  #filter(type == "transposable_element_gene") %>%
  select(gene_id, Derives_from) %>%
  dplyr::rename(Transposon_Name = "Derives_from")

TE_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt",
                   sep = "\t") %>%
  select(Transposon_Name, Transposon_Family, Transposon_Super_Family)

RNA_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
  filter(gene_model_type == "transposable_element_gene") %>%
  dplyr::rename(gene_id = "locus_tag") %>%
  .[,1:4] %>% filter(padj < 0.05)
  


merged = merge.data.frame(RNA_file, ann_file, by = "gene_id") %>%
  merge.data.frame(., TE_file, by = "Transposon_Name")

######################################
grouped_families <- function(superFamily=NULL, Family=NULL) {
  
  # for 'annotate and significant' data frame
  if (!is.null(superFamily)) {
    xx = merged %>% filter(Transposon_Super_Family == superFamily)
    
    final_df = data.frame(families = unique(xx$Transposon_Family), annotate = NA, significant = NA)
    
    for (i.fam in 1:length(final_df$families)) {
      final_df$annotate[i.fam] = xx %>% filter(Transposon_Family == final_df$families[i.fam]) %>% nrow()
      final_df$significant[i.fam] = xx %>% filter(Transposon_Family == final_df$families[i.fam]) %>% filter(pValue < 0.05) %>% nrow()
    }
    
    x = final_df %>% arrange(-annotate)
  }

  # for family results
  if (!is.null(Family)) {
  x = merged %>% filter(Transposon_Family == Family) %>% select(-padj) %>% arrange(pValue)
  }
  
  return(x)
}
######################################

Copia = grouped_families(superFamily = "LTR/Copia")
Gypsy = grouped_families(superFamily = "LTR/Gypsy")
LINE = grouped_families(superFamily = "LINE/L1")

ATHILA2 = grouped_families(Family = "ATHILA2")
ATLINE1_6 = grouped_families(Family = "ATLINE1_6")
TA11 = grouped_families(Family = "TA11")
ATHILA2 = grouped_families(Family = "ATHILA2")
