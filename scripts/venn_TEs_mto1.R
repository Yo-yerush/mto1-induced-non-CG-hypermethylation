library(dplyr)
library(VennDiagram)

CG = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/Transposable_Elements_CG_genom_annotations.csv")
CHG = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/Transposable_Elements_CHG_genom_annotations.csv")
CHH = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/Transposable_Elements_CHH_genom_annotations.csv")

TE_sf = "LTR/Gypsy"

cntx_sf <- function(x) {
  x = x %>%
    filter(Transposon_Super_Family == TE_sf) %>%
    mutate(loc = paste(seqnames,start,end, sep = "_")) %>% # to merge by position in the Venn-Diagram
    select(Transposon_Name, Transposon_Super_Family, loc)
}

CG_sf = cntx_sf(CG)
CHG_sf = cntx_sf(CHG)
CHH_sf = cntx_sf(CHH)




CG_C = merge.data.frame(CG_sf, CHG_sf, by = "loc")
merged = merge.data.frame(CHG_sf, CHH, by = "loc")
merged = merge.data.frame(merged, CHH, by = "loc")


#fill = c("#440154ff", '#21908dff', '#fde725ff')
venn_colors = c("red", "blue", "green")

venn.diagram(
  x = list(CG_sf$loc, CHG_sf$loc, CHH_sf$loc),
  category.names = c("CG", "CHG", "CHH"),
  filename = "C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEs_superFamily_VennDiagram.png",
  disable.logging = T,
  output = TRUE ,
  imagetype="png" ,
  height = 400 , 
  width = 400 , 
  resolution = 250,
  compression = "lzw",
  lwd = 1,
  fill = venn_colors,
  #alpha = rep(0.3, length(x)),
  cex = 0.5,
  fontfamily = "serif",
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.pos = c(0,0,0),
  cat.fontface = 2,
  cat.fontfamily = "sans"
  #    cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  #    col=venn_colors,
  #    rotation = 1
)
