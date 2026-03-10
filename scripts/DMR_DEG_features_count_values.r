library(dplyr)

##################### DEGs IDs
DEGs <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(padj < 0.05) %>%
    rename(gene_id = locus_tag) %>%
    distinct(gene_id)

##################### DMRs IDs
read_DMRs <- function(f) {
    x <- NULL
    for (c in c("CG", "CHG", "CHH")) {
        x <- rbind(x, read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", c, "/", f, "_", c, "_genom_annotations.csv")))
    }
    xx = distinct(x, gene_id)
    n.total =  nrow(xx)
    n.feature = merge.data.frame(xx, DEGs, by = "gene_id") %>% distinct(gene_id) %>% nrow()

    xxx = round((n.feature / n.total) * 100, 1)
    return(xxx)
}

for (feature in c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs")) {
   message("Percentage of overlapped ", feature, ": ", read_DMRs(feature))
}


