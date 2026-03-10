library(GenomicRanges)
library(dplyr)

edit_TE_fam_names <- function(x) {
    x$Transposon_Super_Family <- gsub("LTR/", "", x$Transposon_Super_Family)
    x$Transposon_Super_Family <- gsub("RC/", "", x$Transposon_Super_Family)
    x$Transposon_Super_Family <- gsub("LINE.*", "LINE", x$Transposon_Super_Family)
    x$Transposon_Super_Family <- gsub("RathE.*", "SINE", x$Transposon_Super_Family)
    x$Transposon_Super_Family <- gsub("DNA.*", "TIR", x$Transposon_Super_Family)
    x <- x %>% filter(Transposon_Super_Family != "Unassigned")
    return(x)
}

TE_Super_Family_Frequency = function(context) {
  
  DMRsReplicates_TE_file.0 = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/",context,"/Transposable_Elements_",context,"_genom_annotations.csv")
  
  if (file.exists(DMRsReplicates_TE_file.0)) {
    
    ### total, hyper- and hypo- DMRs df
    DMRsReplicates_TE_file <- read.csv(DMRsReplicates_TE_file.0) %>% edit_TE_fam_names()
    DMRsReplicates_TE_up <- DMRsReplicates_TE_file[DMRsReplicates_TE_file$regionType == "gain", ]
    DMRsReplicates_TE_down <- DMRsReplicates_TE_file[DMRsReplicates_TE_file$regionType == "loss", ]
    
    
    ### frequency of TE super family overlapped with DMRs
    TE_Freq = as.data.frame(table(DMRsReplicates_TE_file$Transposon_Super_Family)) %>%
      setNames(c("Transposon_Super_Family", "total_DMRs"))
    
    TE_Freq_up = as.data.frame(table(DMRsReplicates_TE_up$Transposon_Super_Family)) %>%
      setNames(c("Transposon_Super_Family", "hyper_DMRs"))
    
    TE_Freq_down = as.data.frame(table(DMRsReplicates_TE_down$Transposon_Super_Family)) %>%
      setNames(c("Transposon_Super_Family", "hypo_DMRs"))
    
    
    ### frequency of TE super family unique IDs
    superFamilies = unique(DMRsReplicates_TE_file$Transposon_Super_Family)
    TE_uniqueID = data.frame(Transposon_Super_Family = NA, unique_IDs = NA)
    
    for (sf.i in 1:length(superFamilies)) {
      tryCatch({
        sf.unique = DMRsReplicates_TE_file[DMRsReplicates_TE_file$Transposon_Super_Family == superFamilies[sf.i],]
        TE_uniqueID[sf.i,1] = superFamilies[sf.i]
        TE_uniqueID[sf.i,2] = length(unique(sf.unique$gene_id))
        
      }, error = function(cond) {
        TE_uniqueID[sf.i,1] = superFamilies[sf.i]
        TE_uniqueID[sf.i,2] = 0
      })
    }
    
    
    TE_Freq_df = merge(TE_uniqueID, TE_Freq, by = "Transposon_Super_Family", all = T)
    TE_Freq_df = merge(TE_Freq_df, TE_Freq_up, by = "Transposon_Super_Family", all = T)
    TE_Freq_df = merge(TE_Freq_df, TE_Freq_down, by = "Transposon_Super_Family", all = T) %>%
      arrange(desc(total_DMRs))
    TE_Freq_df[is.na(TE_Freq_df)] = 0
    
    write.csv(TE_Freq_df, paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/",context,"_TE_Super_Family_Freq_new.csv"), row.names = F)
  }
}

TE_Super_Family_Frequency("CG")
TE_Super_Family_Frequency("CHG")
TE_Super_Family_Frequency("CHH")
