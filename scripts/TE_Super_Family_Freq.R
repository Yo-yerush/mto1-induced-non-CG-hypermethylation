library(GenomicRanges)
library(dplyr)

TE_Super_Family_Frequency = function(context, TE.gr) {
  
  DMRsReplicates_TE_file.0 = paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/",context,"/Transposable_Elements_",context,"_genom_annotations.csv")
  
  if (file.exists(DMRsReplicates_TE_file.0)) {
    
    ### total, hyper- and hypo- DMRs df
    DMRsReplicates_TE_file = read.csv(DMRsReplicates_TE_file.0)
    DMRsReplicates_TE_up = DMRsReplicates_TE_file[DMRsReplicates_TE_file$regionType == "gain",]
    DMRsReplicates_TE_down = DMRsReplicates_TE_file[DMRsReplicates_TE_file$regionType == "loss",]
    
    
    ### frequency of TE super family overlapped with DMRs
    TE_Freq = as.data.frame(table(DMRsReplicates_TE_file$Transposon_Super_Family)) %>%
      setNames(c("Transposon_Super_Family", "total_DMRs"))
    
    TE_Freq_up = as.data.frame(table(DMRsReplicates_TE_up$Transposon_Super_Family)) %>%
      setNames(c("Transposon_Super_Family", "hyper_DMRs"))
    
    TE_Freq_down = as.data.frame(table(DMRsReplicates_TE_down$Transposon_Super_Family)) %>%
      setNames(c("Transposon_Super_Family", "hypo_DMRs"))
    
    
    ### frequency of TE super family unique IDs
    superFamilies = unique(TE.gr$Transposon_Super_Family)
    TE_uniqueID = data.frame(Transposon_Super_Family = NA, unique_IDs = NA)
    
    for (sf.i in 1:length(superFamilies)) {
      tryCatch({
        sf.unique = DMRsReplicates_TE_file[DMRsReplicates_TE_file$Transposon_Super_Family == superFamilies[sf.i],]
        TE_uniqueID[sf.i,1] = superFamilies[sf.i]
        TE_uniqueID[sf.i,2] = length(unique(sf.unique$Transposon_Name))
        
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
    
    return(print(TE_Freq_df))
    #write.csv(TE_Freq_df, paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/",context,"_TE_Super_Family_Freq.csv"), row.names = F)
  }
}

edit_TE_file <- function(TE_df) {
  
  TE_df <- TE_df %>%
    mutate(seqnames = NA) %>%  # Add a new column with NA values
    dplyr::select(seqnames,Transposon_min_Start,Transposon_max_End,orientation_is_5prime, everything())
  
  for (i in 1:5) {
    TE_df$seqnames[grep(paste0("AT",i,"TE"),TE_df$Transposon_Name)] = paste0("Chr",i)
  }
  TE_df$orientation_is_5prime = gsub("true","+",TE_df$orientation_is_5prime)
  TE_df$orientation_is_5prime = gsub("false","-",TE_df$orientation_is_5prime)
  
  names(TE_df)[1:4] = c("seqnames","start","end","strand")
  TE_gr = makeGRangesFromDataFrame(TE_df, keep.extra.columns = T)
  
  return(TE_gr)
}

TE_file = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")
TE_file = edit_TE_file(TE_file)

TE_Super_Family_Frequency("CG", TE_file)
TE_Super_Family_Frequency("CHG", TE_file)
TE_Super_Family_Frequency("CHH", TE_file)
