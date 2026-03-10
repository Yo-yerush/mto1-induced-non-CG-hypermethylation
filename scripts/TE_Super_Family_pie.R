library(dplyr)
library(RColorBrewer)

edit_TE_fam_names <- function(x) {
    x$Transposon_Super_Family <- gsub("LTR/", "", x$Transposon_Super_Family)
    x$Transposon_Super_Family <- gsub("RC/", "", x$Transposon_Super_Family)
    x$Transposon_Super_Family <- gsub("LINE.*", "LINE", x$Transposon_Super_Family)
    x$Transposon_Super_Family <- gsub("RathE.*", "SINE", x$Transposon_Super_Family)
    x$Transposon_Super_Family <- gsub("DNA.*", "TIR", x$Transposon_Super_Family)
    x <- x %>% filter(Transposon_Super_Family != "Unassigned")
    return(x)
}

# color index
color_vec <- c(brewer.pal(n = 6, name = "Set2"))
color_vec <- paste0(color_vec, "90")
col_indx <- data.frame(
    SF_name = c("Copia", "Gypsy", "LINE", "SINE", "Helitron", "TIR"),
    col = color_vec
)


for (context in c("CG", "CHG", "CHH")) {
    DMRsReplicates_TE_file <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Transposable_Elements_", context, "_genom_annotations.csv")) %>%
        edit_TE_fam_names()

    ######################### TE pie plot
    TE_Freq <- as.data.frame(table(DMRsReplicates_TE_file$Transposon_Super_Family))
    names(TE_Freq)[1] <- "SF_name"
    TE_Freq <- merge.data.frame(TE_Freq, col_indx, by = "SF_name")

    # display top TE-SuperFamily labels
    top_n <- 6
    labels <- rep("", nrow(TE_Freq))
    indices_top_n <- order(TE_Freq$Freq, decreasing = TRUE)[1:top_n]
    labels[indices_top_n] <- as.character(TE_Freq$SF_name[indices_top_n])

    if (context != "CHH") {
        labels <- gsub("^SINE$", "", labels)
    }


    svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/", context, "_TE_Super_Family_pie.svg"), width = 2.65, height = 2.65, family = "serif")
    par(mar = c(1, 1, 1, 1))
    # par(fig = c(0, 6, 0, 10) / 10)
    par(lwd = 2)
    # pie plot
    pie(TE_Freq$Freq,
        labels = labels,
        # main = paste0("Transposon Super Family:\nDMRs in ",context," context"),
        main = "",
        border = "white",
        col = TE_Freq$col
    )
    # outer circle
    symbols(0, 0,
        circles = 0.8, inches = FALSE, add = TRUE,
        fg = "#575652"
    )
    # title
    # title(
    #    main = paste0("Transposon Super Family:\nDMRs in ", context, " context"),
    #    line = -3, cex.main = 1.25
    # )

    dev.off()
}

# legend
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/legend_TE_Super_Family_pie.svg"), width = 6, height = 10, family = "serif")
plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE)
par(mar = rep(0, 4))
par(lwd = 1.25)
legend("top", legend = as.character(TE_Freq$SF_name), fill = TE_Freq$col, bty = "n")
dev.off()




######################### all TEs distribution (TAIR10)
### upload and edit TE file
TE_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt",
    sep = "\t"
) %>%
    mutate(seqnames = NA) %>% # Add a new column with NA values
    dplyr::select(seqnames, Transposon_min_Start, Transposon_max_End, orientation_is_5prime, everything()) %>%
    dplyr::rename(Derives_from = Transposon_Name) %>%
    edit_TE_fam_names()

for (i in 1:5) {
    TE_file$seqnames[grep(paste0("AT", i, "TE"), TE_file$Derives_from)] <- paste0("Chr", i)
}
TE_file$orientation_is_5prime <- gsub("true", "+", TE_file$orientation_is_5prime)
TE_file$orientation_is_5prime <- gsub("false", "-", TE_file$orientation_is_5prime)
names(TE_file)[1:4] <- c("seqnames", "start", "end", "strand")

### TE pie plot
TE_Freq <- as.data.frame(table(TE_file$Transposon_Super_Family))
names(TE_Freq)[1] <- "SF_name"
TE_Freq <- merge.data.frame(TE_Freq, col_indx, by = "SF_name")

# display top TE-SuperFamily labels
top_n <- 6
labels <- rep("", nrow(TE_Freq))
indices_top_n <- order(TE_Freq$Freq, decreasing = TRUE)[1:top_n]
labels[indices_top_n] <- as.character(TE_Freq$SF_name[indices_top_n])
labels <- gsub("^SINE$", "\nSINE", labels)

svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/total_TE_Super_Family_pie.svg"), width = 2.65, height = 2.65, family = "serif")
par(mar = c(1, 1, 1, 1))
par(lwd = 2)
pie(TE_Freq$Freq,
    labels = labels,
    main = "",
    border = "white",
    col = TE_Freq$col
)
# outer circle
symbols(0, 0,
    circles = 0.8, inches = FALSE, add = TRUE,
    fg = "#575652"
)

dev.off()
