library(ggplot2)
library(dplyr)
library(GenomicRanges)

output_dir <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/mutants_figs"

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/windowSize_for_GRanges_mcol.r")
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/ChrPlots_CX_yo.R")

read_mut_file <- function(mut_name) {
    cat(paste0("read ", mut_name, " file..."))
    x <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/mutants_figs/", mut_name, "_delta_df.csv.gz"))

    return(list(
        cg = x[x$context == "CG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta"),
        chg = x[x$context == "CHG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta"),
        chh = x[x$context == "CHH", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta")
    ))
    cat(" done\n")
}

mto1 <- read_mut_file("mto1")
met1 <- read_mut_file("met1")
cmt2 <- read_mut_file("cmt2")
cmt3 <- read_mut_file("cmt3")
ddm1 <- read_mut_file("ddm1")

tiff(paste0(output_dir, "/ChrPlot_test_stroud_all.tif"), width = 7, height = 4, units = "in", res = 600, family = "serif")
ChrPlots_CX_all(
    meth_var_list = list(mto1,met1,cmt2,cmt3,ddm1), # "test_stroud_all",
    meth_names = c("mto1","met1","cmt2","cmt3","ddm1"),
    y_max_cg = 0,
    y_max_chg = 0.2,
    y_max_chh = 0.05,
    y_mid_cg = NULL,
    y_mid_chg = 0,
    y_mid_chh = 0,
    y_min_cg = -1,
    y_min_chg = -0.5,
    y_min_chh = -0.15,
    italic_legend_names = TRUE,
    ylab_suffix = "(Δ)",
    y_title_cex = 1.1,
    TE_as_gr = "tair10"
    )
dev.off()


