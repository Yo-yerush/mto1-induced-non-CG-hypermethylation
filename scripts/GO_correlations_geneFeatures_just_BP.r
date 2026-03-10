library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(GenomicRanges)
library(topGO)

##################
just_BP = TRUE
#################



RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
    dplyr::select(gene_id, DEG_log2FC, padj, pValue)

feature_file_fun <- function(context) {
    feature_file <- data.frame()
    for (ann in c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs", "TEG")) {
        # DMR results file
        ann_DMRs <- read.csv(paste0(
            "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/",
            context, "/", ann, "_", context, "_genom_annotations.csv"
        )) %>%
            dplyr::select(
                gene_id, log2FC, context, type, Symbol, Computational_description
            ) %>%
            dplyr::rename(DMR_log2FC = log2FC)

        # correlation results file
        ann_corr <- read.csv(paste0(
            "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/Gene_feature/mto1/", context, "/", ann, ".corr.", context, ".mto1.csv"
        )) %>%
            dplyr::select(-padj) %>%
            dplyr::rename(cor_pValue = pval)

        ann_merged <- merge.data.frame(ann_corr, ann_DMRs, by = "gene_id", all.y = T)

        feature_file <- rbind(feature_file, ann_merged)
    }
    return(feature_file)
}

################
remove_dup_DMR <- function(y) {
    y <- as.character(unique(unlist(strsplit(y, ","))))
    paste(y, collapse = ",")
}
################

################
feture_res_table <- rbind(
    feature_file_fun("CG"),
    feature_file_fun("CHG"),
    feature_file_fun("CHH")
) %>%
    mutate(
        DMR_log2FC = round(DMR_log2FC, 2),
        CG_DMRs = ifelse(grepl("CG", context), DMR_log2FC, NA),
        CHG_DMRs = ifelse(grepl("CHG", context), DMR_log2FC, NA),
        CHH_DMRs = ifelse(grepl("CHH", context), DMR_log2FC, NA)
    ) %>%
    merge.data.frame(RNAseq, ., by = "gene_id", all.y = TRUE) %>%
    mutate(tmp = paste(gene_id, type, sep = "_")) %>%
    group_by(tmp) %>%
    summarise(
        across(contains("CG_DMRs") | contains("CHG_DMRs") | contains("CHH_DMRs"), ~ remove_dup_DMR(paste(., collapse = ","))), # apply remove_dup_DMR function for DMR columns
        across(!contains("CG_DMRs") | contains("CHG_DMRs") | contains("CHH_DMRs"), dplyr::first) # for other columns
    ) %>%
    as.data.frame() %>%
    mutate(across(contains("_DMRs"), ~ gsub("NA", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",,", ",", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub("^,", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",$", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",", ", ", .))) %>%
    dplyr::relocate(CG_DMRs, CHG_DMRs, CHH_DMRs, .before = context) %>%
    dplyr::relocate(type, cor, cor_pValue, .after = gene_id) %>%
    dplyr::select(-context, -DMR_log2FC, -tmp) %>%
    # filter(pValue < 0.05) %>%
    arrange(pValue) %>%
    arrange(type) %>%
    mutate(across(contains("padj") | contains("pValue"), ~ gsub(" NA", NA, .)),
        pValue = as.numeric(formatC(.$pValue, format = "e", digits = 2)),
        padj = as.numeric(formatC(.$padj, format = "e", digits = 2)),
        DEG_log2FC = round(DEG_log2FC, 3),
        cor = round(cor, 3),
        cor_pValue = as.numeric(formatC(.$cor_pValue, format = "e", digits = 2))
    ) %>%
    ### add to this script:
    dplyr::select(gene_id, pValue, cor_pValue, DEG_log2FC) %>%
    distinct(gene_id, .keep_all = T)


feture_res_table$pValue[is.na(feture_res_table$pValue)] <- 0.999
feture_res_table$pValue[feture_res_table$pValue == 0] <- 1e-300



direction_fun <- function(x, gainORloss) {
    if (gainORloss == "gain") {
        geneList <- ifelse(x$pValue < 0.05 & x$cor_pValue < 0.05 & x$DEG_log2FC > 0, 1, 0) # filter by pValue, cor_pValue and up/down regulated
        names(geneList) <- x$gene_id
    } else if (gainORloss == "loss") {
        geneList <- ifelse(x$pValue < 0.05 & x$cor_pValue < 0.05 & x$DEG_log2FC < 0, 1, 0) # filter by pValue, cor_pValue and up/down regulated
        names(geneList) <- x$gene_id
    }

    res_list <- list()
    for (GO_type_loop in c("BP", "MF", "CC")) {
        myGOdata <- new("topGOdata",
            ontology = GO_type_loop,
            allGenes = geneList,
            geneSelectionFun = function(x) (x == 1),
            # description = "Test",
            annot = annFUN.org,
            # nodeSize = 5,
            mapping = "org.At.tair.db"
        )

        sg <- sigGenes(myGOdata)
        str(sg)
        numSigGenes(myGOdata)

        resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")


        allRes <- GenTable(myGOdata,
            Fisher = resultFisher,
            orderBy = "Fisher", ranksOf = "Fisher", topNodes = length(resultFisher@score)
        )
        allRes$Fisher <- as.numeric(allRes$Fisher)
        allRes$Term <- gsub(",", ";", allRes$Term)
        allRes$type <- GO_type_loop

        res_list[[GO_type_loop]] <- allRes[allRes$Fisher <= 0.01, ]
    }

    res_bind <- rbind(res_list[["BP"]], res_list[["CC"]], res_list[["MF"]])
    res_bind <- res_bind[!grepl("cellular_component|biological_process|molecular_function", res_bind$Term), ]

    return(res_bind)
}
### plot
gain_bind <- direction_fun(feture_res_table, "gain")
loss_bind <- direction_fun(feture_res_table, "loss")

################
### just BP
if (just_BP) {
    gain_bind = gain_bind %>% filter(type == "BP") %>% mutate(direction = "Up-regulated")
    loss_bind = loss_bind %>% filter(type == "BP") %>% mutate(direction = "Down-regulated")

    gain_bind = gain_bind %>% filter(Term != "biosynthetic process")
}
file_name_addition = ifelse(just_BP, "_just_BP", "")
################
# dot color normalization
gain_col <- "#cf534c"
loss_col <- "#6397eb"

pVal_max = 0.01
pVal_min = min(c(gain_bind$Fisher, loss_bind$Fisher))

gain_bind <- gain_bind %>%
    arrange(Fisher) %>%
    mutate(
        color = {
            scaled <- (Fisher - pVal_min) / (pVal_max - pVal_min)
            scaled[scaled < 0] <- 0
            scaled[scaled > 1] <- 1
            rgb(colorRamp(c(gain_col, "black"))(scaled), maxColorValue = 255)
        }
    )

loss_bind <- loss_bind %>%
    arrange(Fisher) %>%
    mutate(
        color = {
            scaled <- (Fisher - pVal_min) / (pVal_max - pVal_min)
            scaled[scaled < 0] <- 0
            scaled[scaled > 1] <- 1
            rgb(colorRamp(c(loss_col, "black"))(scaled), maxColorValue = 255)
        }
    )

##########
### dot size normalization
normalize_fun <- function(x, min_val, max_val) {
    return((x - min(x)) / (max(x) - min(x)) * (max_val - min_val) + min_val)
}

annotate_max = max(c(gain_bind$Annotated, loss_bind$Annotated))
annotate_min = min(c(gain_bind$Annotated, loss_bind$Annotated))

ann_size_index = data.frame(
    original_row = c(1:nrow(gain_bind), 1:nrow(loss_bind)),
    direction = c(gain_bind$direction, loss_bind$direction),
    Annotated = c(gain_bind$Annotated, loss_bind$Annotated),
    ann_indx = normalize_fun(c(gain_bind$Annotated, loss_bind$Annotated), 1, 6)
)

up_size = ann_size_index[ann_size_index$direction == "Up-regulated", "ann_indx"]
down_size = ann_size_index[ann_size_index$direction == "Down-regulated", "ann_indx"]
#########

bubble_up <- gain_bind %>%
    ggplot(aes(Significant, reorder(Term, Significant))) + #size = Annotated
    #scale_color_gradient("p.value", low = gain_col, high = "black") +
    labs(x = "Significant", y = "") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_point(size = up_size, color = gain_bind$color) +
    facet_grid(rows = vars(direction), scales = "free_y", space = "free_y") +
    guides(color = guide_colorbar(order = 1, barheight = 4))

bubble_down <- loss_bind %>%
    ggplot(aes(Significant, reorder(Term, Significant))) + #size = Annotated
    #scale_color_gradient("p.value", low = loss_col, high = "black") + # theme_classic() +
    labs(x = "Significant", y = "") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_point(size = down_size, color = loss_bind$color) +
    facet_grid(rows = vars(direction), scales = "free_y", space = "free_y") +
    guides(color = guide_colorbar(order = 1, barheight = 4))


Height = 2.5
Width = 9.5
source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/multiplot_ggplot2.R")

# upregulated plot
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_correlation_up", file_name_addition, "_geneFeature.svg"), width = 4.25, height = 1.65, family = "serif")
bubble_up
dev.off()

# downregulated plot
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_correlation_down", file_name_addition, "_geneFeature.svg"), width = 4.25, height = 2.5, family = "serif")
bubble_down
dev.off()

# annotated legend
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_correlation_Annotated_legend.svg"), width = 1, height = 1.5, family = "serif")
legend_plot <- ggplot(ann_size_index, aes(x = 1, y = 1, size = Annotated)) +
    geom_point() +
    scale_size(
        range = c(1, 6),
        name = " ",
        breaks = c(min(ann_size_index$Annotated), 50, 100, max(ann_size_index$Annotated))
    ) +
    theme_void() +
    theme(legend.position = "right")

legend_only <- cowplot::get_legend(legend_plot)
grid.newpage()
grid.draw(legend_only)
dev.off()

## pValue legend
#legend_up <- ggplot(data.frame(x = c(pVal_min, 0.01)), aes(x = x, y = x, color = x)) +
#    geom_point(size = 5) +
#    scale_color_gradient(
#        limits = c(pVal_min, 0.01),
#        breaks = c(pVal_min, 0.01),
#        labels = c("",""),
#        low = gain_col,
#        high = "black",
#        name = "p-Value"
#    ) +
#    theme_void() +
#    theme(legend.position = "right")
#
#legend_down <- ggplot(data.frame(x = c(pVal_min, 0.01)), aes(x = x, y = x, color = x)) +
#    geom_point(size = 5) +
#    scale_color_gradient(
#        limits = c(pVal_min, 0.01),
#        breaks = c(pVal_min, 0.01),
#        labels = c(pVal_min, "0.01"),
#        low = loss_col,
#        high = "black",
#        name = ""
#    ) +
#    theme_void() +
#    theme(legend.position = "right")
#cowplot::get_legend(legend_down)
#svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_correlation_pValue_up_legend.svg"), width = 1.25, height = 1.75, family = "serif")
#grid.draw(cowplot::get_legend(legend_up))
#dev.off()
#
#svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_correlation_pValue_down_legend.svg"), width = 1.25, height = 1.75, family = "serif")
#grid.draw(cowplot::get_legend(legend_down))
#dev.off()
#
#
#colfunc_up <- colorRampPalette(c("black", gain_col))
#colfunc_down <- colorRampPalette(c("black", loss_col))
#
#svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_correlation_pValue_legend.svg"), width = 5, height = 5, family = "serif")
#par(mar = c(1, 1, 1, 1))
#plot.new()
#legend_image_up <- as.raster(matrix(colfunc_up(50), ncol = 1))
#legend_image_down <- as.raster(matrix(colfunc_down(50), ncol = 1))
#
#rasterImage(legend_image_up, 0, 0.86, 0.025, 0.98) # xleft, ybottom, xright, ytop
#rect(0, 0.86, 0.025, 0.98, border = "black") # box for scale
#
#rasterImage(legend_image_down, 0.035, 0.86, 0.06, 0.98) # xleft, ybottom, xright, ytop
#rect(0.035, 0.86, 0.06, 0.98, border = "black") # box for scale
#
#text(x = 0.035, y = 1, labels = expression(bold("log2FC")), cex = 0.75)
#text(x = 0.07, y = c(0.88, 0.965), labels = c(expression(bold("-1")), expression(bold(" 1"))), cex = rep(0.75, 2))
#
#dev.off()
#