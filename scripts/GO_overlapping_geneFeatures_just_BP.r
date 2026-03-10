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
    dplyr::select(gene_id, DEG_log2FC, padj) %>%
    #filter(padj < 0.05) %>%
    dplyr::rename(pValue = padj)

feature_file_fun <- function(context.f, ann.f) {
    feature_file <- data.frame()

    for (context in context.f) {
        for (ann in ann.f) {
            # DMR results file
            ann_DMRs <- read.csv(paste0(
                "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/",
                context, "/", ann, "_", context, "_genom_annotations.csv"
            )) %>%
                dplyr::select(
                    gene_id, log2FC, context, type
                ) %>%
                dplyr::rename(DMR_log2FC = log2FC)

            feature_file <- rbind(feature_file, ann_DMRs)
        }
    }
    return(feature_file)
}

################
merged_list = list(
all_res = feature_file_fun(c("CG", "CHG", "CHH"), c("Genes", "Promoters")) %>%
    merge(RNAseq, ., by = "gene_id", all.x = T) %>%
    distinct(gene_id, .keep_all = T),

CG_genes = feature_file_fun("CG", "Genes") %>%
    merge(RNAseq, ., by = "gene_id", all.x = T) %>%
    distinct(gene_id, .keep_all = T),

non_CG_promoters = feature_file_fun(c("CHG", "CHH"), "Promoters") %>%
    merge(RNAseq, ., by = "gene_id", all.x = T) %>%
    distinct(gene_id, .keep_all = T)
)


for(ii in names(merged_list)) {

dir = ii
feture_res_table = merged_list[[ii]]

feture_res_table$pValue[is.na(feture_res_table$pValue)] <- 0.999
feture_res_table$pValue[feture_res_table$pValue == 0] <- 1e-300



direction_fun <- function(x, gainORloss) {
    if (gainORloss == "gain") {
        geneList <- ifelse(x$pValue < 0.05 & x$DEG_log2FC > 0, 1, 0) # filter by pValue, cor_pValue and up/down regulated
        names(geneList) <- x$gene_id
    } else if (gainORloss == "loss") {
        geneList <- ifelse(x$pValue < 0.05 & x$DEG_log2FC < 0, 1, 0) # filter by pValue, cor_pValue and up/down regulated
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


Height_up <- max(nrow(gain_bind)) / 6.25
Height_down <- max(nrow(loss_bind)) / 6.25
if (Height_up < 1.5) {Height_up = 1.5}
if (Height_down < 1.5) {Height_down = 1.5}

source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/multiplot_ggplot2.R")

dir.create(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap/", dir), showWarnings = F)

write.csv(rbind(gain_bind, loss_bind), (paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap/", dir, "/GO_overlapping_",dir,".csv")), row.names = F)

# upregulated plot
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap/", dir, "/GO_overlapping_up", file_name_addition, "_geneFeature.svg"), width = 4.25, height = Height_up, family = "serif")
print(bubble_up)
dev.off()

# downregulated plot
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap/", dir, "/GO_overlapping_down", file_name_addition, "_geneFeature.svg"), width = 4.25, height = Height_down, family = "serif")
print(bubble_down)
dev.off()

# annotated legend
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap/", dir, "/GO_overlapping_Annotated_legend.svg"), width = 1, height = 1.5, family = "serif")
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

}

