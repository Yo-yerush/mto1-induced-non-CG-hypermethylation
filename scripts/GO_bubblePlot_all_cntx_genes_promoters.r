library(ggplot2)
library(dplyr)
library(grid)

file_name_addition = "all_cntx_gene_promoter"

gain_bind <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap/all_cntx_genes_promoters_0.01/GO_all_cntx_genes_promoters_upregulated_DEGs_overlap.csv")

loss_bind <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap/all_cntx_genes_promoters_0.05/GO_all_cntx_genes_promoters_downregulated_DEGs_overlap.csv")


### keep only top10 by significant count (not fisher)
loss_bind <- loss_bind %>%
    mutate(tmp = row_number()) %>%
    arrange(desc(Significant)) %>%
    .[1:10, ] %>%
    arrange(tmp) %>%
    select(-tmp)


### edit Term names
gain_bind$Term[gain_bind$GO.ID == "GO:0045227"] = "capsule polysaccharide biosynthetic process"
gain_bind$Term[gain_bind$GO.ID == "GO:0033499"] = "galactose catabolic process via UDP-galactose"
gain_bind$Term[gain_bind$GO.ID == "GO:0009793"] = "embryo development ending in seed dormancy"

loss_bind$Term[loss_bind$GO.ID == "GO:0006086"] = "pyruvate decarboxylation to acetyl-CoA"
loss_bind$Term[loss_bind$GO.ID == "GO:0006556"] = "S-adenosylmethionine biosynthetic process"
loss_bind$Term[loss_bind$GO.ID == "GO:1904143"] = "positive regulation of carotenoid biosynthetic process"
loss_bind$Term[loss_bind$GO.ID == "GO:0009089"] = "lysine biosynthetic process via diaminopimelate"
loss_bind$Term[loss_bind$GO.ID == "GO:0030705"] = "cytoskeleton-dependent intracellular transport"
loss_bind$Term[loss_bind$GO.ID == "GO:0043069"] = "negative regulation of programmed cell death"



gain_col <- "#cf534c"
loss_col <- "#6397eb"

pVal_max <- 0.01
pVal_min <- min(c(gain_bind$Fisher, loss_bind$Fisher))

gain_bind <- gain_bind %>%
    mutate(direction = "Up-regulated") %>%
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
    mutate(direction = "Down-regulated") %>%
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

annotate_max <- max(c(gain_bind$Annotated, loss_bind$Annotated))
annotate_min <- min(c(gain_bind$Annotated, loss_bind$Annotated))

ann_size_index <- data.frame(
    original_row = c(1:nrow(gain_bind), 1:nrow(loss_bind)),
    direction = c(gain_bind$direction, loss_bind$direction),
    Annotated = c(gain_bind$Annotated, loss_bind$Annotated),
    ann_indx = normalize_fun(c(gain_bind$Annotated, loss_bind$Annotated), 1, 6)
)

up_size <- ann_size_index[ann_size_index$direction == "Up-regulated", "ann_indx"]
down_size <- ann_size_index[ann_size_index$direction == "Down-regulated", "ann_indx"]
#########

bubble_up <- gain_bind %>%
    ggplot(aes(Significant, reorder(Term, Significant))) + # size = Annotated
    # scale_color_gradient("p.value", low = gain_col, high = "black") +
    labs(x = "Significant", y = "") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_point(size = up_size, color = gain_bind$color) +
    facet_grid(rows = vars(direction), scales = "free_y", space = "free_y") +
    guides(color = guide_colorbar(order = 1, barheight = 4))

bubble_down <- loss_bind %>%
    ggplot(aes(Significant, reorder(Term, Significant))) + # size = Annotated
    # scale_color_gradient("p.value", low = loss_col, high = "black") + # theme_classic() +
    labs(x = "Significant", y = "") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_point(size = down_size, color = loss_bind$color) +
    facet_grid(rows = vars(direction), scales = "free_y", space = "free_y") +
    guides(color = guide_colorbar(order = 1, barheight = 4))


# upregulated plot
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap_up_", file_name_addition, "_geneFeature.svg"), width = 4.25, height = 2, family = "serif")
bubble_up
dev.off()

# downregulated plot
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap_down_", file_name_addition, "_geneFeature.svg"), width = 4.25, height = 2.25, family = "serif")
bubble_down
dev.off()

# annotated legend
svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_overlap_Annotated_legend_", file_name_addition, ".svg"), width = 1, height = 1.5, family = "serif")
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