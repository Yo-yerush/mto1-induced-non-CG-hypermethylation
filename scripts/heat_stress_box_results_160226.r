library(ggplot2)
library(ggbreak)
library(multcompView)
library(dplyr)
library(tidyr)
library(grid)

# Read the CSV file
heat_res <- read.csv("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/heat_sress_40_data_extracted_from_rons_file.csv") %>%
pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
filter(!grepl("mto3", variable)) %>%
mutate(
    variable = gsub("WT", "wt", variable),
    variable = gsub("\\.$", "", variable),
    variable = gsub("\\.+rt", "_c", variable),
    variable = gsub("\\.+heat", "_h", variable)
) %>%
as.data.frame() %>%
rename(name = variable)


# Set factor levels for desired order
heat_res$name <- factor(heat_res$name, levels = c("wt_c", "wt_h", "mto1_c", "mto1_h"))
heat_res$genotype <- factor(ifelse(grepl("wt", heat_res$name), "WT", "mto1"), levels = c("WT", "mto1"))

### use Ron's stat ### # ANOVA + Tukey post-hoc test
### use Ron's stat ### aov_res <- aov(value ~ name * genotype, data = heat_res)
### use Ron's stat ### tukey_res <- TukeyHSD(aov_res)
### use Ron's stat ### print(summary(aov_res))
### use Ron's stat ### print(tukey_res)
### use Ron's stat ### 
### use Ron's stat ### # Extract compact letter display (CLD)
### use Ron's stat ### tukey_pvals <- tukey_res$name[, "p adj"]
### use Ron's stat ### cld <- multcompLetters(tukey_pvals)$Letters
### use Ron's stat ### 
### use Ron's stat ### # Create label dataframe with y positions for each group
### use Ron's stat ### # (adjust y_pos values to position letters above each box)
### use Ron's stat ### label_df <- data.frame(
### use Ron's stat ###     name = factor(names(cld), levels = c("wt_c", "wt_h", "mto1_c", "mto1_h")),
### use Ron's stat ###     letter = toupper(as.character(cld)),
### use Ron's stat ###     y_pos = c(
### use Ron's stat ###         "wt_c" = max(heat_res[heat_res$name == "wt_c", "value"], na.rm = T)*1.75,
### use Ron's stat ###         "wt_h" = max(heat_res[heat_res$name == "wt_h", "value"], na.rm = T)*1.2,
### use Ron's stat ###         "mto1_c" = max(heat_res[heat_res$name == "mto1_c", "value"], na.rm = T)*1.25,
### use Ron's stat ###         "mto1_h" = max(heat_res[heat_res$name == "mto1_h", "value"], na.rm = T)*1.1
### use Ron's stat ###     )[names(cld)]
### use Ron's stat ### )
label_df <- data.frame(
    name = factor(c("wt_c", "wt_h", "mto1_c", "mto1_h"), levels = c("wt_c", "wt_h", "mto1_c", "mto1_h")),
    letter = c("B", "A", "B", "C"),
    y_pos = c(
        "wt_c" = max(heat_res[heat_res$name == "wt_c", "value"], na.rm = T)*1.85,
        "wt_h" = max(heat_res[heat_res$name == "wt_h", "value"], na.rm = T)*1.225,
        "mto1_c" = max(heat_res[heat_res$name == "mto1_c", "value"], na.rm = T)*1.35,
        "mto1_h" = max(heat_res[heat_res$name == "mto1_h", "value"], na.rm = T)*1.125
    )[c("wt_c", "wt_h", "mto1_c", "mto1_h")]
)


{
svg("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/heat_40_from_ron_results_160226.svg", width = 2.85, height = 3, family = "serif")
p <- ggplot(heat_res, aes(x = name, y = value, fill = name)) +
    geom_boxplot(data = heat_res,
                 alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = genotype), width = 0.15, size = 2, alpha = 0.8) +
    geom_vline(xintercept = 2.5, color = "gray", linetype = "dashed", linewidth = 0.5) +
    geom_text(data = label_df, aes(x = name, y = y_pos, label = letter),
              inherit.aes = FALSE, size = 4.5, fontface = "bold", color = "gray35") +
    # scale_fill_manual(values = c("wt_c" = "gray50", "wt_h" = "gray80",
    #                               "mto1_c" = "#eed3aa", "mto1_h" = "#eed3aa")) +
    scale_fill_manual(values = c("wt_c" = "white", "wt_h" = "white",
                                  "mto1_c" = "white", "mto1_h" = "white")) +
    scale_color_manual(values = c("WT" = "gray20", "mto1" = "#b6801d"),
                       labels = c("WT" = "WT", "mto1" = expression(italic("mto1")))) +
    guides(fill = "none", color = guide_legend(title = NULL, override.aes = list(size = 3))) +
    scale_x_discrete(labels = c("wt_c" = "RT", "wt_h" = "40°C",
                                "mto1_c" = "RT", "mto1_h" = "40°C")) +
    labs(x = NULL, y = "Total radiant efficiency (%)") +
    # labs(x = NULL, y = "Total radiant efficiency (% of control)") +
    scale_y_continuous(limits = c(0, 2200)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 13),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          plot.title = element_blank(),
          plot.margin = margin(5, 5, 20, 5))
print(p)

# Add legend inside the plot (top-right) using grid
grid.points(x = unit(0.32, "npc"), y = unit(0.9, "npc"),
            pch = 19, size = unit(0.55, "char"),
            gp = gpar(col = "gray35"))
grid.text("WT", x = unit(0.36, "npc"), y = unit(0.9, "npc"),
          just = "left", gp = gpar(fontsize = 12, fontfamily = "serif"))
grid.points(x = unit(0.32, "npc"), y = unit(0.85, "npc"),
            pch = 19, size = unit(0.55, "char"),
            gp = gpar(col = "#b6801d"))
grid.text(expression(italic("mto1")), x = unit(0.36, "npc"), y = unit(0.85, "npc"),
          just = "left", gp = gpar(fontsize = 12, fontfamily = "serif"))
dev.off()
}
