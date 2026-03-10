library(ggplot2)
library(ggbreak)
library(multcompView)
library(dplyr)
library(tidyr)
library(grid)

# Read the CSV file
EL_res <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/EL_sress_data_extracted_from_rons_file.csv") %>%
pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
filter(!grepl("mto3", variable)) %>%
mutate(
    variable = gsub("WT", "wt", variable),
    variable = gsub("\\.$", "", variable),
    variable = gsub("\\.+control", "_c", variable),
    variable = gsub("\\.+el", "_el", variable)
) %>%
as.data.frame() %>%
rename(name = variable)


# Set factor levels for desired order
EL_res$name <- factor(EL_res$name, levels = c("wt_c", "wt_el", "mto1_c", "mto1_el"))
EL_res$genotype <- factor(ifelse(grepl("wt", EL_res$name), "WT", "mto1"), levels = c("WT", "mto1"))

# ANOVA + Tukey post-hoc test
aov_res <- aov(value ~ name * genotype, data = EL_res)
tukey_res <- TukeyHSD(aov_res)
print(summary(aov_res))
print(tukey_res)

# Extract compact letter display (CLD)
tukey_pvals <- tukey_res$name[, "p adj"]
cld <- multcompLetters(tukey_pvals)$Letters

# Create label dataframe with y positions for each group
# (adjust y_pos values to position letters above each box)
label_df <- data.frame(
    name = factor(names(cld), levels = c("wt_c", "wt_el", "mto1_c", "mto1_el")),
    letter = toupper(as.character(cld)),
    y_pos = c(
        "wt_c" = max(EL_res[EL_res$name == "wt_c", "value"], na.rm = T)*1.3,
        "wt_el" = max(EL_res[EL_res$name == "wt_el", "value"])*1.125,
        "mto1_c" = max(EL_res[EL_res$name == "mto1_c", "value"], na.rm = T)*1.3,
        "mto1_el" = max(EL_res[EL_res$name == "mto1_el", "value"], na.rm = T)*1.1
    )[names(cld)]
)

{
svg("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/EL_from_ron_results_160226.svg", width = 2.85, height = 3, family = "serif")
p <- ggplot(EL_res, aes(x = name, y = value, fill = name)) +
    geom_boxplot(data = EL_res,
                 alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = genotype), width = 0.15, size = 2, alpha = 0.8) +
    geom_vline(xintercept = 2.5, color = "gray", linetype = "dashed", linewidth = 0.5) +
    geom_text(data = label_df, aes(x = name, y = y_pos, label = letter),
              inherit.aes = FALSE, size = 4.5, fontface = "bold", color = "gray35") +
    # scale_fill_manual(values = c("wt_c" = "gray50", "wt_el" = "gray80",
    #                               "mto1_c" = "#eed3aa", "mto1_el" = "#eed3aa")) +
    scale_fill_manual(values = c("wt_c" = "white", "wt_el" = "white",
                                  "mto1_c" = "white", "mto1_el" = "white")) +
    scale_color_manual(values = c("WT" = "gray20", "mto1" = "#b6801d"),
                       labels = c("WT" = "WT", "mto1" = expression(italic("mto1")))) +
    guides(fill = "none", color = guide_legend(title = NULL, override.aes = list(size = 3))) +
    scale_x_discrete(labels = c("wt_c" = "Ctrl", "wt_el" = "EL",
                                "mto1_c" = "Ctrl", "mto1_el" = "EL")) +
    labs(x = NULL, y = "Total radiant efficiency (%)") +
    # labs(x = NULL, y = "Total radiant efficiency (% of control)") +
    scale_y_continuous(breaks = c(0, 250, 500, 750, 1000), labels = c("0", "250", "500", "750", "  1000")) +
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
