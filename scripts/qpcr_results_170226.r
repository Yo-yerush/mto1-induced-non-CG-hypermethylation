pcr_res <- data.frame(
    name = c(
        "wt_c",
        "wt_c",
        "wt_c",
        "wt_c",
        "wt_h",
        "wt_h",
        "wt_h",
        "mto1_c",
        "mto1_c",
        "mto1_c",
        "mto1_c",
        "mto1_h",
        "mto1_h",
        "mto1_h",
        "mto1_h"
    ),
    value = c(
        0.254517588,
        12.35097389,
        1.992312853,
        0.15967005,
        577.7357104,
        3047.866339,
        472.0515931,
        2.814328118,
        1.676886539,
        7.394422182,
        5.861153635,
        61.16714657,
        175.8704458,
        114.4769462,
        163.8271283
    ),
    x = ""
)

pcr_res$x[pcr_res$value == 3047.866339] <- "out"
pcr_res$value <- pcr_res$value / mean(pcr_res[pcr_res$name == "wt_c", 2])

library(ggplot2)
library(ggbreak)
library(multcompView)

# Set factor levels for desired order
pcr_res$name <- factor(pcr_res$name, levels = c("wt_c", "wt_h", "mto1_c", "mto1_h"))
pcr_res$genotype <- factor(ifelse(grepl("wt", pcr_res$name), "WT", "mto1"), levels = c("WT", "mto1"))

# ANOVA + Tukey post-hoc test
aov_res <- aov(value ~ name, data = pcr_res[pcr_res$x != "out", ])
tukey_res <- TukeyHSD(aov_res)
print(summary(aov_res))
print(tukey_res)

# Extract compact letter display (CLD)
tukey_pvals <- tukey_res$name[, "p adj"]
cld <- multcompLetters(tukey_pvals)$Letters

# Create label dataframe with y positions for each group
# (adjust y_pos values to position letters above each box)
label_df <- data.frame(
    name = factor(names(cld), levels = c("wt_c", "wt_h", "mto1_c", "mto1_h")),
    letter = toupper(as.character(cld)),
    y_pos = c(
        "wt_c" = 25, # just above wt_c max
        "wt_h" = 182.5, # above wt_h (in upper panel)
        "mto1_c" = 23, # just above mto1_c max
        "mto1_h" = 75 # just above mto1_h max
    )[names(cld)]
)

svg("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/qpcr_results_100326.svg", width = 2.85, height = 3, family = "serif")
p <- ggplot(pcr_res, aes(x = name, y = value, fill = name)) +
    geom_boxplot(
        data = pcr_res[pcr_res$x != "out", ],
        alpha = 0.7, outlier.shape = NA
    ) +
    geom_jitter(aes(color = genotype), width = 0.15, size = 2, alpha = 0.8) +
    geom_vline(xintercept = 2.5, color = "gray", linetype = "dashed", linewidth = 0.5) +
    geom_text(
        data = label_df, aes(x = name, y = y_pos, label = letter),
        inherit.aes = FALSE, size = 4.5, fontface = "bold", color = "gray35"
    ) +
    # scale_fill_manual(values = c("wt_c" = "gray50", "wt_h" = "gray80",
    #                               "mto1_c" = "#eed3aa", "mto1_h" = "#eed3aa")) +
    scale_fill_manual(values = c(
        "wt_c" = "white", "wt_h" = "white",
        "mto1_c" = "white", "mto1_h" = "white"
    )) +
    scale_color_manual(
        values = c("WT" = "gray20", "mto1" = "#b6801d"),
        labels = c("WT" = "WT", "mto1" = expression(italic("mto1")))
    ) +
    guides(fill = "none", color = guide_legend(title = NULL, override.aes = list(size = 3))) +
    scale_x_discrete(labels = c(
        "wt_c" = "RT", "wt_h" = "37°C",
        "mto1_c" = "RT", "mto1_h" = "37°C"
    )) +
    labs(x = NULL, y = bquote("Relative Expression (" * italic("Copia78") * ")")) +
    scale_y_break(c(182.5, 800), space = 0.5) + # , scales = 150, space = 250
    scale_y_continuous(breaks = c(0, 50, 150, 800), limits = c(0, 840)) + # , labels  = c(0, 50, 150, "", 825, "")
    theme_classic() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(5, 5, 20, 5)
    )
print(p)

# Add legend inside the plot (top-right) using grid
library(grid)
grid.points(
    x = unit(0.78, "npc"), y = unit(0.875, "npc"),
    pch = 19, size = unit(0.55, "char"),
    gp = gpar(col = "gray35")
)
grid.text("WT",
    x = unit(0.82, "npc"), y = unit(0.875, "npc"),
    just = "left", gp = gpar(fontsize = 12, fontfamily = "serif")
)
grid.points(
    x = unit(0.78, "npc"), y = unit(0.825, "npc"),
    pch = 19, size = unit(0.55, "char"),
    gp = gpar(col = "#b6801d")
)
grid.text(expression(italic("mto1")),
    x = unit(0.82, "npc"), y = unit(0.825, "npc"),
    just = "left", gp = gpar(fontsize = 12, fontfamily = "serif")
)
dev.off()


a <- pcr_res[pcr_res$name == "wt_c", 2]
b <- pcr_res[pcr_res$name == "wt_h", 2]
b2 <- pcr_res[pcr_res$x != "out" & pcr_res$name == "wt_h", 2]
c <- pcr_res[pcr_res$name == "mto1_c", 2]
d <- pcr_res[pcr_res$name == "mto1_h", 2]

mean(b) / mean(a)
mean(b2) / mean(a)

mean(d) / mean(c)
