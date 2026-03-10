library(ggplot2)
library(dplyr)

df_x <- as.data.frame(readxl::read_xlsx("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_genes.xlsx",
    sheet = "mto1_vs_wt"
)) %>%
    select(gene_id, Symbol, RNA_log2FC, RNA_pvalue) %>%
    distinct(gene_id, .keep_all = T) %>%
    mutate(Symbol = if_else(is.na(Symbol), gene_id, Symbol)) %>%
    mutate(Symbol = factor(Symbol, levels = c("DML1", "DML2", "DML3", "DME"))) %>%
    mutate(
        color = ifelse(RNA_log2FC > 0, "red", "blue"),
        significance = case_when(
            RNA_pvalue < 0.001 ~ "***",
            RNA_pvalue < 0.01 ~ "**",
            RNA_pvalue < 0.05 ~ "*",
            TRUE ~ ""
        ),
        text_pos = ifelse(RNA_log2FC > 0, RNA_log2FC + 0.05, 0.05)
    ) %>%
    filter(grepl("DML|DME", Symbol)) # keep only thise mentioned in the paper


y_max_plot <- max(df_x$RNA_log2FC)
y_min_plot <- min(df_x$RNA_log2FC)

adding_max <- (y_max_plot + abs(y_min_plot)) * 0.1
adding_min <- (y_max_plot + abs(y_min_plot)) * 0.05

g1 <- ggplot(data = df_x, aes(x = Symbol, y = RNA_log2FC, fill = color)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.75, colour = "black") +
    geom_hline(yintercept = 0, color = "black") +  # Added this line
    geom_text(aes(label = significance, y = text_pos), vjust = 0, color = "black", size = 6) +
    scale_fill_manual(values = c("red" = "#d96c6c", "blue" = "#6c96d9")) + #"gray60")) +
    theme_classic() +
    theme( # plot.margin = unit(c(1, 1, 4, 1), "lines"),
        # title = element_text(size = 9, face="bold"),
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12, face = "bold", angle = 90),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.line.x = element_line(linewidth = 1.1),
        axis.line.y = element_line(linewidth = 1.1),
        axis.ticks = element_line(linewidth = 1.1),
        axis.ticks.length = unit(0.01, "cm"),
        legend.position = "none"
    ) +
    labs(y = "Log2(fold change)") +
    scale_y_continuous(expand = c(0, 0), limits = c(y_min_plot - adding_min, y_max_plot + adding_max))

h <- 2.25
w <- 2
svg("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/deMTs_barPplot_RNAseq.svg", width = w, height = h, family = "serif")
plot(g1)
dev.off()


