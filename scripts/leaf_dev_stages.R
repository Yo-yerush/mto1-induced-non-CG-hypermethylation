library(dplyr)

grow <- read.csv("C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/grow.txt.gz", sep = "\t", header = F)
names(grow) <- c("gene_id", "locus_info", "gene_symbol", "relationship", "developmental_stage", "stage_ontology", "stage_id", "qualifier", "evidence_code", "evidence_description", "with_from", "reference", "assigned_by", "date")

groups_2_keep <- c(
    Seedling = "seedling development stage",
    "LP-2" = "two leaves visible stage",
    "LP-6" = "six leaves visible stage",
    "LP-8" = "eight leaves visible stage",
    "LP-10" = "ten leaves visible stage",
    "LP-12" = "twelve leaves visible stage",
    # "fourteen leaves visible stage",
    # "early rosette growth stage",
    # Rosette = "rosette growth complete stage",
    Senescence = "4 leaf senescence stage"
)

grow <- grow %>%
    arrange(stage_ontology) %>%
    arrange(stage_id) %>%
    # select(gene_id, relationship, developmental_stage, stage_ontology, stage_id, evidence_description) %>%
    select(gene_id, relationship, developmental_stage, evidence_description) %>%
    filter(grepl(paste(groups_2_keep, collapse = "|"), developmental_stage)) %>%
    mutate(tmp = row_number())


rnaseq_res <- read.csv("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all_genes_results_mto1_vs_wt.csv") %>%
    filter(padj < 0.05) %>%
    # filter(pValue < 0.05) %>%
    # filter(gene_model_type == "transposable_element_gene") %>%
    # select(gene_id, log2FoldChange, padj, pValue, Symbol, gene_model_type, short_description)
    select(gene_id, log2FoldChange, pValue, padj)


rnaseq_res %>%
    filter(log2FoldChange > 0) %>%
    nrow()

merged_df <- merge(rnaseq_res, grow, by = "gene_id", all.y = T) %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(tmp) %>%
    select(-tmp)

#########################################

merged_expression_dtage <- merged_df[, grep("gene_id|developmental_stage|log2FoldChange", names(merged_df))]

for (lp in 1:length(groups_2_keep)) {
    merged_expression_dtage$developmental_stage[grep(groups_2_keep[lp], merged_expression_dtage$developmental_stage)] <- names(groups_2_keep)[lp]
}

#########################################

library(ggplot2)

#########################################

# 4. Violin plot with summary statistics
pdf("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/leaf_dev_violinPlot.pdf", width = 4.75, height = 3.5, family = "serif")
ggplot(merged_expression_dtage, aes(x = factor(developmental_stage, levels = names(groups_2_keep)), y = log2FoldChange)) +
    geom_violin(aes(fill = developmental_stage), alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    ggthemes::theme_base() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 10)
    ) +
    labs(
        title = "Violin Plot of log2FoldChange by Developmental Stage",
        x = "Developmental Stage",
        y = "log2FoldChange"
    ) +
    guides(fill = "none")
dev.off()

# 5. Count of genes per stage
# unique genes for each stage
unique_genes_per_stage <- merged_expression_dtage %>%
    group_by(gene_id) %>%
    mutate(n_stages = n()) %>%
    filter(n_stages == 1) %>%
    ungroup() %>%
    count(developmental_stage, name = "unique_gene_count") %>%
    arrange(developmental_stage)

# Genes shared in stages
genes_shared_in_all <- merged_expression_dtage %>%
    group_by(gene_id) %>%
    summarise(n_stages = n_distinct(developmental_stage)) %>%
    filter(n_stages == length(groups_2_keep))

genes_shared_in_LP <- merged_expression_dtage %>%
    filter(developmental_stage != "Seedling") %>%
    filter(developmental_stage != "Senescence") %>%
    group_by(gene_id) %>%
    summarise(n_stages = n_distinct(developmental_stage)) %>%
    filter(n_stages == length(groups_2_keep) - 2)

stage_counts <- merged_expression_dtage %>%
    count(developmental_stage, name = "gene_count") %>%
    arrange(developmental_stage) %>%
    mutate(unique_gene_count = unique_genes_per_stage$unique_gene_count) %>%
    mutate(shared_gene_count = gene_count - unique_gene_count) %>%
    mutate(shared_in_all_count = nrow(genes_shared_in_all)) %>%
    mutate(shared_in_LP_count = nrow(genes_shared_in_LP))
stage_counts$shared_in_LP_count[grep("Seedling|Senescence", stage_counts$developmental_stage)] <- 0

print("Number of unique genes per developmental stage:")
print(unique_genes_per_stage)

pdf("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/leaf_dev_count_barPlot.pdf", width = 5, height = 3, family = "serif")
library(ggbreak)

ggplot(stage_counts, aes(x = factor(developmental_stage, levels = names(groups_2_keep)[length(groups_2_keep):1]))) +
    geom_col(aes(y = shared_gene_count, fill = "Shared"), alpha = 0.8) +
    # geom_col(aes(y = shared_in_LP_count, fill = "shared - LP"), alpha = 0.8) +
    geom_col(aes(y = unique_gene_count, fill = "Unique"), alpha = 0.8) +
    geom_text(aes(y = gene_count, label = gene_count), hjust = -0.1, size = 3) +
    coord_flip() +
    # scale_y_cut(breaks = 60,   scales = "fixed", space = 0.5) +
    ggthemes::theme_base() +
    labs(
        title = "Number of Genes per Developmental Stage",
        x = "Developmental Stage",
        y = "Number of Genes"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(stage_counts$gene_count) * 1.15)) +
    scale_fill_manual(values = c(
        "Shared" = "steelblue",
        # "shared - LP" = "#65aa65",
        "Unique" = "#be0000"
    )) +
    theme(
        plot.title = element_text(size = 11),
        axis.text.y = element_text(size = 12), # , angle = 45),
        axis.text.x = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)
    ) #+ guides(fill = guide_legend(breaks = c(shared_leg, "unique\ngenes")))
dev.off()
print("Gene counts by stage:")
print(stage_counts)

#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################

# 1. Distribution of log2FoldChange by developmental stage
ggplot(merged_expression_dtage, aes(x = developmental_stage, y = log2FoldChange)) +
    geom_boxplot(aes(fill = developmental_stage)) +
    geom_jitter(alpha = 0.3, width = 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
        title = "Distribution of log2FoldChange by Developmental Stage",
        x = "Developmental Stage",
        y = "log2FoldChange"
    ) +
    guides(fill = "none")

# 2. Density plot of log2FoldChange for each stage
ggplot(merged_expression_dtage, aes(x = log2FoldChange, fill = developmental_stage)) +
    geom_density(alpha = 0.7) +
    facet_wrap(~developmental_stage, scales = "free_y") +
    theme_minimal() +
    labs(
        title = "Density Distribution of log2FoldChange by Stage",
        x = "log2FoldChange",
        y = "Density"
    ) +
    guides(fill = "none")

# 3. Summary statistics table
summary_stats <- merged_expression_dtage %>%
    group_by(developmental_stage) %>%
    summarise(
        n_genes = n(),
        mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
        median_log2FC = median(log2FoldChange, na.rm = TRUE),
        sd_log2FC = sd(log2FoldChange, na.rm = TRUE),
        min_log2FC = min(log2FoldChange, na.rm = TRUE),
        max_log2FC = max(log2FoldChange, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(developmental_stage)

print("Summary statistics by developmental stage:")
print(summary_stats)

#########################################
# Additional plots and gene overlap analysis

library(VennDiagram)
library(UpSetR)
library(ComplexHeatmap)
library(tidyr)
library(tibble)

gene_stage_matrix <- merged_df %>%
    select(gene_id, log2FoldChange, developmental_stage) %>%
    # distinct() %>%
    # Create a wide format: genes as rows, stages as columns
    pivot_wider(
        names_from = developmental_stage,
        values_from = log2FoldChange,
        values_fill = NA
    ) %>%
    column_to_rownames("gene_id")

# Check the matrix
print("Gene-stage matrix dimensions:")
print(dim(gene_stage_matrix))
print("First few rows and columns:")
print(gene_stage_matrix[1:5, 1:min(3, ncol(gene_stage_matrix))])


# 6. Scatter plot matrix of log2FoldChange between stages
if (ncol(gene_stage_matrix) >= 2) {
    # Create pairwise scatter plots for first few stages
    stages_subset <- colnames(gene_stage_matrix)[1:min(4, ncol(gene_stage_matrix))]
    pairs(gene_stage_matrix[, stages_subset],
        main = "Pairwise log2FoldChange correlations",
        pch = 16, alpha = 0.6
    )
}

# 7. Heatmap of actual expression values (top varying genes)
# Select top 50 most variable genes across stages
gene_vars <- apply(gene_stage_matrix, 1, function(x) var(x, na.rm = TRUE))
top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(50, length(gene_vars))]

gene_stage_matrix_hm <- gene_stage_matrix[top_var_genes, ] %>%
    filter(if_all(-1, is.finite))

pheatmap(
    as.matrix(gene_stage_matrix_hm),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    # scale = "row",
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = "Top 50 most variable genes across stages",
    show_rownames = FALSE
)

#########################################
# Gene overlap analysis

# Create gene sets for each developmental stage
gene_sets <- merged_df %>%
    group_by(developmental_stage) %>%
    summarise(genes = list(unique(gene_id)), .groups = "drop") %>%
    deframe()

# 8. Gene overlap statistics
overlap_stats <- data.frame(
    stage = names(gene_sets),
    total_genes = sapply(gene_sets, length),
    stringsAsFactors = FALSE
)

print("Gene counts per stage:")
print(overlap_stats)

# 9. Identify unique vs shared genes
all_genes <- unique(merged_df$gene_id)
gene_stage_membership <- merged_df %>%
    select(gene_id, developmental_stage) %>%
    distinct() %>%
    group_by(gene_id) %>%
    summarise(
        n_stages = n(),
        stages = paste(developmental_stage, collapse = ", "),
        .groups = "drop"
    )

# Summary of gene sharing
gene_sharing_summary <- gene_stage_membership %>%
    count(n_stages, name = "n_genes") %>%
    mutate(percentage = round(n_genes / sum(n_genes) * 100, 1))

print("Gene sharing across stages:")
print(gene_sharing_summary)

# 10. Upset plot for gene overlaps (if you have UpSetR installed)
if (requireNamespace("UpSetR", quietly = TRUE)) {
    # Create binary matrix for UpSetR
    binary_matrix <- merged_df %>%
        select(gene_id, developmental_stage) %>%
        distinct() %>%
        mutate(present = 1) %>%
        tidyr::pivot_wider(
            names_from = developmental_stage,
            values_from = present,
            values_fill = 0
        ) %>%
        column_to_rownames("gene_id")

    UpSetR::upset(
        as.data.frame(binary_matrix),
        nsets = ncol(binary_matrix),
        order.by = "freq",
        main.bar.color = "steelblue",
        sets.bar.color = "darkgreen"
    )
}

# 11. Bar plot of gene sharing patterns
ggplot(gene_sharing_summary, aes(x = as.factor(n_stages), y = n_genes)) +
    geom_col(fill = "coral", alpha = 0.7) +
    geom_text(aes(label = paste0(n_genes, "\n(", percentage, "%)")),
        vjust = -0.3, size = 3
    ) +
    theme_minimal() +
    labs(
        title = "Distribution of Genes by Number of Stages",
        subtitle = "How many genes appear in 1, 2, 3, ... stages",
        x = "Number of Stages",
        y = "Number of Genes"
    )

# 12. Expression patterns of genes in multiple stages
multi_stage_genes <- gene_stage_membership %>%
    filter(n_stages > 1) %>%
    pull(gene_id)

if (length(multi_stage_genes) > 0) {
    multi_stage_data <- merged_df %>%
        filter(gene_id %in% multi_stage_genes[1:min(20, length(multi_stage_genes))]) %>%
        mutate(gene_id = factor(gene_id))

    ggplot(multi_stage_data, aes(x = developmental_stage, y = log2FoldChange, group = gene_id)) +
        geom_line(aes(color = gene_id), alpha = 0.7) +
        geom_point(aes(color = gene_id), size = 2) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
        ) +
        labs(
            title = "Expression Patterns of Multi-Stage Genes (top 20)",
            x = "Developmental Stage",
            y = "log2FoldChange"
        )
}

# 13. Stage-specific vs shared genes analysis
stage_specific_genes <- gene_stage_membership %>%
    filter(n_stages == 1) %>%
    separate_rows(stages, sep = ", ") %>%
    count(stages, name = "unique_genes") %>%
    rename(developmental_stage = stages)

shared_genes_per_stage <- merged_df %>%
    select(gene_id, developmental_stage) %>%
    distinct() %>%
    left_join(gene_stage_membership %>% select(gene_id, n_stages), by = "gene_id") %>%
    filter(n_stages > 1) %>%
    count(developmental_stage, name = "shared_genes")

stage_gene_summary <- merge(
    overlap_stats %>% select(stage, total_genes) %>% rename(developmental_stage = stage),
    stage_specific_genes,
    by = "developmental_stage", all.x = TRUE
) %>%
    merge(shared_genes_per_stage, by = "developmental_stage", all.x = TRUE) %>%
    mutate(
        unique_genes = ifelse(is.na(unique_genes), 0, unique_genes),
        shared_genes = ifelse(is.na(shared_genes), 0, shared_genes)
    ) %>%
    tidyr::pivot_longer(
        cols = c(unique_genes, shared_genes),
        names_to = "gene_type",
        values_to = "count"
    )

ggplot(stage_gene_summary, aes(x = developmental_stage, y = count, fill = gene_type)) +
    geom_col(position = "stack") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
        title = "Stage-specific vs Shared Genes",
        x = "Developmental Stage",
        y = "Number of Genes",
        fill = "Gene Type"
    ) +
    scale_fill_manual(values = c("unique_genes" = "lightcoral", "shared_genes" = "lightblue"))

print("Top 10 genes appearing in most stages:")
print(gene_stage_membership %>% arrange(desc(n_stages)) %>% head(10))

#########################################


unique(grow$developmental_stage)
unique(c(
    "flower development stage",
    "flower meristem visible stage",
    "flower meristem notched stage",
    "flower organ development stage",
    "anther development stage",
    "formation of primary parietal and sporogenous cells stage",
    "locules established stage",
    "pollen development stage",
    "pollen mother cell meiosis stage",
    "tetrad stage",
    "dry seed stage",
    "plant embryo stage",
    "proembryo stage",
    "globular stage",
    "bilateral stage",
    "expanded cotyledon stage",
    "mature embryo stage",
    "endosperm development stage",
    "chalazal and micropylar domain establishment stage",
    "primary endosperm cell stage",
    "fertilized ovule stage",
    "seed maturation stage",
    "sporophyte vegetative stage",
    "seed germination stage",
    "seed imbibition stage",
    "seedling development stage",
    "two leaves visible stage",
    "four leaves visible stage",
    "six leaves visible stage",
    "eight leaves visible stage",
    "ten leaves visible stage",
    "twelve leaves visible stage",
    "fourteen leaves visible stage",
    "early rosette growth stage",
    "rosette growth complete stage",
    "sporophyte reproductive stage",
    "inflorescence just visible stage",
    "ripening stage",
    "sporophyte senescent stage",
    "leaf trichome development stage",
    "sporophyte development stage"
)) %>% length()
