source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/multiplot_ggplot2.R")
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/yo_theme_base_ggplot2.r")

average_meth_results_directory <- paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/average.meth.genes.levels/mto1_vs_wt/")
DMRs_results_directory <- paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/")
RNAseq_results_directory <- paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/")
main_output_directory <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/linear_plot_TEGs_width_one_plot/"

var1_name <- "wt"
var2_name <- "mto1"
var1_rnaseq_names <- c("met20", "met22")
var2_rnaseq_names <- c("met14", "met15", "met16")
genes_2_keep <- "filtered_by_DEGs"
additional_plots <- TRUE
pValues_table <- TRUE
var1_col <- "#7F7F7F"
var2_col <- "#bf6828"

library(plyr)
library(dplyr)
library(corrplot)
library(ggpubr)
library(ggplot2)
library(lmerTest)
library(tidyr)
library(MASS)
library(scales)

# pValues function
pval_fun <- function(x, is.inte=T) {
    if (!is.inte) {
        "n.e."
    } else {
        ifelse(x <= 0.001, "***",
            ifelse(x <= 0.01, "**  ",
                ifelse(x <= 0.05, "*    ",
                    "ns  "
                )
            )
        )
    }
}

### filter by 'genes_2_keep'
rna_2_filter <- read.csv(paste0(RNAseq_results_directory, "all_genes_results_", var2_name, "_vs_", var1_name, ".csv"))
gene_list_name <- genes_2_keep

if (genes_2_keep == "filtered_by_DEGs") {
    gene_list <- filter(rna_2_filter, padj < 0.05) %>% distinct(gene_id)
} else if (genes_2_keep == "filtered_by_DMRs" | genes_2_keep == "all_genes") {
    gene_list <- distinct(rna_2_filter, gene_id)
} else {
    gene_list <- filter(rna_2_filter, gene_id %in% genes_2_keep) %>% distinct(gene_id)
    gene_list_name <- "selected_genes"
}

# open new plot for some reason..........
ggplot() +
    theme_void()

annotation_type_names <- c("transposable_element_gene") # "gene",

annotation_type_length <- length(annotation_type_names)

#### addit annotation names ####
edit_ann_names <- function(x) {
    x <- gsub("^promoter$", "Promoters", x)
    # x <- gsub("^gene$", "Genes", x)
    x <- gsub("^intron$", "Introns", x)
    x <- gsub("^five_prime_UTR$", "5'UTRs", x)
    x <- gsub("^three_prime_UTR$", "3'UTRs", x)
    x <- gsub("^transposable_element_gene$", "TEGs", x)
    return(x)
}

# set list for all context plots
all_cntx_points <- setNames(vector("list", 3), c("CG", "CHG", "CHH"))
all_cntx_regression <- setNames(vector("list", 3), c("CG", "CHG", "CHH"))
all_cntx_residuals <- setNames(vector("list", 3), c("CG", "CHG", "CHH"))


#########################
#### theta, R^2 and p values (to plot after all)
theta_df <- data.frame(annotation = rep(edit_ann_names(annotation_type_names), 3), theta_value = NA, theta_se = NA, context = c(rep("CG", annotation_type_length), rep("CHG", annotation_type_length), rep("CHH", annotation_type_length)))
Rsqr_df <- data.frame(annotation = rep(edit_ann_names(annotation_type_names), 3), Rsqr_value = NA, context = c(rep("CG", annotation_type_length), rep("CHG", annotation_type_length), rep("CHH", annotation_type_length)))
ps_Rsqr_df <- data.frame(annotation = rep(edit_ann_names(annotation_type_names), 3), Rsqr_value = NA, context = c(rep("CG", annotation_type_length), rep("CHG", annotation_type_length), rep("CHH", annotation_type_length)))
p_df <- data.frame(annotation = rep(edit_ann_names(annotation_type_names), 3), Meth = NA, genotype = NA, interaction = NA, context = c(rep("CG", annotation_type_length), rep("CHG", annotation_type_length), rep("CHH", annotation_type_length)))

combined_plot_df <- data.frame()

i_perc <- 1
total_perc <- 18
cat(paste0(var2_name, "_vs_", var1_name, ":\n"))

for (context in c("CG", "CHG", "CHH")) {
    # set list for all annotation plots
    # reset it for each context
    annotation_point_plots_list <- setNames(vector("list", length(annotation_type_names)), annotation_type_names)
    annotation_regression_plot_list <- setNames(vector("list", length(annotation_type_names)), annotation_type_names)
    residuals_plots_list <- setNames(vector("list", length(annotation_type_names)), annotation_type_names)

    for (annotation_type in annotation_type_names) {
        path_2_save.0 <- main_output_directory
        path_2_save.1 <- paste0(path_2_save.0, gene_list_name, "/")
        path_2_save.2 <- paste0(path_2_save.1, var2_name, "/")
        dir.create(path_2_save.0, showWarnings = F)
        dir.create(path_2_save.1, showWarnings = F)
        dir.create(path_2_save.2, showWarnings = F)

        # Load the methylation files of the promoters of interestL
        meth_matrix <- read.csv(paste0(average_meth_results_directory, context, "/meth.", annotation_type, ".", context, ".", var2_name, "_vs_", var1_name, ".csv"))

        # filter by DMRs-gene_list (if exist)
        if (gene_list_name == "filtered_by_DMRs") {
            meth_matrix <- filter(meth_matrix, gene_id %in% gene_list$gene_id)
        }

        # Load the RNAseq '*.gene.results' files (output from RSEM pipeline)
        RNA <- read.csv(paste0(RNAseq_results_directory, "norm_counts_", var2_name, "_vs_", var1_name, ".csv"))

        # filter by DMRs-gene_list (if exist)
        if (gene_list_name == "filtered_by_DEGs") {
            RNA <- filter(RNA, gene_id %in% gene_list$gene_id)
        }

        #### Gene expression ####
        # edit table and names

        var1_names_pos <- grep(paste(var1_rnaseq_names, collapse = "|"), names(RNA))
        var2_names_pos <- grep(paste(var2_rnaseq_names, collapse = "|"), names(RNA))

        if (length(var1_names_pos) == 0 | length(var2_names_pos) == 0) {
            stop(paste0("<<", paste(var1_rnaseq_names, collapse = ", "), ">> or <<", paste(var2_rnaseq_names, collapse = ", "), ">> not found in RNAseq results. Please check the names in the RNAseq results file."))
        }

        names(RNA)[var1_names_pos] <- paste0(var1_name, ".", seq(var1_names_pos), "_RNA")
        names(RNA)[var2_names_pos] <- paste0(var2_name, ".", seq(var2_names_pos), "_RNA")

        RNA <- RNA[, c(1, var1_names_pos, var2_names_pos)]


        #### Methylation ####
        meth_matrix <- meth_matrix[, grep(paste("gene_id", var2_name, var1_name, sep = "|"), names(meth_matrix))]
        names(meth_matrix)[-1] <- paste0(names(meth_matrix)[-1], "_meth")
        meth_matrix <- na.omit(meth_matrix) # remove rows contain 'NA'

        #### Merge methylation and RNA dataframes ####
        meth_matrix_RNA <- merge(meth_matrix, RNA, by = "gene_id", all = FALSE)

        # filter by other gene_list (if exist)
        if (gene_list_name != "all_genes" & gene_list_name != "filtered_by_DMRs" & gene_list_name != "filtered_by_DEGs") {
            meth_matrix_RNA <- filter(meth_matrix_RNA, gene_id %in% gene_list$gene_id)
        }

        meth_new <- meth_matrix_RNA[, c(grep("gene_id", names(meth_matrix_RNA)), grep("meth", names(meth_matrix_RNA)))]
        RNA_new <- meth_matrix_RNA[, c(grep("gene_id", names(meth_matrix_RNA)), grep("RNA", names(meth_matrix_RNA)))]

        names(meth_new) <- gsub("_meth", "", names(meth_new))
        names(RNA_new) <- gsub("_RNA", "", names(RNA_new))

        # move the columns/rows:
        fin_meth <- meth_new %>%
            pivot_longer(!c(gene_id), names_to = "Sample", values_to = "Meth")

        fin_RNA <- RNA_new %>%
            pivot_longer(!c(gene_id), names_to = "Sample", values_to = "normCounts")

        fin_meth_RNA <- cbind(fin_meth, normCounts = fin_RNA$normCounts)
        # fin_meth_RNA <- fin_meth_RNA %>%
        #  dplyr::select(gene_id, everything())
        fin_meth_RNA$Meth[is.na(fin_meth_RNA$Meth)] <- 0 # i dont know if its the right way

        #### final DF for plotting ####
        plot_df <- fin_meth_RNA %>%
            mutate(genotype = ifelse(grepl(var1_name, Sample), var1_name, var2_name))

        #####################################################

        tryCatch(
            {
                #### linear model ####
                plot_df$normCounts <- as.numeric(plot_df$normCounts)
                plot_df$Meth <- as.numeric(plot_df$Meth)
                plot_df$genotype <- as.factor(plot_df$genotype)

                # do it with/without interaction, and initialize lm_model results
                lm_model_0 <- NULL
                is.interaction <- FALSE
                tryCatch(
                    {
                        lm_model_0 <<- glm.nb(normCounts ~ Meth * genotype, data = plot_df)
                        is.interaction <<- TRUE
                    },
                    error = function(e) {
                        lm_model_0 <<- glm.nb(normCounts ~ Meth + genotype, data = plot_df)
                        is.interaction <<- FALSE
                    }
                )
                lm_model <- summary(lm_model_0)

                # fitted vs. residuals
                plot_df$residuals <- residuals(lm_model_0)
                plot_df$fitted <- fitted(lm_model_0)

                # Meth and genotype rows
                meth_row <- grep("^Meth$", rownames(lm_model$coefficients))
                genotype_row <- grep("^genotypewt$", rownames(lm_model$coefficients))
                if (is.interaction) {interaction_row <- grep("Meth:genotypewt", rownames(lm_model$coefficients))}

                #### linear model p-Value ####
                M_pval <- pval_fun(lm_model$coefficients[meth_row, 4])
                G_pval <- pval_fun(lm_model$coefficients[genotype_row, 4])
                MnG_pval <- pval_fun(lm_model$coefficients[interaction_row, 4], is.interaction)

                #####################################################

                ann_cntx_row_number <- which(Rsqr_df$annotation == edit_ann_names(annotation_type) & Rsqr_df$context == context)

                #### p value ####
                p_df[ann_cntx_row_number, ]$Meth <- lm_model$coefficients[meth_row, 4]
                p_df[ann_cntx_row_number, ]$genotype <- lm_model$coefficients[genotype_row, 4]
                p_df[ann_cntx_row_number, ]$interaction <- ifelse(!is.interaction, "n.e.", lm_model$coefficients[interaction_row, 4])

                #### R^2 value ####
                res_dev <- lm_model_0$deviance
                null_dev <- lm_model_0$null.deviance
                calc_R2 <- 1 - (res_dev / null_dev)
                Rsqr_df[ann_cntx_row_number, ]$Rsqr_value <- as.numeric(calc_R2)
                Rsqr_df[ann_cntx_row_number, ]$context <- context

                #### Pseudo-R^2 value ####
                invisible(capture.output({ # run it quietly
                    ps_Rsqr <- pscl::pR2(lm_model_0)[4]
                }))
                ps_Rsqr_df[ann_cntx_row_number, ]$Rsqr_value <- as.numeric(ps_Rsqr)
                ps_Rsqr_df[ann_cntx_row_number, ]$context <- context

                #### theta value ####
                theta_value <- lm_model_0$theta
                theta_se <- lm_model_0$SE.theta
                theta_df[ann_cntx_row_number, ]$theta_value <- theta_value
                theta_df[ann_cntx_row_number, ]$theta_se <- theta_se
                theta_df[ann_cntx_row_number, ]$context <- context

                #####################################################

                # accumulate data for combined plot (all contexts in one)
                plot_df$context <- context
                combined_plot_df <- rbind(combined_plot_df, plot_df)


                #### residuals plot ####
                residuals_plots_list[[annotation_type]] <- ggplot(plot_df, aes(x = fitted, y = residuals, color = genotype)) +
                    # ggrastr::geom_point_rast(alpha = 0.075, size = 0.25) +
                    geom_point(alpha = 0.1, size = 0.25) +
                    # geom_smooth() +
                    scale_color_manual(values = c(var2_col, var1_col)) +
                    yo_theme_base(base_rect_size = 1) +
                    theme(legend.position = "none") +
                    labs(
                        title = edit_ann_names(annotation_type),
                        x = "Fitted",
                        y = "Residuals"
                    )

                # message(paste0("created plot to - '", annotation_type, "' annotations in '", context, "' context successfully"))
            },
            error = function(e) {
                message(paste0("\ncan't use 'glm.nb' and create plot to - ", annotation_type, " annotations in ", context, " context\n\n"))
                annotation_point_plots_list[[annotation_type]] <- ggplot() +
                    theme_void()
                annotation_regression_plot_list[[annotation_type]] <- ggplot() +
                    theme_void()
                residuals_plots_list[[annotation_type]] <- ggplot() +
                    theme_void()
            }
        )
        cat(paste0("\rcreate main plot: ", round((i_perc / total_perc) * 100, 0), "%  "))
        i_perc <- i_perc + 1
    }

    ### list of all plots from contexts ###
    all_cntx_points[[context]] <- annotation_point_plots_list
    # all_cntx_regression[[context]] <- annotation_regression_plot_list
    all_cntx_residuals[[context]] <- residuals_plots_list
}
cat("\n")

total_plots <- ifelse(additional_plots, 5, 2)
cat(paste0("\rsave plots: [", 0, "/", total_plots, "]"))






#### Combined plot - all contexts, shapes for genotype, colors for context ####
context_colors <- c("CG" = "#8C8C8C", "CHG" = "#5A7D9A", "CHH" = "#9C7A8E")
genotype_shapes <- setNames(c(16, 17), c(var1_name, var2_name))
genotype_linetypes <- setNames(c("solid", "dashed"), c(var1_name, var2_name))

combined_plot_df$context <- factor(combined_plot_df$context, levels = c("CG", "CHG", "CHH"))
combined_plot_df$genotype <- factor(combined_plot_df$genotype, levels = c(var1_name, var2_name))

combined_plot <- ggplot(combined_plot_df, aes(x = Meth, y = normCounts, color = context, shape = genotype)) +
    geom_point(alpha = 0.25, size = 1) +
    stat_smooth(aes(linetype = genotype), method = "glm.nb", formula = y ~ x, linewidth = 0.95, se = F) +
    scale_y_continuous(
        trans = scales::log1p_trans(),
        breaks = scales::rescale(c(0, 0.01, 0.25, 1), to = range(plot_df$normCounts)),
        labels = function(x) {
            ifelse(x < 10000, label_number()(x), label_scientific(digits = 0)(x))
        }
    ) +
    scale_color_manual(values = context_colors) +
    scale_shape_manual(values = genotype_shapes) +
    scale_linetype_manual(values = genotype_linetypes) +
    yo_theme_base(base_rect_size = 1) +
    theme(
        legend.position = "right",
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12, color = "#4b4b4b"),
        axis.title = element_text(size = 16)
    ) +
    labs(
        x = "Methylation level",
        y = "Norm. Counts",
        color = "Context",
        shape = "Genotype",
        linetype = "Genotype"
    )

### save combined plot
tiff(paste0(path_2_save.2, "lm.stats.plot.", var2_name, ".tif"), width = 6, height = 3, units = "in", res = 600, family = "serif")
print(combined_plot)
dev.off()
cat(paste0("\rsave plots: [", 1, "/", total_plots, "]"))























### ### residuals plot
### tiff(paste0(path_2_save.2, "lm.residuals.plot.", var2_name, ".tif"), width = 6, height = 2, units = "in", res = 600, family = "serif")
### do.call(multiplot, c(
###     c(all_cntx_residuals$CG, all_cntx_residuals$CHG, all_cntx_residuals$CHH),
###     list(
###         layout = matrix(1:6, nrow = 1, ncol = 3, byrow = TRUE),
###         heights = c(1, 1, 1),
###         widths = c(1, 0.8, 0.8)
###     )
### ))
### dev.off()
### cat(paste0("\rsave plots: [", 2, "/", total_plots, "]"))
### 
### ################################################################
### 
### #### legend ####
### legend_plot_df <- data.frame(
###     Meth = runif(100, min = 0, max = 1),
###     normCounts = rnorm(100, mean = 10, sd = 2),
###     genotype = sample(c(var1_name, var2_name), 100, replace = TRUE)
### )
### p_lg <- ggplot(legend_plot_df, aes(x = Meth, y = normCounts, fill = genotype)) +
###     geom_point(alpha = 1, size = 6.5, shape = 21, stroke = 0.5, color = "black") +
###     scale_fill_manual(
###         values = c(var2_col, var1_col),
###         breaks = c(var2_name, var1_name) # , labels = c(expression(italic(var1_name)), var1_name)
###     ) +
###     ggthemes::theme_base() +
###     theme(
###         legend.text = element_text(size = 16),
###         legend.position = "left"
###     ) + # Increase legend text size and position it at the top left
###     labs(fill = "")
### 
### legend_p_0 <- ggplotGrob(p_lg)$grobs[[which(sapply(ggplotGrob(p_lg)$grobs, function(x) x$name) == "guide-box")]]
### legend_p <- ggplot() +
###     theme_void() +
###     annotation_custom(legend_p_0, xmin = 0.1, xmax = 0.65, ymin = 0.5, ymax = 1)
### # legend_p = grid::grid.draw(legend_p)
### 
### svg(paste0(path_2_save.2, "legend_lm_plots_", var2_name, ".svg"), width = 2, height = 1, family = "serif")
### print(legend_p)
### dev.off()
### 
### ################################################################
### if (additional_plots) {
###     ### theta bar plot
###     theta_df$annotation <- factor(theta_df$annotation, levels = unique(theta_df$annotation))
###     theta_bar_plot <- ggplot(theta_df, aes(x = annotation, y = theta_value, fill = context)) +
###         geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "black") +
###         geom_errorbar(aes(ymin = theta_value - theta_se, ymax = theta_value + theta_se), width = 0.2, position = position_dodge(0.8)) +
###         scale_fill_brewer(palette = "Set3") +
###         theme_bw() +
###         labs(
###             x = "",
###             y = "ϴ",
###             fill = "Context"
###         ) +
###         theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
### 
###     svg(paste0(path_2_save.2, "theta_values_", var2_name, ".svg"), width = 3, height = 3, family = "serif")
###     print(theta_bar_plot)
###     dev.off()
###     cat(paste0("\rsave plots: [", 3, "/", total_plots, "]"))
### 
###     ### R^2 bar plot
###     Rsqr_df$annotation <- factor(Rsqr_df$annotation, levels = unique(Rsqr_df$annotation))
###     end_brk_1 <- max(Rsqr_df[Rsqr_df$annotation != "TEGs", ]$Rsqr_value)
###     start_brk_2 <- min(c(Rsqr_df[Rsqr_df$annotation == "TEGs", ]$Rsqr_value))
### 
###     Rsqr_bar_plot <- ggplot(Rsqr_df, aes(x = annotation, y = Rsqr_value, fill = context)) +
###         geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "black") +
###         scale_fill_brewer(palette = "Set3") +
###         theme_bw() +
###         labs(
###             x = "",
###             y = "Deviance R²",
###             fill = "Context"
###         ) # +
###         # theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
###         # ggbreak::scale_y_break(
###         #     c(end_brk_1 * 1.075, start_brk_2 * 0.925),
###         #     space = 0.055,
###         #     scales = 0.5
###         # )
### 
###     svg(paste0(path_2_save.2, "Rsqr_values_", var2_name, ".svg"), width = 3, height = 3, family = "serif")
###     print(Rsqr_bar_plot)
###     dev.off()
###     cat(paste0("\rsave plots: [", 4, "/", total_plots, "]"))
### 
###     ### Pseudo-R^2 bar plot
###     ps_Rsqr_df$annotation <- factor(ps_Rsqr_df$annotation, levels = unique(ps_Rsqr_df$annotation))
###     end_brk_1 <- max(ps_Rsqr_df[ps_Rsqr_df$annotation != "TEGs", ]$Rsqr_value)
###     start_brk_2 <- min(c(ps_Rsqr_df[ps_Rsqr_df$annotation == "TEGs", ]$Rsqr_value))
### 
###     ps_Rsqr_bar_plot <- ggplot(ps_Rsqr_df, aes(x = annotation, y = Rsqr_value, fill = context)) +
###         geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "black") +
###         scale_fill_brewer(palette = "Set3") +
###         theme_bw() +
###         labs(
###             x = "",
###             y = "Pseudo-R²",
###             fill = "Context"
###         ) # +
###         # theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
###         # ggbreak::scale_y_break(
###         #     c(end_brk_1 * 1.075, start_brk_2 * 0.925),
###         #     space = 0.055,
###         #     scales = 0.5
###         # )
### 
###     svg(paste0(path_2_save.2, "Pseudo_Rsqr_values_", var2_name, ".svg"), width = 3, height = 3, family = "serif")
###     print(ps_Rsqr_bar_plot)
###     dev.off()
###     cat(paste0("\rsave plots: [", 5, "/", total_plots, "]"))
###     cat("\n")
### }
### 
### if (pValues_table) {
###     ### p values data frame
###     write.csv(p_df, paste0(path_2_save.2, "p_values_", var2_name, ".csv"), row.names = FALSE)
### }
### cat("\n\n")
