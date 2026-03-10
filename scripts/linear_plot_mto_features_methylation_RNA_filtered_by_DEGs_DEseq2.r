library(plyr)
library(dplyr)
library(corrplot)
library(ggpubr)
library(ggplot2)
library(lmerTest)
library(tidyr)
library(MASS)
library(scales)
library(hexbin)
source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/multiplot_ggplot2.R")
source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/yo_theme_base_ggplot2.R")


# # # # # # # # # # # # # # # #

make_it_light <- FALSE # if TRUE, will reduce the number of points in the plot

# # # # # # # # # # # # # # # #

treatment <- "mto1"
mto_rnaseq_names <- c("met14", "met15", "met16")
wt_rnaseq_names <- c("met20", "met22")

average_meth_results_directory <- paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/average.meth.genes.levels/", treatment, "_vs_wt/")

DMRs_results_directory <- paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/", treatment, "_vs_wt/genome_annotation/")

RNAseq_results_directory <- paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/", treatment, "_vs_wt/")

main_output_directory <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/Linear_correlation/"

### filtered by differentially expressed genes
DEGs_genes_list <- list(
  gene_list_name = "filtered_by_DEGs",
  gene_list = read.csv(paste0(RNAseq_results_directory, "all_genes_results_", treatment, "_vs_wt.csv")) %>%
    filter(padj < 0.05) %>%
    distinct(gene_id)
)

# open new plot for some reason..........
ggplot() +
  theme_void()

is.transcript <- F

annotation_type_names <- c("promoter", "CDS", "intron", "five_prime_UTR", "three_prime_UTR", "transposable_element_gene") # "gene",

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

# set list for all annotation plots
annotation_point_plots_list <- setNames(vector("list", length(annotation_type_names)), annotation_type_names)
annotation_regression_plot_list <- setNames(vector("list", length(annotation_type_names)), annotation_type_names)
residuals_plots_list <- setNames(vector("list", length(annotation_type_names)), annotation_type_names)

# set list for all context plots
all_cntx_points <- setNames(vector("list", 3), c("CG", "CHG", "CHH"))
all_cntx_regression <- setNames(vector("list", 3), c("CG", "CHG", "CHH"))
all_cntx_residuals <- setNames(vector("list", 3), c("CG", "CHG", "CHH"))

###############################

gene_list <- DEGs_genes_list$gene_list
gene_list_name <- DEGs_genes_list$gene_list_name


#########################
#### theta, R^2 and p values (to plot after all)
theta_df <- data.frame(annotation = rep(edit_ann_names(annotation_type_names), 3), theta_value = NA, theta_se = NA, context = c(rep("CG", annotation_type_length), rep("CHG", annotation_type_length), rep("CHH", annotation_type_length)))
Rsqr_df <- data.frame(annotation = rep(edit_ann_names(annotation_type_names), 3), Rsqr_value = NA, context = c(rep("CG", annotation_type_length), rep("CHG", annotation_type_length), rep("CHH", annotation_type_length)))
ps_Rsqr_df <- data.frame(annotation = rep(edit_ann_names(annotation_type_names), 3), Rsqr_value = NA, context = c(rep("CG", annotation_type_length), rep("CHG", annotation_type_length), rep("CHH", annotation_type_length)))
p_df <- data.frame(annotation = rep(edit_ann_names(annotation_type_names), 3), Meth = NA, genotype = NA, context = c(rep("CG", annotation_type_length), rep("CHG", annotation_type_length), rep("CHH", annotation_type_length)))

i_perc <- 1
total_perc <- 18

for (context in c("CG", "CHG", "CHH")) {
  for (annotation_type in annotation_type_names) {
    path_2_save.0 <- main_output_directory
    path_2_save.1 <- paste0(path_2_save.0, gene_list_name, "/")
    path_2_save.2 <- paste0(path_2_save.1, treatment, "/")
    dir.create(path_2_save.0, showWarnings = F)
    dir.create(path_2_save.1, showWarnings = F)
    dir.create(path_2_save.2, showWarnings = F)

    # Load the methylation files of the promoters of interestL
    meth_matrix <- read.csv(paste0(average_meth_results_directory, context, "/meth.", annotation_type, ".", context, ".", treatment, "_vs_wt.csv"))

    # filter by DMRs-gene_list (if exist)
    if (gene_list_name == "filtered_by_DMRs") {
      meth_matrix <- filter(meth_matrix, gene_id %in% gene_list$gene_id)
    }

    # Load the RNAseq '*.gene.results' files (output from RSEM pipeline)
    RNA <- read.csv(paste0(RNAseq_results_directory, "norm_counts_", treatment, "_vs_wt.csv"))

    # filter by DMRs-gene_list (if exist)
    if (gene_list_name == "filtered_by_DEGs") {
      RNA <- filter(RNA, gene_id %in% gene_list$gene_id)
    }

    #### Gene expression ####
    # edit table and names
    mto_names_pos <- grep(paste(mto_rnaseq_names, collapse = "|"), names(RNA))
    names(RNA)[mto_names_pos] <- paste0(treatment, ".", 1:3, "_RNA")

    wt_names_pos <- grep(paste(wt_rnaseq_names, collapse = "|"), names(RNA))
    names(RNA)[wt_names_pos] <- paste0("wt.", 1:2, "_RNA")

    RNA <- RNA[, c(1, wt_names_pos, mto_names_pos)]


    #### Methylation ####
    meth_matrix <- meth_matrix[, grep("gene_id|mto|wt", names(meth_matrix))]
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
      mutate(genotype = ifelse(grepl("^wt", Sample), "wt", treatment))

    #####################################################

    tryCatch(
      {
        #### linear model ####
        plot_df$Meth <- as.numeric(plot_df$Meth)
        plot_df$genotype <- as.factor(plot_df$genotype)
        lm_model_0 <- glm.nb(normCounts ~ Meth + genotype, data = plot_df) #  + Meth:genotype
        lm_model <- summary(lm_model_0)

        # fitted vs. residuals
        plot_df$residuals <- residuals(lm_model_0)
        plot_df$fitted <- fitted(lm_model_0)

        #### linear model p-Value ####
        pval_fun <- function(x) {
          ifelse(x <= 0.001, "***",
            ifelse(x <= 0.01, "**",
              ifelse(x <= 0.05, "*",
                "nf"
              )
            )
          )
        }
        M_pval <- pval_fun(lm_model$coefficients["Meth", 4])
        G_pval <- pval_fun(lm_model$coefficients["genotypewt", 4])
        # MnG_pval <- pval_fun(lm_model$coefficients["Meth:genotypewt", 4])

        #####################################################

        ann_cntx_row_number <- which(Rsqr_df$annotation == edit_ann_names(annotation_type) & Rsqr_df$context == context)

        #### p value ####
        p_df[ann_cntx_row_number, ]$Meth <- lm_model$coefficients["Meth", 4]
        p_df[ann_cntx_row_number, ]$genotype <- lm_model$coefficients["genotypewt", 4]

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

        #### reduce points
        # dense_keep <- sample(nrow(plot_df), size = nrow(plot_df) * 0.75)
        # plot_df_sub <- plot_df[dense_keep, ]
        #
        #
        # hb <- hexbin(plot_df$Meth, plot_df$normCounts, xbins = 200, IDs = T, )
        # plot_df$bin_id <- hb@cID
        # max_per_bin <- 200
        # plot_df_sub <- plot_df %>%
        #  group_by(bin_id) %>%
        #  slice_sample(n = min(nrow(plot_df), max_per_bin)) %>%
        #  ungroup()

        if (make_it_light) {
          n_bins <- 1000
          max_per_bin <- 100
          plot_df$x_bin <- cut(plot_df$Meth, breaks = n_bins)
          plot_df <- plot_df %>%
            group_by(x_bin) %>%
            slice_sample(n = min(nrow(plot_df), max_per_bin)) %>%
            ungroup()
        }

        if (context == "CG") {
          title_var <- element_text(size = 16)
        } else {
          title_var <- element_blank()
        }

        if (context == "CHH") {
          xaxis_var <- element_text(size = 16)
        } else {
          xaxis_var <- element_blank()
        }

        if (annotation_type == "transposable_element_gene") {
          yaxis_var <- element_blank()
          yaxis_txt_var <- element_text(size = 13)
        } else {
          yaxis_var <- element_blank()
          yaxis_txt_var <- element_blank()
        }

        if (annotation_type == "promoter") {
          yaxis_var <- element_text(size = 16)
          yaxis_txt_var <- element_text(size = 13)
          context_label <- context
        } else {
          yaxis_var <- element_blank()
          context_label <- ""
        }


        #### plot ####
        p <- ggplot(plot_df, aes(x = Meth, y = normCounts, color = genotype))

        p_points <- p + geom_point(alpha = 0.6, size = 0.3) +
          stat_smooth(method = "glm.nb", formula = y ~ x) +
          scale_y_continuous(
            trans = scales::log1p_trans(),
            breaks = scales::rescale(c(0, 0.01, 0.25, 1), to = range(plot_df$normCounts)),
            labels = function(x) {
              ifelse(x < 10000, label_number()(x), label_scientific(digits = 0)(x))
            }
          ) +
          scale_color_manual(values = c("#bf6828", "gray50")) +
          annotate("text", x = Inf, y = Inf, label = context_label, hjust = 1.25, vjust = 1.5, size = 6) +
          #### updated
          # annotate("text", x = Inf, y = Inf, label = "M", hjust = 4, vjust = 2, size = 3.5) +
          # annotate("text", x = Inf, y = Inf, label = "G", hjust = 4.75, vjust = 4, size = 3.5) +
          # annotate("text", x = Inf, y = Inf, label = paste0("'", M_pval, "'"), hjust = 1.15, vjust = 2, size = 3.5) +
          # annotate("text", x = Inf, y = Inf, label = paste0("'", G_pval, "'"), hjust = 1.15, vjust = 4, size = 3.5) +
          ####
          yo_theme_base() +
          theme(
            legend.position = "none",
            #### updated
            title = title_var,
            axis.text.x = element_text(size = 13),
            axis.text.y = yaxis_txt_var,
            axis.title.x = xaxis_var,
            axis.title.y = yaxis_var
            ####
          ) +
          labs(
            title = edit_ann_names(annotation_type),
            x = "Methylation level",
            y = "Norm. Counts"
          )

        p_regression <- p + stat_smooth(method = "glm.nb", formula = y ~ x) +
          scale_y_continuous(trans = scales::log1p_trans()) +
          scale_color_manual(values = c("#bf6828", "gray50")) +
          # 'M' and 'G' text
          # annotate("text", x = Inf, y = Inf, label = "M", hjust = 4, vjust = 2, size = 3.5) +
          # annotate("text", x = Inf, y = Inf, label = "G", hjust = 4.75, vjust = 4, size = 3.5) +
          # annotate("text", x = Inf, y = Inf, label = paste0("'", M_pval, "'"), hjust = 1.15, vjust = 2, size = 3.5) +
          # annotate("text", x = Inf, y = Inf, label = paste0("'", G_pval, "'"), hjust = 1.15, vjust = 4, size = 3.5) +
          yo_theme_base() +
          theme(
            legend.position = "none",
            #### updated
            # plot.title = element_text(size = 14),
            title = element_blank(),
            axis.text = element_text(size = 12),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
            # axis.title.y = element_text(size = 14)
            ####
          ) +
          labs(
            title = edit_ann_names(annotation_type),
            x = paste0(context, " Methylation levels"),
            y = "Norm. Counts"
          )

        annotation_point_plots_list[[annotation_type]] <- p_points
        annotation_regression_plot_list[[annotation_type]] <- p_regression


        #### residuals plot ####
        residuals_plots_list[[annotation_type]] <- ggplot(plot_df, aes(x = fitted, y = residuals, color = genotype)) +
          # ggrastr::geom_point_rast(alpha = 0.075, size = 0.25) +
          geom_point(alpha = 0.1, size = 0.25) +
          # geom_smooth() +
          scale_color_manual(values = c("#bf6828", "gray50")) +
          yo_theme_base() +
          theme(legend.position = "none") +
          labs(
            title = edit_ann_names(annotation_type),
            x = "Fitted",
            y = "Residuals"
          )

        # message(paste0("created plot to - '", annotation_type, "' annotations in '", context, "' context successfully"))
      },
      error = function(e) {
        message(paste0("\ncan't use 'glm.nb' and create plot to - ", annotation_type, " annotations in ", context, " context\n\n", e))
        annotation_point_plots_list[[annotation_type]] <- ggplot() +
          theme_void()
        annotation_regression_plot_list[[annotation_type]] <- ggplot() +
          theme_void()
        residuals_plots_list[[annotation_type]] <- ggplot() +
          theme_void()
      }
    )
    cat(paste0("\r", treatment, ": ", round((i_perc / total_perc) * 100, 0), "%  "))
    i_perc <- i_perc + 1
  }

  #### legend ####
  legend_plot_df <- data.frame(
    Meth = runif(100, min = 0, max = 1),
    normCounts = rnorm(100, mean = 10, sd = 2),
    genotype = sample(c("wt", treatment), 100, replace = TRUE)
  )
  p <- ggplot(legend_plot_df, aes(x = Meth, y = normCounts, fill = genotype)) +
    geom_point(alpha = 1, size = 6.5, shape = 21, stroke = 0.5, color = "black") +
    scale_fill_manual(
      values = c("mto1" = "#bf6828", "wt" = "gray50"),
      breaks = c(treatment, "wt") # , labels = c(expression(italic("mto1")), "wt")
    ) +
    ggthemes::theme_base() +
    theme(
      legend.text = element_text(size = 16),
      legend.position = "left"
    ) + # Increase legend text size and position it at the top left
    labs(fill = "")

  legend_p_0 <- ggplotGrob(p)$grobs[[which(sapply(ggplotGrob(p)$grobs, function(x) x$name) == "guide-box")]]
  legend_p <- ggplot() +
    theme_void() +
    annotation_custom(legend_p_0, xmin = 0, xmax = 0.5, ymin = 0.5, ymax = 1)
  # legend_p = grid::grid.draw(legend_p)

  ### list of all plots from contexts ###
  all_cntx_points[[context]] <- annotation_point_plots_list
  all_cntx_regression[[context]] <- annotation_regression_plot_list
  all_cntx_residuals[[context]] <- residuals_plots_list
}
cat("\n")

#### save all plot as one, for each context ####
### main plot - with points
png(paste0(path_2_save.2, "lm.stats.plot.", treatment, ".points.png"), width = 18, height = 7.5, units = "in", res = 300, family = "serif")
do.call(multiplot, c(
  c(all_cntx_points$CG, all_cntx_points$CHG, all_cntx_points$CHH),
  list(
    layout = matrix(1:18, nrow = 3, ncol = 6, byrow = TRUE),
    heights = c(1, 0.9, 1),
    widths = c(1, 0.8, 0.8, 0.8, 0.8, 0.9)
  )
))
dev.off()


################################################################
################################################################

### p values data frame
write.csv(p_df, paste0(path_2_save.2, "p_values_", treatment, ".csv"), row.names = FALSE)

### theta bar plot
theta_df$annotation <- factor(theta_df$annotation, levels = unique(theta_df$annotation))
theta_bar_plot <- ggplot(theta_df, aes(x = annotation, y = theta_value, fill = context)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = theta_value - theta_se, ymax = theta_value + theta_se), width = 0.2, position = position_dodge(0.8)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  labs(
    x = "",
    y = "ϴ",
    fill = "Context"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

svg(paste0(path_2_save.2, "theta_values_", treatment, ".svg"), width = 5, height = 3, family = "serif")
print(theta_bar_plot)
dev.off()

### R^2 bar plot
Rsqr_df$annotation <- factor(Rsqr_df$annotation, levels = unique(Rsqr_df$annotation))
end_brk_1 <- max(Rsqr_df[Rsqr_df$annotation != "TEGs", ]$Rsqr_value)
start_brk_2 <- min(c(Rsqr_df[Rsqr_df$annotation == "TEGs", ]$Rsqr_value))

Rsqr_bar_plot <- ggplot(Rsqr_df, aes(x = annotation, y = Rsqr_value, fill = context)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "black") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  labs(
    x = "",
    y = "R²",
    fill = "Context"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  ggbreak::scale_y_break(
    c(end_brk_1 * 1.075, start_brk_2 * 0.925),
    space = 0.055,
    scales = 0.5
  )

svg(paste0(path_2_save.2, "Rsqr_values_", treatment, ".svg"), width = 5, height = 3, family = "serif")
print(Rsqr_bar_plot)
dev.off()

### Pseudo-R^2 bar plot
ps_Rsqr_df$annotation <- factor(ps_Rsqr_df$annotation, levels = unique(ps_Rsqr_df$annotation))
end_brk_1 <- max(ps_Rsqr_df[ps_Rsqr_df$annotation != "TEGs", ]$Rsqr_value)
start_brk_2 <- min(c(ps_Rsqr_df[ps_Rsqr_df$annotation == "TEGs", ]$Rsqr_value))

ps_Rsqr_bar_plot <- ggplot(ps_Rsqr_df, aes(x = annotation, y = Rsqr_value, fill = context)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "black") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  labs(
    x = "",
    y = "Pseudo-R²",
    fill = "Context"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  ggbreak::scale_y_break(
    c(end_brk_1 * 1.075, start_brk_2 * 0.925),
    space = 0.055,
    scales = 0.5
  )

svg(paste0(path_2_save.2, "Pseudo_Rsqr_values_", treatment, ".svg"), width = 5, height = 3, family = "serif")
print(ps_Rsqr_bar_plot)
dev.off()


################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################









#### main plot - with points
# svg(paste0(path_2_save.3, "lm.stats.plot.", context, ".", treatment, ".points.svg"), width = 16, height = 2.5, family = "serif")
# multiplot(
#  all_cntx_points$CG,
#  all_cntx_points$CHG,
#  all_cntx_points$CHH,
#  # legend_p,
#  cols = 6,
#  rows = 3
# )
# dev.off()
#
#### main plot - regression
# svg(paste0(path_2_save.3, "lm.stats.plot.", context, ".", treatment, ".regression.svg"), width = 16, height = 2.5, family = "serif") # width = 12, height = 5
# multiplot(
#  annotation_regression_plot_list[[1]],
#  annotation_regression_plot_list[[2]],
#  annotation_regression_plot_list[[3]],
#  annotation_regression_plot_list[[4]],
#  annotation_regression_plot_list[[5]],
#  annotation_regression_plot_list[[6]],
#  # annotation_regression_plot_list[[7]],
#  # legend_p,
#  cols = 6
# )
# dev.off()
#
#### residuals plot
# svg(paste0(path_2_save.3, "lm.residuals.plot.", context, ".", treatment, ".svg"), width = 16, height = 2.5, family = "serif")
# multiplot(
#  residuals_plots_list[[1]],
#  residuals_plots_list[[2]],
#  residuals_plots_list[[3]],
#  residuals_plots_list[[4]],
#  residuals_plots_list[[5]],
#  residuals_plots_list[[6]],
#  # residuals_plots_list[[7]],
#  # legend_p,
#  cols = 6
# )
# dev.off()
#

