# run it on darwin (linux server)
library(dplyr)
library(ggplot2)
library(DMRcaller)
library(org.At.tair.db)
library(GenomicFeatures)
library(plyranges)
library(parallel)
library(data.table)

source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/Genes_metaPlot_fun.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/load_replicates.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/edit_TE_file.R")

TE_file_path <- "https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/TAIR10_Transposable_Elements.txt"
samples_path_df <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/Methylome.At/samples_table/samples_table_mto1.txt"
n.cores <- 30

#############################################

# load replicates raw 'CX' data
var_table <- read.csv(samples_path_df, header = F, sep = "\t")
vars_vector <- unique(var_table[, 1])
var1_path <- var_table[grep(vars_vector[1], var_table[, 1]), 2]
var2_path <- var_table[grep(vars_vector[2], var_table[, 1]), 2]
var_args <- list(
    list(path = var1_path, name = "wt"),
    list(path = var2_path, name = "mto1")
)
load_vars <- mclapply(var_args, function(x) load_replicates(x$path, n.cores, x$name, T, "CX_report"), mc.cores = 5)
meth_wt <- trimm_and_rename(load_vars[[1]])
meth_mto1 <- trimm_and_rename(load_vars[[2]])

TE_df <- read.csv(TE_file_path, sep = "\t")
te_width <- data.frame(
    te_id = TE_df$Transposon_Name,
    width = (TE_df$Transposon_max_End - TE_df$Transposon_min_Start),
    superfamily = (TE_df$Transposon_Super_Family)
)
long_tes <- te_width[te_width$width >= 4000, 1]
short_tes <- te_width[te_width$width <= 500, 1] # %>% sample(10000) # random IDs

TE_gr <- edit_TE_file(TE_df)
filtered_te_width <- te_width %>% filter(width < 20000)

#############################################

#############################################
## TE methylation levels and size
scatter_te_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots/scatter_plots"
dir.create(scatter_te_path, showWarnings = F)
setwd(scatter_te_path)

# Function to calculate average methylation levels for TEs
calculate_te_methylation <- function(meth_data, TE_gr, context, is.delta = F) {
    meth_data <- meth_data[meth_data$context == context]
    if (!is.delta) {
        meth_data$Proportion <- meth_data$readsM / meth_data$readsN
    }
    # Find overlaps between methylation data and TEs
    overlaps <- findOverlaps(TE_gr, meth_data)
    overlaps_df <- data.frame(
        te_id = queryHits(overlaps),
        meth_idx = subjectHits(overlaps)
    )

    # Get TE IDs and methylation levels for overlapping regions
    overlaps_df$te_id <- TE_gr$gene_id[overlaps_df$te_id]
    overlaps_df$meth_level <- meth_data$Proportion[overlaps_df$meth_idx]

    # Calculate average methylation level for each TE
    te_avg_meth <- overlaps_df %>%
        group_by(te_id) %>%
        summarise(avg_meth = mean(meth_level, na.rm = TRUE), .groups = "drop")
    te_avg_meth$context <- context
    as.data.frame(te_avg_meth) %>% merge(., filtered_te_width, by = "te_id")
}

te_size_plot <- function(x, is.delta = F) {
    if (!is.delta) {
        smooth_p <- geom_smooth(
            method = lm, formula = y ~ splines::bs(x, 3), # "gam",
            se = TRUE,
            linewidth = 0.5
        )
        y_dashed_line <- NULL
        x_dashed_line <- NULL
        color_values <- c("WT" = "#7F7F7F", "mto1" = "#bf6828")
        y_scale <- scale_y_continuous(limits = c(0, 1), expand = c(0,0), breaks = c(0.005, 0.5, 0.995), labels = c(0, 0.5, 1))
    } else {
        smooth_p <- geom_smooth(
            method = lm, formula = y ~ splines::bs(x, 3),
            se = TRUE,
            linewidth = 0.5,
            color = "red4"
        )
        y_dashed_line <- geom_hline(yintercept = 0, linetype = "dashed", color = "black") # , alpha = 0.6)
        x_dashed_line <- geom_vline(xintercept = 4000, linetype = "dashed", color = "gray")
        color_values <- c("delta" = "gray20")
        y_scale <- scale_y_continuous(limits = c(-0.3, 0.3), expand = c(0,0), breaks = c(-0.296, 0, 0.296), labels = c(-0.3, 0, 0.3))
    }

    ggplot(x, aes(x = width, y = avg_meth, color = sample)) +
        x_dashed_line +
        geom_point(alpha = ifelse(is.delta, 0.6, 0.4), size = 0.3, shape = 20) +
        y_dashed_line +
        smooth_p +
        scale_color_manual(values = color_values) +
        labs(
            title = paste(cntx, "context"),
            x = "TE size (Kbp)",
            y = paste0(ifelse(is.delta, "Δ ", ""), "Methylation") # ,
            # color = "Sample"
        ) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            plot.title = element_text(hjust = 0.5, size = 10),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 9)
        ) +
        scale_x_continuous(
            breaks = c(65, 5000, 10000, 15000, 19935),
            labels = seq(0, 20, by = 5),
            limits = c(0, 20000), expand = c(0, 0)
        ) +
        y_scale
}

# Calculate average methylation for both samples
meth_delta <- meth_wt[, 1]
meth_delta$Proportion <- (meth_mto1$readsM / meth_mto1$readsN) - (meth_wt$readsM / meth_wt$readsN)
meth_delta$Proportion[is.nan(meth_delta$Proportion)] <- 0

# plots
for (cntx in c("CG", "CHG", "CHH")) {
    # Combine the data for both samples
    te_meth_wt <- calculate_te_methylation(meth_wt, TE_gr, cntx)
    te_meth_mto1 <- calculate_te_methylation(meth_mto1, TE_gr, cntx)
    combined_data <- rbind(
        data.frame(te_meth_wt, sample = "WT"),
        data.frame(te_meth_mto1, sample = "mto1")
    )

    # delta
    te_meth_delta <- calculate_te_methylation(meth_delta, TE_gr, cntx, T)
    te_meth_delta$sample <- "delta"

    png(paste0("TE_size_scatter_", cntx, ".png"), width = 3500, height = 2500, res = 1200, family = "serif")
    print(te_size_plot(combined_data))
    dev.off()

    png(paste0("TE_size_delta_scatter_", cntx, ".png"), width = 3500, height = 2500, res = 1200, family = "serif")
    print(te_size_plot(te_meth_delta, is.delta = T))
    dev.off()
}

#############################################
## TE methylation levels and distance from centromer plots
distance_te_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots/distance_from_centromere"
dir.create(distance_te_path, showWarnings = F)
setwd(distance_te_path)

# centromers positions
cen_pos <- c(14.845, 3.44, 13.855, 3.13, 11.795) * 1e6
te_distance <- data.frame(
    chr = as.character(seqnames(TE_gr)),
    pos = (as.numeric(start(TE_gr)) + as.numeric(end(TE_gr))) / 2,
    te_id = TE_gr$gene_id,
    centromere = NA
)
for (cen_i in 1:5) {
    te_distance$centromere[te_distance$chr == paste0("Chr", cen_i)] <- cen_pos[cen_i]
}
te_distance$distance <- abs(te_distance$centromere - te_distance$pos) / 1e6
all_cx_dis <- rbind(
    calculate_te_methylation(meth_delta, TE_gr, "CG", T),
    calculate_te_methylation(meth_delta, TE_gr, "CHG", T),
    calculate_te_methylation(meth_delta, TE_gr, "CHH", T)
)
keep_col <- grep("te_id|distance", names(te_distance))
te_distance_merged <- merge(te_distance[, keep_col], all_cx_dis, by = "te_id")

te_distance_cntx <- rbind(
    te_distance_merged %>%
        filter(context == "CG") %>%
        mutate(window = floor(distance)) %>% # floor(distance / 1e6) * 1e6) %>%
        group_by(window, context) %>%
        summarise(
            avg_meth = mean(avg_meth, na.rm = TRUE),
            distance = mean(distance, na.rm = TRUE),
            .groups = "drop"
        ),
    te_distance_merged %>%
        filter(context == "CHG") %>%
        mutate(window = floor(distance)) %>% # floor(distance / 1e6) * 1e6) %>%
        group_by(window, context) %>%
        summarise(
            avg_meth = mean(avg_meth, na.rm = TRUE),
            distance = mean(distance, na.rm = TRUE),
            .groups = "drop"
        ),
    te_distance_merged %>%
        filter(context == "CHH") %>%
        mutate(window = floor(distance)) %>% # floor(distance / 1e6) * 1e6) %>%
        group_by(window, context) %>%
        summarise(
            avg_meth = mean(avg_meth, na.rm = TRUE),
            distance = mean(distance, na.rm = TRUE),
            .groups = "drop"
        )
) %>%
    as.data.frame() %>%
    filter(avg_meth != max(avg_meth, na.rm = TRUE))

plot_colors <- c("#3d53b4", "#3b8f3e", "#bb4949")

te_distance_plot <- ggplot(data = te_distance_cntx, aes(x = distance, y = avg_meth, color = context, group = context)) +
    geom_line(linewidth = 0.85) +
    theme_bw() + # theme_classic() +
    labs(
        title = " ",
        x = "Distance from centromer (Mbp)",
        y = "Δ Methylation"
    ) +
    theme(
        legend.position = "none",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 9)
    ) +
    scale_x_continuous(
        limits = c(0, 15), # max(-te_distance_cntx$distance)),
        breaks = c(0.06, 5, 10, 14.94),
        labels = seq(0, 15, by = 5),
        expand = c(0, 0)
        ) +
    scale_y_continuous(
        limits = c(-0.001, max(te_distance_cntx$avg_meth)),
        breaks = c(0, 0.02, 0.04)
        ) +
    annotate("text",
        x = 12.75, # 12.25,
        y = max(te_distance_cntx$avg_meth) * 0.98,
        label = c("CG", "\nCHG", "\n\nCHH"),
        hjust = 0, vjust = 0.75, size = 3.25,
        color = plot_colors, fontface = "bold"
    )

png(paste0("centromere_distance.png"), width = 3500, height = 2500, res = 1200, family = "serif")
print(te_distance_plot)
dev.off()


#############################################

#############################################
## run long and short TEs
metaPlot_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots"
long_te_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots/long_TEs"
short_te_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots/short_TEs"
dir.create(metaPlot_path, showWarnings = F)
dir.create(long_te_path, showWarnings = F)
dir.create(short_te_path, showWarnings = F)

setwd(long_te_path)
Genes_metaPlot(meth_wt, meth_mto1, "wt", "mto1", TE_gr, long_tes, 6, n.cores, is_TE = T)

setwd(short_te_path)
Genes_metaPlot(meth_wt, meth_mto1, "wt", "mto1", TE_gr, short_tes, 6, n.cores, is_TE = T)

########################
## run all TEs
all_tes <- te_width[, 1]
all_te_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots/all_TEs"
dir.create(all_te_path, showWarnings = F)

setwd(all_te_path)
Genes_metaPlot(meth_wt, meth_mto1, "wt", "mto1", TE_gr, all_tes, 6, n.cores, is_TE = T)

########################
## run each super-families group
SF_list <- list(
    Gypsy = te_width %>% filter(grepl("Gypsy", superfamily)) %>% .[, 1],
    Copia = te_width %>% filter(grepl("Copia", superfamily)) %>% .[, 1],
    LINE = te_width %>% filter(grepl("LINE", superfamily)) %>% .[, 1],
    Helitron = te_width %>% filter(grepl("Helitron", superfamily)) %>% .[, 1],
    TIR = te_width %>% filter(grepl("DNA", superfamily)) %>% .[, 1],
    SINE = te_width %>% filter(grepl("SINE|Rath", superfamily)) %>% .[, 1],
    Unassigned = te_width %>% filter(grepl("Unassigned", superfamily)) %>% .[, 1]
)
superfamilies_te_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots/superfamilies_TEs/"
dir.create(superfamilies_te_path, showWarnings = F)

for (i_te_list in seq(length(SF_list))) {
    new_dir_SF <- paste0(superfamilies_te_path, names(SF_list)[i_te_list])
    dir.create(new_dir_SF, showWarnings = F)
    setwd(new_dir_SF)
    Genes_metaPlot(meth_wt, meth_mto1, "wt", "mto1", TE_gr, SF_list[[i_te_list]], 6, n.cores, is_TE = T)
}

########################
setwd(metaPlot_path)

# Write summary information to a text file - thanks GPT
summary_info <- paste0(
    "Transposable Elements (TEs) Methylation Analysis Summary\n",
    "======================================================\n\n",
    "Analysis Date: ", Sys.Date(), "\n",
    "Analysis Time: ", Sys.time(), "\n\n",
    "Data Sources:\n",
    "- TE annotation file: ", TE_file_path, "\n",
    "- Samples table: ", samples_path_df, "\n\n",
    "Sample Information:\n",
    "- Wild type (wt): ", length(var1_path), " replicates\n",
    "- mto1 mutant: ", length(var2_path), " replicates\n\n",
    "TE Categories Analyzed:\n",
    "- Total TEs in dataset: ", nrow(TE_df), "\n",
    "- Long TEs (≥4000 bp): ", length(long_tes), " TEs\n",
    "- Short TEs (≤500 bp): ", length(short_tes), " TEs\n",
    "- All TEs analyzed: ", length(all_tes), " TEs\n\n",
    "Superfamilies Analyzed:\n",
    "- Gypsy: ", nrow(SF_list[[1]]), " TEs\n",
    "- Copia: ", nrow(SF_list[[2]]), " TEs\n",
    "- LINE: ", nrow(SF_list[[3]]), " TEs\n",
    "- Helitron: ", nrow(SF_list[[4]]), " TEs\n",
    "- TIR (DNA): ", nrow(SF_list[[5]]), " TEs\n",
    "- SINE: ", nrow(SF_list[[6]]), " TEs\n",
    "- Unassigned: ", nrow(SF_list[[7]]), " TEs\n\n",
    "Analysis Parameters:\n",
    "- Number of cores used: ", n.cores, "\n",
    "- Bin size for metaplots: 20\n",
    "- Minimum reads per cytosine: 6\n",
    "- Context analyzed: CX (all cytosines)\n\n",
    "Output Directories:\n",
    "- Main output: ", metaPlot_path, "\n",
    "- Long TEs: ", long_te_path, "\n",
    "- Short TEs: ", short_te_path, "\n",
    "- All TEs: ", all_te_path, "\n",
    "- Superfamilies: ", superfamilies_te_path, "\n"
)
writeLines(summary_info, file.path(metaPlot_path, "analysis_summary.txt"))
#
