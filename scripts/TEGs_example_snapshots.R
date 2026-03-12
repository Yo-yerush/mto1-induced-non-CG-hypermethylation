#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(Gviz)
  library(grid)
})

root_dir <- normalizePath(".", mustWork = TRUE)
output_dir <- file.path(root_dir, "TEGs_snapshots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

targets <- c("AT5G19097", "AT4G15590")
flank_size <- 4000L
grid_step <- 25L
smooth_radius <- 150L

cx_files <- c(
  wt_1 = "methylome_CX_reports/wt_1_bismark_pe.CX_report.txt.gz",
  wt_2 = "methylome_CX_reports/wt_2_bismark_pe.CX_report.txt.gz",
  mto1_1 = "methylome_CX_reports/mto1_1_bismark_pe.CX_report.txt.gz",
  mto1_2 = "methylome_CX_reports/mto1_2_bismark_pe.CX_report.txt.gz",
  mto1_3 = "methylome_CX_reports/mto1_3_bismark_pe.CX_report.txt.gz"
)

rna_bam_files <- c(
  wt_1 = "RNAseq_sorted_bam_files/wt_1_sorted.bam",
  wt_2 = "RNAseq_sorted_bam_files/wt_2_sorted.bam",
  wt_3 = "RNAseq_sorted_bam_files/wt_3_sorted.bam",
  mto1_1 = "RNAseq_sorted_bam_files/mto1_1_sorted.bam",
  mto1_2 = "RNAseq_sorted_bam_files/mto1_2_sorted.bam",
  mto1_3 = "RNAseq_sorted_bam_files/mto1_3_sorted.bam"
)

ctx_colors <- c(CG = "#33a02c", CHG = "#4a9744", CHH = "#50754d")
sample_colors <- c(
  wt_1 = "#565656",
  wt_2 = "#717171",
  mto1_1 = "#d95f02",
  mto1_2 = "#fc8d62",
  mto1_3 = "#b15928"
)
group_colors <- c(wt = "#6c6c6c", mto1 = "#d95f02")
context_ymax <- c(CG = 100, CHG = 80, CHH = 60)

load_gff_table <- function(gff_path) {
  dt <- fread(
    cmd = sprintf("grep -v '^#' %s", shQuote(gff_path)),
    sep = "\t",
    header = FALSE,
    quote = "",
    col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  )
  dt[, ID := fifelse(grepl("(^|;)ID=", attributes), sub(".*(?:^|;)ID=([^;]+).*", "\\1", attributes, perl = TRUE), NA_character_)]
  dt[, Name := fifelse(grepl("(^|;)Name=", attributes), sub(".*(?:^|;)Name=([^;]+).*", "\\1", attributes, perl = TRUE), NA_character_)]
  dt[, Parent := fifelse(grepl("(^|;)Parent=", attributes), sub(".*(?:^|;)Parent=([^;]+).*", "\\1", attributes, perl = TRUE), NA_character_)]
  dt
}

get_target_model <- function(gff_dt, target_id) {
  gene_row <- gff_dt[
    type %in% c("gene", "transposable_element_gene") &
      (ID == target_id | Name == target_id)
  ]
  if (nrow(gene_row) != 1L) {
    stop(sprintf("Could not resolve a unique gene record for %s", target_id))
  }
  tx_row <- gff_dt[
    type %in% c("mRNA", "mRNA_TE_gene") &
      Parent == target_id
  ]
  if (nrow(tx_row) < 1L) {
    stop(sprintf("Could not resolve a transcript model for %s", target_id))
  }
  list(
    chromosome = gene_row$seqid[1],
    gene_start = as.integer(gene_row$start[1]),
    gene_end = as.integer(gene_row$end[1]),
    strand = gene_row$strand[1],
    transcript_id = tx_row$ID[1]
  )
}

read_dmrs <- function(path, chromosome, start_pos, end_pos) {
  dt <- fread(path)
  dt[
    seqnames == chromosome &
      end >= start_pos &
      start <= end_pos
  ]
}

extract_cx_windows <- function(path, windows_dt) {
  clauses <- apply(windows_dt, 1, function(row) {
    sprintf("($1==\"%s\" && $2>=%d && $2<=%d)", row[["chromosome"]], as.integer(row[["window_start"]]), as.integer(row[["window_end"]]))
  })
  cmd <- sprintf(
    "gzip -cd %s | awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"} %s {print}'",
    shQuote(path),
    paste(clauses, collapse = " || ")
  )
  dt <- fread(
    cmd = cmd,
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "pos", "strand", "meth", "unmeth", "context", "trinuc")
  )
  if (!nrow(dt)) {
    dt[, coverage := numeric()]
    dt[, pct := numeric()]
    return(dt)
  }
  dt[, coverage := meth + unmeth]
  dt[, pct := fifelse(coverage > 0, 100 * meth / coverage, NA_real_)]
  dt
}

smooth_methylation <- function(dt, context_label, start_pos, end_pos, step = grid_step, radius = smooth_radius) {
  grid <- data.table(pos = seq.int(start_pos, end_pos, by = step))
  subset_dt <- dt[context == context_label & coverage > 0]
  if (!nrow(subset_dt)) {
    grid[, score := NA_real_]
    return(grid)
  }

  scores <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    distances <- abs(subset_dt$pos - grid$pos[i])
    keep <- distances <= radius
    if (!any(keep)) {
      scores[i] <- NA_real_
      next
    }
    kernel <- exp(-0.5 * (distances[keep] / (radius / 2))^2)
    total_weight <- kernel * subset_dt$coverage[keep]
    scores[i] <- sum(total_weight * subset_dt$pct[keep]) / sum(total_weight)
  }
  grid[, score := scores]
  grid
}

make_signal_track <- function(position_dt, chromosome, score_col, signal_color, track_name, ylim_vals, track_type = "l", line_width = 0.85) {
  gr <- GRanges(seqnames = chromosome, ranges = IRanges(position_dt$pos, width = 1L))
  DataTrack(
    range = gr,
    data = position_dt[[score_col]],
    genome = "TAIR10",
    chromosome = chromosome,
    type = track_type,
    ylim = ylim_vals,
    col = signal_color,
    fill = signal_color,
    lwd = line_width,
    name = track_name,
    na.rm = TRUE,
    cex.title = 0.72,
    cex.axis = 0.5,
    col.axis = "#444444",
    background.title = "#f3f3f3"
  )
}

make_overlay_track <- function(track_list, track_name) {
  ov <- OverlayTrack(trackList = track_list, name = track_name)
  displayPars(ov) <- list(background.title = "#f3f3f3", cex.title = 0.72)
  ov
}

make_dmr_track <- function(dt, chromosome, context_label) {
  if (!nrow(dt)) {
    return(
      AnnotationTrack(
        range = GRanges(chromosome, IRanges(1, 1))[0],
        genome = "TAIR10",
        chromosome = chromosome,
        name = sprintf("%s DMR", context_label),
        fill = ctx_colors[[context_label]],
        col = ctx_colors[[context_label]],
        background.title = "#f3f3f3",
        cex.title = 0.72,
        just.group = "left"
      )
    )
  }
  gr <- GRanges(
    seqnames = dt$seqnames,
    ranges = IRanges(dt$start, dt$end),
    direction = ifelse(dt$direction > 0, "gain", "loss")
  )
  AnnotationTrack(
    range = gr,
    genome = "TAIR10",
    chromosome = chromosome,
    name = sprintf("%s DMR", context_label),
    stacking = "squish",
    fill = ctx_colors[[context_label]],
    col = ctx_colors[[context_label]],
    background.title = "#f3f3f3",
    cex.title = 0.72,
    just.group = "left",
    group = mcols(gr)$direction
  )
}

draw_gene_direction <- function(target_row) {
  panel_left <- 0.18
  panel_right <- 0.995
  gene_x0 <- panel_left + ((target_row$gene_start - target_row$window_start) / (target_row$window_end - target_row$window_start)) * (panel_right - panel_left)
  gene_x1 <- panel_left + ((target_row$gene_end - target_row$window_start) / (target_row$window_end - target_row$window_start)) * (panel_right - panel_left)
  arrow_label <- if (target_row$strand == "-") "<" else ">"
  arrow_count <- max(5L, min(9L, floor((gene_x1 - gene_x0) / 0.05)))
  x_vals <- seq(gene_x0 + 0.012, gene_x1 - 0.012, length.out = arrow_count)
  grid.text(
    label = rep(arrow_label, length(x_vals)),
    x = unit(x_vals, "npc"),
    y = unit(0.712, "npc"),
    gp = gpar(col = "white", fontsize = 10, fontface = "bold")
  )
}

get_mapped_reads <- function(bam_path) {
  idx_dt <- fread(
    text = system2("samtools", c("idxstats", bam_path), stdout = TRUE),
    sep = "\t",
    header = FALSE,
    col.names = c("ref", "length", "mapped", "unmapped")
  )
  idx_dt[ref != "*", sum(mapped)]
}

extract_rna_profile <- function(bam_path, transcript_id, strand_value, gene_start, gene_end, window_start, window_end, mapped_reads) {
  depth_lines <- system2("samtools", c("depth", "-aa", "-r", transcript_id, bam_path), stdout = TRUE)
  grid <- data.table(pos = seq.int(window_start, window_end, by = grid_step))
  if (!length(depth_lines)) {
    grid[, score := 0]
    return(grid)
  }
  depth_dt <- fread(
    text = depth_lines,
    sep = "\t",
    header = FALSE,
    col.names = c("ref", "tx_pos", "depth")
  )
  if (strand_value == "-") {
    depth_dt[, genomic_pos := gene_end - tx_pos + 1L]
  } else {
    depth_dt[, genomic_pos := gene_start + tx_pos - 1L]
  }
  depth_dt[, norm_depth := depth * 1e6 / mapped_reads]
  depth_dt <- depth_dt[genomic_pos >= window_start & genomic_pos <= window_end]
  depth_dt[, pos := window_start + ((genomic_pos - window_start) %/% grid_step) * grid_step]
  binned <- depth_dt[, .(score = mean(norm_depth, na.rm = TRUE)), by = pos]
  merged <- merge(grid, binned, by = "pos", all.x = TRUE)
  merged[is.na(score), score := 0]
  merged
}

build_methyl_context_track <- function(target_id, chromosome, window_start, window_end, context_label, cx_cache) {
  track_list <- lapply(names(cx_files), function(sample_name) {
    smoothed <- smooth_methylation(cx_cache[[sample_name]][[target_id]], context_label, window_start, window_end)
    make_signal_track(
      smoothed,
      chromosome,
      "score",
      sample_colors[[sample_name]],
      sprintf("\t%s  \n\t\tmeth. (%%)", context_label),
      c(0, context_ymax[[context_label]])
    )
  })
  make_overlay_track(track_list, sprintf("\t%s  \n\t\tmeth. (%%)", context_label))
}

plot_target <- function(target_id, dmr_paths, windows_dt, cx_cache, rna_cache) {
  target_name <- target_id
  target_row <- windows_dt[windows_dt[["target_id"]] == target_name, ][1]
  chromosome <- target_row$chromosome
  gene_start <- target_row$gene_start
  gene_end <- target_row$gene_end
  window_start <- target_row$window_start
  window_end <- target_row$window_end

  axis_track <- GenomeAxisTrack(littleTicks = TRUE, cex = 0.7, col = "#444444")
  target_track <- AnnotationTrack(
    range = GRanges(chromosome, IRanges(gene_start, gene_end)),
    genome = "TAIR10",
    chromosome = chromosome,
    name = "Genes",
    fill = "#111111",
    col = "#111111",
    background.title = "#f3f3f3",
    shape = "box",
    just.group = "left"
  )

  dmr_tracks <- lapply(names(dmr_paths), function(ctx) {
    make_dmr_track(read_dmrs(dmr_paths[[ctx]], chromosome, window_start, window_end), chromosome, ctx)
  })
  names(dmr_tracks) <- names(dmr_paths)

  rna_ylim <- c(0, 0.605) # c(0, max(rna_cache[[target_id]][, .(wt, mto1)], na.rm = TRUE) * 1.05 + 1e-06)
  rna_track <- make_overlay_track(
    list(
      make_signal_track(rna_cache[[target_id]], chromosome, "wt", grDevices::adjustcolor(group_colors[["wt"]], alpha.f = 0.7), "\tRNA  \n\t\tcov. (RPM)", rna_ylim, track_type = "h", line_width = 1.1),
      make_signal_track(rna_cache[[target_id]], chromosome, "mto1", grDevices::adjustcolor(group_colors[["mto1"]], alpha.f = 0.7), "\tRNA  \n\t\tcov. (RPM)", rna_ylim, track_type = "h", line_width = 1.1)
    ),
    "\tRNA  \n\t\tcov. (RPM)"
  )

  methyl_tracks <- lapply(c("CG", "CHG", "CHH"), function(ctx) {
    build_methyl_context_track(target_id, chromosome, window_start, window_end, ctx, cx_cache)
  })

  track_list <- c(list(axis_track, target_track, rna_track), dmr_tracks, methyl_tracks)
  sizes <- c(0.28, 0.042, 0.3, 0.06, 0.06, 0.06, 0.22, 0.22, 0.22)

  svg_path <- file.path(output_dir, sprintf("%s_results_snapshot.svg", target_id))
  png_path <- file.path(output_dir, sprintf("%s_results_snapshot.png", target_id))

  for (out_path in c(svg_path, png_path)) {
    if (grepl("\\.svg$", out_path)) {
      svg(out_path, width = 10.08, height = 7.92, family = "serif")
    } else {
      png(out_path, width = 1500, height = 1180, res = 180, family = "serif")
    }
    plotTracks(
      track_list,
      from = window_start,
      to = window_end,
      sizes = sizes,
      main = sprintf(
        "%s plus 4kb flanks (%s:%s-%s)",
        target_id,
        chromosome,
        format(window_start, big.mark = ","),
        format(window_end, big.mark = ",")
      ),
      cex.main = 1,
      cex.title = 0.72,
      cex.axis = 0.5,
      fontfamily = "serif",
      fontfamily.title = "serif",
      fontcolor.title = "#111111",
      background.panel = "#ffffff",
      background.title = "#f3f3f3",
      col.title = "#111111",
      title.width = 1.75,
      rotation.title = 0
    )
    draw_gene_direction(target_row)
    dev.off()
  }

  data.table(
    target_id = target_id,
    chromosome = chromosome,
    gene_start = gene_start,
    gene_end = gene_end,
    window_start = window_start,
    window_end = window_end,
    svg = svg_path,
    png = png_path
  )
}

gff_path <- file.path(root_dir, "TAIR10", "TAIR10_GFF3_genes.gff")
dmr_paths <- c(
  CG = file.path(root_dir, "DMRs_files", "DMRs_CG_mto1_vs_wt.csv"),
  CHG = file.path(root_dir, "DMRs_files", "DMRs_CHG_mto1_vs_wt.csv"),
  CHH = file.path(root_dir, "DMRs_files", "DMRs_CHH_mto1_vs_wt.csv")
)

gff_dt <- load_gff_table(gff_path)

windows_dt <- rbindlist(lapply(targets, function(target_id) {
  model <- get_target_model(gff_dt, target_id)
  data.table(
    target_id = target_id,
    chromosome = model$chromosome,
    gene_start = model$gene_start,
    gene_end = model$gene_end,
    strand = model$strand,
    transcript_id = model$transcript_id,
    window_start = max(1L, model$gene_start - flank_size),
    window_end = model$gene_end + flank_size
  )
}))

cx_cache <- list()
for (sample_name in names(cx_files)) {
  sample_dt <- extract_cx_windows(file.path(root_dir, cx_files[[sample_name]]), windows_dt)
  cx_cache[[sample_name]] <- lapply(targets, function(target_id) {
    target_name <- target_id
    target_row <- windows_dt[windows_dt[["target_id"]] == target_name, ][1]
    sample_dt[
      chr == target_row$chromosome &
        pos >= target_row$window_start &
        pos <= target_row$window_end
    ]
  })
  names(cx_cache[[sample_name]]) <- targets
}

rna_norm_factors <- setNames(
  vapply(file.path(root_dir, rna_bam_files), get_mapped_reads, numeric(1)),
  names(rna_bam_files)
)

rna_cache <- list()
for (target_id in targets) {
  target_name <- target_id
  target_row <- windows_dt[windows_dt[["target_id"]] == target_name, ][1]
  sample_profiles <- lapply(names(rna_bam_files), function(sample_name) {
    prof <- extract_rna_profile(
      bam_path = file.path(root_dir, rna_bam_files[[sample_name]]),
      transcript_id = target_row$transcript_id,
      strand_value = target_row$strand,
      gene_start = target_row$gene_start,
      gene_end = target_row$gene_end,
      window_start = target_row$window_start,
      window_end = target_row$window_end,
      mapped_reads = rna_norm_factors[[sample_name]]
    )
    setnames(prof, "score", sample_name)
    prof
  })
  merged_profiles <- Reduce(function(x, y) merge(x, y, by = "pos", all = TRUE), sample_profiles)
  merged_profiles[is.na(merged_profiles)] <- 0
  merged_profiles[, wt := rowMeans(.SD), .SDcols = c("wt_1", "wt_2", "wt_3")]
  merged_profiles[, mto1 := rowMeans(.SD), .SDcols = c("mto1_1", "mto1_2", "mto1_3")]
  rna_cache[[target_id]] <- merged_profiles[, .(pos, wt, mto1)]
}

summaries <- rbindlist(lapply(targets, function(target_id) {
  plot_target(target_id, dmr_paths, windows_dt, cx_cache, rna_cache)
}))

fwrite(summaries, file.path(output_dir, "results_snapshot_summary.tsv"), sep = "\t")
