##############
# discrete
##############

setwd("/ix1/ctseng/hac377/cnv_calling/letter")
df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)

plot_bin_heatmap <- function(plot_df,
                             method_cols,
                             method_display_names = NULL,
                             sample_name = NULL,
                             output_plot_heatmap = "bin_heatmap.png",
                             width = 20, height = 7, dpi = 300, text.prop = 1) {
  
  stopifnot(all(c("chr", "start_position") %in% colnames(plot_df)))
  stopifnot(all(method_cols %in% colnames(plot_df)))
  
  df <- plot_df[, c("chr", "start_position", method_cols), drop = FALSE]
  df <- as.data.frame(df)
  
  chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
  chr_levels <- chr_levels[chr_levels %in% unique(as.character(df$chr))]
  df$chr <- factor(as.character(df$chr), levels = chr_levels)
  
  df <- df %>% arrange(chr, start_position)
  df$counted_pos <- seq_len(nrow(df))
  
  chr_boundaries <- df %>%
    group_by(chr) %>%
    summarize(
      start_chr = min(counted_pos),
      mean_chr = mean(counted_pos),
      .groups = "drop"
    )
  
  if (is.null(method_display_names)) {
    method_display_names <- setNames(method_cols, method_cols)
  }
  
  stopifnot(all(method_cols %in% names(method_display_names)))
  
  plot_data <- reshape2::melt(
    df,
    id.vars = c("chr", "start_position", "counted_pos"),
    measure.vars = method_cols
  )
  
  plot_data$variable <- method_display_names[as.character(plot_data$variable)]
  plot_data$variable <- factor(
    plot_data$variable,
    levels = method_display_names[method_cols]
  )
  
  ## robust 3-state conversion
  plot_data$value <- as.character(plot_data$value)
  plot_data$value <- trimws(tolower(plot_data$value))
  
  plot_data$value <- dplyr::recode(
    plot_data$value,
    "-1" = "loss",
    "0"  = "base",
    "1"  = "gain",
    "loss" = "loss",
    "base" = "base",
    "neutral" = "base",
    "gain" = "gain",
    .default = NA_character_
  )
  
  plot_data$value <- factor(plot_data$value, levels = c("loss", "base", "gain"))
  
  g <- ggplot(plot_data, aes(x = counted_pos, y = variable, fill = value)) +
    geom_tile() +
    theme_bw() +
    scale_fill_manual(
      name = NULL,
      values = c(
        loss = "#BDAAD9",
        base = "white",
        gain = "#E1AB9E"
      ),
      breaks = c("loss", "base", "gain"),
      labels = c("loss", "base", "gain"),
      na.value = "grey80",
      guide = guide_legend(
        direction = "horizontal",
        title = NULL,
        label.position = "bottom"
      )
    ) +
    labs(title = sample_name) +
    xlab("Chromosome position") +
    ylab("Method") +
    geom_vline(xintercept = chr_boundaries$start_chr) +
    scale_x_continuous(
      breaks = chr_boundaries$mean_chr,
      labels = chr_boundaries$chr
    ) +
    scale_y_discrete(limits = rev) +
    coord_cartesian(xlim = c(1, max(plot_data$counted_pos)), expand = FALSE) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 16 * text.prop),
      legend.box = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.spacing.x = unit(0.6, "cm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      text = element_text(size = 22 * text.prop),
      plot.title = element_text(hjust = 0, face = "bold", size = 22 * text.prop)
    )
  
  ggsave(output_plot_heatmap, plot = g, width = width, height = height, dpi = dpi)
  return(g)
}


########################
## methods
########################
method_cols <- c(
  "aneufinder_original_state",
  "scevan_original_state",
  "copykat_original_state",
  "numbat_original_state",
  "infercnv_original_state",
  "conicsmat_original_state",
  "casper_original_state"
)

method_display_names <- c(
  aneufinder_original_state = "Aneufinder_Schmid_R",
  scevan_original_state = "SCEVAN (Expr)",
  copykat_original_state = "copyKat",
  numbat_original_state = "Numbat(CNV)",
  infercnv_original_state = "InferCNV(Expr)",
  conicsmat_original_state = "CONICsmat",
  casper_original_state = "CaSpER"
)


for (sampleID in df$sample) {
  cat("Plotting sample:", sampleID, "\n")
  
  plot_df <- readRDS(
    paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor2_3/",
      sampleID, "_combined_methods.rds"
    )
  )
  
  plot_bin_heatmap(
    plot_df = plot_df,
    method_cols = method_cols,
    method_display_names = method_display_names,
    sample_name = sampleID,
    output_plot_heatmap = paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor2_3/Karyogram/",
      sampleID, "_Karyogram_state.png"
    )
  )
}


############################
# ratio 
############################
setwd("/ix1/ctseng/hac377/cnv_calling/letter")
df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)

plot_bin_heatmap <- function(plot_df,
                             method_cols,
                             method_display_names = NULL,
                             sample_name = NULL,
                             output_plot_heatmap = "bin_heatmap.png",
                             cap_val = 5,
                             width = 20, height = 8, dpi = 300, text.prop = 1) {
  ## plot_df must already contain:
  ## chr, start_position, and method columns
  
  stopifnot(all(c("chr", "start_position") %in% colnames(plot_df)))
  stopifnot(all(method_cols %in% colnames(plot_df)))
  
  df <- plot_df[, c("chr", "start_position", method_cols), drop = FALSE]
  df <- as.data.frame(df)
  
  if (!is.factor(df$chr)) {
    chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
    chr_levels <- chr_levels[chr_levels %in% unique(as.character(df$chr))]
    df$chr <- factor(df$chr, levels = chr_levels)
  }
  
  df <- df %>% arrange(chr, start_position)
  df$counted_pos <- seq_len(nrow(df))
  
  chr_boundaries <- df %>%
    group_by(chr) %>%
    summarize(
      start_chr = min(counted_pos),
      mean_chr = mean(counted_pos),
      .groups = "drop"
    )
  
  if (is.null(method_display_names)) {
    method_display_names <- setNames(method_cols, method_cols)
  }
  
  plot_data <- reshape2::melt(
    df,
    id.vars = c("chr", "start_position", "counted_pos"),
    measure.vars = method_cols
  )
  
  plot_data$variable <- method_display_names[as.character(plot_data$variable)]
  plot_data$variable <- factor(
    plot_data$variable,
    levels = method_display_names[method_cols]
  )
  
  g <- ggplot(plot_data, aes(x = counted_pos, y = variable, fill = value)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradient2(
      name = NULL,
      low = "darkblue",
      mid = "white",
      high = "darkred",
      midpoint = 0,
      limits = c(-cap_val, cap_val),
      breaks = c(-cap_val, 0, cap_val),
      labels = c("loss", "base", "gain"),
      guide = guide_colorbar(
        direction = "horizontal",
        title = NULL,
        label.position = "bottom",
        barwidth = unit(5, "cm"),
        barheight = unit(0.5, "cm"),
        ticks = FALSE,
        frame.colour = NA
      )
    ) +
    labs(title = sample_name) +
    xlab("Chromosome position") +
    ylab("Method") +
    geom_vline(xintercept = chr_boundaries$start_chr) +
    scale_x_continuous(
      breaks = chr_boundaries$mean_chr,
      labels = chr_boundaries$chr
    ) +
    scale_y_discrete(limits = rev) +
    coord_cartesian(xlim = c(1, max(plot_data$counted_pos)), expand = FALSE) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 16 * text.prop),
      legend.box = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.spacing.x = unit(0.6, "cm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      text = element_text(size = 21 * text.prop),
      plot.title = element_text(hjust = 0, face = "bold", size = 22 * text.prop)
    )
  
  ggsave(output_plot_heatmap, plot = g, width = width, height = height, dpi = dpi)
  return(g)
}


########################
## methods
########################
# method_cols <- c(
#   "aneufinder_original",
#   "aneufinder_ploidy",
#   "infercnv",
#   "copykat",
#   "numbat",
#   "scevan",
#   "casper",
#   "conicsmat"
# )
# 
# method_display_names <- c(
#   aneufinder_original = "Schmid",
#   aneufinder_ploidy = "AneuFinder_ploidy",
#   infercnv  = "inferCNV",
#   copykat   = "CopyKAT",
#   numbat    = "Numbat",
#   scevan    = "SCEVAN",
#   casper    = "CaSpER",
#   conicsmat = "CONICSmat"
# )


method_cols <- c(
  "aneufinder_original",
  "infercnv",
  "copykat",
  "numbat",
  "scevan",
  "casper",
  "conicsmat"
)

method_display_names <- c(
  aneufinder_original = "Aneufinder_Schmid_R",
  infercnv  = "inferCNV (Expr)",
  copykat   = "copyKat",
  numbat    = "Numbat (CNV)",
  scevan    = "SCEVAN (Expr)",
  casper    = "CaSpER",
  conicsmat = "CONICSmat"
)


cap_val <- 5

for (sampleID in df$sample) {
  cat("Plotting sample:", sampleID, "\n")
  
  plot_df <- readRDS(
    paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/ratio/",
      sampleID, "_combined_methods_ratio.rds"
    )
  )
  
  plot_df[, method_cols] <- lapply(plot_df[, method_cols], function(x) {
    pmax(pmin(x, cap_val), -cap_val)
  })
  
  plot_bin_heatmap(
    plot_df = plot_df,
    method_cols = method_cols,
    method_display_names = method_display_names,
    sample_name = sampleID,
    output_plot_heatmap = paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/ratio/Karyogram/",
      sampleID, "_Karyogram_ratio.png"
    ),
    cap_val = cap_val
  )
}



########################
## gt aneufinder only
########################
method_cols_gt <- c(
  "aneufinder_original",
  "aneufinder_ploidy"
)

method_display_names_gt <- c(
  aneufinder_original = "Aneufinder_Schmid_R",
  aneufinder_ploidy = "Aneufinder_Karyotyping"
)

cap_val <- 5

for (sampleID in df$sample) {
  cat("Plotting GT-only sample:", sampleID, "\n")
  
  plot_df <- readRDS(
    paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/Karyogram_gt/",
      sampleID, "_combined_gt.rds"
    )
  )
  
  plot_df[, method_cols_gt] <- lapply(plot_df[, method_cols_gt], function(x) {
    pmax(pmin(x, cap_val), -cap_val)
  })
  
  plot_bin_heatmap(
    plot_df = plot_df,
    method_cols = method_cols_gt,
    method_display_names = method_display_names_gt,
    sample_name = sampleID,
    output_plot_heatmap = paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/Karyogram_gt/aneufinder_only/",
      sampleID, "_Karyogram_gt.png"
    ),
    cap_val = cap_val,
    width = 15, height = 3, text.prop = 0.8
  )
}








# ########################
# ## gt
# ########################
# method_cols_gt <- c(
#   "aneufinder_original",
#   "aneufinder_ploidy",
#   "copykit_state",
#   "copykit_intensity"
# )
# 
# method_display_names_gt <- c(
#   aneufinder_original = "Schmid",
#   aneufinder_ploidy = "AneuFinder_ploidy",
#   copykit_state = "CopyKit_state",
#   copykit_intensity   = "CopyKit_intensity"
# )
# 
# cap_val <- 5
# 
# for (sampleID in df$sample) {
#   cat("Plotting GT-only sample:", sampleID, "\n")
#   
#   plot_df <- readRDS(
#     paste0(
#       "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/Karyogram_gt/",
#       sampleID, "_combined_gt.rds"
#     )
#   )
#   
#   plot_df[, method_cols_gt] <- lapply(plot_df[, method_cols_gt], function(x) {
#     pmax(pmin(x, cap_val), -cap_val)
#   })
#   
#   plot_bin_heatmap(
#     plot_df = plot_df,
#     method_cols = method_cols_gt,
#     method_display_names = method_display_names_gt,
#     sample_name = sampleID,
#     output_plot_heatmap = paste0(
#       "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/Karyogram_gt/",
#       sampleID, "_Karyogram_gt.png"
#     ),
#     cap_val = cap_val,
#     width = 15, height = 4, text.prop = 0.8
#   )
# }
# 
# 
# 
# 
# ########################
# ## state only
# ########################
# method_cols_gt <- c(
#   "aneufinder_original",
#   "aneufinder_ploidy",
#   "copykit_state"
# )
# 
# method_display_names_gt <- c(
#   aneufinder_original = "Schmid",
#   aneufinder_ploidy = "AneuFinder_ploidy",
#   copykit_state = "CopyKit_state"
# )
# 
# cap_val <- 5
# 
# for (sampleID in df$sample) {
#   cat("Plotting GT-only sample:", sampleID, "\n")
# 
#   plot_df <- readRDS(
#     paste0(
#       "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/Karyogram_gt/",
#       sampleID, "_combined_gt.rds"
#     )
#   )
# 
#   plot_df[, method_cols_gt] <- lapply(plot_df[, method_cols_gt], function(x) {
#     pmax(pmin(x, cap_val), -cap_val)
#   })
# 
#   plot_bin_heatmap(
#     plot_df = plot_df,
#     method_cols = method_cols_gt,
#     method_display_names = method_display_names_gt,
#     sample_name = sampleID,
#     output_plot_heatmap = paste0(
#       "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/Karyogram_gt/state_only/",
#       sampleID, "_Karyogram_gt.png"
#     ),
#     cap_val = cap_val,
#     width = 13, height = 3, text.prop = 0.8
#   )
# }
# 
