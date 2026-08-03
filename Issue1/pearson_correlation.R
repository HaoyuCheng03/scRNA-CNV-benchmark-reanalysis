############################
## Pearson barplot
############################

library(dplyr)
library(ggplot2)
library(reshape2)

plot_pearson_matrix <- function(plot_df,
                                gt_cols,
                                method_cols,
                                gt_display_names = NULL,
                                method_display_names = NULL,
                                output_plot = "pearson_matrix.png",
                                output_cor = "pearson_matrix.txt") {
  stopifnot(all(c(gt_cols, method_cols) %in% colnames(plot_df)))

  if (is.null(gt_display_names)) {
    gt_display_names <- setNames(gt_cols, gt_cols)
  }
  if (is.null(method_display_names)) {
    method_display_names <- setNames(method_cols, method_cols)
  }

  cor_mat <- matrix(
    NA_real_,
    nrow = length(method_cols),
    ncol = length(gt_cols),
    dimnames = list(method_cols, gt_cols)
  )

  for (i in seq_along(method_cols)) {
    for (j in seq_along(gt_cols)) {
      x <- plot_df[[method_cols[i]]]
      y <- plot_df[[gt_cols[j]]]
      keep <- complete.cases(x, y)

      if (sum(keep) >= 3) {
        cor_mat[i, j] <- cor(x[keep], y[keep], method = "pearson")
      }
    }
  }

  ## save correlation matrix, ordered by aneufinder_original (high to low)
  cor_mat_out <- as.data.frame(cor_mat)
  cor_mat_out$Method <- rownames(cor_mat_out)
  cor_mat_out <- cor_mat_out[order(cor_mat_out$aneufinder_original, decreasing = TRUE), ]
  cor_mat_out <- cor_mat_out[, c("Method", gt_cols)]
  write.table(
    cor_mat_out,
    file = output_cor,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  ## use the same order in the plot
  ordered_methods <- cor_mat_out$Method

  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(cor_df) <- c("Method", "GroundTruth", "Pearson")

  cor_df$Method <- factor(
    cor_df$Method,
    levels = rev(ordered_methods),
    labels = rev(method_display_names[ordered_methods])
  )

  cor_df$GroundTruth <- factor(
    cor_df$GroundTruth,
    levels = gt_cols,
    labels = gt_display_names[gt_cols]
  )

  cor_df$label <- ifelse(is.na(cor_df$Pearson), "NA", sprintf("%.2f", cor_df$Pearson))

  g <- ggplot(cor_df, aes(x = GroundTruth, y = Method, fill = Pearson)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = 8) +   # keep correlation value size
    scale_fill_gradient2(
      name = "Pearson r",
      low = "darkblue",
      mid = "white",
      high = "darkred",
      midpoint = 0.5,
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    theme_bw(base_size = 22) +
    xlab("Ground truth") +
    ylab("Method") +
    theme(
      panel.grid = element_blank(),
      
      axis.text.x = element_text(size = 22, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 22),
      axis.title.x = element_text(size = 26),
      axis.title.y = element_text(size = 26),
      
      legend.title = element_text(size = 24),
      legend.text = element_text(size = 22),
      
      plot.title = element_text(size = 26, hjust = 0.5),
      strip.text = element_text(size = 24)
    )

  ggsave(output_plot, plot = g, width = 10, height = 9, dpi = 300)
  return(g)
}


setwd("/ix1/ctseng/hac377/cnv_calling/letter")
df <- read.table("samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## 2 ground truths (columns)
gt_cols <- c(
  "aneufinder_original",
  "aneufinder_ploidy"
)

## 6 methods (rows)
method_cols <- c(
  "infercnv",
  "copykat",
  "numbat",
  "scevan",
  "casper",
  "conicsmat"
)

gt_display_names <- c(
  aneufinder_original = "Aneufinder_Replicate",
  aneufinder_ploidy = "Aneufinder_Karyotyping"
)

method_display_names <- c(
  infercnv = "inferCNV",
  copykat = "CopyKAT",
  numbat = "Numbat",
  scevan = "SCEVAN",
  casper = "CaSpER",
  conicsmat = "CONICSmat"
)

for (sampleID in df$sample) {
  cat("Plotting sample:", sampleID, "\n")

  plot_df <- readRDS(
    paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs/ratio/",
      sampleID, "_combined_methods_ratio.rds"
    )
  )

  p <- plot_pearson_matrix(
    plot_df = plot_df,
    gt_cols = gt_cols,
    method_cols = method_cols,
    gt_display_names = gt_display_names,
    method_display_names = method_display_names,
    output_plot = paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs/ratio/Pearson/",
      sampleID, "_PearsonMatrix_ratio.png"
    ),
    output_cor = paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs/ratio/Pearson/",
      sampleID, "_PearsonMatrix_ratio.txt"
    )
  )

  print(p)
}







####################
## paired line plot
####################
library(dplyr)
library(ggplot2)
library(reshape2)


plot_pearson_slope <- function(plot_df,
                               gt_cols,
                               method_cols,
                               gt_display_names = NULL,
                               method_display_names = NULL) {
  stopifnot(all(c(gt_cols, method_cols) %in% colnames(plot_df)))
  
  if (is.null(gt_display_names)) {
    gt_display_names <- setNames(gt_cols, gt_cols)
  }
  
  if (is.null(method_display_names)) {
    method_display_names <- setNames(method_cols, method_cols)
  }
  
  ## Calculate Pearson correlation matrix
  cor_mat <- matrix(
    NA_real_,
    nrow = length(method_cols),
    ncol = length(gt_cols),
    dimnames = list(method_cols, gt_cols)
  )
  
  for (i in seq_along(method_cols)) {
    for (j in seq_along(gt_cols)) {
      x <- plot_df[[method_cols[i]]]
      y <- plot_df[[gt_cols[j]]]
      keep <- complete.cases(x, y)
      
      if (sum(keep) >= 3) {
        cor_mat[i, j] <- cor(x[keep], y[keep], method = "pearson")
      }
    }
  }
  
  ## Order methods by first ground truth
  cor_mat_out <- as.data.frame(cor_mat)
  cor_mat_out$Method <- rownames(cor_mat_out)
  cor_mat_out <- cor_mat_out[order(cor_mat_out[[gt_cols[1]]], decreasing = TRUE), ]
  
  ordered_methods <- cor_mat_out$Method
  
  ## Convert to long format
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(cor_df) <- c("Method", "GroundTruth", "Pearson")
  
  cor_df$Method <- factor(
    cor_df$Method,
    levels = ordered_methods,
    labels = method_display_names[ordered_methods]
  )
  
  cor_df$GroundTruth <- factor(
    cor_df$GroundTruth,
    levels = gt_cols,
    labels = gt_display_names[gt_cols]
  )
  
  cor_df$label <- ifelse(
    is.na(cor_df$Pearson),
    "NA",
    sprintf("%.2f", cor_df$Pearson)
  )
  
  ## Left/right label alignment by group
  cor_df$hjust <- ifelse(
    cor_df$GroundTruth == levels(cor_df$GroundTruth)[1],
    1.2,
    -0.2
  )
  
  ## Adjust label positions within each group to reduce overlap
  cor_df <- cor_df %>%
    group_by(GroundTruth) %>%
    arrange(Pearson) %>%
    mutate(
      label_y = Pearson,
      min_gap = 0.05
    )
  
  for (grp in unique(cor_df$GroundTruth)) {
    idx <- which(cor_df$GroundTruth == grp)
    for (k in 2:length(idx)) {
      if ((cor_df$label_y[idx[k]] - cor_df$label_y[idx[k - 1]]) < cor_df$min_gap[idx[k]]) {
        cor_df$label_y[idx[k]] <- cor_df$label_y[idx[k - 1]] + cor_df$min_gap[idx[k]]
      }
    }
  }
  
  cor_df$label_y <- pmin(cor_df$label_y, 1)
  
  ## Plot slope chart
  g <- ggplot(
    cor_df,
    aes(
      x = GroundTruth,
      y = Pearson,
      group = Method,
      color = Method
    )
  ) +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +
    geom_text(
      aes(
        y = label_y,
        label = label,
        hjust = hjust
      ),
      vjust = 0.5,
      size = 7,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25)
    ) +
    theme_bw(base_size = 24) +
    xlab("Ground truth") +
    ylab("Pearson r") +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 22),
      axis.text.y = element_text(size = 24),
      axis.title.x = element_text(size = 26),
      axis.title.y = element_text(size = 26),
      legend.title = element_blank(),
      legend.text = element_text(size = 24),
      legend.position = "right"
    )
  
  return(g)
}


setwd("/ix1/ctseng/hac377/cnv_calling/letter")

df <- read.table(
  "samples.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

## Ground truths
gt_cols <- c(
  "aneufinder_original",
  "aneufinder_ploidy"
)

## Methods
method_cols <- c(
  "infercnv",
  "copykat",
  "numbat",
  "scevan",
  "casper",
  "conicsmat"
)

## Display names
gt_display_names <- c(
  aneufinder_original = "Aneufinder_Schmid_R",
  aneufinder_ploidy = "Aneufinder_Karyotyping"
)

method_display_names <- c(
  infercnv = "inferCNV",
  copykat = "CopyKAT",
  numbat = "Numbat",
  scevan = "SCEVAN",
  casper = "CaSpER",
  conicsmat = "CONICSmat"
)

for (sampleID in df$sample) {
  
  cat("Plotting sample:", sampleID, "\n")
  
  plot_df <- readRDS(
    paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/ratio/",
      sampleID,
      "_combined_methods_ratio.rds"
    )
  )
  
  p <- plot_pearson_slope(
    plot_df = plot_df,
    gt_cols = gt_cols,
    method_cols = method_cols,
    gt_display_names = gt_display_names,
    method_display_names = method_display_names
  )
  
  print(p)
  
  
  ggsave(
    filename = paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/ratio/Pearson/Pearson_Slope/",
      sampleID,
      "_PearsonSlope_ratio.png"
    ),
    plot = p,
    width = 13,
    height = 8,
    dpi = 500
  )
}




############################
## Wilcoxon signed-rank test
############################

library(dplyr)

group_col <- "group"

all_delta_results <- list()

for (sampleID in df$sample) {
  
  cat("Processing:", sampleID, "\n")
  
  plot_df <- readRDS(
    paste0(
      "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor1/ratio/",
      sampleID,
      "_combined_methods_ratio.rds"
    )
  )
  
  sample_results <- data.frame(
    sample = character(),
    method = character(),
    pearson_original = numeric(),
    pearson_ploidy = numeric(),
    delta = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (method in method_cols) {
    
    x <- plot_df[[method]]
    
    y_original <- plot_df[["aneufinder_original"]]
    y_ploidy <- plot_df[["aneufinder_ploidy"]]
    
    keep1 <- complete.cases(x, y_original)
    keep2 <- complete.cases(x, y_ploidy)
    
    r_original <- cor(
      x[keep1],
      y_original[keep1],
      method = "pearson"
    )
    
    r_ploidy <- cor(
      x[keep2],
      y_ploidy[keep2],
      method = "pearson"
    )
    
    sample_results <- rbind(
      sample_results,
      data.frame(
        sample = sampleID,
        method = method,
        pearson_original = r_original,
        pearson_ploidy = r_ploidy,
        delta = r_ploidy - r_original
      )
    )
  }
  
  all_delta_results[[sampleID]] <- sample_results
}

all_delta_results <- bind_rows(all_delta_results)

## add group info
all_delta_results <- all_delta_results %>%
  left_join(df[, c("sample", group_col)], by = "sample")


####################################################
## Wilcoxon signed-rank test within each group
####################################################

group_names <- unique(all_delta_results[[group_col]])

for (grp in group_names) {
  
  cat("\n============================\n")
  cat("Group:", grp, "\n")
  
  delta_values <- all_delta_results %>%
    filter(.data[[group_col]] == grp) %>%
    pull(delta)
  
  print(delta_values)
  
  test_res <- wilcox.test(
    delta_values,
    mu = 0,
    alternative = "greater",
    exact = FALSE
  )
  
  print(test_res)
}
