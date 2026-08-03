library(pROC)

calc_partial_auc <- function(truth, method_vec, sampleID, method_name, group_label) {
  
  if (length(truth) != length(method_vec)) {
    stop("truth and method_vec must have same length")
  }
  
  keep <- !is.na(truth) & !is.na(method_vec)
  truth <- truth[keep]
  method_vec <- method_vec[keep]
  
  valid_states <- c("base", "loss", "gain")
  if (!all(truth %in% valid_states)) {
    stop("truth must only contain: 'base', 'loss', 'gain'")
  }
  
  ########################################
  ## Hardcoded plot directory
  ########################################
  plot_dir <- paste0(
    "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor2_3/roc/",
    group_label,
    "/",
    sampleID
  )
  
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  ########################################
  ## Gain ROC
  ########################################
  gain_truth <- truth == "gain"
  
  if (sum(gain_truth) == 0 || sum(!gain_truth) == 0) {
    auc_gain <- NA_real_
    pauc_gain <- NA_real_
    s_max_gain <- NA_real_
    roc_gain <- NULL
    
  } else {
    roc_gain <- roc(
      response = gain_truth,
      predictor = method_vec,
      quiet = TRUE,
      direction = "<"
    )
    
    auc_gain <- as.numeric(auc(roc_gain))
    
    s_max_gain <- sum(method_vec > 0 & gain_truth) / sum(gain_truth)
    
    if (s_max_gain <= 0) {
      pauc_gain <- NA_real_
    } else {
      pauc_gain <- as.numeric(
        auc(
          roc_gain,
          partial.auc = c(0, s_max_gain),
          partial.auc.focus = "sensitivity",
          partial.auc.correct = FALSE
        )
      )
    }
  }
  
  ########################################
  ## Loss ROC
  ########################################
  loss_truth <- truth == "loss"
  
  if (sum(loss_truth) == 0 || sum(!loss_truth) == 0) {
    auc_loss <- NA_real_
    pauc_loss <- NA_real_
    s_max_loss <- NA_real_
    roc_loss <- NULL
    
  } else {
    loss_scores <- method_vec * -1
    
    roc_loss <- roc(
      response = loss_truth,
      predictor = loss_scores,
      quiet = TRUE,
      direction = "<"
    )
    
    auc_loss <- as.numeric(auc(roc_loss))
    
    s_max_loss <- sum(loss_scores > 0 & loss_truth) / sum(loss_truth)
    
    if (s_max_loss <= 0) {
      pauc_loss <- NA_real_
    } else {
      pauc_loss <- as.numeric(
        auc(
          roc_loss,
          partial.auc = c(0, s_max_loss),
          partial.auc.focus = "sensitivity",
          partial.auc.correct = FALSE
        )
      )
    }
  }
  
  ########################################
  ## Plot ROC curves
  ########################################
  png(
    filename = paste0(
      plot_dir,
      "/",
      method_name,
      "_ROC.png"
    ),
    width = 2600,
    height = 1200,
    res = 200
  )
  
  par(mfrow = c(1, 2))
  
  ########################################
  ## Gain plot
  ########################################
  if (!is.null(roc_gain)) {
    
    plot(
      roc_gain,
      main = paste0(method_name, " Gain ROC (", group_label, ")"),
      legacy.axes = TRUE
    )
    
    abline(h = s_max_gain, lty = 2)
    
    text(
      x = 0.6,
      y = s_max_gain,
      labels = paste0("s_max = ", round(s_max_gain, 3)),
      pos = 3
    )
    
  } else {
    plot.new()
    title("Gain ROC unavailable")
  }
  
  ########################################
  ## Loss plot
  ########################################
  if (!is.null(roc_loss)) {
    
    plot(
      roc_loss,
      main = paste0(method_name, " Loss ROC (", group_label, ")"),
      legacy.axes = TRUE
    )
    
    abline(h = s_max_loss, lty = 2)
    
    text(
      x = 0.6,
      y = s_max_loss,
      labels = paste0("s_max = ", round(s_max_loss, 3)),
      pos = 3
    )
    
  } else {
    plot.new()
    title("Loss ROC unavailable")
  }
  
  dev.off()
  
  ########################################
  ## Return
  ########################################
  list(
    gain = list(
      auc = auc_gain,
      partial_auc = pauc_gain,
      s_max = s_max_gain
    ),
    loss = list(
      auc = auc_loss,
      partial_auc = pauc_loss,
      s_max = s_max_loss
    )
  )
}

calc_f1 <- function(pred, truth, class_name) {
  tp <- sum(pred == class_name & truth == class_name)
  fp <- sum(pred == class_name & truth != class_name)
  fn <- sum(pred != class_name & truth == class_name)
  
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  
  ifelse(
    precision + recall == 0,
    0,
    2 * precision * recall / (precision + recall)
  )
}


calc_F1_ratio <- function(truth, method_vec, max_thresholds = 200) {
  
  stopifnot(length(truth) == length(method_vec))
  
  keep <- !is.na(truth) & !is.na(method_vec)
  truth <- truth[keep]
  method_vec <- method_vec[keep]
  
  valid_states <- c("base", "loss", "gain")
  if (!all(truth %in% valid_states)) {
    stop("truth must only contain: base, loss, gain")
  }
  
  scores <- sort(unique(method_vec))
  
  loss_thresholds <- scores[scores < 0]
  gain_thresholds <- scores[scores > 0]
  
  if (length(loss_thresholds) > max_thresholds) {
    loss_thresholds <- loss_thresholds[
      round(seq(1, length(loss_thresholds), length.out = max_thresholds))
    ]
  } else if (length(loss_thresholds) == 0) {
    loss_thresholds <- -1
  }
  
  if (length(gain_thresholds) > max_thresholds) {
    gain_thresholds <- gain_thresholds[
      round(seq(1, length(gain_thresholds), length.out = max_thresholds))
    ]
  } else if (length(gain_thresholds) == 0) {
    gain_thresholds <- 1
  }
  
  ## also test baseline threshold directly
  loss_thresholds <- sort(unique(c(loss_thresholds, 0)))
  gain_thresholds <- sort(unique(c(0, gain_thresholds)))
  
  grid <- expand.grid(
    loss_threshold = loss_thresholds,
    gain_threshold = gain_thresholds
  )

  
  calc_multi_f1 <- function(loss_t, gain_t) {
    pred <- rep("base", length(method_vec))
    pred[method_vec < loss_t] <- "loss"
    pred[method_vec > gain_t] <- "gain"
    
    f1_loss <- calc_f1(pred, truth, "loss")
    f1_base <- calc_f1(pred, truth, "base")
    f1_gain <- calc_f1(pred, truth, "gain")
    
    mean(c(f1_loss, f1_base, f1_gain))
  }
  
  grid$f1 <- mapply(
    calc_multi_f1,
    grid$loss_threshold,
    grid$gain_threshold
  )
  
  best <- grid[which.max(grid$f1), ]
  
  pred <- rep("base", length(method_vec))
  pred[method_vec < best$loss_threshold] <- "loss"
  pred[method_vec > best$gain_threshold] <- "gain"
  
  
  list(
    loss_threshold = best$loss_threshold,
    gain_threshold = best$gain_threshold,
    maximal_F1 = best$f1,
    prediction = pred
    )
}


calc_F1_state <- function(truth, pred) {
  
  if (length(truth) != length(pred)) {
    stop("truth and pred must have same length")
  }
  
  # remove NA
  keep <- !is.na(truth) & !is.na(pred)
  truth <- truth[keep]
  pred <- pred[keep]
  
  # convert
  truth_num <- rep(NA_integer_, length(truth))
  truth_num[truth == "loss"] <- -1L
  truth_num[truth == "base"] <- 0L
  truth_num[truth == "gain"] <- 1L
  
  if (any(is.na(truth_num))) {
    stop("truth must only contain: 'loss', 'base', 'gain'")
  }
  
  # check pred
  if (!all(pred %in% c(-1, 0, 1))) {
    stop("pred must only contain: -1, 0, 1")
  }
  
  
  # get per-class F1
  f1_loss <- calc_f1(pred, truth_num, -1)
  f1_base <- calc_f1(pred, truth_num, 0)
  f1_gain <- calc_f1(pred, truth_num, 1)
  
  # get multi-class F1
  multi_f1 <- mean(c(f1_loss, f1_base, f1_gain))

  return(multi_f1)
}


calc_pearson <- function(gt_vec, method_vec) {
  
  if (length(gt_vec) != length(method_vec)) {
    stop("gt_vec and method_vec must have same length")
  }
  
  keep <- !is.na(gt_vec) & !is.na(method_vec)
  
  cor(
    gt_vec[keep],
    method_vec[keep],
    method = "pearson"
  )
}



method_ratio_cols <- c("infercnv_ratio", "copykat_ratio", "numbat_ratio", 
                       "scevan_ratio", "casper_ratio", "conicsmat_ratio")
method_original_state_cols <- c("infercnv_original_state", "copykat_original_state", "numbat_original_state", 
                                "scevan_original_state", "casper_original_state", "conicsmat_original_state")
method_true_state_cols <- c("infercnv_true_state", "copykat_true_state", "numbat_true_state", 
                            "scevan_true_state", "casper_true_state", "conicsmat_true_state")

method_names <- c("inferCNV", "CopyKAT", "Numbat", "SCEVAN", "CaSpER", "CONICSmat")
names(method_ratio_cols) <- method_names
names(method_original_state_cols) <- method_names
names(method_true_state_cols) <- method_names

sample_df <- read.table("/ix1/ctseng/hac377/cnv_calling/letter/samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# auc_debug <- sapply(method_ratio_cols, function(x) {
#   c(
#     auc_gain_smax = mean(df[[x]][df$aneufinder_original_state == "gain"] > 0, na.rm = TRUE),
#     auc_loss_smax = mean(df[[x]][df$aneufinder_original_state == "loss"] < 0, na.rm = TRUE)
#   )
# })
# 
# auc_debug



outDir <- "/ix1/ctseng/hac377/cnv_calling/letter/figs_factor2_3/results/"
if (!dir.exists(outDir)) {dir.create(outDir)}

for (sampleID in sample_df$sample) {
  
  cat("Processing sample:", sampleID, "\n")
  
  infile <- paste0("/ix1/ctseng/hac377/cnv_calling/letter/figs_factor2_3/",
                   sampleID, "_combined_methods.rds")
  
  sample_outDir <- paste0(outDir, sampleID)
  if (!dir.exists(sample_outDir)) {dir.create(sample_outDir)}
  
  df <- readRDS(infile)
  
  ################
  ## 0, 0
  ## AUC/maximal F1: gt_original_state -- method_ratio
  ## Pearson: gt_original_ratio -- method_ratio
  ################
  outfile_00 <- paste0(sample_outDir, "/00_res.rds")
  
  results_df <- data.frame(
    pearson = rep(NA_real_, length(method_names)),
    auc_gain = rep(NA_real_, length(method_names)),
    auc_loss = rep(NA_real_, length(method_names)),
    pauc_gain = rep(NA_real_, length(method_names)),
    pauc_loss = rep(NA_real_, length(method_names)),
    F1 = rep(NA_real_, length(method_names)),
    row.names = method_names,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(method_ratio_cols)) {
    
    method_col <- method_ratio_cols[i]
    method_name <- names(method_ratio_cols)[i]
    
    # AUC
    # auc_res <- calc_partial_auc(
    #   df$aneufinder_original_state,
    #   df[[method_col]]
    # )
    
    auc_res <- calc_partial_auc(
      truth = df$aneufinder_original_state,
      method_vec = df[[method_col]],
      sampleID = sampleID,
      method_name = method_name,
      group_label = "00"
    )
    
    gain_auc <- auc_res$gain$auc
    loss_auc <- auc_res$loss$auc
    gain_partial_auc <- auc_res$gain$partial_auc
    loss_partial_auc <- auc_res$loss$partial_auc
    
    # F1
    F1_res <- calc_F1_ratio(
      df$aneufinder_original_state,
      df[[method_col]]
    )
    
    maximal_F1 <- F1_res$maximal_F1
    
    # Pearson corr
    pearson_corr <- calc_pearson(
      df$aneufinder_original_ratio,
      df[[method_col]]
    )
    
    # Fill dataframe
    results_df[method_name, ] <- c(
      pearson_corr,
      gain_auc,
      loss_auc,
      gain_partial_auc,
      loss_partial_auc,
      maximal_F1
    )
  }
  
  saveRDS(results_df, outfile_00)
  
  
  
  
  ################
  ## 0, 1
  ## AUC/multi F1: gt_original_state -- method_original_state
  ## Pearson: gt_original_ratio -- method_ratio
  ################
  outfile_01 <- paste0(sample_outDir, "/01_res.rds")
  
  results_df <- data.frame(
    pearson = rep(NA_real_, length(method_names)),
    auc_gain = rep(NA_real_, length(method_names)),
    auc_loss = rep(NA_real_, length(method_names)),
    pauc_gain = rep(NA_real_, length(method_names)),
    pauc_loss = rep(NA_real_, length(method_names)),
    F1 = rep(NA_real_, length(method_names)),
    row.names = method_names,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(method_ratio_cols)) {
    
    method_ratio_col <- method_ratio_cols[i]
    method_original_state_col <- method_original_state_cols[i]
    method_name <- names(method_ratio_cols)[i]
    
    # AUC
    # auc_res <- calc_partial_auc(
    #   df$aneufinder_original_state,
    #   df[[method_original_state_col]]
    # )
    
    auc_res <- calc_partial_auc(
      truth = df$aneufinder_original_state,
      method_vec = df[[method_original_state_col]],
      sampleID = sampleID,
      method_name = method_name,
      group_label = "01"
    )
    
    
    gain_auc <- auc_res$gain$auc
    loss_auc <- auc_res$loss$auc
    gain_partial_auc <- auc_res$gain$partial_auc
    loss_partial_auc <- auc_res$loss$partial_auc
    
    # F1
    maximal_F1 <- calc_F1_state(
      df$aneufinder_original_state,
      df[[method_original_state_col]]
    )
    
    # Pearson
    pearson_corr <- calc_pearson(
      df$aneufinder_original_ratio,
      df[[method_ratio_col]]
    )
    
    results_df[method_name, ] <- c(
      pearson_corr,
      gain_auc,
      loss_auc,
      gain_partial_auc,
      loss_partial_auc,
      maximal_F1
    )
  }
  
  saveRDS(results_df, outfile_01)
  
  
  
  ################
  ## 1, 0
  ## AUC/maximal F1: gt_true_state -- method_ratio
  ## Pearson: gt_true_ratio -- method_ratio
  ################
  outfile_10 <- paste0(sample_outDir, "/10_res.rds")
  
  results_df <- data.frame(
    pearson = rep(NA_real_, length(method_names)),
    auc_gain = rep(NA_real_, length(method_names)),
    auc_loss = rep(NA_real_, length(method_names)),
    pauc_gain = rep(NA_real_, length(method_names)),
    pauc_loss = rep(NA_real_, length(method_names)),
    F1 = rep(NA_real_, length(method_names)),
    row.names = method_names,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(method_ratio_cols)) {
    
    method_ratio_col <- method_ratio_cols[i]
    method_name <- names(method_ratio_cols)[i]
    
    # AUC
    # auc_res <- calc_partial_auc(
    #   df$aneufinder_ploidy_state,
    #   df[[method_ratio_col]]
    # )
    
    auc_res <- calc_partial_auc(
      truth = df$aneufinder_ploidy_state,
      method_vec = df[[method_ratio_col]],
      sampleID = sampleID,
      method_name = method_name,
      group_label = "10"
    )
    
    gain_auc <- auc_res$gain$auc
    loss_auc <- auc_res$loss$auc
    gain_partial_auc <- auc_res$gain$partial_auc
    loss_partial_auc <- auc_res$loss$partial_auc
    
    # F1
    maximal_F1 <- calc_F1_ratio(
      df$aneufinder_ploidy_state,
      df[[method_ratio_col]]
    )$maximal_F1
    
    # Pearson
    pearson_corr <- calc_pearson(
      df$aneufinder_true_ratio,
      df[[method_ratio_col]]
    )
    
    results_df[method_name, ] <- c(
      pearson_corr,
      gain_auc,
      loss_auc,
      gain_partial_auc,
      loss_partial_auc,
      maximal_F1
    )
  }
  
  saveRDS(results_df, outfile_10)
  
  
  
  ################
  ## 1, 1
  ## AUC/multi F1: gt_true_state -- method_true_state
  ## Pearson: gt_true_ratio -- method_ratio
  ################
  
  outfile_11 <- paste0(sample_outDir, "/11_res.rds")
  
  results_df <- data.frame(
    pearson = rep(NA_real_, length(method_names)),
    auc_gain = rep(NA_real_, length(method_names)),
    auc_loss = rep(NA_real_, length(method_names)),
    pauc_gain = rep(NA_real_, length(method_names)),
    pauc_loss = rep(NA_real_, length(method_names)),
    F1 = rep(NA_real_, length(method_names)),
    row.names = method_names,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(method_ratio_cols)) {
    
    method_ratio_col <- method_ratio_cols[i]
    method_true_state_col <- method_true_state_cols[i]
    method_name <- names(method_ratio_cols)[i]
    
    # AUC
    # auc_res <- calc_partial_auc(
    #   df$aneufinder_ploidy_state,
    #   df[[method_true_state_col]]
    # )
    # 
    auc_res <- calc_partial_auc(
      truth = df$aneufinder_ploidy_state,
      method_vec = df[[method_true_state_col]],
      sampleID = sampleID,
      method_name = method_name,
      group_label = "11"
    )
    
    gain_auc <- auc_res$gain$auc
    loss_auc <- auc_res$loss$auc
    gain_partial_auc <- auc_res$gain$partial_auc
    loss_partial_auc <- auc_res$loss$partial_auc
    
    # F1
    maximal_F1 <- calc_F1_state(
      df$aneufinder_ploidy_state,
      df[[method_true_state_col]]
    )
    
    # Pearson
    pearson_corr <- calc_pearson(
      df$aneufinder_true_ratio,
      df[[method_ratio_col]]
    )
    
    results_df[method_name, ] <- c(
      pearson_corr,
      gain_auc,
      loss_auc,
      gain_partial_auc,
      loss_partial_auc,
      maximal_F1
    )
  }
  
  saveRDS(results_df, outfile_11)
  
}



