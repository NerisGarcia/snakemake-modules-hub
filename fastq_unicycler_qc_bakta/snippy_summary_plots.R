suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

snippy_report <- function(
  snippy_txt,
  snippy_vcf,
  outdir = "data/2_mapping/snippy_report",
  window_size = 5000,
  mask_window_size = 2000
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # ---------- 1) Read Snippy TXT summary ----------
  df <- fread(snippy_txt, sep = "\t", header = TRUE, data.table = FALSE)

  # Remove synthetic reference row if present
  df <- df %>% filter(ID != "Reference")

  # Strip trailing _snippy from sample IDs
  df$ID <- sub("_snippy$", "", df$ID)

  # Derived metrics
  df <- df %>%
    mutate(
      COVERAGE_PCT = 100 * ALIGNED / LENGTH,
      UNALIGNED_PCT = 100 * UNALIGNED / LENGTH,
      LOWCOV_PCT = 100 * LOWCOV / LENGTH
    )

  metrics <- c("UNALIGNED", "LOWCOV", "HET", "VARIANT", "COVERAGE_PCT")

  get_outliers <- function(x) {
    q1 <- as.numeric(quantile(x, 0.25, na.rm = TRUE))
    q3 <- as.numeric(quantile(x, 0.75, na.rm = TRUE))
    iqr <- q3 - q1
    c(low = q1 - 1.5 * iqr, high = q3 + 1.5 * iqr)
  }

  long_df <- df %>%
    select(ID, all_of(metrics)) %>%
    pivot_longer(cols = -ID, names_to = "metric", values_to = "value")

  bounds <- long_df %>%
    group_by(metric) %>%
    summarise(
      low = get_outliers(value)[["low"]],
      high = get_outliers(value)[["high"]],
      .groups = "drop"
    )

  outliers <- long_df %>%
    left_join(bounds, by = "metric") %>%
    filter(value < low | value > high)

  # ---------- QC flag table ----------
  # Direction: COVERAGE_PCT → flag only LOW outliers (more coverage = better)
  # All others (UNALIGNED, LOWCOV, HET, VARIANT) → flag only HIGH outliers
  flag_df <- long_df %>%
    left_join(bounds, by = "metric") %>%
    mutate(
      outlier_flag = case_when(
        metric == "COVERAGE_PCT" & value < low  ~ "LOW_COVERAGE",
        metric != "COVERAGE_PCT" & value > high ~ paste0("HIGH_", metric),
        TRUE                                    ~ NA_character_
      )
    ) %>%
    filter(!is.na(outlier_flag)) %>%
    select(ID, metric, outlier_flag)

  # Wide table: one row per sample, one column per metric + flag columns
  qc_wide <- df %>%
    select(ID, ALIGNED, UNALIGNED, LOWCOV, HET, VARIANT,
           COVERAGE_PCT, UNALIGNED_PCT, LOWCOV_PCT) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))

  # Combine flags per sample into a single string
  flag_summary <- flag_df %>%
    group_by(ID) %>%
    summarise(FLAGS = paste(outlier_flag, collapse = "; "), .groups = "drop")

  qc_report <- qc_wide %>%
    left_join(flag_summary, by = "ID") %>%
    mutate(
      FLAGS   = ifelse(is.na(FLAGS), "", FLAGS),
      PASS_QC = ifelse(FLAGS == "", "PASS", "FLAG")
    ) %>%
    arrange(PASS_QC, COVERAGE_PCT)

  write.csv(df, file.path(outdir, "snippy_summary_with_derived_metrics.csv"), row.names = FALSE)
  write.csv(outliers, file.path(outdir, "snippy_metric_outliers.csv"), row.names = FALSE)
  write.csv(qc_report, file.path(outdir, "snippy_qc_report.csv"), row.names = FALSE)

  # Metric labels for display
  metric_labels <- c(
    UNALIGNED   = "Unaligned (bp)",
    LOWCOV      = "Low coverage (bp)",
    HET         = "Heterozygous sites",
    VARIANT     = "Variant sites",
    COVERAGE_PCT = "Coverage (%)"
  )
  long_df$metric_label <- metric_labels[long_df$metric]
  outliers$metric_label <- metric_labels[outliers$metric]

  p_box <- ggplot(long_df, aes(x = "", y = value)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey30") +
    geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
    geom_point(data = outliers, color = "red", size = 1.8) +
    ggrepel::geom_text_repel(
      data = outliers,
      aes(label = ID),
      color = "red",
      size = 2.2,
      max.overlaps = 20
    ) +
    facet_wrap(~metric_label, scales = "free_y", nrow = 1) +
    theme_bw() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text   = element_text(size = 9)
    ) +
    labs(title = "Snippy per-sample metrics", x = NULL, y = "Value")

  # Attach metric_label to bounds for the vline faceting
  bounds$metric_label <- metric_labels[bounds$metric]
  p_hist <- ggplot(long_df, aes(x = value)) +
    geom_histogram(bins = 40, fill = "steelblue", color = "white", linewidth = 0.2) +
    geom_vline(
      data = bounds,
      aes(xintercept = high),
      linetype = "dashed", color = "red", linewidth = 0.5
    ) +
    facet_wrap(~metric_label, scales = "free", nrow = 1) +
    theme_bw() +
    theme(strip.text = element_blank()) +
    labs(x = "Value", y = "Count")

  p_combined <- p_box / p_hist + plot_layout(heights = c(2, 1))

  ggsave(file.path(outdir, "snippy_metrics_boxplot_outliers.png"), p_combined, width = 14, height = 7, dpi = 300)

  # ---------- 2) Read VCF and compute diversity accumulation ----------
  vcf_cmd <- if (grepl("\\.gz$", snippy_vcf)) {
    paste0("zcat ", shQuote(snippy_vcf), " | grep -v '^##'")
  } else {
    paste0("grep -v '^##' ", shQuote(snippy_vcf))
  }

  vcf <- fread(cmd = vcf_cmd, sep = "\t", header = TRUE, data.table = FALSE, quote = "", fill = TRUE)

  stopifnot(all(c("#CHROM", "POS", "REF", "ALT", "FORMAT") %in% colnames(vcf)))

  sample_cols <- colnames(vcf)[10:ncol(vcf)]
  if (length(sample_cols) == 0) stop("No sample genotype columns found in VCF.")

  # Strip _snippy suffix from VCF sample names
  sample_cols_clean <- sub("_snippy$", "", sample_cols)

  gt_mat <- as.matrix(vcf[, sample_cols, drop = FALSE])
  gt_mat <- apply(gt_mat, 2, function(x) sub(":.*$", "", x))
  if (is.null(dim(gt_mat))) {
    gt_mat <- matrix(gt_mat, ncol = 1)
    colnames(gt_mat) <- sample_cols[1]
  }

  is_var <- !(gt_mat %in% c("0/0", "0|0", "./.", ".|.", ".", "0"))
  if (is.null(dim(is_var))) {
    is_var <- matrix(is_var, ncol = ncol(gt_mat))
    colnames(is_var) <- colnames(gt_mat)
  }
  storage.mode(is_var) <- "integer"
  colnames(is_var) <- sample_cols_clean

  pos <- as.integer(vcf$POS)

  # ---- Diversity line plot (window_size, default 5000 bp) ----
  bins <- ((pos - 1L) %/% window_size) * window_size + 1L

  var_site <- rowSums(is_var) > 0
  site_df <- data.frame(bin = bins, var_site = var_site) %>%
    group_by(bin) %>%
    summarise(n_variant_sites = sum(var_site), .groups = "drop") %>%
    arrange(bin)

  write.csv(site_df, file.path(outdir, "vcf_variant_density_by_window.csv"), row.names = FALSE)

  p_density <- ggplot(site_df, aes(x = bin, y = n_variant_sites)) +
    geom_line(color = "#3182bd", linewidth = 0.7) +
    geom_point(color = "#3182bd", size = 1.2) +
    theme_bw() +
    labs(
      title = "Accumulated diversity along genome",
      subtitle = paste0("Window size: ", format(window_size, big.mark = ","), " bp"),
      x = "Genome position (bp)",
      y = "Number of variant sites"
    )

  ggsave(file.path(outdir, "vcf_variant_density_by_window.png"), p_density, width = 12, height = 4.5, dpi = 300)

  # ---- Samples-per-window plot (mask_window_size, default 2000 bp) ----
  # Counts how many SAMPLES have at least one variant in each window.
  # High-sample-count windows flag candidate regions to mask (phage, MGE, duplications).
  mbins <- ((pos - 1L) %/% mask_window_size) * mask_window_size + 1L

  n_samples_with_var <- vapply(
    split(seq_len(nrow(is_var)), mbins),
    function(idx) {
      sub_mat <- is_var[idx, , drop = FALSE]
      sum(colSums(sub_mat) > 0)
    },
    integer(1)
  )

  mask_df <- data.frame(
    bin           = as.integer(names(n_samples_with_var)),
    n_samples     = n_samples_with_var
  ) %>%
    filter(n_samples > 0) %>%
    arrange(bin)

  write.csv(mask_df, file.path(outdir, "vcf_samples_per_window_masking.csv"), row.names = FALSE)

  n_total_samples <- ncol(is_var)
  p_mask <- ggplot(mask_df, aes(x = bin, y = n_samples)) +
    geom_col(fill = "#e6550d", width = mask_window_size * 0.85) +
    geom_hline(
      yintercept = n_total_samples * 0.9,
      linetype = "dashed", color = "black", linewidth = 0.5
    ) +
    annotate(
      "text", x = max(mask_df$bin), y = n_total_samples * 0.9,
      label = "90% of samples", vjust = -0.4, hjust = 1, size = 3
    ) +
    theme_bw() +
    labs(
      title = "Samples with a variant per genome window",
      subtitle = paste0(
        "Window: ", format(mask_window_size, big.mark = ","),
        " bp | Windows with zero variants removed | n samples = ", n_total_samples
      ),
      x = "Genome position (bp)",
      y = "Number of samples with a variant"
    )

  ggsave(file.path(outdir, "vcf_samples_per_window_masking.png"), p_mask, width = 12, height = 4.5, dpi = 300)

  invisible(list(
    summary_table  = df,
    outliers       = outliers,
    qc_report      = qc_report,
    window_density = site_df,
    mask_windows   = mask_df,
    boxplot        = p_box,
    density_plot   = p_density,
    mask_plot      = p_mask
  ))
}
