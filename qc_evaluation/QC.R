# =============================================================================
# Genomic Description of the Isolate Sequences and Genomes
# REFACTORED: ALL ANALYSIS FIRST, THEN ALL PLOTS AT END
# =============================================================================

# Libraries -------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(readODS)
library(patchwork)
library(phangorn)

# personal libraries
library(Itools)
library(customR)

write_ods_file <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_ods(x, path)
}


# ── Paths ──────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

# ——————————————————————————————————————————————————————————————————————---------
# SECTION 5: ASSEMBLY ANALYSIS
# ——————————————————————————————————————————————————————————————————————---------

# Read sample list
dataset_name <- "Lcrispatus.v1"
#dataset_name <-  args[1]

sample_list =  "data/1_datasets/Lcrispatus.v1.accesions.txt.csv"
#sample_list <- args[2]

qc_input_dir <-  "data/0_input_data/3_QC"
#qc_input_dir <- args[3]

dataset_output_dir <- paste0("data/1_datasets", "/", dataset_name)
#dataset_output_dir <- args[4]

figures_output_dir <-  paste0("figures/1_dataset_description", "/", dataset_name)
#figures_output_dir <- args[5]


#figure_log <- args[6]

sample_ids <- readLines(sample_list, warn = FALSE) 

quast_report_files <- 
list.files(paste0(qc_input_dir, "/QUAST"), 
pattern = "quast.transposed_report.tsv", 
full.names = TRUE, 
  recursive = TRUE)


if (length(sample_ids) > 0) {
  sample_pattern <- str_c(str_replace_all(sample_ids, "([.|()\\^{}+$*?]|\\[|\\]|\\\\)", "\\\\\\1"), collapse = "|")
  quast_report_files <- quast_report_files[str_detect(basename(quast_report_files), sample_pattern)]
}

quast_stats <- 
tibble(file = quast_report_files) %>%
  mutate(data = map(file, ~ read.csv(.x, sep = "\t", header = TRUE))) %>%
  unnest(data) %>%
  select(-file) %>%
  rename_with(~ str_squish(str_replace_all(str_remove(.x, "^X"), "\\.+", " "))) %>%
  mutate(
    across(
      starts_with("predicted genes"),
      list(
        complete = ~ as.numeric(str_extract(.x, "^\\s*\\d+")),
        partial  = ~ as.numeric(str_match(.x, "\\+\\s*(\\d+)\\s*part")[, 2])
      ),
      .names = "{.col} {.fn}"
    )
  ) %>%
  select(!starts_with("predicted genes") | ends_with("complete") | ends_with("partial")) %>%
  mutate(across(-Assembly, ~ (as.numeric(.x))))

# Busco results
busco_summary_files <-
list.files(paste0(qc_input_dir, "/BUSCO"),
pattern = "_busco_full_table.tsv",
full.names = TRUE,
  recursive = TRUE)


if (length(sample_ids) > 0) {
  busco_pattern <- str_c(str_replace_all(sample_ids, "([.|()\\^{}+$*?]|\\[|\\]|\\\\)", "\\\\\\1"), collapse = "|")
  busco_summary_files <- busco_summary_files[str_detect(basename(busco_summary_files), busco_pattern)]
}

busco_stats <-
tibble(file = busco_summary_files) %>%
  mutate(data = map(file, ~ read.csv(.x, sep = "\t", header = TRUE, skip=2))) %>%
  unnest(data) %>%
  mutate(Assembly = str_remove(basename(file), "_busco_full_table.tsv")) %>%
  select(-file)   %>% 
  select(Assembly, everything())   %>% 
  group_by(Assembly) %>%  
  count(Status)

busco_stats_wide <-
busco_stats %>%
  pivot_wider(names_from = Status, values_from = n, values_fill = 0, names_prefix = "BUSCO_") 

joined_stats <-
left_join(quast_stats, busco_stats_wide, by = "Assembly")

# Metrics used for assembly QC filtering and plots
assembly_metrics <- c(
  "contigs",
  "Total length",
  "GC",
  "predicted genes unique complete",
  "predicted genes 0 bp complete",
  "predicted genes 0 bp partial", 
  "BUSCO_Complete", 
  "BUSCO_Fragmented", 
  "BUSCO_Missing"
)

metric_filtering_stats <- 
joined_stats %>%
  select(Assembly, all_of(assembly_metrics)) %>%
  pivot_longer(cols = -Assembly, names_to = "metric", values_to = "value") %>%
  filter(metric %in% assembly_metrics) %>%
  filter(!is.na(value)) %>%
  group_by(metric) %>%
  mutate(
    median_v = median(value, na.rm = TRUE),
    mad_v = mad(value, na.rm = TRUE),
    lower_threshold = median_v - (3.5 * mad_v),
    upper_threshold = median_v + (3.5 * mad_v),
  ) %>%
  mutate(lower_threshold = case_when(
    str_starts(metric, regex("contigs", ignore_case = TRUE)) ~ 0, 
    metric == "BUSCO_Fragmented" ~ 0,
    metric == "BUSCO_Missing" ~ 0,
    metric == "BUSCO_Duplicated" ~ 0,
    str_starts(metric, regex("partial", ignore_case = TRUE)) ~ 0, 
    metric == "BUSCO_Complete" ~ median_v - (10 * mad_v),
    TRUE ~ lower_threshold
  ),
    upper_threshold = case_when(
      metric %in% c("BUSCO_Fragmented", "BUSCO_Missing") ~ 15,
      metric == "Largest contig" ~ Inf,      
      TRUE ~ upper_threshold
    )) 

metric_filtering_stats_out <- 
metric_filtering_stats  %>% 
select(metric, median_v, mad_v, lower_threshold, upper_threshold)  %>% 
distinct()  

write_ods_file(metric_filtering_stats_out, paste0(dataset_output_dir, "/", dataset_name, "_assembly_qc_filtering_thresholds.ods"))


joined_flags <-
metric_filtering_stats  %>% 
  mutate(is_outlier = case_when(
    metric == "BUSCO_Complete" ~ value < lower_threshold,
    metric %in% c("BUSCO_Fragmented", "BUSCO_Missing") ~ value > upper_threshold,
    str_starts(metric, "BUSCO_") ~ FALSE,
    TRUE ~ value < lower_threshold | value > upper_threshold
  )) %>%
  ungroup() %>%
  select(Assembly, metric, value, median_v, mad_v, lower_threshold, upper_threshold, is_outlier) %>%
  filter(is_outlier) %>%
  group_by(Assembly) %>%
  mutate(
    assembly_flags = str_c(sort(unique(metric)), collapse = "; ")
  ) %>%
  ungroup()  %>% 
  select(Assembly, assembly_flags)


joined_stats_flags <- left_join(joined_stats, joined_flags, by = "Assembly")  %>%  ungroup()  %>%  distinct() %>%
  mutate(assembly_flags = replace_na(assembly_flags, "PASS"))

write_ods_file(joined_stats_flags, paste0(dataset_output_dir, "/", dataset_name, "_assembly_qc_summary.ods"))



# PLOTS


# ———— Assembly Plots ——————————————————————————————————————————————
# Prepare assembly data (ALL FILTERS from qc.metadata.final, exclude reference samples)

assembly_plot_data <-
  metric_filtering_stats  %>% 
  mutate(is_outlier = case_when(
    metric == "BUSCO_Complete" ~ value < lower_threshold,
    metric %in% c("BUSCO_Fragmented", "BUSCO_Missing") ~ value > upper_threshold,
    str_starts(metric, "BUSCO_") ~ FALSE,
    TRUE ~ value < lower_threshold | value > upper_threshold
  )) %>%
  ungroup() %>%
    filter(metric %in% assembly_metrics)  %>% 
      filter(is_outlier, metric != "BUSCO_Complete") %>%
     group_by(Assembly) %>%
  mutate(
    assembly_flags = str_c(sort(unique(metric)), collapse = "; ")
  ) %>%
  ungroup()  %>% 
  select(Assembly, assembly_flags)  %>% 
  right_join(., joined_stats, by = "Assembly") %>%
  select(Assembly, all_of(assembly_metrics), assembly_flags) %>%
  mutate(has_flag = ifelse(is.na(assembly_flags), "PASS", "OUTLIER")) %>%
  pivot_longer(
    cols = all_of(assembly_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  # left_join(assembly_outlier_samples, by = c("CultureID"))  %>%
  mutate(is_outlier = ifelse(is.na(assembly_flags), FALSE, TRUE))
  



flag_levels <- sort(unique(assembly_plot_data$assembly_flags))
non_pass_flags <- setdiff(flag_levels, "PASS")
non_pass_colors <- if (length(non_pass_flags) > 0) {
  setNames(grDevices::hcl.colors(length(non_pass_flags), palette = "Dark 3"), non_pass_flags)
} else {
  character(0)
}

flag_color_key <- c("PASS" = "lightgrey", non_pass_colors)


quast_rank_data <-
  assembly_plot_data %>%
  group_by(metric) %>%
  arrange(desc(value), .by_group = TRUE) %>%
  mutate(sample_order = row_number()) %>%
  ungroup()


quast_box_data <- quast_rank_data %>%
  group_by(metric) %>%
  mutate(box_x = max(sample_order, na.rm = TRUE) * 0.5) %>%
  ungroup()

# Create assembly plots
# Ranked QUAST metrics (ordered high to low, colored and shaped by FILTER - no faceting)





# Threshold lines for each metric facet
threshold_lines <- metric_filtering_stats_out %>%
  filter(metric %in% assembly_metrics) %>%
  pivot_longer(
    cols = c(lower_threshold, upper_threshold),
    names_to = "threshold_type",
    values_to = "threshold"
  ) %>%
  mutate(
    threshold_type = recode(
      threshold_type,
      lower_threshold = "Lower threshold",
      upper_threshold = "Upper threshold"
    )
  )

p_quast_ranked <-
  ggplot(quast_rank_data, aes(x = sample_order, y = value)) +
  geom_violin(fill = "#A4DDEF", alpha = 0.5) +
  geom_hline(
    data = threshold_lines,
    aes(yintercept = threshold),
    inherit.aes = FALSE,
    colour = "#660b0b",
    linewidth = 0.2
  ) +
  geom_point(
    data = dplyr::filter(quast_rank_data, has_flag == "PASS"),
    aes(colour = has_flag, size = has_flag, shape = has_flag)
  ) +
  geom_point(
    data = dplyr::filter(quast_rank_data, has_flag != "PASS"),
    aes(y = value, fill = has_flag, colour = assembly_flags, shape = has_flag)
  ) +
  scale_size_manual(values = c("PASS" = 2, "OUTLIER" = 3), guide = "none") +
  scale_shape_manual(values = c("PASS" = 1, "OUTLIER" = 21)) +
  scale_color_manual(values = flag_color_key) +
  scale_fill_manual(values = flag_color_key, guide = "none") +
  facet_wrap(~factor(metric, levels = assembly_metrics), scales = "free_y", ncol = 3) +
  labs(
    title = "Assembly QC per Sample",
    x = "Samples ordered high to low",
    y = "Metric value",
    caption = paste0("Dataset: ", dataset_name)  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )


print(p_quast_ranked)
 

saveplot(p_quast_ranked, file.path(figures_output_dir, paste0(dataset_name, "_qc_assembly.violin_plot")), 
  layout = "custom", height_cm = 15, width_cm = 12, mode = "editing", tag = "filtered")



# LOG ————————————————————————————————————————————————————————————————————
# write log
filelog <- file(figure_log)
writeLines(paste("Figure done in", Sys.Date(), ":", figure_log), filelog)
close(filelog)
