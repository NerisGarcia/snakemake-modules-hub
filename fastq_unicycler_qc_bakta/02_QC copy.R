# =============================================================================
# Genomic Description of the Isolate Sequences and Genomes
# REFACTORED: ALL ANALYSIS FIRST, THEN ALL PLOTS AT END
# =============================================================================

# Libraries -------------------------------------------------------------------
library(dplyr)
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

source("code/design.R")

# ── Paths ──────────────────────────────────────────────────────────────
plot_dir <- "figures/1_QC_Assembly_mapping"
itol_dir <- "figures/itol"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
#dir.create(itol_dir, recursive = TRUE, showWarnings = FALSE)

# ————————————————————————————————————————————————————————————————————————
# SECTION 1: DATA LOADING AND METADATA
# ————————————————————————————————————————————————————————————————————————

# Load manifest data
# isoaltes_metadata <- read.csv("data/20260423.01.ZB_Lcrispatus_manifest_parsed.csv")

args <- commandArgs(trailingOnly = TRUE)

isoaltes_metadata.file <- args[1]
raw_metadata_out <- if (length(args) >= 7 && nzchar(args[7])) args[7] else "data/Lcrispatus417_NCBI.raw_reads_QC_metadata_combined.ods"
qc_metadata_out <- if (length(args) >= 8 && nzchar(args[8])) args[8] else "data/Lcrispatus417.metadata.QC_stats_combined.ods"
qc_criteria_out <- if (length(args) >= 9 && nzchar(args[9])) args[9] else "data/Lcrispatus417.metadata.QC_outlier_criteria_summary.ods"

write_ods_file <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_ods(x, path)
}

isoaltes_metadata <- read.csv(isoaltes_metadata.file)


isoaltes_metadata %>% head()
backuptable(isoaltes_metadata.file)

# ——————————————————————————————————————————————————————————————————————---------
# SECTION 2: RAW AND SUBSAMPLED READS QC ANALYSIS
# ——————————————————————————————————————————————————————————————————————---------

# Load raw reads QC (for reference - not used for filtering)
# qc_raw_path <- "data/0_raw_reads_QC/ZB_Lcrispatus417.reads.seqkitstats.csv"

qc_raw_path <- args[2]
backuptable(qc_raw_path)

qc_raw <- read.table(
  qc_raw_path,
  sep = "",
  header = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE
) %>%
  filter(!is.na(V1) & V1 != "file") %>%
  transmute(
    CultureID = V1 %>% basename() %>% str_remove("\\.R[12]\\.fq\\.gz$"),
    format = V2,
    type = V3,
    raw_num_seqs = as.numeric(str_remove_all(V4, ",")),
    raw_sum_len = as.numeric(str_remove_all(V5, ",")),
    raw_min_len = as.numeric(V6),
    raw_avg_len = as.numeric(V7),
    raw_max_len = as.numeric(V8),
    raw_total_reads = raw_num_seqs * 2,
    raw_avg_depth = (raw_total_reads * raw_avg_len) / 2416053
  ) %>%
  select(c(CultureID, raw_total_reads, raw_avg_len, raw_avg_depth))

head(qc_raw)

# Merge raw reads with manifest
qc.metadata.raw_temp <- isoaltes_metadata %>%
  left_join(qc_raw, by = "CultureID")

qc.metadata.raw_temp %>% head()

# Load subsampled reads QC (fastp processed, ~3M reads - actual data used in pipeline)
qc_sub_path <- args[3]

# qc_sub_path <- "data/0_raw_reads_QC/ZB_Lcrispatus417.reads.sub3M.seqkitstats.csv"

backuptable(qc_sub_path)

qc_subsampled <-
  read.table(
    qc_sub_path,
    sep = "",
    header = FALSE,
    stringsAsFactors = FALSE,
    fill = TRUE
  ) %>%
  filter(!is.na(V1) & V1 != "file") %>%
  transmute(
    CultureID = V1 %>% basename() %>% str_remove("\\.sub3M\\.R[12]\\.fq\\.gz$"),
    format = V2,
    type = V3,
    num_seqs = as.numeric(str_remove_all(V4, ",")),
    sum_len = as.numeric(str_remove_all(V5, ",")),
    min_len = as.numeric(V6),
    fastp_avg_len = as.numeric(V7),
    max_len = as.numeric(V8),
    fastp_total_reads = num_seqs * 2,
    fastp_avg_depth = (fastp_total_reads * fastp_avg_len) / 2416053
  ) %>%
  select(-c(format, type, min_len, max_len, num_seqs, sum_len))

head(qc_subsampled)

# Final metadata: raw reads (for reference) + subsampled reads (for analysis/filtering)
qc.metadata.raw <- qc.metadata.raw_temp %>%
  left_join(qc_subsampled, by = "CultureID") %>%
  select(
    CultureID, SampleID, Cohort, ega_del_week, preterm, subCST, SamplePrefix, ColonyType, ColonyNum,
    PTB_category, PTB_category_det, raw_avg_len, raw_total_reads, raw_avg_depth, fastp_avg_len, fastp_total_reads, fastp_avg_depth
  )

head(qc.metadata.raw)
dim(qc.metadata.raw)
# 417

# Subsampled read flags ——————————————————————————————————————————
# Outlier detection using 3 SD threshold on subsampled reads (the data actually used in pipeline)
subsampled_metrics <- c("fastp_avg_len", "fastp_total_reads", "fastp_avg_depth")

subsampled_reads_outlier_per_metric <- qc.metadata.raw %>%
  select(CultureID, all_of(subsampled_metrics)) %>%
  pivot_longer(cols = all_of(subsampled_metrics), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(metric) %>%
  mutate(
    mean_v = mean(value, na.rm = TRUE),
    sd_v = sd(value, na.rm = TRUE),
    lower_threshold = mean_v - 3 * sd_v,
    upper_threshold = mean_v + 3 * sd_v,
    is_outlier = value < lower_threshold | value > upper_threshold
  ) %>%
  ungroup() %>%
  select(CultureID, metric, value, mean_v, sd_v, lower_threshold, upper_threshold, is_outlier)

subsampled_reads_outlier_per_metric %>%
  select(metric, lower_threshold, upper_threshold, mean_v, sd_v) %>%
  distinct()

# metric            lower_threshold upper_threshold   mean_v      sd_v
# fastp_avg_len               128.             152.     140.      3.96
# fastp_total_reads        332275.         4741604. 2536940. 734888.
# fastp_avg_depth              20.9            273.     147.     41.9

subsampled_outlier_samples <- subsampled_reads_outlier_per_metric %>%
  filter(is_outlier) %>%
  group_by(CultureID) %>%
  summarise(
    subsampled_reads_outlier_flags = str_c(sort(unique(metric)), collapse = "; "),
    .groups = "drop"
  )

subsampled_outlier_samples %>%
  count(subsampled_reads_outlier_flags)

# subsampled_reads_outlier_flags         n
# fastp_avg_depth; fastp_total_reads    10
# fastp_avg_len                          5

temp <- subsampled_outlier_samples %>%
  filter(subsampled_reads_outlier_flags == "fastp_avg_len") %>%
  left_join(qc.metadata.raw, by = "CultureID")

rm(temp)

raw_qc.metrics_to_filter <- c("fastp_avg_depth", "fastp_total_reads")

# Create metadata table with subsampled reads quality flags
qc.metadata.flagged.rawreads <- qc.metadata.raw %>%
  left_join(subsampled_outlier_samples, by = "CultureID") %>%
  mutate(
    FILTER = case_when(
      grepl(paste(raw_qc.metrics_to_filter, collapse = "|"), subsampled_reads_outlier_flags)
      ~ subsampled_reads_outlier_flags,
      TRUE ~ NA_character_
    ),
    filter_by = case_when(
      grepl(paste(raw_qc.metrics_to_filter, collapse = "|"), subsampled_reads_outlier_flags)
      ~ "fastp",
      TRUE ~ NA_character_
    )
  ) %>%
  rename(fastp_flags = subsampled_reads_outlier_flags) %>%
  mutate(FLAGS = fastp_flags)

qc.metadata.flagged.rawreads %>% head()

qc.metadata.flagged.rawreads %>% count(fastp_flags)

qc.metadata.flagged.rawreads %>% count(FILTER)


# ——————————————————————————————————————————————————————————————————————---------
# SECTION 3: KRAKEN2 ANALYSIS
# ——————————————————————————————————————————————————————————————————————---------

species_file <- args[4]
#species_file <- "data/0_raw_reads_QC/3_Kraken2/Lcrispatus402.Kraken.all.summary"


backuptable(species_file)

kraken_raw <- read.table(species_file,
  sep = "\t", header = TRUE,
  quote = "", fill = TRUE, stringsAsFactors = FALSE
)
head(kraken_raw)

kraken_sp <- kraken_raw %>%
  filter(rank == "S") %>%
  mutate(name = trimws(name), percent = as.numeric(percent))

kraken_raw %>%
  group_by(sample) %>%
  filter(cur_group_id() == 1) %>%
  group_by(rank) %>%
  summarise(tr = sum(percent))
head


# Master Kraken table: per-sample LC%, other species%, unclassified%
kraken_master <- kraken_raw %>%
  mutate(
    percent = as.numeric(percent),
    taxid = as.numeric(taxid),
    rank = trimws(rank)
  ) %>%
  group_by(sample) %>%
  summarise(
    lc_pct = sum(percent[rank == "S" & taxid == 47770], na.rm = TRUE),
    other_lactobacillus_pct = sum(percent[rank == "S" & taxid != 47770 & grepl("^Lactobacillus", name)], na.rm = TRUE),
    other_species_pct = sum(percent[rank == "S" & taxid != 47770], na.rm = TRUE),
    l_genus_pct = sum(percent[rank == "G" & taxid == 1578], na.rm = TRUE),
    other_l_genus_pct = sum(percent[rank == "G" & taxid != 1578], na.rm = TRUE),
    unclassified_pct = sum(percent[rank == "U"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lc_mean = mean(lc_pct, na.rm = TRUE),
    lc_sd = sd(lc_pct, na.rm = TRUE),
    os_mean = mean(other_species_pct, na.rm = TRUE),
    os_sd = sd(other_species_pct, na.rm = TRUE),
    lg_mean = mean(l_genus_pct, na.rm = TRUE),
    lg_sd = sd(l_genus_pct, na.rm = TRUE),
    olg_mean = mean(other_l_genus_pct, na.rm = TRUE),
    olg_sd = sd(other_l_genus_pct, na.rm = TRUE),
    uc_mean = mean(unclassified_pct, na.rm = TRUE),
    uc_sd = sd(unclassified_pct, na.rm = TRUE),
    lc_low_fence = lc_mean - 2 * lc_sd,
    os_high_fence = os_mean + 2 * os_sd,
    l_genus_low_fence = lg_mean - 2 * lg_sd,
    other_l_genus_high_fence = olg_mean + 1 * olg_sd,
    uc_high_fence = uc_mean + 2 * uc_sd,
    low_lc_outlier = lc_pct < lc_low_fence,
    high_os_outlier = other_species_pct > os_high_fence,
    low_l_genus_outlier = l_genus_pct < l_genus_low_fence,
    high_other_l_genus_outlier = other_l_genus_pct > other_l_genus_high_fence,
    high_uc_outlier = unclassified_pct > uc_high_fence,
    kraken_flag = str_c(
      if_else(low_lc_outlier, "LOW Lcrispatus PCT; ", ""),
      if_else(high_os_outlier, "HIGH Other Species PCT; ", ""),
      if_else(low_l_genus_outlier, "LOW L Genus PCT; ", ""),
      if_else(high_other_l_genus_outlier, "HIGH Other Genus PCT; ", ""),
      if_else(high_uc_outlier, "HIGH Unclassified PCT; ", "")
    ),
    kraken_flag = str_trim(str_replace(kraken_flag, ";\\s*$", "")),
    kraken_flag = na_if(kraken_flag, "")
  )

kraken_master %>% head()

# Extract metric-level fences and per-sample values for Kraken metrics
kraken_metric_fences <- kraken_master %>%
  summarise(
    lc_mean = first(lc_mean),
    lc_sd = first(lc_sd),
    os_mean = first(os_mean),
    os_sd = first(os_sd),
    lg_mean = first(lg_mean),
    lg_sd = first(lg_sd),
    olg_mean = first(olg_mean),
    olg_sd = first(olg_sd),
    uc_mean = first(uc_mean),
    uc_sd = first(uc_sd),
    lc_low_fence = first(lc_low_fence),
    os_high_fence = first(os_high_fence),
    l_genus_low_fence = first(l_genus_low_fence),
    other_l_genus_high_fence = first(other_l_genus_high_fence),
    uc_high_fence = first(uc_high_fence)
  ) %>%
  pivot_longer(cols = everything(), names_to = "stat", values_to = "value")

kraken_metric_fences %>% arrange(stat)

# A tibble: 15 × 2#stat                      value#<chr>                     <dbl>
# l_genus_low_fence        82.6
# lc_low_fence             76.4
# lc_mean                  88.1
# lc_sd                     5.85
# lg_mean                  94.1
# lg_sd                     5.78
# olg_mean                  0.296
# olg_sd                    1.74
# os_high_fence             4.30 *
# os_mean                   0.793
# os_sd                     1.75
# other_l_genus_high_fence  2.03*
# uc_high_fence            14.9
# uc_mean                   3.98
# uc_sd                     5.46

kraken_master %>%
  filter(grepl("HIGH Other Genus PCT", kraken_flag)) %>%
  select(sample, other_l_genus_pct)

kraken_master %>%
  arrange(desc(other_l_genus_pct)) %>%
  select(sample, other_l_genus_pct, kraken_flag)

kraken_master  %>% 
  mutate(
    lc_category = case_when(
      as.numeric(lc_pct) >= 90 ~ "above 90",
      as.numeric(lc_pct) >= 80 ~ "above 80",
      TRUE ~ "else"
    )
  )  %>% 
  group_by(lc_category)  %>% 
  summarise(
    count = n(),
    lc_categry_percent = n() / nrow(kraken_master) * 100
  )

#  lc_category  count lc_categry_percent
#  above 80       125              30.0
#  above 90       254              60.9
#  else           38               9.11

kraken_master %>%
  mutate(
    os_category = case_when(
      as.numeric(other_species_pct) >= 4.3 ~ "above 4.3 - os_high_fence",
      as.numeric(other_species_pct) >= 2 ~ "above 1",
      TRUE ~ "else"
    )
  ) %>%
  group_by(os_category) %>%
  summarise(
    count = n(),
    os_category_percent = n() / nrow(kraken_master) * 100
  )
#  os_category               count os_category_percent
# above 1                      36                8.63
# above 4.3 - os_high_fence     6                1.44
# else                        375               89.9

kraken_master  %>% 
  mutate(
    olg_category = case_when(
      as.numeric(other_l_genus_pct) >= 5 ~ "above 5",
      as.numeric(other_l_genus_pct) >= 2 ~ "above 2",
      TRUE ~ "else"
    )
  )  %>% 
  group_by(olg_category)  %>% 
  summarise(
    count = n(),
    olg_category_percent = n() / nrow(kraken_master) * 100
  )

# olg_category count olg_category_percent
# <chr>        <int>                <dbl>
# above 2          2                0.480
# above 5          5                1.20
# else           410               97.8

kraken.metric_to_filter <- c("HIGH Other Genus PCT")

# Add kraken flags to metadata
qc.metadata.flagged.kraken <- qc.metadata.flagged.rawreads %>%
  left_join(
    kraken_master %>%
      select(
        sample, lc_pct, other_species_pct, l_genus_pct, other_l_genus_pct,
        unclassified_pct, kraken_flag
      ),
    by = c("CultureID" = "sample")
  ) %>%
  mutate(
    FLAGS = paste(FLAGS, kraken_flag, sep = "; "),
    FILTER = case_when(
      is.na(FILTER) & grepl(kraken.metric_to_filter, kraken_flag) ~ kraken_flag,
      TRUE ~ FILTER
    ),
    filter_by = case_when(
      is.na(FILTER) & grepl(kraken.metric_to_filter, kraken_flag) ~ "KRAKEN",
      TRUE ~ FILTER
    )
  )
qc.metadata.flagged.kraken %>% count(FILTER)

write_ods_file(qc.metadata.flagged.kraken, raw_metadata_out)


# ——————————————————————————————————————————————————————————————————————---------
# SECTION 4: MAPPING ANALYSIS
# ——————————————————————————————————————————————————————————————————————---------

#mapping_stats_file <- "data/2_mapping/2_snippy_cores/Lcrispatus_417/Lcrispatus_417.txt"

mapping_stats_file <- args[5]

backuptable(mapping_stats_file)

mapping_stats <- read.csv(mapping_stats_file, sep = "\t") %>%
  mutate(ID = str_remove(ID, "_snippy")) %>%
  mutate(
    COV_PCT = ALIGNED / LENGTH * 100,
    UNALIGNED_PCT = 100 * UNALIGNED / LENGTH,
    LOWCOV_PCT = 100 * LOWCOV / LENGTH
  )

mapping_stats %>% head()
mapping_stats %>% dim()

# Mapping flags
mapping_metrics <- c("COV_PCT", "UNALIGNED_PCT", "LOWCOV_PCT", "VARIANT", "HET")

mapping_outlier_per_metric <- mapping_stats %>%
  select(ID, all_of(mapping_metrics)) %>%
  pivot_longer(cols = all_of(mapping_metrics), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(metric) %>%
  mutate(
    mean_v = mean(value, na.rm = TRUE),
    sd_v = sd(value, na.rm = TRUE),
    lower_threshold = mean_v - 3 * sd_v,
    upper_threshold = mean_v + 3 * sd_v,
    is_outlier = value < lower_threshold | value > upper_threshold
  ) %>%
  ungroup() %>%
  select(ID, metric, value, mean_v, sd_v, lower_threshold, upper_threshold, is_outlier)

mapping_outlier_per_metric %>%
  filter(is_outlier) %>%
  count(metric, is_outlier)

mapping_outlier_samples <- mapping_outlier_per_metric %>%
  filter(is_outlier) %>%
  group_by(ID) %>%
  summarise(
    mapping_flags = str_c(sort(unique(metric)), collapse = "; "),
    .groups = "drop"
  )

mapping_outlier_samples %>% count(mapping_flags)

mapping_outlier_per_metric %>%
  select(metric, lower_threshold, upper_threshold, mean_v, sd_v) %>%
  distinct()
#  metric        lower_threshold upper_threshold  mean_v    sd_v
#  COV_PCT                 40.6            113.    76.6    12.0
# UNALIGNED_PCT           -1.83            44.0   21.1     7.64
# LOWCOV_PCT             -25.6             30.1    2.23    9.27
# VARIANT              -5749.            8297.  1274.   2341.
# HET                  -4928.            6861.   967.   1965.

qc.metadata.mapping.raw <-
  full_join(mapping_stats, mapping_outlier_samples, by = c("ID" = "ID"))

qc.metadata.mapping.raw %>% head()

mapping.metrics_to_filter <- c("HET")

# Mapping flags

qc.metadata.flagged.mapping <- qc.metadata.flagged.kraken %>%
  left_join(qc.metadata.mapping.raw, by = c("CultureID" = "ID")) %>%
  mutate(
    FLAGS = paste(FLAGS, mapping_flags, sep = "; "),
    FILTER = case_when(
      is.na(FILTER) & grepl(paste0(mapping.metrics_to_filter, collapse = "|"), mapping_flags) ~ mapping_flags,
      TRUE ~ FILTER
    ),
    filter_by = case_when(
      is.na(FILTER) & grepl(paste0(mapping.metrics_to_filter, collapse = "|"), mapping_flags) ~ "MAPPING",
      TRUE ~ FILTER
    )
  )

qc.metadata.flagged.mapping %>%
  arrange(desc(FILTER)) %>%
  head(20)

qc.metadata.flagged.mapping %>%
  arrange(desc(FILTER)) %>%
  count(FILTER)


# ——————————————————————————————————————————————————————————————————————---------
# SECTION 5: ASSEMBLY ANALYSIS
# ——————————————————————————————————————————————————————————————————————---------

#quast_summary_file <- "data/1_Assembly/2_QC/1_Quast/quast_summary.tsv"

quast_summary_file <- args[6]

quast_stats <- read.csv(quast_summary_file, sep = "\t") %>%
  mutate(CultureID = str_remove(Assembly, "\\.contigs$")) %>%
  rename(
    `Contigs_1000bp` = `X..contigs.....1000.bp.`,
    `Contigs_10000bp` = `X..contigs.....10000.bp.`,
    `Contigs_50000bp` = `X..contigs.....50000.bp.`,
    `Total_Length_1000bp` = `Total.length.....1000.bp.`,
    `Total_Length_10000bp` = `Total.length.....10000.bp.`,
    `Total_Length_50000bp` = `Total.length.....50000.bp.`,
    `Largest_Contig` = `Largest.contig`,
    `GC_Content` = `GC....`,
    `N90` = `N90`,
    `L50` = `L50`,
    `N_per_100kbp` = `X..N.s.per.100.kbp`,
    `Predicted_Genes_0bp` = `X..predicted.genes.....0.bp.`,
    `Predicted_Genes_1500bp` = `X..predicted.genes.....1500.bp.`,
    `Contigs_0bp` = `X..contigs.....0.bp.`,
    `Contigs_5000bp` = `X..contigs.....5000.bp.`,
    `Contigs_25000bp` = `X..contigs.....25000.bp.`,
    `Total_Length_0bp` = `Total.length.....0.bp.`,
    `Total_Length_5000bp` = `Total.length.....5000.bp.`,
    `Total_Length_25000bp` = `Total.length.....25000.bp.`,
    `Contigs` = `X..contigs`,
    `Total_Length` = `Total.length`,
    `N50` = `N50`,
    `auN` = `auN`,
    `L90` = `L90`,
    `Predicted_Genes_Unique` = `X..predicted.genes..unique.`,
    `Predicted_Genes_300bp` = `X..predicted.genes.....300.bp.`,
    `Predicted_Genes_3000bp` = `X..predicted.genes.....3000.bp.`
  ) %>%
  mutate(
    Predicted_Genes_0bp_Complete = as.numeric(str_extract(Predicted_Genes_0bp, "^\\s*\\d+")),
    Predicted_Genes_0bp_Partial = as.numeric(str_match(Predicted_Genes_0bp, "\\+\\s*(\\d+)\\s*part")[, 2]),
    Predicted_Genes_300bp_Complete = as.numeric(str_extract(Predicted_Genes_300bp, "^\\s*\\d+")),
    Predicted_Genes_300bp_Partial = as.numeric(str_match(Predicted_Genes_300bp, "\\+\\s*(\\d+)\\s*part")[, 2]),
    Predicted_Genes_3000bp_Complete = as.numeric(str_extract(Predicted_Genes_3000bp, "^\\s*\\d+")),
    Predicted_Genes_3000bp_Partial = as.numeric(str_match(Predicted_Genes_3000bp, "\\+\\s*(\\d+)\\s*part")[, 2])
  )

quast_stats %>% head()

assembly_metrics <- c(
  "Contigs",
  "Total_Length",
  "GC_Content",
  "Predicted_Genes_Unique",
  "Predicted_Genes_0bp_Complete",
  "Predicted_Genes_0bp_Partial"
)

# Assembly flags
assembly_outlier_per_metric <- quast_stats %>%
  select(CultureID, all_of(assembly_metrics)) %>%
  pivot_longer(cols = all_of(assembly_metrics), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(metric) %>%
  mutate(
    mean_v = mean(value, na.rm = TRUE),
    sd_v = sd(value, na.rm = TRUE),
    lower_threshold = mean_v - 2 * sd_v,
    upper_threshold = mean_v + 2 * sd_v,
  ) %>%
  mutate(lower_threshold = case_when(
    metric == "Contigs" ~ 0
  )) %>%
  mutate(is_outlier = value < lower_threshold | value > upper_threshold) %>%
  ungroup() %>%
  select(CultureID, metric, value, mean_v, sd_v, lower_threshold, upper_threshold, is_outlier) %>%
  filter(is_outlier) %>%
  group_by(CultureID) %>%
  mutate(
    assembly_flags = str_c(sort(unique(metric)), collapse = "; "),
    .groups = "drop"
  ) %>%
  ungroup()  %>% 
  select(CultureID, assembly_flags)

assembly_outlier_per_metric %>%
  select(metric, lower_threshold, upper_threshold, mean_v, sd_v) %>%
  distinct()
# A tibble: 6 × 5
#  metric                       lower_threshold upper_threshold    mean_v    sd_v
# Contigs                                    0          439.      2.51e2 9.39e+1
# Total_Length                              NA      3019198.      2.24e6 3.89e+5
# GC_Content                                NA           37.3     3.68e1 2.46e-1
# Predicted_Genes_Unique                    NA         1356.      9.81e2 1.88e+2
# Predicted_Genes_0bp_Complete              NA         1351.      9.78e2 1.87e+2
# Predicted_Genes_0bp_Partial               NA            9.53    2.61e0 3.46e+0

ass.metrics_to_filter <- c(
  "Contigs",
  "Total_length",
  "GC_content",
  "Predicted_genes_unique",
  "Predicted_genes_0bp_complete",
  "Predicted_genes_0bp_partial"
)

quast_stats %>%
  select(CultureID, all_of(assembly_metrics)) %>% 
  left_join(assembly_outlier_per_metric, by = "CultureID") %>%
  count(assembly_flags)

# FINAL COMBINED TABLE with all FILTER information
qc.metadata.ass <-
quast_stats %>%
  select(CultureID, all_of(assembly_metrics)) %>% 
  left_join(assembly_outlier_per_metric, by = "CultureID")

qc.metadata.ass  %>%  head

qc.metadata.ass %>%
  count(assembly_flags)

qc.metadata.final <-
  qc.metadata.flagged.mapping %>%
  left_join(qc.metadata.ass, by = "CultureID") %>%
  relocate(FILTER, .after = last_col()) %>%
  mutate(
    FLAGS = paste(FLAGS, assembly_flags, sep = "; "),
    FILTER = case_when(
      is.na(FILTER) & grepl(paste0(ass.metrics_to_filter, collapse = "|"), assembly_flags) ~ assembly_flags,
      TRUE ~ FILTER
    ),
    filter_by = case_when(
      is.na(FILTER) & grepl(paste0(ass.metrics_to_filter, collapse = "|"), assembly_flags) ~ "ASSEMBLY",
      TRUE ~ FILTER
    )
  ) %>%
  mutate(FILTER = ifelse(is.na(FILTER), "PASS", FILTER))  %>% 
  ungroup()  %>% 
  distinct()  

qc.metadata.final %>% head()

qc.metadata.final %>%
  count(FILTER)

write_ods_file(qc.metadata.final, qc_metadata_out)

#backuptable("data/*.Lcrispatus417.metadata.QC_stats_combined.ods")

qc.metadata.final %>%
  count(FILTER)


# ——————————————————————————————————————————————————————————————————————---------
# SECTION 6: OUTLIER CRITERIA SUMMARY TABLE
# ——————————————————————————————————————————————————————————————————————---------

# fastp (subsampled reads) ——————————————————————————————————————————————
fastp_criteria <- subsampled_reads_outlier_per_metric %>%
  select(metric, mean_v, sd_v, lower_threshold, upper_threshold) %>%
  distinct() %>%
  mutate(
    software       = "fastp",
    direction      = "both",
    n_sd           = (upper_threshold - mean_v) / sd_v,
    filter_applied = metric %in% raw_qc.metrics_to_filter
  )

# Kraken2 ———————————————————————————————————————————————————————————————
# Fences use asymmetric SD multipliers: lc/os/lg/uc = 2 SD, other_l_genus = 1 SD
kraken_ref <- kraken_master %>% slice(1)

kraken_criteria <- tibble(
  software        = "Kraken2",
  metric          = c("lc_pct", "other_species_pct", "l_genus_pct", "other_l_genus_pct", "unclassified_pct"),
  mean_v          = c(kraken_ref$lc_mean, kraken_ref$os_mean, kraken_ref$lg_mean, kraken_ref$olg_mean, kraken_ref$uc_mean),
  sd_v            = c(kraken_ref$lc_sd, kraken_ref$os_sd, kraken_ref$lg_sd, kraken_ref$olg_sd, kraken_ref$uc_sd),
  lower_threshold = c(kraken_ref$lc_low_fence, NA_real_, kraken_ref$l_genus_low_fence, NA_real_, NA_real_),
  upper_threshold = c(NA_real_, kraken_ref$os_high_fence, NA_real_, kraken_ref$other_l_genus_high_fence, kraken_ref$uc_high_fence)
) %>%
  mutate(
    direction = case_when(
      !is.na(lower_threshold) & !is.na(upper_threshold) ~ "both",
      !is.na(lower_threshold) ~ "low",
      TRUE ~ "high"
    ),
    n_sd = abs(coalesce(lower_threshold, upper_threshold) - mean_v) / sd_v,
    filter_applied = metric == "other_l_genus_pct"
  )

# snippy (mapping) —————————————————————————————————————————————————————
# Note: current filter logic uses all mapping_metrics (not mapping.metrics_to_filter)
mapping_criteria <- mapping_outlier_per_metric %>%
  select(metric, mean_v, sd_v, lower_threshold, upper_threshold) %>%
  distinct() %>%
  mutate(
    software       = "snippy",
    direction      = "both",
    n_sd           = (upper_threshold - mean_v) / sd_v,
    filter_applied = metric %in% mapping_metrics
  )

# QUAST (assembly) —————————————————————————————————————————————————————
# lower_threshold is NA for all non-Contigs metrics (case_when without TRUE clause)
# Contigs lower is set to 0 (practically never triggered); all effectively high-only
assembly_criteria <- assembly_outlier_per_metric %>%
  select(metric, mean_v, sd_v, lower_threshold, upper_threshold) %>%
  distinct() %>%
  mutate(
    software       = "QUAST",
    direction      = "high",
    n_sd           = (upper_threshold - mean_v) / sd_v,
    filter_applied = metric %in% ass.metrics_to_filter
  )

# Combined outlier criteria summary ————————————————————————————————————
outlier_criteria_summary <- bind_rows(
  fastp_criteria,
  kraken_criteria,
  mapping_criteria,
  assembly_criteria
) %>%
  select(
    software, metric, mean_v, sd_v, n_sd, direction,
    lower_threshold, upper_threshold, filter_applied
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

outlier_criteria_summary

write_ods_file(outlier_criteria_summary, qc_criteria_out)

#backuptable("data/20260518.01.Lcrispatus417.metadata.QC_outlier_criteria_summary.ods")

qc.metadata.final  %>%  distinct  %>% group_by(CultureID)  %>%  count()  %>%  arrange(desc(n))
