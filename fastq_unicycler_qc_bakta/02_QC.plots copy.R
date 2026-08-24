# Libraries -------------------------------------------------------------------
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(readODS)
library(patchwork)
library(phangorn)
library(cowplot)

# personal libraries
library(Itools)
library(customR)

source("code/design.R")

plot_dir <- "figures/1_QC_Assembly_mapping"
# =============================================================================
# SECTION 6: ALL PLOTTING (using qc.metadata.final with complete FILTER info)
# =============================================================================

# mash distances —————————————————————————————————————————————————————————


args <- commandArgs(trailingOnly = TRUE)
mash_dist <- args[1]
# Create symbolic link in data/input_files


# mash_dist <- "data/1_Assembly/4_mashtree/Lcrispatus_417.mash_distances.tsv"
# backuptable(mash_dist)

mash_dist.df <-
  read.csv(mash_dist, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  rename(colony1 = ".") %>%
  pivot_longer(
    cols = -colony1,
    names_to = "colony2",
    values_to = "mash_distance"
  ) %>%
  mutate(sample1 = str_remove(colony1, "[A-Z][0-9].contigs"), sample2 = str_remove(colony2, "[A-Z][0-9].contigs")) %>%
  # deal with dup files
  mutate(
    sample1 = str_replace(colony1, "[A-Z][0-9](-dup)?\\.contigs$", "\\1"),
    sample2 = str_replace(colony2, "[A-Z][0-9](-dup)?\\.contigs$", "\\1")
  ) %>%
  mutate(type = case_when(
    colony1 == colony2 ~ "self",
    sample1 == sample2 ~ "intra-sample",
    sample1 != sample2 ~ "inter-sample"
  ))


p_mash1 <-
  mash_dist.df %>%
  dplyr::slice_sample(prop = 0.1) %>%
  ggplot() +
  geom_histogram(aes(x = mash_distance, fill = type),
    binwidth = 0.0001,
    boundary = 0,
    closed = "left",
    pad = TRUE,
    color = "black"
  ) +
  labs(
    title = "Assembly mash distances",
    x = NULL,
    y = "Comparisons",
    caption = "417 colonies. In A, only 10% of comparisons shown"
  ) +
  coord_cartesian(xlim = c(0, 0.03), ylim = c(0, NA), expand = FALSE) +
  geom_vline(xintercept = 0.0002, linetype = "dashed", color = "black", linewidth = 0.2) +
  geom_text(x = 0.0002, y = 500, label = "0.0002", hjust = -0.1, vjust = -1, color = "black", size = 3) +
  theme(legend.position = "none")

# p_mash1

p_mash2 <-
  mash_dist.df %>%
  filter(mash_distance <= 0.003) %>%
  ggplot() +
  geom_histogram(aes(x = mash_distance, fill = type),
    binwidth = 0.0001,
    boundary = 0,
    closed = "left",
    pad = TRUE,
    color = "black"
  ) +
  labs(
    x = "mash distance",
    y = "Comparisons"
  ) +
  coord_cartesian(xlim = c(0, 0.003), ylim = c(0, NA), expand = FALSE) +
  geom_vline(xintercept = 0.0002, linetype = "dashed", color = "black", linewidth = 0.2) +
  geom_text(x = 0.0002, y = 500, label = "0.0002", hjust = -0.1, vjust = -1, color = "black", size = 3)

# p_mash2

p_mash <- p_mash1 / p_mash2

# p_mash

saveplot(
  p_mash,
  file.path(plot_dir, "mash_distances_histogram"),
  layout = "barplot",
  height_cm = 20,
  mode = "editing"
)

# QC!  ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# qc.metadata.final <- read_ods(Sys.glob("data/Lcrispatus417.metadata.QC_stats_combined.ods"))

qc.metadata.final <-
  read_ods(Sys.glob(args[2]))

# qc.filters <- read_ods(Sys.glob("data/*.Lcrispatus417.metadata.QC_outlier_criteria_summary.ods"))


qc.filters <- read_ods(args[3])

# ———— Table —————————————————————————————————————————————————————————————

qc.metadata.final %>%
  count(FILTER)  %>% 
  arrange(desc(n))

# FILTER                                                                      n
# Contigs; GC_Content; Predicted_Genes_0bp_Complete; Predicted_Genes_0bp…     1
# Contigs                                                                     2
# Contigs; Predicted_Genes_0bp_Partial                                        2
# HET                                                                         5
# HIGH Other Genus PCT                                                        2
# HIGH Other Species PCT; HIGH Other Genus PCT                                1
# LOW Lcrispatus PCT; HIGH Other Species PCT; LOW L Genus PCT; HIGH Othe…     4
# Predicted_Genes_0bp_Partial                                                 1
# Total_Length                                                                1
# fastp_avg_depth; fastp_total_reads                                         10
# PASS                                                                      388


qc.metadata.final %>%
  group_by(SampleID) %>%
  mutate(ncolonies = n()) %>%
  filter(FILTER != "PASS") %>%
  mutate(nremoved = n()) %>%
  select(SampleID, ncolonies, nremoved) %>%
  distinct() %>%
  arrange(desc(nremoved)) %>%
  filter(ncolonies == nremoved)


qc.metadata.final %>%
  filter(FILTER == "PASS") %>%
  count()
# 388

qc.metadata.final %>%
  filter(FILTER == "PASS") %>%
  distinct(SampleID) %>%
  count()
# 88

# ———— filter colors and shapes ——————————————————————————————————————————
# Plot 1: Total reads barplot

# Shared FILTER style vectors (used across plots)
filter_levels_global <- qc.metadata.final %>%
  mutate(FILTER = ifelse(is.na(FILTER), "PASS", FILTER)) %>%
  distinct(FILTER) %>%
  pull(FILTER)

filter_shapes <- setNames(
  ifelse(filter_levels_global == "PASS", 1, 23),
  filter_levels_global
)

filter_sizes <- setNames(
  ifelse(filter_levels_global == "PASS", 1, 2),
  filter_levels_global
)

filter_colors_bw <- setNames(
  ifelse(filter_levels_global == "PASS", "lightblue", "black"),
  filter_levels_global
)
filter_fill_bw <- filter_colors_bw

filter_levels_nonpass <- filter_levels_global[filter_levels_global != "PASS"]
filter_palette_nonpass <- grDevices::hcl.colors(max(1, length(filter_levels_nonpass)), palette = "Dark 3")

filter_colors_type <- c(
  "PASS" = "lightblue",
  setNames(filter_palette_nonpass[seq_along(filter_levels_nonpass)], filter_levels_nonpass)
)
filter_colors_type <- filter_colors_type[filter_levels_global]
filter_fill_type <- filter_colors_type

# Plot: FILTER-to-color mapping key
filter_color_key <- tibble(
  FILTER = names(filter_colors_type),
  value = 1
) %>%
  mutate(
    FILTER = stringr::str_wrap(as.character(FILTER), width = 35),
    FILTER = factor(FILTER, levels = unique(FILTER))
  )

filter_colors_type_wrapped <- setNames(
  unname(filter_colors_type),
  stringr::str_wrap(names(filter_colors_type), width = 35)
)

p_filter_color_type_key <- ggplot(filter_color_key, aes(x = FILTER, y = value, fill = FILTER)) +
  geom_tile(width = 0.8, colour = "black") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = filter_colors_type_wrapped, guide = "none") +
  scale_y_continuous(limits = c(0.5, 2)) +
  labs(
    title = "FILTER color key",
    x = NULL,
    y = NULL
  )

print(p_filter_color_type_key)

saveplot(
  p_filter_color_type_key,
  file.path(plot_dir, "FILTER.color_key"),
  layout = "medium",
  mode = "editing"
)

# Plot histogram of FILTER types
pfilter <-
  qc.metadata.final %>%
  mutate(FILTER = ifelse(is.na(FILTER), "PASS", FILTER)) %>%
  filter(FILTER != "PASS") %>%
  count(FILTER) %>%
  mutate(FILTER = factor(FILTER, levels = names(filter_colors_type))) %>%
  ggplot(aes(x = FILTER, y = n, fill = FILTER)) +
  coord_flip(clip = "off") +
  geom_col(color = "black") +
  geom_text(aes(label = n), hjust = -0.6) +
  scale_fill_manual(values = filter_colors_type) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(as.character(x), width = 35)) +
  labs(
    title = "Number of colonies per FILTER category",
    x = "FILTER category",
    y = "Number of colonies"
  ) +
  theme(legend.position = "none", plot.margin = margin(r = 10))
pfilter

saveplot(
  pfilter,
  file.path(plot_dir, "FILTER.histogram"),
  layout = "barplot",
  mode = "editing"
)



# ———— TRee —————————————————————————————————————————————————————

# tree_file <- "data/3_pangenomes/Lcrispatus.417/core_gene_alignment_filtered.treefile"

tree_file <- args[4]


tree <- ape::read.tree(tree_file)

tree_tip_metadata <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(
    qc.metadata.final %>%
      filter(!is.na(Cohort)) %>%
      select(CultureID, FILTER) %>%
      distinct(CultureID, .keep_all = TRUE) %>%
      rename(label = CultureID),
    by = "label"
  ) %>%
  mutate(
    is_filtered = !is.na(FILTER),
    filter_type = ifelse(is.na(FILTER), "PASS", FILTER)
  ) %>%
  filter(!is.na(is_filtered))

tree_tip_metadata

# tree
p_tree <-
  ggtree::`%<+%`(ggtree::ggtree(midpoint(tree), layout = "circular", size = 0.4), tree_tip_metadata) +
  ggtree::geom_tippoint(
    aes(color = filter_type, shape = filter_type, fill = filter_type, size = filter_type),
    stroke = 0.3
  ) +
  ggtree::geom_tiplab(size = 3, aes(color = filter_type)) +
  scale_color_manual(values = filter_colors_type, guide = "none") +
  scale_fill_manual(values = filter_fill_type, guide = "none") +
  scale_shape_manual(values = filter_shapes, guide = "none") +
  scale_size_manual(values = filter_sizes, guide = "none")

p_tree

saveplot(
  p_tree,
  file.path(plot_dir, "ML_tree.QC_filtered_tips"),
  layout = "large",
  mode = "editing"
)


# REads ————————————————————————————————————————
# Prepare QC plot data (exclude reference samples)
qc_reads <-
  qc.metadata.final %>%
  filter(!is.na(Cohort)) %>%
  arrange(desc(fastp_total_reads)) %>%
  mutate(CultureID = factor(CultureID, levels = CultureID))

qc_reads %>% count(Cohort)

# Plot 1: Total reads barplot (by Cohort)
p_total_reads <- ggplot(qc_reads, aes(x = CultureID, y = fastp_total_reads, fill = Cohort)) +
  geom_col() +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = cohort_colors) +
  labs(
    title = "Total Reads per CultureID",
    x = NULL
  ) +
  geom_hline(yintercept = 3000000, linetype = "dashed", color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 332275, linetype = "dashed", color = "black", linewidth = 0.2) +
  annotate("text", x = Inf, y = 3000000, label = "~200x", hjust = 0, vjust = 2, color = "black", size = 6) +
  annotate("text", x = Inf, y = 332275, label = "~21x", hjust = 0, vjust = 1, color = "black", size = 6) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 40, 5.5, 5.5),
    legend.position = "bottom"
  )


p_total_reads_ptb <- ggplot(qc_reads, aes(x = CultureID, y = fastp_total_reads, fill = PTB_category_det)) +
  geom_col() +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = ptb_colors) +
  labs(
    title = "Total Reads per CultureID",
    x = "CultureID"
  ) +
  geom_hline(yintercept = 3000000, linetype = "dashed", color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 332275, linetype = "dashed", color = "black", linewidth = 0.2) +
  annotate("text", x = Inf, y = 3000000, label = "~200x", hjust = .05, vjust = 2, color = "black", size = 6) +
  annotate("text", x = Inf, y = 332275, label = "~21x", hjust = .05, vjust = 1, color = "black", size = 6) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 40, 5.5, 5.5),
    legend.position = "bottom"
  )

p_total_reads_combined <- p_total_reads_ptb | p_total_reads

print(p_total_reads_combined)

saveplot(p_total_reads_combined, file.path(plot_dir, "fastp.total_reads_by_cohort_and_PTBcat"), layout = "large", mode = "editing")

# Plot 2: QC metric distributions
qc_long <- qc.metadata.final %>%
  filter(!is.na(Cohort)) %>%
  select(CultureID, fastp_total_reads) %>%
  pivot_longer(
    cols = -CultureID,
    names_to = "metric",
    values_to = "value"
  )

p_distribution <- ggplot(qc_long, aes(x = value)) +
  geom_histogram(bins = 50, fill = "#1f77b4", color = "white") +
  facet_wrap(~metric, scales = "free", ncol = 2) +
  geom_vline(xintercept = 332275, linetype = "dashed", color = "black", linewidth = 0.2) +
  annotate("text", x = 332275, y = Inf, label = "~21x", hjust = -0.1, vjust = 1.5, color = "black", size = 3) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Distribution of QC Metrics Across Samples",
    x = "Number of reads",
    y = "Isolates"
  )

print(p_distribution)
saveplot(p_distribution, file.path(plot_dir, "fastp.total_reads_histogram"), layout = "medium", mode = "editing")


# ———— Kraken Plots ——————————————————————————————————————————————

# lc_pct plotting table (backward-compatible with downstream plots)
lc_pct <-
  qc.metadata.final %>%
  mutate(FILTER = ifelse(is.na(FILTER), "PASS", FILTER))


lc_low_fence <-
  qc.filters %>%
  filter(software == "Kraken2", metric == "lc_pct") %>%
  pull(lower_threshold)

# Plot 1: Kraken distribution (histogram)
p_kraken_dist <-
  ggplot(lc_pct, aes(x = lc_pct)) +
  geom_histogram(aes(fill = FILTER),
    binwidth = 1, boundary = 0,
    colour = "black", alpha = 0.8
  ) +
  geom_vline(
    xintercept = lc_low_fence, linetype = "dashed",
    colour = "black", linewidth = 0.2
  ) +
  annotate("text",
    x = lc_low_fence - 0.3, y = Inf,
    label = "low fence", hjust = 1, vjust = 1.4,
    colour = "black", size = 3
  ) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Distribution of % of reads classified as L. crispatus by Kraken",
    x = "L. crispatus (% of reads)", y = "Number of colonies"
  ) +
  scale_fill_manual(values = filter_colors_type, name = "FILTER")

print(p_kraken_dist)

ggplot(lc_pct, aes(x = lc_pct)) +
  geom_density(
    binwidth = 5, boundary = 0,
    colour = "black", alpha = 0.8
  ) +
  geom_vline(
    xintercept = lc_low_fence, linetype = "dashed",
    colour = "black", linewidth = 0.2
  ) +
  annotate("text",
    x = lc_low_fence - 0.3, y = Inf,
    label = "low fence", hjust = 1, vjust = 1.4,
    colour = "black", size = 3
  ) +
  labs(
    title = "Distribution of L. crispatus classification %",
    x = "L. crispatus (% of reads)", y = "Number of samples"
  )

print(p_kraken_dist)

saveplot(p_kraken_dist, file.path(plot_dir, "Kraken.distribution_LCreads_by_QCfilters"), layout = "medium", mode = "editing")

# Plot 2: Kraken per-sample sorted dot plot (ALL FILTERS)
lc_sorted <- lc_pct %>%
  arrange(lc_pct) %>%
  mutate(rank_order = seq_len(n()))

p_kraken_boxplot <- ggplot(lc_sorted, aes(x = rank_order, y = lc_pct)) +
  geom_boxplot(fill = "#62a1cf", alpha = 0.5, width = 2) +
  geom_point(
    data = dplyr::filter(lc_sorted, FILTER == "PASS"),
    aes(colour = FILTER, size = FILTER, shape = FILTER),
    alpha = 1
  ) +
  geom_point(
    data = dplyr::filter(lc_sorted, FILTER != "PASS"),
    aes(fill = FILTER, size = FILTER, shape = FILTER),
    alpha = 1
  ) +
  geom_text(
    data = dplyr::filter(lc_sorted, FILTER != "PASS"),
    aes(x = rank_order, y = 95, label = "*"),
    inherit.aes = FALSE,
    colour = "grey40",
    size = 2
  ) +
  geom_hline(
    yintercept = lc_low_fence, linetype = "dashed",
    colour = "lightgrey", linewidth = 0.7
  ) +
  # ggrepel::geom_text_repel(
  #  data = dplyr::filter(lc_sorted, FILTER != "PASS"),
  #  aes(label = CultureID),
  #  size = 2,
  #  max.overlaps = 30,
  #  na.rm = TRUE,
  #  colour = "black"
  # ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
  scale_size_manual(values = filter_sizes, guide = "none") +
  scale_shape_manual(values = filter_shapes) +
  scale_color_manual(values = filter_colors_bw, guide = "none") +
  scale_fill_manual(values = filter_fill_bw, guide = "none") +
  labs(
    title = "L. crispatus % per sample (sorted) - ALL FILTERS",
    x = "Sample",
    y = "L. crispatus (% of reads)"
  ) +
  theme(legend.position = "none")

print(p_kraken_boxplot)

saveplot(p_kraken_boxplot,
  file.path(plot_dir, "Kraken.boxplot_perct_Lcrispatus_persample"),
  layout = "barplot",
  height_cm = 20,
  mode = "editing"
)

# Plot 3: Kraken species composition (ALL FILTERS)
sample_order <- lc_pct %>%
  arrange(desc(lc_pct)) %>%
  pull(CultureID)

#species_file <- "data/0_raw_reads_QC/3_Kraken2/Lcrispatus402.Kraken.all.summary"

species_file <- args[5]

kraken_raw <- read.table(species_file,
  sep = "\t", header = TRUE,
  quote = "", fill = TRUE, stringsAsFactors = FALSE
)
head(kraken_raw)

kraken_sp <- kraken_raw %>%
  filter(rank == "S") %>%
  mutate(name = trimws(name), percent = as.numeric(percent))

kraken_all_sp <-
  kraken_sp %>%
  mutate(name = str_replace_all(trimws(name), " ", "_")) %>%
  mutate(
    species_label = ifelse(percent < 1, "Other (<1%)", name),
    sample = factor(sample, levels = sample_order)
  ) %>%
  left_join(
    qc.metadata.final %>%
      select(CultureID, FILTER, Cohort),
    by = c("sample" = "CultureID")
  ) %>%
  mutate(FILTER = ifelse(is.na(FILTER), "PASS", FILTER)) %>%
  filter(!is.na(Cohort)) %>%
  filter(species_label != "unclassified")

kraken_all_sp %>% head()

species_freq <- kraken_all_sp %>%
  group_by(species_label) %>%
  summarise(mean_pct = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_pct) %>%
  pull(species_label)

kraken_all_sp <- kraken_all_sp %>%
  mutate(species_label = factor(species_label, levels = species_freq))

main_species <- species_freq[species_freq != "Other (<1%)"]
species_colors <- c(setNames(rep("#cccccc", 1), "Other (<1%)"))
species_colors <- c(setNames(rep("#8f8d8d", 1), "unclassified"))

if (length(main_species) > 0) {
  fallback_colors <- grDevices::rgb(runif(length(main_species)), runif(length(main_species)), runif(length(main_species)))
  matched_colors <- unname(taxa_colors[main_species])
  final_colors <- ifelse(is.na(matched_colors), fallback_colors, matched_colors)
  species_colors <- c(species_colors, setNames(final_colors, main_species))
}

species_colors

p_kraken_all_sp <- ggplot(
  kraken_all_sp,
  aes(x = sample, y = (percent), fill = (species_label))
) +
  geom_col() +
  geom_hline(yintercept = 100, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "black", linewidth = 0.2) +
  # facet_grid(. ~ FILTER, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = species_colors) +
  scale_x_discrete(limits = sample_order, drop = TRUE) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
  geom_text(
    data = dplyr::filter(kraken_all_sp, FILTER != "PASS"),
    aes(x = sample, y = 95, label = "*"),
    inherit.aes = FALSE,
    colour = "grey40",
    size = 2
  ) +
  labs(
    title = "Kraken species composition - common species + Other",
    x = "Sample", y = "% of reads", fill = "Species"
  ) +
  theme(
    legend.text = element_text(face = "italic", size = 5),
    legend.title = element_text(size = 6),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

print(p_kraken_all_sp)

saveplot(p_kraken_all_sp, file.path(plot_dir, "Kraken.composition"), layout = "medium", mode = "editing")

aligned_kraken <- cowplot::plot_grid(
  p_kraken_boxplot,
  p_kraken_all_sp,
  nrow = 1,
  rel_widths = c(0.5, 1.4),
  labels = c("A", "B"),
  align = "hv",
  axis = "tblr"
)

print(aligned_kraken)

saveplot(
  aligned_kraken,
  file.path(plot_dir, "Kraken.panel_boxplot_and_species_composition"),
  layout = "large",
  mode = "editing"
)

# Plot 4: Kraken tree with bars
# Prepare Kraken data structures for plotting (exclude reference samples)

# Prepare all tree data structures for tree plot


# Reuse master table (computed in Section 3) — rename for tree compatibility
kraken_tree_stats <- qc.metadata.final %>%
  select(
    CultureID,
    unclassified_pct,
    non_lcrispatus_species_pct = other_species_pct,
    other_genus_pct = other_l_genus_pct
  )

kraken_tree_stats

tree_species_df <- kraken_sp %>%
  mutate(
    name = str_replace_all(name, " ", "_"),
    species_label = ifelse(percent < 1, "Other (<1%)", name),
    sample = factor(sample, levels = tree$tip.label)
  ) %>%
  filter(!is.na(sample))

tree_species_df

species_freq_tree <- tree_species_df %>%
  group_by(species_label) %>%
  summarise(mean_pct = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_pct) %>%
  pull(species_label)

species_freq_tree

tree_species_df <- tree_species_df %>%
  mutate(species_label = factor(species_label, levels = species_freq_tree))


tree_cst_df <- data.frame(sample = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(
    qc.metadata.final %>%
      filter(!is.na(Cohort)) %>%
      select(CultureID, subCST) %>%
      distinct(CultureID, .keep_all = TRUE) %>%
      transmute(sample = CultureID, subCST = ifelse(is.na(subCST), "NA", subCST)),
    by = "sample"
  ) %>%
  filter(!is.na(subCST)) %>%
  mutate(
    sample = factor(sample, levels = tree$tip.label),
    panel = "CST"
  )

tree_ptb_df <- data.frame(sample = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(
    qc.metadata.final %>%
      filter(!is.na(Cohort)) %>%
      select(CultureID, PTB_category_det) %>%
      distinct(CultureID, .keep_all = TRUE) %>%
      transmute(sample = CultureID, PTB_category_det = ifelse(is.na(PTB_category_det), "NA", PTB_category_det)),
    by = "sample"
  ) %>%
  filter(!is.na(PTB_category_det)) %>%
  mutate(
    sample = factor(sample, levels = tree$tip.label),
    panel = "PTB"
  )

tree_bar_df <- data.frame(sample = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(kraken_tree_stats, by = c("sample" = "CultureID")) %>%
  mutate(
    unclassified_pct = replace(unclassified_pct, is.na(unclassified_pct), 0),
    non_lcrispatus_species_pct = replace(non_lcrispatus_species_pct, is.na(non_lcrispatus_species_pct), 0),
    other_genus_pct = replace(other_genus_pct, is.na(other_genus_pct), 0),
    sample = factor(sample, levels = tree$tip.label)
  )

p_kraken_tree_bars <-
  ggtree::`%<+%`(ggtree::ggtree(midpoint(tree), size = 0.4), tree_tip_metadata) +
  ggtree::geom_tippoint(
    aes(color = filter_type, shape = filter_type, fill = filter_type, size = filter_type),
    stroke = 0.3
  ) +
  scale_color_manual(values = filter_colors_type, guide = "none") +
  scale_fill_manual(values = filter_fill_type, guide = "none") +
  scale_shape_manual(values = filter_shapes, guide = "none") +
  scale_size_manual(values = filter_sizes, guide = "none") +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    data = tree_cst_df,
    geom = geom_tile,
    mapping = aes(y = sample, fill = subCST),
    pwidth = 0.0001,
    show.legend = FALSE,
    offset = 0.01
  ) +
  scale_fill_manual(values = c(subCST_colors, "NA" = "grey85"), guide = "none") +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    data = tree_ptb_df,
    geom = geom_tile,
    mapping = aes(y = sample, fill = PTB_category_det),
    pwidth = 0.0001,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(ptb_colors, "NA" = "grey85"), guide = "none") +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    data = tree_species_df,
    geom = geom_bar,
    mapping = aes(y = sample, x = percent, fill = species_label),
    stat = "identity",
    orientation = "y",
    pwidth = 0.45,
    axis.params = list(
      axis = "x",
      title = "Species composition (%)",
      text.size = 2,
      title.angle = 0,
      vjust = -0.5,
      title.height = 0.01
    ),
    show.legend = FALSE
  ) +
  ggtreeExtra::geom_fruit(
    data = tree_bar_df,
    geom = geom_bar,
    mapping = aes(y = sample, x = non_lcrispatus_species_pct),
    stat = "identity",
    orientation = "y",
    fill = "#8DB859",
    pwidth = 0.45,
    axis.params = list(
      axis = "x",
      title = "Non LC Species (%)",
      text.size = 2,
      title.angle = 0,
      vjust = -0.5,
      title.height = 0.01
    ),
    show.legend = FALSE
  ) +
  ggtreeExtra::geom_fruit(
    data = tree_bar_df,
    geom = geom_bar,
    mapping = aes(y = sample, x = other_genus_pct),
    stat = "identity",
    orientation = "y",
    fill = "#C97B63",
    pwidth = 0.45,
    axis.params = list(
      axis = "x",
      title = "Other Genus (%)",
      text.size = 2,
      title.angle = 0,
      vjust = -0.5,
      title.height = 0.01
    ),
    show.legend = FALSE
  ) +
  ggtreeExtra::geom_fruit(
    data = tree_bar_df,
    geom = geom_bar,
    mapping = aes(y = sample, x = unclassified_pct),
    stat = "identity",
    orientation = "y",
    fill = "grey55",
    pwidth = 0.45,
    axis.params = list(
      axis = "x",
      title = "Unclassified (%)",
      text.size = 2,
      title.angle = 0,
      vjust = -0.5,
      title.height = 0.01
    ),
    show.legend = FALSE
  ) +
  labs(y = NULL) +
  scale_fill_manual(values = species_colors) +
  guides(
    fill = guide_legend(position = "bottom")
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6)
  )

print(p_kraken_tree_bars)

saveplot(
  p_kraken_tree_bars,
  file.path(plot_dir, "Kraken.tree_composition_nonclassified"),
  layout = "large",
  mode = "editing",
  height_cm = 30
)

# ———— Mapping Plots ——————————————————————————————————————————————

# Prepare data (ALL FILTERS from qc.metadata.final, exclude reference samples)
mapping_plot <- qc.metadata.final %>%
  filter(!is.na(Cohort), !is.na(PTB_category_det)) %>%
  mutate(PTB_category_det = factor(PTB_category_det, levels = ptb_levels)) %>%
  select(CultureID, SampleID, subCST, PTB_category_det, COV_PCT, VARIANT, UNALIGNED_PCT, HET, FILTER, mapping_flags)

mapping_plot %>% head()

# Faceted ranked mapping plot (styled like p_kraken_boxplot)
mapping_plot_ranked <- mapping_plot %>%
  mutate(FILTER = ifelse(is.na(FILTER), "PASS", FILTER)) %>%
  pivot_longer(
    cols = c(COV_PCT, VARIANT, UNALIGNED_PCT, HET),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  group_by(metric) %>%
  arrange(value, .by_group = TRUE) %>%
  mutate(rank_order = row_number()) %>%
  ungroup()

p_mapping_samples_faceted <- ggplot(mapping_plot_ranked, aes(x = rank_order, y = value)) +
  geom_boxplot(fill = "#62a1cf", alpha = 0.5) +
  geom_point(
    data = dplyr::filter(mapping_plot_ranked, FILTER == "PASS"),
    aes(colour = FILTER, size = FILTER, shape = FILTER),
    alpha = 1
  ) +
  geom_point(
    data = dplyr::filter(mapping_plot_ranked, FILTER != "PASS"),
    aes(fill = FILTER, shape = FILTER, size = FILTER),
    alpha = 1,
    colour = "black"
  ) +
  #ggrepel::geom_text_repel(
  #  data = dplyr::filter(mapping_plot_ranked, FILTER != "PASS"),
  #  aes(label = CultureID),
  #  size = 2.2,
  #  max.overlaps = 30,
  #  na.rm = TRUE,
  #  colour = "black"
  #) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_size_manual(values = filter_sizes, guide = "none") +
  scale_shape_manual(values = filter_shapes, guide = "none") +
  scale_color_manual(values = filter_colors_type, guide = "none") +
  scale_fill_manual(
    values = filter_fill_type,
    guide = guide_legend(override.aes = list(shape = 23, size = 3, colour = "black"), ncol = 2)
  ) +
  labs(
    title = "Mapping metrics per sample (sorted) - ALL FILTERS",
    x = "Sample rank (low to high)",
    y = "Metric value"
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical"
  )

print(p_mapping_samples_faceted)

saveplot(p_mapping_samples_faceted,
  file.path(plot_dir, "mapping.boxplot_faceted_by_metric"),
  layout = "large", mode = "editing"
)


# Alternative faceted mapping plot with FILTER-type colors
p_mapping_samples_faceted_type <- ggplot(mapping_plot_ranked, aes(x = rank_order, y = value)) +
  geom_boxplot(fill = "#62a1cf", alpha = 0.5) +
  geom_point(
    data = dplyr::filter(mapping_plot_ranked, FILTER == "PASS"),
    aes(colour = FILTER, size = FILTER, shape = FILTER),
    alpha = 1
  ) +
  geom_point(
    data = dplyr::filter(mapping_plot_ranked, FILTER != "PASS"),
    aes(fill = FILTER, colour = FILTER, size = FILTER, shape = FILTER),
    alpha = 1
  ) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(mapping_plot_ranked, FILTER != "PASS"),
    aes(label = CultureID),
    size = 2.2,
    max.overlaps = 30,
    na.rm = TRUE,
    colour = "black"
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_size_manual(values = filter_sizes, guide = "none") +
  scale_shape_manual(values = filter_shapes) +
  scale_color_manual(values = filter_colors_type) +
  scale_fill_manual(values = filter_fill_type, guide = "none") +
  labs(
    title = "Mapping metrics per sample (sorted) - FILTER type colors",
    x = "Sample rank (low to high)",
    y = "Metric value"
  ) +
  theme(legend.position = "bottom")

print(p_mapping_samples_faceted_type)

saveplot(p_mapping_samples_faceted_type,
  file.path(plot_dir, "mapping.boxplot_faceted_by_metric_byfilter2"),
  layout = "large", mode = "editing"
)

# Create mapping plots
# Coverage, Variants, HET by PTB Category and subCST (grouped)
p_cov_ptb_subcst_grouped <- ggplot(mapping_plot, aes(x = PTB_category_det, y = COV_PCT, color = subCST)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge()) +
  facet_wrap(~FILTER, scales = "free_x") +
  labs(
    title = "Coverage by PTB Category and subCST (All FILTER groups)",
    x = NULL,
    y = "Coverage"
  ) +
  scale_color_manual(values = subCST_colors, guide = "none")

print(p_cov_ptb_subcst_grouped)

p_variant_ptb_subcst_grouped <- ggplot(mapping_plot, aes(x = PTB_category_det, y = VARIANT, color = subCST)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge()) +
  facet_wrap(~FILTER, scales = "free_x") +
  labs(
    title = "Variants by PTB Category and subCST (All FILTER groups)",
    y = "Variant Count",
    x = NULL
  ) +
  scale_color_manual(values = subCST_colors) +
  theme(legend.position = "right")

print(p_variant_ptb_subcst_grouped)

p_het_ptb_subcst_grouped <- ggplot(mapping_plot, aes(x = PTB_category_det, y = HET, color = subCST)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge()) +
  facet_wrap(~FILTER, scales = "free_x") +
  labs(
    title = "Heterozygotes by PTB Category and subCST (All FILTER groups)",
    y = "HET Count",
    x = NULL
  ) +
  scale_color_manual(values = subCST_colors, guide = "none")

print(p_het_ptb_subcst_grouped)

p_cov_variant_het_combined_ptb_subcst_grouped <- (p_cov_ptb_subcst_grouped / p_variant_ptb_subcst_grouped / p_het_ptb_subcst_grouped)

p_cov_variant_het_combined_ptb_subcst_grouped

saveplot(p_cov_variant_het_combined_ptb_subcst_grouped, file.path(plot_dir, "mapping.stats_cov_variant_het_by_ptb_subcst"), layout = "large", mode = "editing", tag = "filtered")


# Coverage, Variants, HET by subCST (with PTB category overlay)
p_cov_subcst_ptb <- ggplot(mapping_plot, aes(x = subCST, y = COV_PCT)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = PTB_category_det), position = position_jitterdodge()) +
  facet_wrap(~FILTER, scales = "free_x") +
  labs(
    title = "Coverage by subCST (All FILTER groups)",
    x = NULL,
    y = NULL
  ) +
  scale_color_manual(values = ptb_colors, guide = "none")

p_variant_subcst_ptb <- ggplot(mapping_plot, aes(x = subCST, y = VARIANT)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = PTB_category_det), position = position_jitterdodge()) +
  facet_wrap(~FILTER, scales = "free_x") +
  labs(
    title = "Variants by subCST (All FILTER groups)",
    y = NULL,
    x = NULL
  ) +
  scale_color_manual(values = ptb_colors, guide = "none")

p_het_subcst_ptb <- ggplot(mapping_plot, aes(x = subCST, y = HET)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = PTB_category_det), position = position_jitterdodge()) +
  facet_wrap(~FILTER, scales = "free_x") +
  labs(
    title = "Heterozygotes by subCST (All FILTER groups)",
    y = NULL,
    x = NULL
  ) +
  scale_color_manual(values = ptb_colors) +
  theme(legend.position = "bottom")

p_cov_variant_het_combined_subcst_ptb <- (p_cov_subcst_ptb / p_variant_subcst_ptb / p_het_subcst_ptb)

p_cov_variant_het_combined_subcst_ptb

# Coverage, Variants, HET by PTB Category (with subCST overlay)
p_cov_ptb_subcst_overlay <- ggplot(mapping_plot, aes(x = PTB_category_det, y = COV_PCT)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = subCST), position = position_jitterdodge()) +
  facet_wrap(~FILTER, scales = "free_x") +
  labs(
    title = "Coverage by PTB Category and subCST (All FILTER groups)",
    x = NULL,
    y = "Coverage"
  ) +
  scale_color_manual(values = subCST_colors, guide = "none")

p_variant_ptb_subcst_overlay <- ggplot(mapping_plot, aes(x = PTB_category_det, y = VARIANT)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = subCST), position = position_jitterdodge()) +
  facet_wrap(~FILTER, scales = "free_x") +
  labs(
    title = "Variants by PTB Category and subCST (All FILTER groups)",
    y = "Variant Count",
    x = NULL
  ) +
  scale_color_manual(values = subCST_colors, guide = "none")

p_het_ptb_subcst_overlay <- ggplot(mapping_plot, aes(x = PTB_category_det, y = HET)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = subCST), position = position_jitterdodge()) +
  facet_wrap(~FILTER, scales = "free_x") +
  labs(
    title = "Heterozygotes by PTB Category and subCST (All FILTER groups)",
    y = "HET Count",
    x = NULL
  ) +
  scale_color_manual(values = subCST_colors) +
  theme(legend.position = "bottom")

p_cov_variant_het_combined_ptb_subcst_overlay <- (p_cov_ptb_subcst_overlay / p_variant_ptb_subcst_overlay / p_het_ptb_subcst_overlay)

p_cov_variant_het_combined_all <- p_cov_variant_het_combined_ptb_subcst_overlay | p_cov_variant_het_combined_subcst_ptb

print(p_cov_variant_het_combined_all)

saveplot(p_cov_variant_het_combined_all, file.path(plot_dir, "mapping.cov_variant_het.CST_nd_PTB"), layout = "large", mode = "editing", tag = "filtered")

# ———— Assembly Plots ——————————————————————————————————————————————
assembly_metrics <- c(
  "Contigs",
  "Total_Length",
  "GC_Content",
  "Predicted_Genes_Unique",
  "Predicted_Genes_0bp_Complete",
  "Predicted_Genes_0bp_Partial"
)

# Prepare assembly data (ALL FILTERS from qc.metadata.final, exclude reference samples)

assembly_plot_data <-
  qc.metadata.final %>%
  filter(!is.na(Cohort)) %>%
  select(CultureID, subCST, PTB_category_det, FILTER, all_of(assembly_metrics), assembly_flags) %>%
  mutate(FILTER = ifelse(is.na(FILTER), "PASS", FILTER)) %>%
  pivot_longer(
    cols = all_of(assembly_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  # left_join(assembly_outlier_samples, by = c("CultureID"))  %>%
  mutate(is_outlier = ifelse(is.na(assembly_flags), FALSE, TRUE))

assembly_plot_data %>% head()

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
p_quast_ranked <- 
  ggplot(quast_rank_data, aes(x = sample_order, y = value)) +
  geom_boxplot(fill = "#62a1cf", alpha = 0.5) +
  geom_point(
    data = dplyr::filter(quast_rank_data, FILTER == "PASS"),
    aes(colour = FILTER, size = FILTER, shape = FILTER)
  ) +
  geom_point(
    data = dplyr::filter(quast_rank_data, FILTER != "PASS"),
    aes(y = value, fill = FILTER, colour = FILTER, shape = FILTER)
  ) +
  scale_size_manual(values = filter_sizes, guide = "none") +
  scale_shape_manual(values = filter_shapes) +
  scale_color_manual(values = filter_colors_type) +
  scale_fill_manual(values = filter_fill_type, guide = "none") +
  # ggrepel::geom_text_repel(
  #  data = dplyr::filter(quast_rank_data, !is.na(assembly_flags)),
  #  aes(label = CultureID),
  #  size = 2.2,
  #  max.overlaps = 60,
  #  na.rm = TRUE,
  #  colour = "black"
  # ) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  labs(
    title = "QUAST Metrics per Sample (Ordered High to Low) - FILTER type colors",
    x = "Samples ordered high to low",
    y = "Metric value"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

print(p_quast_ranked)

saveplot(p_quast_ranked, file.path(plot_dir, "quast.metrics_ranked_per_sample_v1"), layout = "large", mode = "editing", tag = "filtered")


# Only clean ———————————————————————

# Create assembly plots
# Ranked QUAST metrics (ordered high to low, colored and shaped by FILTER - no faceting)
p_quast_ranked <- ggplot(quast_rank_data, aes(x = sample_order, y = value)) +
  geom_point(
    data = dplyr::filter(quast_rank_data, FILTER == "PASS"),
    aes(colour = FILTER, size = FILTER, shape = FILTER, alpha = is_outlier)
  ) +
  scale_alpha_manual(values = c("FALSE" = 0.6, "TRUE" = 1), guide = "none") +
  scale_size_manual(values = filter_sizes, guide = "none") +
  scale_shape_manual(values = filter_shapes) +
  scale_color_manual(values = filter_colors_type) +
  scale_fill_manual(values = filter_fill_type, guide = "none") +
  geom_boxplot(data = quast_box_data, aes(x = box_x, y = value), inherit.aes = FALSE, width = 30, outlier.shape = NA, color = "grey50", alpha = 0.2) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  labs(
    title = "QUAST Metrics per Sample (Ordered High to Low) - FILTER type colors",
    x = "Samples ordered high to low",
    y = "Metric value"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

print(p_quast_ranked)

saveplot(p_quast_ranked, file.path(plot_dir, "quast.metrics_ranked_per_sample_v2"), layout = "large", mode = "editing", tag = "filtered")
