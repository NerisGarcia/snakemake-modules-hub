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

initialize_saveplot_theme( font="helvetica_neue")
source("code/design.R")

args <- commandArgs(trailingOnly = TRUE)

qc_stats_file <- "data/20260529.09.Lcrispatus417.metadata.QC_stats_combined.ods"

output_culture_list <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "data/Lcrispatus388.CultureID.txt"
plot_dir <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "figures/2_Lcrispatus_v2.388colonies/"

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
# =============================================================================
# SECTION 6: ALL PLOTTING (using qc.metadata.final with complete FILTER info)
# =============================================================================

qc.metadata.final <- 
read_ods(qc_stats_file)  %>% 
filter(FILTER== "PASS") 



qc.metadata.final  %>%  select(CultureID)  %>% count 

write.csv(file = output_culture_list, row.names = FALSE, quote = FALSE)


# by samples
label_df <- 
qc.metadata.final %>% 
group_by(SampleID) %>%
filter(row_number() == 1) %>%
ungroup() %>%
  count(Cohort, PTB_category_det) 


p_cst_samp <- 
qc.metadata.final %>% 
group_by(SampleID) %>%
filter(row_number() == 1) %>%
ungroup() %>%
  count(Cohort, PTB_category_det, subCST) %>%
  ggplot(aes(x = factor(PTB_category_det, levels = ptb_levels), y = n, fill = subCST)) +
  geom_col(colour = "black", position = "stack", width = 0.7) +
geom_text(
  data = label_df,
  aes(x = factor(PTB_category_det, levels = ptb_levels), y = n + 5, label = n),
  inherit.aes = FALSE
)+  scale_fill_manual(values = subCST_colors, name = "subCST")  +
  #scale_y_continuous(labels = comma, limits = c(0, 100)) +
  facet_wrap(~Cohort, ncol=1,  strip.position = "left") +
  theme(axis.text.x  = element_text(angle = 35, hjust = 1),
        strip.background = element_rect(fill = "grey92"), 
        legend.position = "none"
      ) +
  labs(title = "Samples",
      x = NULL, y = "Number of samples",
      caption = paste0("v2 - n=88 samples")) 

p_cst_samp


# by colonies

# by samples
label_df_col <- 
qc.metadata.final %>% 
  count(Cohort, PTB_category_det) 


p_cohort <-
qc.metadata.final  %>%
  count(Cohort, PTB_category_det) %>%
  ggplot(aes(x = PTB_category_det, y = n, fill = Cohort)) +
  geom_col(colour = "black", width = 0.7) +
  geom_text(aes(label = n),  size = 3.5) +
  scale_fill_manual(values = cohort_colors, name = NULL) +
  labs(title = "Colonies",
      x = NULL, y = "Number of colonies",
      caption = paste0("v2 dataset colonies")) +
  theme(axis.text.x  = element_text(angle = 35, hjust = 1),
        legend.position = "none",
        strip.background = element_rect(fill = "grey92"))

p_cst <- 
qc.metadata.final %>% 
  count(Cohort, PTB_category_det, subCST) %>%
  ggplot(aes(x = factor(PTB_category_det, levels = ptb_levels), y = n, fill = subCST)) +
  geom_col(colour = "black", position = "stack", width = 0.7) +
#  geom_text(aes(label = n)) +
  scale_fill_manual(values = subCST_colors, name = "subCST")  +
  #scale_y_continuous(labels = comma, limits = c(0, 100)) +
  geom_text(
  data = label_df_col,
  aes(x = factor(PTB_category_det, levels = ptb_levels), y = n + 5, label = n),
  inherit.aes = FALSE
)+ 
  facet_wrap(~Cohort, ncol=1) +
  theme(axis.text.x  = element_text(angle = 35, hjust = 1),
        legend.position = "right",
         strip.text = element_blank(),       # Removes the text labels
    strip.background = element_blank()  # Removes the background boxes

     ) +
  labs(title = "Colonies",
      x = NULL, y = "Number of colonies",
      caption = paste0("v2 - n=388 colonies")) 

p_cst

plot <- p_cst_samp + p_cst +
plot_layout( ncol = 2)

plot

saveplot(plot, 
filename = paste0(plot_dir, "colony_counts_per_subCST_and_cohort"), 
layout = "small", 
mode="editing")
