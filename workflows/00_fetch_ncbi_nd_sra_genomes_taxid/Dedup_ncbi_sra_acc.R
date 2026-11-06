library(readxl)
library(dplyr)
library(stringr)
library(dplyr)
library(tidyr)


#  Database genomes ________________________________________________________
# ncbi and sra genomes dowloaded the 5 of november of 2026.
# MAGS are excluded

ncbi_genomes_raw <- read.csv("data/0_input_data/0_genomes_05Nov2026/1578_ncbi_genomes_summary.tsv", sep = "\t", quote = "", check.names = TRUE) 

ncbi_genomes_raw %>% 
colnames()

ncbi_genomes_raw %>%
select( Assembly.BioSample.Accession,  Organism.Name  , Assembly.BioSample.Description.Organism.Taxonomic.ID) %>% 
distinct()  %>% 
count(Organism.Name)   %>%
arrange(desc(n)) %>%
  head(20)

#                                 Organism.Name    n
#1                      Lactobacillus crispatus 1159 *
#2  Lactobacillus delbrueckii subsp. bulgaricus  481
#3                        Lactobacillus gasseri  334 *
#4                     Lactobacillus helveticus  236
#5      Lactobacillus delbrueckii subsp. lactis  199
#6                          Lactobacillus iners  175 *
#7                       Lactobacillus jensenii  161 *
#8                    Lactobacillus acidophilus  148
#9                    Lactobacillus paragasseri  125 *
#10                     Lactobacillus johnsonii  121
#11                   Lactobacillus delbrueckii  106
#12                    Lactobacillus amylovorus   96
#13                      Lactobacillus mulieris   59 *
#14 Lactobacillus delbrueckii subsp. allosunkii   48
#15                   Lactobacillus taiwanensis   31
#16                    Lactobacillus gallinarum   13
#17                 Lactobacillus acetotolerans   11
#18              Lactobacillus helsingborgensis   11
#19                  Lactobacillus intestinalis   11
#20                          Lactobacillus apis   10

#extract unqieu biosamples and 1 genome per biosample
# assemblies are duplicated (GCF + GCA) we only keep fnal GCF version or GCA. 
# for each biosample, select the GCF (preferred) with the largest version;
# if no GCF is available for that biosample, select the GCA with the largest version
ncbi_genomes_dedup <- ncbi_genomes_raw %>%
  select(Assembly.Accession, Assembly.BioSample.Accession, Organism.Name, Assembly.BioProject.Accession, Assembly.BioSample.Description.Organism.Taxonomic.ID, Assembly.Sequencing.Tech) %>%
  distinct() %>%
  mutate(
    Assembly_Prefix  = str_extract(Assembly.Accession, "^GC[AF]"),
    Assembly_Number  = as.numeric(str_extract(Assembly.Accession, "(?<=GC[AF]_)[0-9]+")),
    Assembly_Version = as.numeric(str_extract(Assembly.Accession, "(?<=\\.)[0-9]+$"))
  ) %>%
  group_by(Assembly.BioSample.Accession) %>%
  arrange(
    desc(Assembly_Prefix == "GCF"),
    desc(Assembly_Version),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

ncbi_genomes_dedup %>%
  head(20)

ncbi_genomes_dedup  %>%  colnames()

# SRA genomes _____________________________________________________________
# Downloaded the 5 of november of 2026.
# Only RANDOM WGS Illumina paired-end reads were considered (filtered)

sra_genomes <- read.csv("data/0_input_data/0_genomes_05Nov2026/1578_sra_runexperiments_filtered.metadata.tsv", sep = ",", quote = "", check.names = TRUE)

sra_genomes %>% head()

sra_genomes %>% colnames()

sra_genomes %>%  # corrected variable name
distinct()  %>% 
count(TaxID)   %>%
arrange(desc(n)) %>%
  head(20)
#     TaxID   n
#1    47770 387 *crispatus
#2     1596 323
#3    33959 276
#4       NA 212 ?
#5     1579 206
#6   109790 144
#7     1587 133
#8   324831 129
#9     1584  69
#10 2107999  65
#11    1604  49
#12    1591  47
#13  655183  42
#14  147802  36 *iners
#15    1585  34
#16 2508708  31
#17    1578  22
#18    1600  14
#19   52242  13
#20  151781  12


#keep only biosamples not in ncbi_genomes_dedup
sra_metadata_filtered <- 
sra_genomes %>%
select(Run, Experiment, BioProject, BioSample, TaxID)  %>% 
  filter(!BioSample %in% ncbi_genomes_dedup$Assembly.BioSample.Accession) 


sra_metadata_filtered  %>% 
distinct()  %>% 
count(TaxID)   %>%
arrange(desc(n)) %>%
  head(20)
#     TaxID   n
#1    47770 346 *crispatus
#2     1596 285 *gasseri
#3    33959 250
#4       NA 212
#5     1579 190
#6   324831 128
#7     1587 122
#8   109790 106
#9     1591  47
#10    1584  46
#11  655183  42
#12 2107999  32
#13    1585  30
#14  147802  26 *iners
#15    1578  22
#16 2508708  21
#17    1604  19
#18  151781   8
#19  508451   4
#20  883160   4

# merge db genomes: ncbi + sra
# select biosample, bioproject, and accession (Assembly.Accession for ncbi, Run for sra)
# also keep track of the source (ncbi or sra)
# and taxid

ncbi_selected <- ncbi_genomes_dedup %>%
  select(Assembly.BioSample.Accession, Assembly.BioProject.Accession, Assembly.Accession, Assembly.BioSample.Description.Organism.Taxonomic.ID)

sra_selected <- sra_metadata_filtered %>%
  select(BioSample, BioProject, Run, TaxID)

merged_genomes <- bind_rows(
  ncbi_selected %>% rename(BioSample = Assembly.BioSample.Accession,
                           BioProject = Assembly.BioProject.Accession,
                           Accession = Assembly.Accession,
                           TaxID = Assembly.BioSample.Description.Organism.Taxonomic.ID) %>% mutate(Source = "ncbi"),
  sra_selected %>% rename(Accession = Run) %>% mutate(Source = "sra") %>% mutate(TaxID = TaxID)
)

merged_genomes  %>%  head()



merged_genomes  %>% 
count(TaxID, Source)  %>% 
arrange(desc(n)) %>%
pivot_wider(names_from = Source, values_from = n)  %>% 
#add total column
mutate(Total = rowSums(across(c(ncbi, sra)), na.rm = TRUE)) %>%
head(20)

#     TaxID  ncbi   sra Total
#     <int> <int> <int> <dbl>
# 1   47770  1159   346  1505  *crispatus
# 2    1585   481    30   511
# 3    1596   337   285   622 *gasseri
# 4   33959   121   250   371
# 5    1587   236   122   358
# 6      NA    NA   212   212
# 7   29397   199     3   202
# 8    1579   148   190   338
# 9  147802   175    26   201 *iners
#10  109790   161   106   267 *jensenii
#11  324831     5   128   133
#12 2107999   121    32   153 *paragasseri
#13    1584   106    46   152
#14    1604    95    19   114
#15 2508708    59    21    80 *mulieris
#16 1050107    48    NA    48
#17    1591    NA    47    47
#18  655183    NA    42    42
#19  508451    31     4    35
#20    1578    NA    22    22




#export merged_genomes in two table
write.csv(merged_genomes %>% filter(Source == "ncbi")  %>%  select(Accession), "data/0_input_data/0_genomes_05Nov2026/1578.ncbi_acc_todown.csv", row.names = FALSE, quote = FALSE)

write.csv(merged_genomes %>% filter(Source == "sra") %>% select(Accession), "data/0_input_data/0_genomes_05Nov2026/1578.sra_acc_todown.csv", row.names = FALSE, quote = FALSE)


