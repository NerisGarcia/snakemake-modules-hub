# VISTA

VISTA is an open-source classifier that assigns metagenomic subspecies (mgSs) and metagenomic community state types (mgCSTs) from VIRGO2-derived vaginal metagenomic gene profiles.

to use it, first you have to install VISTA manually as there is no conda or pip.

``` bash
git clone https://github.com/JHolm-Lab/VISTA.git
cd VISTA

#add VISTA script path to vista snakemake config file
```

Run snakemake

``` bash
snakemake --snakefile /Users/ngarcia-gonzalez/Desktop/RESEARCH/snakemake-modules-hub/modules/10_metagenomic_taxonomy/vista/Snakefile --configfile /Users/ngarcia-gonzalez/Desktop/RESEARCH/snakemake-modules-hub/modules/10_metagenomic_taxonomy/vista/config.yaml --use-conda  --conda-frontend conda --cores 8 -n
```