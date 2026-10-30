# Valencia v1.2

- v1.2 is an adaptation for VIRGO2-style **relative-abundance** data and taxon naming, while the older script is for count tables. The `ref_ave` and `ref_std` values used in the optional plot are the same in both.

- `Valencia.py` is the simpler version. It expects `sampleID`, `read_count`, and taxon read counts, then divides counts by `read_count` before scoring.

- `Valencia_v1.2.py` adds a `--metagenome` option for VIRGO2 MG/MT output: `sampleID` plus relative-abundance columns. Without that flag, it follows the count-normalization path.

- `Valencia_v1.2.py`**Taxon harmonization:** It combines *L. mulieris* with *L. jensenii*, *L. paragasseri* with *L. gasseri*, and groups several genera together; it also renames the Gardnerella aggregate. This is intended to align input taxa with the CST reference features (`combine_taxa`).

# Run snakemake rule

``` bash
conda activate snakemake
# first do a dry run
snakemake --snakefile /Users/ngarcia-gonzalez/Desktop/RESEARCH/snakemake-modules-hub/modules/10_metagenomic_taxonomy/virgo2valencia/Snakefile --configfile /Users/ngarcia-gonzalez/Desktop/RESEARCH/snakemake-modules-hub/modules/10_metagenomic_taxonomy/virgo2valencia/config.yaml -n

# if everything went ok, run
```

# Taxa names adaptation from VIRGO2 to VALENCIA

To use Virgo2 results from metagenomis in Valencia, there is some adjustment we have to make on taxonomy. As Valencia was based on 16S, the centroids were made for genus, not species level, whith only the exception of lactobacillus species. So taxa should be accumualted in genus as follows. Also some Lacto and Gardenlla need reassigments due to lastest taxonomic changes not cntemplated in VALENCIA.

For lactobacillus: Keep all lactobacillus as their species excep for: Lactobacillus mulieris has to be combinend with jensenii as L.jensenii Lactobacillus paragasseri jhas to be combines with gasseri as gasseri

For ther genuses: #combining together species where genus level typifies the CST combined_data\['g_Bifidobacterium'\] = combined_data\[\[col for col in combined_data.columns if col.startswith('Bifidobacterium')\]\].sum(axis=1) combined_data\['g_Streptococcus'\] = combined_data\[\[col for col in combined_data.columns if col.startswith('Streptococcus')\]\].sum(axis=1) combined_data\['g_Staphylococus'\] = combined_data\[\[col for col in combined_data.columns if col.startswith('Staphylococcus')\]\].sum(axis=1) combined_data\['g_Enterococcus'\] = combined_data\[\[col for col in combined_data.columns if col.startswith('Enterococcus')\]\].sum(axis=1)

Gardenella is specieal, cause VALENCIA need all gardenella as G.vaginalis. combined_data\['g_Gardnerella'\] = combined_data\[\[col for col in combined_data.columns if col.startswith('Gardnerella')\]\].sum(axis=1) combined_data=combined_data.rename(columns={'g_Gardnerella':'Gardnerella_vaginalis'})

Candidatus_Lachnocurva_vaginae has to be changed Ca_Lachnocurva_vaginae

# CENTROIDS

is works with VALENCIA2_CST_centroids_19Aug2024.csv to download

``` bash
wget https://github.com/ravel-lab/speciateIT/blob/master/VALENCIA2_CST_centroids_19Aug2024.csv
```