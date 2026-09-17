
# [shovill] Need spades version >= 3.14 and < 4.0
conda search -c bioconda spades


conda create -n genomics_env

conda activate genomics_env

conda install -c bioconda fastp shovill==1.1.0 bakta busco 

conda activate genomics_env





checkm2 database --download --path /Users/ngarcia-gonzalez/Desktop/SOFTWARE_NGG/DATABASES/checkm2

conda create -n checkm2_env python=3.7

conda install -c defaults -c conda-forge -c bioconda checkm2

conda activate checkm2_env

 conda deactivate
 conda remove -n checkm2_env --all