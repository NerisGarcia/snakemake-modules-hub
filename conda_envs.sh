
# [shovill] Need spades version >= 3.14 and < 4.0
conda search -c bioconda spades


conda create -n genomics_env

conda activate genomics_env

conda create -n genomics_env -c bioconda fastp shovill==1.1.0 bakta busco 

conda activate genomics_env





checkm2 database --download --path /Users/ngarcia-gonzalez/Desktop/SOFTWARE_NGG/DATABASES/checkm2

conda create -n checkm2_env python=3.7

conda install -c defaults -c conda-forge -c bioconda checkm2

conda activate checkm2_env

 conda deactivate
 conda remove -n checkm2_env --all




rsync -avhP \
  --include='*/' \
  --include='*.gff3' \
  --include='*.faa' \
  --include='*.fna' \
  --exclude='*' \
  /Users/ngarcia-gonzalez/Desktop/RESEARCH/Projects/01_PTB_GWAS/data/0_input_data/3_Annotation/ \
  ngarcia-gonzalez@ravellab.igs.umaryland.edu:/local/scratch/neris/01_PTB_GWAS/data/3_Annotation/



  rsync -avhP  ngarcia-gonzalez@ravellab.igs.umaryland.edu:/local/projects-t3/PGATES/Nova2/01.RawData/Pilar_2* data/test/metagenomics/VIRGO2/raw_data/




snakemake --snakefile workflows/10_fastq2CST_VIRGO2/Snakefile --cores 4 --configfile /Users/ngarcia-gonzalez/Desktop/RESEARCH/snakemake-modules-hub/workflows/10_fastq2CST_VIRGO2/config.yml --use-conda --rerun-triggers mtime

conda create -n virgo2_env

conda activate virgo2_env

conda install -c bioconda fastp

#kraken2 test db

mkdir -p data/dbs/kraken2_test_db

kraken2-build --download-taxonomy --db data/dbs/kraken2_test_db --use-ftp
kraken2-build --download-library bacteria --db data/dbs/kraken2_test_db --use-ftp
kraken2-build --build --db data/dbs/kraken2_test_db --use-ftp


kraken2-build --download-taxonomy --db micro_db --use-ftp


mkdir -p data/test_db/library

cat > data/test_db/library/test.fa << EOF
>taxid|562|ecoli
ATCGATCGATCGATCGATCGATCGATCG
EOF


kraken2-build --add-to-library data/test_db/libraÇry/test.fa --db data/test_db 
kraken2-build --build --db data/test_db


# virgo -
# remember to install from zenodo, github doesnt work
conda activate virgo2_env

ln -s /Users/ngarcia-gonzalez/Desktop/RESEARCH/snakemake-modules-hub/VIRGO2/VIRGO2.py $CONDA_PREFIX/bin/VIRGO2.py
chmod +x /Users/ngarcia-gonzalez/Desktop/RESEARCH/snakemake-modules-hub/VIRGO2/VIRGO2.py