
conda create --name ncbi_download -y 

conda activate ncbi_download

conda install -c conda-forge ncbi-datasets-cli



#install sra-tools
cd ~/Desktop/SOFTWARE_NGG
curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz


#server
curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-alma_linux64.tar.gz
tar -vxzf sratoolkit.tar.gz

#append to conda path

#echo 'export PATH=~/Desktop/SOFTWARE_NGG/sratoolkit.3.0.1-mac64/bin:$PATH' >> ~/.bash_profile

ln -s /home/ngarcia-gonzalez/software/sratoolkit.current-alma_linux64/bin/* $CONDA_PREFIX/bin/


# install entrez 
conda install bioconda::entrez-direct