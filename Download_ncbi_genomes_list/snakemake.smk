configfile: "Download_ncbi_genomes_list/config.yml"

ACCESSIONS_FILE = config.get("accessions_file")
DATASET_DIR = config.get("dataset_out_dir")
DATASET_NAME = config.get("dataset_name")

ZIPFILE = f"{DATASET_DIR}/{DATASET_NAME}.zip"
GENOMES_DIR = f"{DATASET_DIR}/genomes"


DONEFILE = f"{DATASET_DIR}/download.log"


from datetime import date
DATE = date.today().isoformat()  #todays date in YYYY-MM-DD format

rule all:
	input:
		DONEFILE,
		directory(GENOMES_DIR)


rule download_genomes_dehydrated:
	input:
		ACCESSIONS_FILE
	output:
		ZIPFILE
	params:
		dir=DATASET_DIR
	conda:
		"ncbi_download"
	shell:
		"""
		mkdir -p {params.dir}
		datasets download genome accession --inputfile {input} --dehydrated --filename {output}
		"""


rule unzip_dataset:
	input:
		ZIPFILE
	output:
		directory(f"{DATASET_DIR}/{DATASET_NAME}")
	conda:
		"ncbi_download"
	shell:
		"""
		rm -rf {output}
		mkdir -p {output}
		unzip -o {input} -d {output}

		"""

rule rehydrate_dataset:
	input:
		f"{DATASET_DIR}/{DATASET_NAME}"		
	output:
		DONEFILE
	params:
		date=DATE
	conda:
		"ncbi_download"
	shell:
		"""
		datasets rehydrate --directory {input}
		echo "Downloaded genomes and metadata from accession list on {params.date}" > {output}
		"""

rule fix_folders:
	input:
		DONEFILE,
		DATASET_ROOT=f"{DATASET_DIR}/{DATASET_NAME}"
	output:
		directory(GENOMES_DIR)
	shell:
		"""
		mkdir -p {output}
		
		find {input.DATASET_ROOT} \
			-type f -name "*.fna" \
			! -path "{output}/*" \
			-exec mv -f {{}} {output}/ \;
		"""
