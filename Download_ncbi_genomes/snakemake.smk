configfile: "Download_ncbi_genomes/config.yml"

TAXON = config.get("taxon")
DATASET_DIR = config.get("dataset_out_dir")



DATASET_NAME = f"{TAXON}_ncbi_genomes"
ZIPFILE = f"{DATASET_DIR}/{DATASET_NAME}.zip"
SUMMARYFILE = f"{DATASET_DIR}/{DATASET_NAME}_summary.tsv"
GENOMES_DIR = f"{DATASET_DIR}/genomes"


DONEFILE = f"{DATASET_DIR}/download.log"


from datetime import date
DATE = date.today().isoformat()  #todays date in YYYY-MM-DD format

rule all:
	input:
		DONEFILE,
		SUMMARYFILE,
		directory(GENOMES_DIR)


rule download_lactobacillus_summary:
	output:
		SUMMARYFILE
	params:
		taxon=TAXON,
		dir=DATASET_DIR
	conda:
		"ncbi_download"
	shell:
		"""
		mkdir -p {params.dir}
		datasets summary genome taxon {params.taxon} --mag exclude --as-json-lines | dataformat tsv genome > {output}

		"""

rule download_lactobacillus_genomes_dehydrated:
	output:
		ZIPFILE
	params:
		taxon=TAXON,
		dir=DATASET_DIR
	conda:
		"ncbi_download"
	shell:
		"""
		mkdir -p {params.dir}
		datasets download genome taxon "{params.taxon}" --mag exclude --dehydrated --filename {output}
		"""


rule unzip_lactobacillus_dataset:
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

rule rehydrate_lactobacillus_dataset:
	input:
		f"{DATASET_DIR}/{DATASET_NAME}"		
	output:
		DONEFILE
	params:
		taxon=TAXON,
		date=DATE
	conda:
		"ncbi_download"
	shell:
		"""
		datasets rehydrate --directory {input}
		echo "Downloaded {params.taxon} genomes and metadata on {params.date}" > {output}
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