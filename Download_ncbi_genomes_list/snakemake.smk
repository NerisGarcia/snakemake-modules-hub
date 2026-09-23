configfile: "Download_ncbi_genomes_list/config.yml"

ACCESSIONS_FILE = config.get("accessions_file")
DATASET_DIR = config.get("dataset_out_dir")
DATASET_NAME = config.get("dataset_name")

ZIPFILE = f"{DATASET_DIR}/{DATASET_NAME}.zip"
GENOMES_DIR = f"{DATASET_DIR}/{DATASET_NAME}_fna"
CLEAN_FOLDER_DONEFILE = f"{GENOMES_DIR}/.folder_cleaned"
CLEAN_DONEFILE = f"{GENOMES_DIR}/.cleaned"


DONEFILE = f"{DATASET_DIR}/{DATASET_NAME}_download.log"


from datetime import date
DATE = date.today().isoformat()  #todays date in YYYY-MM-DD format

rule all:
	input:
		DONEFILE,
		CLEAN_DONEFILE


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

rule clean_genomes_folder:
	input:
		DONEFILE
	output:
		CLEAN_FOLDER_DONEFILE
	params:
		dir=f"{DATASET_DIR}/{DATASET_NAME}"
	shell:
		"""
		mkdir -p {GENOMES_DIR}

		if [ -d "{DATASET_DIR}/{DATASET_NAME}" ]; then
			find {DATASET_DIR}/{DATASET_NAME} -type f -name "*.fna" -print0 | while IFS= read -r -d '' f; do
				mv -f "$f" "{GENOMES_DIR}/$(basename "$f")"
			done
		fi

		rm -rf {params.dir}
		touch {output}
		"""


rule clean_fna_names:
	input:
		CLEAN_FOLDER_DONEFILE
	output:
		CLEAN_DONEFILE
	shell:
		"""
		rename_fna() {{
			local f="$1"
			base=$(basename "$f")
			accession=$(printf '%s\n' "$base" | sed -E 's/^(GC[AF]_[0-9]+\.[0-9]+).*/\\1/')
			target="{GENOMES_DIR}/$accession.fna"
			if [ "$f" != "$target" ]; then
				mv -f "$f" "$target"
			fi
		}}

		find {GENOMES_DIR} -maxdepth 1 -type f -name "*.fna" -print0 | while IFS= read -r -d '' f; do
			rename_fna "$f"
		done

		touch {output}
		"""
