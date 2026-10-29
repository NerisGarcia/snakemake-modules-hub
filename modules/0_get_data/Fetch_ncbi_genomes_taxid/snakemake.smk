configfile: "Fetch_ncbi_genomes_taxid/config.yml"

TAXID = config.get("taxid")
DATASET_DIR = config.get("dataset_out_dir")

DATASET_NAME = f"{TAXID}_ncbi_genomes"
SUMMARYFILE = f"{DATASET_DIR}/{DATASET_NAME}_summary.tsv"
ACCESSIONS_FILE = f"{DATASET_DIR}/{DATASET_NAME}_accessions.txt"


rule all:
	input:
		SUMMARYFILE,
		ACCESSIONS_FILE


rule download_summary:
	output:
		SUMMARYFILE
	params:
		taxid=TAXID,
		dir=DATASET_DIR
	conda:
		"ncbi_download"
	shell:
		"""
		mkdir -p {params.dir}
		datasets summary genome taxon {params.taxid} --mag exclude --as-json-lines | dataformat tsv genome > {output}

		"""

rule extract_accessions:
	input:
		SUMMARYFILE
	output:
		ACCESSIONS_FILE
	shell:
		"""
		awk -F'\t' 'NR==1 {{ for (i=1;i<=NF;i++) if ($i=="Assembly Accession") col=i; next }} {{ print $col }}' {input} > {output}
		"""
