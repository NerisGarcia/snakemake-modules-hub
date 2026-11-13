# antismash

Detects biosynthetic gene clusters (secondary metabolites) in each genome with antiSMASH.

## Input
- `sample_file`: one sample name per line.
- `input_dir`, `input_extension`: genome read from `input_dir/{sample}{input_extension}`, a nucleotide FASTA or an annotated GenBank file.

## Output
- `out_dir/antismash/{sample}_antismash/`: antiSMASH results (`{sample}.json`, `{sample}.gbk`, HTML report); `logs/` holds the tool log.

## Config
See [config.yaml](config.yaml). Key variables: `out_dir`, `db_dir`, `antismash.taxon`, `antismash.genefinding_tool`, `antismash.extra_args`, `threads`.

## Environment / databases
- Conda env: `antismash.yaml`
- Databases: `db_dir` (`null` = installed with the env; run `download-antismash-databases` once).

## Usage
```bash
conda activate snakemake
snakemake --use-conda --cores 8 --snakefile modules/07_genome_characterization/antismash/Snakefile \
    --configfile modules/07_genome_characterization/antismash/config.yaml
```

## Notes
For a FASTA, antiSMASH predicts genes with prodigal. For a GenBank input (e.g. Bakta `.gbff`) set `antismash.genefinding_tool: none`.
