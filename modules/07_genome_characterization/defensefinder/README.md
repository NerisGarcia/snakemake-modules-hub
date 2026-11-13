# defensefinder

Detects anti-phage defense systems in the predicted proteins of each sample with DefenseFinder.

## Input
- `sample_file`: one sample name per line.
- `input_dir`, `input_extension`: protein FASTA read from `input_dir/{sample}{input_extension}`.

## Output
- `out_dir/defensefinder/{sample}/{sample}_defense_finder_{systems,genes,hmmer}.tsv`

## Config
See [config.yaml](config.yaml). Key variables: `out_dir`, `db_dir`, `defensefinder.extra_args`, `threads`.

## Environment / databases
- Conda env: `defensefinder.yaml`
- Models: `db_dir` (`null` = tool default; download with `defense-finder update`).

## Usage
```bash
conda activate snakemake
snakemake --use-conda --cores 8 --snakefile modules/07_genome_characterization/defensefinder/Snakefile \
    --configfile modules/07_genome_characterization/defensefinder/config.yaml
```
