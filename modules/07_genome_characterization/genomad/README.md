# genomad

Identifies plasmids and viruses (including proviruses) in each assembly with geNomad `end-to-end`.

## Input
- `sample_file`: one sample name per line.
- `input_dir`, `input_extension`: assemblies are read from `input_dir/{sample}{input_extension}`.

## Output
- `out_dir/genomad/{sample}_genomad/`: geNomad results; `{sample}_summary.log` marks completion.

## Config
See [config.yaml](config.yaml). Key variables: `out_dir`, `db_dir`, `genomad.splits`, `genomad.extra_args`, `threads`.

## Environment / databases
- Conda env: `genomad.yaml`
- Database: `db_dir` (required; download with `genomad download-database`).

## Usage
```bash
conda activate snakemake
snakemake --use-conda --cores 8 --snakefile modules/07_genome_characterization/genomad/Snakefile \
    --configfile modules/07_genome_characterization/genomad/config.yaml
```

## Notes
`--cleanup` removes intermediate files, so reruns start from scratch for a sample.
