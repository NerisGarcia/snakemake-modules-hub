# amrfinderplus

Finds antimicrobial resistance, stress and virulence genes in each assembly with NCBI AMRFinderPlus.

## Input
- `sample_file`: one sample name per line.
- `input_dir`, `input_extension`: nucleotide FASTA read from `input_dir/{sample}{input_extension}`.

## Output
- `out_dir/amrfinderplus/{sample}.amrfinderplus.tsv`: AMRFinderPlus table (`logs/` holds the tool log).

## Config
See [config.yaml](config.yaml). Key variables: `out_dir`, `db_dir`, `amrfinderplus.organism`, `amrfinderplus.extra_args`, `threads`.

## Environment / databases
- Conda env: `amrfinderplus.yaml`
- Database: `db_dir` (`null` = the database in the env; update it with `amrfinder --update`).

## Usage
```bash
conda activate snakemake
snakemake --use-conda --cores 8 --snakefile modules/07_genome_characterization/amrfinderplus/Snakefile \
    --configfile modules/07_genome_characterization/amrfinderplus/config.yaml
```

## Notes
Runs in nucleotide-only mode; protein/GFF input (which can improve results) is not used.
