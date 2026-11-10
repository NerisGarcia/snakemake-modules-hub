# gtdb-tk_classify_wf

Taxonomically classifies genome assemblies against the Genome Taxonomy Database (GTDB) with GTDB-Tk `classify_wf`.

## Input

- `genome_list`: a file containing one complete genome path per line. This
  module uses a path list rather than sample names and a common input directory.

## Output

- `out_dir/GTDB-Tk/gtdbtk.bac120.summary.tsv`: GTDB classification per genome (bacteria). This is the file tracked by Snakemake.
- `out_dir/GTDB-Tk/gtdbtk.ar53.summary.tsv`: archaeal classifications, written by GTDB-Tk only if archaeal genomes are present (not tracked).
- `out_dir/GTDB-Tk/batchfile.tsv`: genome path and name table passed to GTDB-Tk.

## Config

See [config.yaml](config.yaml). Key variables: `genome_list`, `out_dir`, `db_dir`, `threads`, and the `gtdbtk` block:

- `skip_ani_tree`: whether to skip ANI tree placement.
- `pplacer_cpus`: CPUs for pplacer, which is RAM-heavy (`1`).
- `extra_args`: any additional `gtdbtk classify_wf` options.

## Environment / databases

- Conda env: `envs/cenv.gtdbtk.yaml`.
- Database: GTDB-Tk reference data (\~100 GB) in `db_dir`, exported as `GTDBTK_DATA_PATH`. Download it with `gtdbtk download` or `download-db.sh`.

## Usage

``` bash
conda activate snakemake

snakemake --use-conda --cores 8 --snakefile modules/07_genome_characterization/gtdb-tk_classify_wf/Snakefile --configfile modules/07_genome_characterization/gtdb-tk_classify_wf/config.yaml
```

## Notes

- pplacer needs a lot of memory (roughly 100 GB or more).
- Cite: Chaumeil et al., *GTDB-Tk v2*, Bioinformatics 2022.