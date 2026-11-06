# gtdb-tk_classify_wf

Taxonomically classifies genome assemblies against the Genome Taxonomy Database (GTDB) with GTDB-Tk `classify_wf`.

## Input

- `sample_file`: file with one sample name per line.
- `raw_dir`: directory with the genomes as `{sample}.{gtdbtk.extension}` (default `.fna`), e.g. from the annotation module.

## Output

- `out_dir/gtdbtk.bac120.summary.tsv`: GTDB classification per genome (bacteria). This is the file tracked by Snakemake.
- `out_dir/gtdbtk.ar53.summary.tsv`: archaeal classifications, written by GTDB-Tk only if archaeal genomes are present (not tracked).
- `out_dir/batchfile.tsv`: genome path and name table passed to GTDB-Tk.

## Config

See [config.yaml](config.yaml). Key variables: `sample_file`, `raw_dir`, `out_dir`, `db_dir`, `threads`, and the `gtdbtk` block:

- `extension`: genome file extension (`fna`).
- `skip_ani_screen`: skip the ANI pre-screen against GTDB genomes (`true`).
- `pplacer_cpus`: CPUs for pplacer, which is RAM-heavy (`1`).
- `extra_args`: any additional `gtdbtk classify_wf` options.

## Environment / databases

- Conda env: `environment.yaml` (`gtdbtk>=2.4`).
- Database: GTDB-Tk reference data (\~100 GB) in `db_dir`, exported as `GTDBTK_DATA_PATH`. Download it with `gtdbtk download` or `download-db.sh`.

## Usage

``` bash
conda activate snakemake

snakemake --use-conda --cores 8 --configfile modules/07_genome_chaacterization/gtdb-tk_classify_wf/config.yaml
```

## Notes

- pplacer needs a lot of memory (roughly 100 GB or more).
- Cite: Chaumeil et al., *GTDB-Tk v2*, Bioinformatics 2022.