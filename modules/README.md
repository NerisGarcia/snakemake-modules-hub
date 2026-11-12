# Module layout

Modules are organized by numbered, snake_case category and tool directories:

``` text
modules/NN_category/tool_name/
  Snakefile
  config.yaml
  environment.yaml
  scripts/
  README.md
```

Each module owns its tool rules and tool-specific postprocessing. Module
settings use snake_case keys and are passed to imported modules through
Snakemake's `module` configuration. `config.yaml`, `environment.yaml`,
`scripts/`, and `README.md` are included when applicable.

## Local module tests

The implemented module `config.yaml` files use existing local files in
`data/test/`, except for the unfinished modules and the external resources
listed below. Data acquisition modules were excluded. From the repository
root, run an individual module with its config and Conda environment, for
example:

```bash
conda run -n snakemake snakemake --cores 4 --use-conda \
  -s modules/02_assembly/quast/Snakefile \
  --configfile modules/02_assembly/quast/config.yaml
```

Replace both paths with those of the module to test. Results are written to
`tests_runs/modules/` so existing inputs and results are not overwritten.
The read files are several gigabytes each, so read-based modules can take
substantial time and storage. The existing reads are metagenomic rather than
isolate reads, so Shovill's test may not produce a useful assembly.
`data/` and `tests_runs/` are ignored by Git; these configs depend on local
data and are not self-contained in a clone.
The alignment and genome path lists in `data/test/module_test/` refer to
existing files in `data/test/genomics/`.

Some tools also require references: Bakta, CheckM2, VIRGO2 and VISTA use
local installation paths in their configs; change these if needed. BUSCO may
download its lineage. Cleanifier needs a matching index, Kraken2 and Bracken
need a built database, and GTDB-Tk needs its reference data; supply valid
paths in their configs before running. PanACoTA needs a formatted database and
its current workflow does not define a producer for its final tree. FastTree
is empty, and trimAl and MetaBAT are unfinished. Their configs alone cannot
make them runnable. HmmCleaner is not tested locally. The gene alignment
cleaning test uses a high variable-site threshold to avoid running HmmCleaner,
so it does not test the cleaning branch.

## Sample-based input contract

The [configuration template](New_module_template/config.yaml) defines the
common names. `sample_file` contains one sample identifier per line, `input_dir`
is the input directory, and `out_dir` is the parent of the module's
tool-specific output directory. Empty lines and comment lines beginning with
`#` are ignored.

- A single per-sample input is `input_dir/{sample}<input_extension>`. For
  example, `".fna"` selects `input_dir/sample.fna` and `".contigs.fa"` selects
  `input_dir/sample.contigs.fa`. Configure one suffix per run; there is no
  extension guessing.
- Paired reads are `input_dir/{sample}<pairend_extension_pattern><input_extension>`.
  The pattern substitutes `{read}` with `1` and `2`: `".fastp.R{read}"` and
  `".fastq.gz"` select `sample.fastp.R1.fastq.gz` and `sample.fastp.R2.fastq.gz`.
- When supported, `unpaired_input_pattern` is a path relative to `input_dir`
  containing `{sample}`; append `input_extension`. Set it to null if no
  unpaired input is available.

Set these same values when importing a module into a workflow. Explicit filenames
allow Snakemake to connect outputs that upstream rules have not created yet.
Standalone targets select the listed samples; workflows can exclude `all`
and choose their own targets.

Some modules have fundamentally different inputs: acquisition modules use
taxonomic IDs or accession lists; GTDB-Tk uses a list of complete genome
paths; and other modules consume shared tables, alignments, or Panaroo
results. These inputs are documented in the relevant module configuration.
PanACoTA was excluded from this migration; a future revision may use
`out_dir/PanACoTA`. FastTree is empty and trimAl and MetaBAT are unfinished;
they are not usable standardized modules yet.

## Migration

The standardized keys and tool subdirectories replace previous module-specific
input/output keys. Existing custom configs and downstream scripts must be
updated: rename the former `raw_dir` key to `input_dir` in configs and module
imports. Consult each module's config and README for its exact output paths.
The former `input_pattern` setting is no longer used: inputs must be in
`input_dir` with the sample name followed by `input_extension`. Shovill exposes
its assemblies as flat links while retaining the original per-sample outputs.
Fastp filtering is now configured using `fastp.extra_args`; its paired outputs
are named `{sample}.fastp.R1{input_extension}` and
`{sample}.fastp.R2{input_extension}`, and it also produces combined unpaired
reads and reports.