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

## Sample-based input contract

The [configuration template](New_module_template/config.yaml) defines the
common names. `sample_file` contains one sample identifier per line, `raw_dir`
is the input directory, and `out_dir` is the parent of the module's
tool-specific output directory. Empty lines and comment lines beginning with
`#` are ignored.

- A single per-sample input is `raw_dir/{sample}<input_extension>`. For
  example, `".fna"` selects `raw_dir/sample.fna` and `".contigs.fa"` selects
  `raw_dir/sample.contigs.fa`. Configure one suffix per run; there is no
  extension guessing.
- Paired reads are `raw_dir/{sample}<pairend_extension_pattern><input_extension>`.
  The pattern substitutes `{read}` with `1` and `2`: `".fastp.R{read}"` and
  `".fastq.gz"` select `sample.fastp.R1.fastq.gz` and `sample.fastp.R2.fastq.gz`.
- When supported, `unpaired_input_pattern` is a path relative to `raw_dir`
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
updated; consult each module's config and README for its exact output paths.
The former `input_pattern` setting is no longer used: inputs must be in
`raw_dir` with the sample name followed by `input_extension`. Shovill exposes
its assemblies as flat links while retaining the original per-sample outputs.
Fastp filtering is now configured using `fastp.extra_args`; its paired outputs
are named `{sample}.fastp.R1{input_extension}` and
`{sample}.fastp.R2{input_extension}`, and it also produces combined unpaired
reads and reports.