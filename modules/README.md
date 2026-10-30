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

Each module owns its tool rules and tool-specific postprocessing. Keep sample lists and pipeline-level reporting in workflows. Module settings use snake_case keys and are passed to imported modules through Snakemake's `module` configuration. `config.yaml`, `environment.yaml`, `scripts/`, and `README.md` are included when applicable.