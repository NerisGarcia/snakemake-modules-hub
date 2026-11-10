# module_name

One-sentence description of what this module does and which tool(s) it wraps.

## Input
- `sample_file`: ...
- `raw_dir`, `input_extension`: ...
- `pairend_extension_pattern` and optional `unpaired_input_pattern`: ...
- `input_file`: ... (for a specific shared file, when applicable)

## Output
- `out_dir/<tool>/...`: ...

## Config
See [config.yaml](config.yaml). Key variables: `dataset_name`, `out_dir`, `threads`, ...

## Environment / databases
- Conda env: `environment.yaml`
- Databases/references: ... (`null` = module default)

## Usage
```bash
conda activate snakemake
snakemake --use-conda --cores 8 --configfile modules/<category>/module_name/config.yaml
```

## Notes
Citation / known limitations (optional).
