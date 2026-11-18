# Genome Characterization Workflow

Genome characterization workflow that runs AMRFinderPlus, DefenseFinder, and geNomad in parallel. antiSMASH is run independently and is not included here.

## Workflow Overview

This workflow processes assembled bacterial genomes to identify:
- **AMRFinderPlus**: Antimicrobial resistance, stress, and virulence genes
- **DefenseFinder**: Anti-phage defense systems
- **geNomad**: Plasmids and viral sequences (including proviruses)

## Input Requirements

1. **Sample list**: One sample identifier per line (text file)
2. **Assembled genomes**: Nucleotide FASTA files
   - Location: `input_dir/{sample}{input_extension}`
   - Format: `.fna` or other FASTA extension
3. **Annotated genomes** (for DefenseFinder): Bakta output with protein predictions
   - Location: `annotation_dir/{sample}.faa`
   - Format: Protein FASTA

## Output Structure

```
results_dir/07_genome_characterization/
├── amrfinderplus/
│   ├── {sample}.amrfinderplus.tsv
│   └── logs/
├── defensefinder/
│   ├── {sample}/
│   │   ├── {sample}_defense_finder_systems.tsv
│   │   ├── {sample}_defense_finder_genes.tsv
│   │   └── {sample}_defense_finder_hmmer.tsv
│   └── logs/
├── genomad/
│   ├── {sample}_genomad/
│   │   ├── {sample}_summary.log
│   │   ├── {sample}_virus_summary.txt
│   │   ├── {sample}_plasmid_summary.txt
│   │   └── other outputs...
│   └── logs/
```

## Configuration

See [config.yaml](config.yaml) for all options:

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sample_file` | Path to sample list file | `data/test/sample.list` |
| `input_dir` | Directory containing assembly FASTA files | `data/test/assembly` |
| `input_extension` | File extension for assemblies | `.fna` |
| `annotation_dir` | Directory containing protein predictions (bakta output) | `data/test/results/02_annotation/bakta` |
| `genomics_results_dir` | Output directory for all results | `data/test/results` |

### Database Configuration

- `amrfinderplus_db`: AMRFinderPlus database (null = env default; update with `amrfinder --update`)
- `defensefinder_db`: DefenseFinder models (null = env default; update with `defense-finder update`)
- `genomad_db`: geNomad database (required; download with `genomad download-database`)

### Tool-Specific Options

```yaml
amrfinderplus:
  organism: null           # e.g., "Escherichia" for species-specific results
  extra_args: ""

defensefinder:
  extra_args: ""

genomad:
  splits: 16              # Number of genomic splits (memory/speed tradeoff)
  extra_args: ""
```

### Thread Allocation

```yaml
threads:
  amrfinderplus: 4
  defensefinder: 4
  genomad: 8
```

## Usage

### 1. Download and Prepare Databases

The databases must be downloaded separately before running the workflow. This can be done in a single session with conda:

```bash
# Activate the conda environment with all tools
conda activate snakemake

# Download/update databases (can take 10-30 minutes depending on tool)
mkdir -p /path/to/databases

# AMRFinderPlus (required; ~200 MB)
amrfinder --update --database /path/to/databases/amrfinderplus_db

# DefenseFinder (required; ~100 MB)
defense-finder update --models-dir /path/to/databases/defensefinder_db

# geNomad (required; ~1 GB)
genomad download-database /path/to/databases/genomad_db
```

### 2. Configure the Workflow

Edit `config.yaml` with your project paths:

```bash
nano workflows/07_amr_smash_defense_genomad/config.yaml
```

Set these paths:

```yaml
sample_file: "projects/my_project/data/samples.list"
input_dir: "projects/my_project/data/assembly"
input_extension: ".fna"
annotation_dir: "projects/my_project/data/results/02_annotation/bakta"
genomics_results_dir: "projects/my_project/data/results"

# Point to the downloaded databases
amrfinderplus_db: "/path/to/databases/amrfinderplus_db"
defensefinder_db: "/path/to/databases/defensefinder_db"
genomad_db: "/path/to/databases/genomad_db"
```

If you omit the `db_dir` settings (keep them as `null`), Snakemake will use the databases installed with the conda environment instead.

### 3. Run the Workflow

```bash
# Dry-run to see what will execute
./workflows/07_amr_smash_defense_genomad/run.sh dry-run

# Execute the workflow
./workflows/07_amr_smash_defense_genomad/run.sh run

# Or manually with custom threads
snakemake --use-conda --cores 32 \
  --snakefile workflows/07_amr_smash_defense_genomad/Snakefile \
  --configfile workflows/07_amr_smash_defense_genomad/config.yaml
```

### 4. Check Results

After completion, results will be in:
```
{genomics_results_dir}/07_genome_characterization/
├── amrfinderplus/
├── defensefinder/
├── genomad/
```

## Basic Run

```bash
cd /local/scratch/neris/snakemake-modules-hub

# Configure config.yaml with your paths and settings
# Then run:
snakemake --use-conda --cores 32 \
  --snakefile workflows/07_amr_smash_defense_genomad/Snakefile \
  --configfile workflows/07_amr_smash_defense_genomad/config.yaml
```

### Advanced Options

```bash
# Dry run to see what will execute
snakemake --use-conda --cores 32 --dry-run \
  --snakefile workflows/07_amr_smash_defense_genomad/Snakefile \
  --configfile workflows/07_amr_smash_defense_genomad/config.yaml

# Run specific tool only
snakemake --use-conda --cores 32 \
  --snakefile workflows/07_amr_smash_defense_genomad/Snakefile \
  --configfile workflows/07_amr_smash_defense_genomad/config.yaml \
  amrfinderplus_all

# Run with job submission (HPC)
snakemake --use-conda --cores 32 --profile slurm \
  --snakefile workflows/07_amr_smash_defense_genomad/Snakefile \
  --configfile workflows/07_amr_smash_defense_genomad/config.yaml
```

## Dependencies

Each tool uses its own conda environment:
- `amrfinderplus.yaml` → AMRFinderPlus
- `defensefinder.yaml` → DefenseFinder
- `genomad.yaml` → geNomad

Snakemake will automatically create these environments with `--use-conda`.

## Notes

### Workflow Characteristics

- **Independent tools**: All three tools run in parallel (no dependencies between them)
- **Input flexibility**: Accepts any nucleotide FASTA format for assemblies
- **Failure tolerance**: Samples whose assembly is missing or empty are skipped and do not stop the workflow
- **Logging**: All tool outputs and logs are saved in `logs/` subdirectories

### DefenseFinder Special Note

DefenseFinder requires protein sequences (`.faa` files) as input. This workflow expects these to come from Bakta annotation. If proteins are unavailable, DefenseFinder will be skipped for those samples (marked as FAIL).

### Database Notes

Each tool requires its database to be downloaded before running the workflow:

1. **Download to a central location** (recommended for sharing across projects):
   ```bash
   mkdir -p ~/databases/genome_characterization
   amrfinder --update --database ~/databases/genome_characterization/amrfinderplus_db
   defense-finder update --models-dir ~/databases/genome_characterization/defensefinder_db
   genomad download-database ~/databases/genome_characterization/genomad_db
   ```

2. **Point config.yaml to these paths** (see Configuration section above)

3. **Alternative: Use conda environment defaults** by setting all db_dir to `null` - tools will use their embedded databases (simpler but slower for first run)

## Database Setup Details

| Tool | Required? | Download Command | Size | Location in Config |
|------|-----------|-----------------|------|-------------------|
| AMRFinderPlus | Yes | `amrfinder --update --database <path>` | ~200 MB | `amrfinderplus_db` |
| DefenseFinder | Yes | `defense-finder update --models-dir <path>` | ~100 MB | `defensefinder_db` |
| geNomad | Yes | `genomad download-database <path>` | ~1 GB | `genomad_db` |

**Total database disk space: ~1-2 GB**

## Troubleshooting

### Common Issues

**KeyError for database paths:**
- Ensure all required database paths are correctly set in config.yaml
- For null values, verify the databases are installed in the conda environment

**Missing protein files for DefenseFinder:**
- Ensure Bakta has been run successfully first
- Protein files should be at `{annotation_dir}/{sample}.faa`

**geNomad database errors:**
- Verify database was downloaded correctly: `ls {genomad_db}/` should show multiple files
- Run: `genomad download-database {genomad_db}` to update or repair

## References

- [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/)
- [DefenseFinder](https://github.com/mdmparis/defense-finder)
- [geNomad](https://github.com/apcamargo/genomad)
