"""Regression tests for module input/output contracts using lightweight fixtures."""

import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = shutil.which("snakemake")


@unittest.skipUnless(SNAKEMAKE, "Snakemake is not available on PATH")
class ModuleContractTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.raw = self.root / "raw"
        self.raw.mkdir()
        self.samples = self.root / "samples.txt"
        self.samples.write_text("alpha\nbeta\n")

    def run_snakemake(self, snakefile, config, *args, env=None):
        configfile = self.root / "config.json"
        configfile.write_text(json.dumps(config))
        return subprocess.run(
            [
                SNAKEMAKE,
                "--snakefile",
                str(ROOT / snakefile),
                "--configfile",
                str(configfile),
                "--directory",
                str(self.root),
                "--cores",
                "1",
                *args,
            ],
            env=env,
            text=True,
            capture_output=True,
            check=True,
        )

    def test_quast_targets_each_sample_with_compound_extension(self):
        for sample in ("alpha", "beta"):
            (self.raw / f"{sample}.contigs.fa").write_text(">contig\nACGT\n")
        result = self.run_snakemake(
            "modules/02_assembly/quast/Snakefile",
            {
                "sample_file": str(self.samples),
                "raw_dir": str(self.raw),
                "input_extension": ".contigs.fa",
                "out_dir": str(self.root / "results"),
                "threads": 1,
            },
            "--dry-run",
            "--printshellcmds",
        )
        output = result.stdout + result.stderr
        for sample in ("alpha", "beta"):
            self.assertIn(f"{sample}.contigs.fa", output)
            self.assertIn(f"QUAST/{sample}.quast.transposed_report.tsv", output)

    def test_fastp_combines_unpaired_reads_and_generates_reports(self):
        for sample in ("alpha", "beta"):
            for read in (1, 2):
                with gzip.open(self.raw / f"{sample}_R{read}.fastq.gz", "wt") as stream:
                    stream.write(f"@{sample}\nACGT\n+\nIIII\n")

        bindir = self.root / "bin"
        bindir.mkdir()
        executable = bindir / "fastp"
        executable.write_text(
            "#!/usr/bin/env python3\n"
            "import gzip, json, sys\n"
            "args = sys.argv[1:]\n"
            "def value(flag): return args[args.index(flag) + 1]\n"
            "for flag in ('--out1', '--out2', '--unpaired1', '--unpaired2'):\n"
            "    with gzip.open(value(flag), 'wt') as handle: handle.write('@read\\nACGT\\n+\\nIIII\\n')\n"
            "with open(value('--json'), 'w') as handle:\n"
            "    json.dump({'summary': {'before_filtering': {'total_reads': 2, 'gc_content': 0.5}, "
            "'after_filtering': {'total_reads': 2, 'gc_content': 0.5}}}, handle)\n"
            "with open(value('--html'), 'w') as handle: handle.write('report')\n"
        )
        executable.chmod(0o755)
        env = os.environ.copy()
        env["PATH"] = f"{bindir}{os.pathsep}{env['PATH']}"
        self.run_snakemake(
            "modules/01_read_preprocessing/fastp/Snakefile",
            {
                "sample_file": str(self.samples),
                "raw_dir": str(self.raw),
                "pairend_extension_pattern": "_R{read}",
                "input_extension": ".fastq.gz",
                "out_dir": str(self.root / "results"),
                "threads": 1,
                "fastp": {"extra_args": "--length_required 4"},
            },
            env=env,
        )
        output_dir = self.root / "results" / "fastp"
        for sample in ("alpha", "beta"):
            with gzip.open(output_dir / f"{sample}.fastp.unpaired.fastq.gz", "rt") as stream:
                self.assertEqual(stream.read().count("@read\n"), 2)
            self.assertTrue((output_dir / f"{sample}.fastp.R1.fastq.gz").is_file())
            self.assertTrue((output_dir / f"{sample}.fastp.R2.fastq.gz").is_file())
            self.assertIn(sample, (output_dir / f"{sample}_report.tsv").read_text())
            self.assertFalse((output_dir / f"{sample}.fastp.unpaired.R1.fastq.gz").exists())
            self.assertFalse((output_dir / f"{sample}.fastp.unpaired.R2.fastq.gz").exists())

    def test_cleanifier_optional_unpaired_reads(self):
        for sample in ("alpha", "beta"):
            for suffix in ("_1", "_2", ".unpaired"):
                with gzip.open(self.raw / f"{sample}{suffix}.fastq.gz", "wt") as stream:
                    stream.write(f"@{sample}\nACGT\n+\nIIII\n")
        bindir = self.root / "bin"
        bindir.mkdir()
        cleanifier = bindir / "cleanifier"
        cleanifier.write_text(
            "#!/usr/bin/env python3\n"
            "import gzip, sys\n"
            "args = sys.argv[1:]\n"
            "prefix = args[args.index('--out') + 1]\n"
            "suffixes = ('_1_keep.fq.gz', '_2_keep.fq.gz') if '--pairs' in args else ('_keep.fq.gz',)\n"
            "for suffix in suffixes:\n"
            "    with gzip.open(prefix + suffix, 'wt') as handle: handle.write('@read\\nACGT\\n+\\nIIII\\n')\n"
        )
        cleanifier.chmod(0o755)
        seqkit = bindir / "seqkit"
        seqkit.write_text("#!/bin/sh\nprintf 'sample statistics\\n'\n")
        seqkit.chmod(0o755)
        env = os.environ.copy()
        env["PATH"] = f"{bindir}{os.pathsep}{env['PATH']}"
        config = {
            "sample_file": str(self.samples),
            "raw_dir": str(self.raw),
            "pairend_extension_pattern": "_{read}",
            "input_extension": ".fastq.gz",
            "unpaired_input_pattern": "{sample}.unpaired",
            "out_dir": str(self.root / "results"),
            "cleanifier_index": "test-index",
            "threads": 1,
        }
        self.run_snakemake("modules/01_read_preprocessing/cleanifier/Snakefile", config, env=env)
        output_dir = self.root / "results" / "cleanifier"
        for sample in ("alpha", "beta"):
            for suffix in ("R1", "R2", "unpaired"):
                self.assertTrue((output_dir / f"{sample}.cleanifier.{suffix}.fastq.gz").is_file())
            self.assertIn("statistics", (output_dir / "qc" / f"{sample}.cleanifier.seqkit_stats.txt").read_text())

    def test_quast_accepts_upstream_generated_files(self):
        workflow = self.root / "Snakefile"
        module = ROOT / "modules/02_assembly/quast/Snakefile"
        results = self.root / "results"
        workflow.write_text(
            f"module quast:\n    snakefile: {str(module)!r}\n"
            f"    config: {{'sample_file': {str(self.samples)!r}, "
            f"'raw_dir': {str(results / 'bakta')!r}, "
            f"'input_extension': '.fna', 'out_dir': {str(results)!r}, 'threads': 1}}\n"
            "use rule * from quast exclude all as quast_*\n"
            "rule all:\n    default_target: True\n"
            f"    input: {str(results / 'QUAST/alpha.quast.transposed_report.tsv')!r}\n"
            "rule upstream:\n"
            f"    output: {str(results / 'bakta/{sample}.fna')!r}\n"
            "    shell: \"touch {output}\"\n"
        )
        result = subprocess.run(
            [SNAKEMAKE, "--snakefile", str(workflow), "--directory", str(self.root), "--cores", "1", "--dry-run"],
            text=True,
            capture_output=True,
            check=True,
        )
        output = result.stdout + result.stderr
        self.assertIn("upstream", output)
        self.assertIn("quast_assembly_qc_sample", output)
        self.assertIn("bakta/alpha.fna", output)

    def test_shovill_links_existing_contigs_for_bakta(self):
        output_dir = self.root / "results" / "shovill"
        for sample in ("alpha", "beta"):
            for read in (1, 2):
                (self.raw / f"{sample}_{read}.fastq.gz").write_text("reads\n")
            nested = output_dir / f"{sample}_shovill"
            nested.mkdir(parents=True)
            (nested / f"{sample}.contigs.fa").write_text(f">{sample}\nACGT\n")
        self.run_snakemake(
            "modules/02_assembly/shovill/Snakefile",
            {
                "sample_file": str(self.samples),
                "raw_dir": str(self.raw),
                "pairend_extension_pattern": "_{read}",
                "input_extension": ".fastq.gz",
                "out_dir": str(self.root / "results"),
                "threads": 1,
            },
        )
        for sample in ("alpha", "beta"):
            link = output_dir / f"{sample}.contigs.fa"
            self.assertTrue(link.is_symlink())
            self.assertEqual(link.read_text(), f">{sample}\nACGT\n")
            self.assertEqual(os.readlink(link), f"{sample}_shovill/{sample}.contigs.fa")
        result = self.run_snakemake(
            "modules/03_gene_annotation/bakta/Snakefile",
            {
                "sample_file": str(self.samples),
                "raw_dir": str(output_dir),
                "input_extension": ".contigs.fa",
                "out_dir": str(self.root / "results"),
                "bakta_db": "test-db",
                "threads": 1,
            },
            "--dry-run",
        )
        for sample in ("alpha", "beta"):
            self.assertIn(f"shovill/{sample}.contigs.fa", result.stdout + result.stderr)

    def test_single_file_sample_modules_select_each_sample(self):
        for sample in ("alpha", "beta"):
            (self.raw / f"{sample}.fasta").write_text(f">{sample}\nACGT\n")
        common = {
            "sample_file": str(self.samples),
            "raw_dir": str(self.raw),
            "input_extension": ".fasta",
            "out_dir": str(self.root / "results"),
            "threads": 1,
        }
        cases = [
            ("modules/03_gene_annotation/bakta/Snakefile", {"bakta_db": "test-db"}, "bakta/{sample}.gff3"),
            ("modules/02_assembly/busco/Snakefile", {"busco_lineage": "bacteria"}, "BUSCO/{sample}_busco_short_summary.txt"),
            ("modules/02_assembly/checkm2/Snakefile", {"checkm2_db": str(self.samples)}, "CheckM2/{sample}/quality_report.tsv"),
        ]
        for snakefile, extra, expected in cases:
            with self.subTest(snakefile=snakefile):
                result = self.run_snakemake(snakefile, {**common, **extra}, "--dry-run")
                output = result.stdout + result.stderr
                for sample in ("alpha", "beta"):
                    self.assertIn(f"{sample}.fasta", output)
                    self.assertIn(expected.format(sample=sample), output)

    def test_duplicate_sample_names_are_rejected(self):
        self.samples.write_text("alpha\nalpha\n")
        config = {
            "sample_file": str(self.samples),
            "raw_dir": str(self.raw),
            "input_extension": ".fasta",
            "out_dir": str(self.root / "results"),
            "threads": 1,
        }
        with self.assertRaises(subprocess.CalledProcessError) as error:
            self.run_snakemake("modules/02_assembly/quast/Snakefile", config, "--dry-run")
        self.assertIn("duplicate identifiers", error.exception.stdout + error.exception.stderr)

    def test_virgo2_combines_configured_paired_and_unpaired_reads(self):
        for suffix in ("_1", "_2", ".unpaired"):
            with gzip.open(self.raw / f"alpha{suffix}.fastq.gz", "wt") as stream:
                stream.write(f"@{suffix}\nACGT\n+\nIIII\n")
        target = self.root / "results/VIRGO2/reads/alpha.se.fastq.gz"
        self.run_snakemake(
            "modules/10_metagenomic_taxonomy/virgo2/Snakefile",
            {
                "dataset_name": "example",
                "sample_file": str(self.samples),
                "raw_dir": str(self.raw),
                "pairend_extension_pattern": "_{read}",
                "unpaired_input_pattern": "{sample}.unpaired",
                "input_extension": ".fastq.gz",
                "out_dir": str(self.root / "results"),
                "virgo2_script": "not-used-for-read-combination",
                "threads": 1,
            },
            str(target),
        )
        with gzip.open(target, "rt") as stream:
            reads = stream.read()
        self.assertIn("@_1\n", reads)
        self.assertIn("@_2\n", reads)
        self.assertIn("@.unpaired\n", reads)


if __name__ == "__main__":
    unittest.main()
