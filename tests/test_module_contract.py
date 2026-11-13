"""Regression tests for module input/output contracts using lightweight fixtures."""

import gzip
import hashlib
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

    def panaroo_inputs(self):
        source = self.root / "previous"
        source.mkdir()
        rtab = source / "presence.tsv"
        rtab.write_text("Gene\talpha\tbeta\nfamily1\t1\t1\n")
        roary = source / "genes.csv"
        roary.write_text("Gene,Non-unique Gene name,alpha,beta\nfamily1,,a,b\n")
        alignments = source / "families"
        alignments.mkdir()
        (alignments / "family1.aln.fas").write_text(">alpha;1\nACGT\n>beta;2\nACGA\n")
        return {"rtab": str(rtab), "roary": str(roary), "gene_sequences": str(alignments)}

    def test_quast_targets_each_sample_with_compound_extension(self):
        for sample in ("alpha", "beta"):
            (self.raw / f"{sample}.contigs.fa").write_text(">contig\nACGT\n")
        result = self.run_snakemake(
            "modules/02_assembly/quast/Snakefile",
            {
                "sample_file": str(self.samples),
                "input_dir": str(self.raw),
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

    def test_amas_summary_uses_sample_file_alignments(self):
        for sample in ("alpha", "beta"):
            (self.raw / f"{sample}.fasta").write_text(f">{sample}\nACGT\n")
        result = self.run_snakemake(
            "modules/06_msa/amas_summary/Snakefile",
            {
                "sample_file": str(self.samples),
                "input_dir": str(self.raw),
                "input_extension": ".fasta",
                "out_dir": str(self.root / "results"),
            },
            "--dry-run",
            "--printshellcmds",
        )
        output = result.stdout + result.stderr
        self.assertIn("alpha.fasta", output)
        self.assertIn("beta.fasta", output)
        self.assertIn("amas_alignment_summary.tsv", output)

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
                "input_dir": str(self.raw),
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
            "input_dir": str(self.raw),
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
            f"'input_dir': {str(results / 'bakta')!r}, "
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
                "input_dir": str(self.raw),
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
                "input_dir": str(output_dir),
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
            "input_dir": str(self.raw),
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
            "input_dir": str(self.raw),
            "input_extension": ".fasta",
            "out_dir": str(self.root / "results"),
            "threads": 1,
        }
        with self.assertRaises(subprocess.CalledProcessError) as error:
            self.run_snakemake("modules/02_assembly/quast/Snakefile", config, "--dry-run")
        self.assertIn("duplicate identifiers", error.exception.stdout + error.exception.stderr)

    def test_core_subset_writes_parameters_for_optional_sample_lists(self):
        subset = self.root / "subset.txt"
        subset.write_text("alpha\n")
        panaroo = self.panaroo_inputs()
        for use_core, use_analysis in ((True, False), (False, True), (True, True)):
            with self.subTest(use_core=use_core, use_analysis=use_analysis):
                subcore = {"core_threshold": 0.8}
                subsample = {}
                if use_core:
                    subcore.update(core_sample_list=str(subset), label="core")
                if use_analysis:
                    subsample.update(analysis_sample_list=str(subset), label="analysis")
                prefix = "_".join(
                    part for part, enabled in (("subcore", use_core), ("subsamp", use_analysis))
                    if enabled
                )
                label = "__".join(
                    part for part, enabled in (("core", use_core), ("analysis", use_analysis))
                    if enabled
                )
                run_name = f"{prefix}.{label}.core_80"
                output = self.root / "results" / "core_subset" / run_name / "params.txt"
                config = {
                    "dataset_name": "example",
                    "out_dir": str(self.root / "results"),
                    "panaroo": panaroo,
                    "subcore": subcore,
                    "subsample": subsample,
                }
                self.run_snakemake("modules/04_pangenome/core_subset/Snakefile", config, str(output))
                params = dict(line.split(": ", 1) for line in output.read_text().splitlines())
                checksum = hashlib.md5(subset.read_bytes()).hexdigest()
                for name, enabled, expected_label in (
                    ("core", use_core, "core"),
                    ("analysis", use_analysis, "analysis"),
                ):
                    label_key = "subcore_label" if name == "core" else "subsample_label"
                    self.assertEqual(params[label_key], expected_label if enabled else "NA")
                    self.assertEqual(params[f"{name}_sample_list"], str(subset) if enabled else "NA")
                    self.assertEqual(params[f"{name}_sample_list_md5"], checksum if enabled else "NA")
                    self.assertEqual(params[f"{name}_samples"], "1" if enabled else "2")
                self.assertEqual(params["panaroo_samples"], "2")
                self.assertEqual(params["panaroo_samples_md5"], hashlib.md5(b"alpha\nbeta\n").hexdigest())
                sample_output = output.parent / "analysis_samples.txt"
                self.run_snakemake("modules/04_pangenome/core_subset/Snakefile", config, str(sample_output))
                expected_analysis = "alpha\n" if use_analysis else "alpha\nbeta\n"
                expected_core = ["alpha"] if use_core else ["alpha", "beta"]
                self.assertEqual(sample_output.read_text(), expected_analysis)
                selections = (output.parent / "sample_selection.tsv").read_text().splitlines()
                core_selections = [
                    line.split("\t", 1)[1] for line in selections
                    if line.startswith("core_definition\t")
                ]
                self.assertEqual(core_selections, expected_core)

    def test_editcore_uses_existing_panaroo_inputs_and_checks_samples(self):
        panaroo = self.panaroo_inputs()
        core = self.root / "core.txt"
        core.write_text("alpha\n")
        analysis = self.root / "analysis.txt"
        analysis.write_text("beta\n")
        results = self.root / "results"
        config = {
            "dataset_name": "example",
            "out_dir": str(results),
            "panaroo": panaroo,
            "subcore": {"core_sample_list": str(core), "core_threshold": 0.8, "label": "core"},
            "subsample": {"analysis_sample_list": str(analysis), "label": "analysis"},
            "threads": {"mafft": 1, "snp_dists": 1, "iqtree": 1},
            "gene_aln_cleaning": {},
            "iqtree": {"model": "GTR", "fast": True},
        }
        snakefile = "workflows/02_panaroo_then_editcore/Snakefile"
        run_dir = results / "core_subset/subcore_subsamp.core__analysis.core_80"
        groups = run_dir / "selected_core_groups.tsv"
        self.run_snakemake(snakefile, config, str(groups))
        self.assertIn("family1", groups.read_text())
        self.assertEqual(
            (run_dir / "selected_gene_fastas/subcore_subsamp.core__analysis.core_80_family1.fasta").read_text(),
            ">beta\nACGA\n",
        )

        for selected in (core, analysis):
            with self.subTest(invalid_list=selected.name):
                original = selected.read_text()
                selected.write_text("unknown\n")
                with self.assertRaises(subprocess.CalledProcessError) as error:
                    self.run_snakemake(snakefile, config, "--dry-run")
                self.assertIn("not in the Panaroo Rtab", error.exception.stdout + error.exception.stderr)
                selected.write_text(original)

        roary = Path(panaroo["roary"])
        roary.write_text("Gene,alpha\nfamily1,a\n")
        with self.assertRaises(subprocess.CalledProcessError) as error:
            self.run_snakemake(snakefile, config, "--dry-run")
        self.assertIn("exactly one column", error.exception.stdout + error.exception.stderr)

        roary.write_text("Gene,alpha,beta\nfamily1,a,b\n")
        Path(panaroo["rtab"]).write_text("Gene\talpha\nfamily1\t1\n")
        with self.assertRaises(subprocess.CalledProcessError) as error:
            self.run_snakemake(snakefile, config, "--dry-run")
        self.assertIn("not in the Panaroo Rtab", error.exception.stdout + error.exception.stderr)

        Path(panaroo["rtab"]).write_text("Gene\talpha\talpha\nfamily1\t1\t1\n")
        with self.assertRaises(subprocess.CalledProcessError) as error:
            self.run_snakemake(snakefile, config, "--dry-run")
        self.assertIn("duplicate sample columns", error.exception.stdout + error.exception.stderr)

        Path(panaroo["rtab"]).unlink()
        with self.assertRaises(subprocess.CalledProcessError) as error:
            self.run_snakemake(snakefile, config, "--dry-run")
        self.assertIn("Missing or empty Panaroo file", error.exception.stdout + error.exception.stderr)

    def test_paired_kraken_and_panaroo_use_input_dir(self):
        for sample in ("alpha", "beta"):
            for suffix in ("_1.fastq", "_2.fastq", ".unpaired.fastq", ".gff3"):
                (self.raw / f"{sample}{suffix}").write_text("fixture\n")
        common = {"sample_file": str(self.samples), "input_dir": str(self.raw),
                  "out_dir": str(self.root / "results"), "threads": 1}
        cases = [
            (
                "modules/01_read_preprocessing/kraken2/Snakefile",
                {**common, "input_extension": ".fastq", "pairend_extension_pattern": "_{read}",
                 "kraken_db": "test-db", "bracken_level": "S", "bracken_read_len": 150, "bracken_threshold": 10},
                ("{sample}_1.fastq", "{sample}_2.fastq"),
            ),
            (
                "modules/01_read_preprocessing/kraken2_paired_end/Snakefile",
                {**common, "input_extension": ".fastq", "pairend_extension_pattern": "_{read}",
                 "unpaired_input_pattern": "{sample}.unpaired", "kraken_db": "test-db"},
                ("{sample}_1.fastq", "{sample}_2.fastq", "{sample}.unpaired.fastq"),
            ),
            (
                "modules/04_pangenome/panaroo/Snakefile",
                {**common, "input_extension": ".gff3", "core_threshold": 0.9,
                 "clean_mode": "strict", "alignment": "pan", "family_threshold": 0.7,
                 "refind_mode": "default", "aligner": "mafft"},
                ("{sample}.gff3",),
            ),
        ]
        for snakefile, config, inputs in cases:
            with self.subTest(snakefile=snakefile):
                result = self.run_snakemake(snakefile, config, "--dry-run")
                output = result.stdout + result.stderr
                for sample in ("alpha", "beta"):
                    for input_name in inputs:
                        self.assertIn(str(self.raw / input_name.format(sample=sample)), output)

    def test_workflows_connect_modules_with_input_dir(self):
        samples = ("alpha", "beta")
        for sample in samples:
            (self.raw / f"{sample}.fna").write_text(f">{sample}\nACGT\n")
            for read in (1, 2):
                (self.raw / f"{sample}_{read}.fastq.gz").write_text("reads\n")
        genomes = self.root / "genomes.txt"
        genomes.write_text("".join(f"{self.raw / (sample + '.fna')}\n" for sample in samples))
        results = self.root / "results"
        common = {
            "sample_file": str(self.samples),
            "input_dir": str(self.raw),
            "input_extension": ".fna",
            "bakta_db": "test-db",
            "checkm2_db": str(self.samples),
            "busco_lineage": "bacteria",
            "threads": dict.fromkeys(("fastp", "kraken2", "shovill", "bakta", "busco", "quast", "checkm2", "gtdbtk"), 1),
        }
        cases = [
            (
                "workflows/01_fasta_qc_bakta/Snakefile",
                {**common, "genomics_results_dir": str(results), "summary_table": str(results / "summary.tsv")},
                ("bakta_bakta", "busco_busco", "quast_assembly_qc_sample", "checkm2_checkm2"),
            ),
            (
                "workflows/gtdbk_bakta_quast/Snakefile",
                {**common, "genome_list": str(genomes), "results_dir": str(results),
                 "gtdbtk_db": "test-db", "gtdbtk": {"skip_ani_tree": True, "pplacer_cpus": 1, "extra_args": ""}},
                ("bakta_bakta", "quast_assembly_qc_sample", "gtdbtk_gtdbtk_classify"),
            ),
            (
                "workflows/01_fastq_shovill_qc_bakta/Snakefile",
                {**common, "sra_acc_list": str(self.samples), "sra_raw": str(self.raw),
                 "clean_dir": str(self.root / "clean"), "genomics_results_dir": str(results),
                 "summary_table": str(results / "summary.tsv"), "fastp": {"extra_args": ""},
                 "kraken_db": "test-db", "bracken_read_len": 150, "bracken_level": "S",
                 "bracken_threshold": 10},
                ("fastp_fastp", "shovill_shovill", "bakta_bakta", "quast_assembly_qc_sample"),
            ),
        ]
        for snakefile, config, expected_rules in cases:
            with self.subTest(snakefile=snakefile):
                result = self.run_snakemake(snakefile, config, "--dry-run")
                output = result.stdout + result.stderr
                for rule in expected_rules:
                    self.assertIn(rule, output)

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
                "input_dir": str(self.raw),
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
