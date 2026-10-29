import csv
import math


def alignment_stats(path):
    """Return (number of sequences, alignment length) from a FASTA alignment."""
    n_seqs, length, in_first = 0, 0, False
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                n_seqs += 1
                in_first = n_seqs == 1
            elif in_first:
                length += len(line.strip())
    return n_seqs, length


def pct(part, whole):
    return f"{100 * part / whole:.2f}" if whole else "NA"


threshold = float(snakemake.params.core_threshold)

with open(snakemake.input.rtab) as handle:
    reader = csv.reader(handle, delimiter="\t")
    samples = next(reader)[1:]
    presence = [sum(int(v) for v in row[1:]) for row in reader if row]

n_samples = len(samples)
total_genes = len(presence)
min_samples = math.ceil(threshold * n_samples - 1e-9)
n_core_all = sum(1 for count in presence if count >= min_samples)
aln_seqs, aln_len = alignment_stats(snakemake.input.alignment)

rows = [
    ("samples_in_pangenome", n_samples),
    ("total_genes_detected", total_genes),
    ("core_threshold", threshold),
    ("min_samples_for_core_gene", f"{min_samples} of {n_samples}"),
    ("core_genes_at_threshold", n_core_all),
    ("core_genes_pct_of_total_genes", pct(n_core_all, total_genes)),
]

core_genes_used = n_core_all
core_groups = getattr(snakemake.input, "core_groups", None)
if core_groups:
    with open(core_groups, newline="") as handle:
        groups = list(csv.DictReader(handle, delimiter="\t"))
    selected = [g for g in groups if g["Status"] == "selected"]
    core_genes_used = len(selected)
    n_core_samples = int(groups[0]["Core_samples_total"]) if groups else 0
    rows += [
        ("core_definition_samples", n_core_samples),
        ("core_definition_samples_pct_of_pangenome_samples", pct(n_core_samples, n_samples)),
        ("analysis_samples", snakemake.params.analysis_samples),
        ("selected_core_genes", core_genes_used),
        ("selected_core_genes_pct_of_total_genes", pct(core_genes_used, total_genes)),
    ]

rows += [
    ("alignment_file", snakemake.input.alignment),
    ("alignment_samples", aln_seqs),
    ("alignment_length_bp", aln_len),
    ("alignment_mean_bp_per_core_gene", f"{aln_len / core_genes_used:.1f}" if core_genes_used else "NA"),
]

with open(snakemake.output.summary, "w") as handle:
    handle.write("metric\tvalue\n")
    for key, value in rows:
        handle.write(f"{key}\t{value}\n")
