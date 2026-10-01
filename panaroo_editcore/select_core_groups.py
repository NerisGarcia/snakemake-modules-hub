import csv
import os
import shutil
from collections import Counter


def read_samples(path):
    """Read unique sample IDs, allowing blank lines and # comments."""
    with open(path) as handle:
        samples = [
            line.strip()
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]
    if not samples:
        raise ValueError(f"Sample list is empty: {path}")
    if len(samples) != len(set(samples)):
        raise ValueError(f"Sample list contains duplicate names: {path}")
    return samples


def read_fasta(path):
    """Load a FASTA into {first header token: sequence} records."""
    records = {}
    name = None
    sequence = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(sequence)
                name = line[1:].split()[0]
                sequence = []
            elif name is None:
                raise ValueError(f"Invalid FASTA; sequence before header in {path}")
            else:
                sequence.append(line)
    if name is not None:
        records[name] = "".join(sequence)
    return records


def write_fasta(path, records):
    """Write sequence records with sample IDs as headers for downstream alignment."""
    with open(path, "w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")


# These three lists have different jobs: core_samples define gene prevalence,
# analysis_samples define retained sample columns/sequences, and all_samples
# identify which Roary columns are samples rather than metadata.
core_samples = read_samples(snakemake.input.core_samples)
analysis_samples = read_samples(snakemake.input.analysis_samples)
all_samples = read_samples(snakemake.input.all_samples)
all_sample_set = set(all_samples)

# The threshold is a fraction, e.g. 0.9 means present in at least 90% of the
# core-definition samples. Presence below is determined from family FASTAs,
# not CDS IDs in the Roary table.
threshold = float(snakemake.params.core_threshold)
if not 0 < threshold <= 1:
    raise ValueError("custom_alignment.core_threshold must be greater than 0 and at most 1")

# Require exactly one Roary column per sample. Repeated sample columns indicate
# duplicate sample IDs were supplied to Panaroo, so stop rather than silently
# choosing or merging their values.
with open(snakemake.input.roary, newline="") as handle:
    reader = csv.reader(handle)
    header = next(reader, None)
    if not header or "Gene" not in header:
        raise ValueError("Roary CSV must have a header containing a Gene column")
    rows = list(reader)

header_counts = Counter(header)
duplicate_sample_columns = sorted(
    sample
    for sample in all_sample_set
    if header_counts[sample] > 1
)
if duplicate_sample_columns:
    raise ValueError(
        "Roary has repeated sample columns; check sample_list for duplicate IDs: "
        + ", ".join(duplicate_sample_columns)
    )

gene_column = header.index("Gene")
sample_columns = {
    sample: [index for index, column in enumerate(header) if column == sample]
    for sample in all_samples
}
missing_roary_columns = sorted(
    sample
    for sample in set(core_samples + analysis_samples)
    if not sample_columns.get(sample)
)
if missing_roary_columns:
    raise ValueError(
        "Samples are not Roary CSV columns: " + ", ".join(missing_roary_columns)
    )

# Selection lists must refer to the full Panaroo dataset, otherwise a missing
# sample could be mistaken for a gene-absence result.
outside_dataset = sorted(set(core_samples + analysis_samples) - all_sample_set)
if outside_dataset:
    raise ValueError(
        "Core/analysis samples are not in sample_list: " + ", ".join(outside_dataset)
    )

    # The source directory contains one FASTA per family, named exactly after the
    # Roary Gene value (for example group_10758.aln.fas or group_10758.fasta). Keep only metadata columns
# and one output column for each selected analysis sample.
metadata_columns = [
    index for index, column in enumerate(header) if column not in all_sample_set
]
os.makedirs(os.path.dirname(snakemake.output.filtered_roary), exist_ok=True)
if os.path.exists(snakemake.output.gene_fastas):
    shutil.rmtree(snakemake.output.gene_fastas)
os.makedirs(snakemake.output.gene_fastas)

manifest = []
filtered_rows = []
source_dir = snakemake.input.gene_sequence_dir
file_prefix = snakemake.params.file_prefix

for row in rows:
    # Pad short rows so every indexed lookup has a predictable result.
    row += [""] * (len(header) - len(row))
    family = row[gene_column].strip()
    if not family:
        continue
    if os.path.basename(family) != family:
        raise ValueError(f"Invalid family name in Roary Gene column: {family}")

    # Panaroo FASTA headers start with the sample ID and then a semicolon plus
    # sequence metadata (e.g. GCF_...;0_8_108). Count records per sample so
    # duplicates are represented as absence for that sample, not as a reason to
    # discard the entire family.
    family_paths = [
        os.path.join(source_dir, f"{family}.aln.fas"),
        os.path.join(source_dir, f"{family}.fasta"),
    ]
    family_path = next((path for path in family_paths if os.path.isfile(path)), None)
    if family_path is None:
        raise FileNotFoundError(
            f"No Panaroo family FASTA for {family}; expected one of: "
            + ", ".join(family_paths)
        )
    family_records = read_fasta(family_path)
    records_by_sample = {}
    for sequence_id, sequence in family_records.items():
        sample = sequence_id.split(";", 1)[0]
        if sample in all_sample_set:
            records_by_sample.setdefault(sample, []).append((sequence_id, sequence))

    # A sample is considered present only if exactly one family sequence is
    # available. Zero copies and multiple copies both count as absent.
    present_core_samples = [
        sample for sample in core_samples if len(records_by_sample.get(sample, [])) == 1
    ]
    core_count = len(present_core_samples)
    is_core = core_count / len(core_samples) >= threshold

    # Construct the filtered Roary row. If a selected sample has zero or
    # multiple family records, emit an empty cell; otherwise preserve its
    # original Roary sequence identifier when available.
    filtered_row = [row[index] for index in metadata_columns]
    for sample in analysis_samples:
        sample_records = records_by_sample.get(sample, [])
        if len(sample_records) == 1:
            original_values = [row[index].strip() for index in sample_columns[sample]]
            cell = next((value for value in original_values if value), sample_records[0][0])
        else:
            cell = ""
        filtered_row.append(cell)
    filtered_rows.append(filtered_row)

    if not is_core:
        continue

    # Keep one FASTA record per analysis sample, using sample IDs as headers so
    # taxa names agree across genes during AMAS concatenation. Duplicate-copy
    # samples are simply omitted here; the finalizer later pads them with gaps.
    alignment_records = {
        sample: records_by_sample[sample][0][1]
        for sample in analysis_samples
        if len(records_by_sample.get(sample, [])) == 1
    }
    duplicate_core_samples = [
        sample for sample in core_samples if len(records_by_sample.get(sample, [])) > 1
    ]
    duplicate_analysis_samples = [
        sample for sample in analysis_samples if len(records_by_sample.get(sample, [])) > 1
    ]

    status = "selected" if alignment_records else "no_single_copy_analysis_sequences"
    reason = (
        ""
        if status == "selected"
        else "no analysis sample has exactly one family sequence"
    )
    file_id = f"{file_prefix}_{family}"
    if status == "selected":
        output_path = os.path.join(
            snakemake.output.gene_fastas, f"{file_id}.fasta"
        )
        write_fasta(output_path, alignment_records)
    manifest.append(
        (
            family,
            file_id,
            core_count,
            len(core_samples),
            len(alignment_records),
            ";".join(duplicate_core_samples),
            ";".join(duplicate_analysis_samples),
            status,
            reason,
        )
    )

# Write one sample column per selected sample even if the original Roary header
# repeated a genome name. This file records duplicate-bearing samples as absent.
with open(snakemake.output.filtered_roary, "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([header[index] for index in metadata_columns] + analysis_samples)
    writer.writerows(filtered_rows)

# The manifest retains core-threshold families that could not be aligned, with
# duplicate samples and reasons made explicit for review.
with open(snakemake.output.core_groups, "w") as handle:
    handle.write(
        "Gene\tFile_ID\tCore_samples_present\tCore_samples_total\t"
        "Analysis_sequences\tDuplicate_core_samples\tDuplicate_analysis_samples\t"
        "Status\tReason\n"
    )
    for row in manifest:
        handle.write("\t".join(map(str, row)) + "\n")

# Save the exact effective sample sets alongside this named derived run.
with open(snakemake.output.sample_log, "w") as handle:
    handle.write("selection\tsample\n")
    for sample in core_samples:
        handle.write(f"core_definition\t{sample}\n")
    for sample in analysis_samples:
        handle.write(f"analysis_subsample\t{sample}\n")