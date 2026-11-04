# This script is run once per selected core gene (rule prepare_gene_alignment)
# after MAFFT has aligned that gene's sequences for the analysis samples. Its
# job is to turn a raw MAFFT output FASTA into a "finalized" per-gene FASTA
# that is safe to concatenate across genes with AMAS: every gene FASTA must
# contain exactly the same set of sample headers, in the same alignment
# length, so that columns from different genes line up sample-for-sample.


def read_samples(path):
    """Read the full analysis sample list (snakemake.input.samples).

    This defines the complete, ordered set of samples that must appear in the
    finalized alignment, including samples that MAFFT's input FASTA did not
    actually contain (see "records.get(sample, ...)" below).
    """
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
    """Parse a FASTA file into {header_id: sequence}.

    Only the first whitespace-delimited token of each header line is kept as
    the ID (mirrors how select_core_groups.py and MAFFT itself identify
    sequences), and sequence lines are concatenated/stripped of newlines.
    """
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
    """Write {header_id: sequence} back out as a standard single-line FASTA."""
    with open(path, "w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")


# samples: the authoritative list of analysis sample IDs/order for this run.
# aligned: the MAFFT-aligned sequences for this one gene, keyed by whatever
#          header MAFFT wrote out (which may be decorated, see below).
samples = read_samples(snakemake.input.samples)
aligned = read_fasta(snakemake.input.alignment)

# MAFFT was run with --adjustdirection, which reverse-complements sequences it
# believes are on the wrong strand and marks them by prefixing "_R_" to the
# header. Strip that marker so the header goes back to being a plain sample
# ID; otherwise a reoriented sequence would silently fail to match any sample
# in `samples` and would be dropped as "unexpected".
normalized = {}
for name, sequence in aligned.items():
    if name.endswith("_R_"):
        name = name[:-3]
    elif name.startswith("_R_"):
        name = name[3:]
    normalized[name] = sequence

# Sanity checks on the normalized MAFFT output before trusting it. These are
# defensive checks: none of them should normally trigger, and each one means
# a different kind of problem happened before this script ran.

# 1) MAFFT must have produced at least one sequence. An empty `normalized`
#    means the input FASTA from select_core_groups.py was itself empty (0
#    analysis samples had this gene) — that gene should not have reached this
#    rule at all, so this points to a bug in the upstream gene-selection step.
if not normalized:
    raise ValueError(f"MAFFT returned no sequences for {snakemake.wildcards.gene}")

# 2) Every header remaining after the _R_ strip must be one of the known
#    analysis sample IDs. MAFFT does not invent or rename sequences, so a
#    header ending up here that isn't in `samples` means the ID round-trip
#    broke somewhere (e.g. the input FASTA used a different ID scheme, or the
#    _R_ stripping above was insufficient for some exotic header). Fail loudly
#    instead of silently writing a bogus sample into the alignment.
if not set(normalized).issubset(samples):
    unexpected = sorted(set(normalized) - set(samples))
    raise ValueError("Unexpected sequence IDs after cleaning: " + ", ".join(unexpected))

# 3) MAFFT's whole job is to make all input sequences the same length by
#    inserting "-" gaps, so a correctly-finished alignment always has exactly
#    one distinct length across all sequences. If `lengths` has more than one
#    value here, MAFFT did not actually align this gene (e.g. it crashed
#    partway, or the input was pre-aligned/malformed), so the output cannot be
#    trusted and must not be passed on to AMAS concatenation.
lengths = {len(sequence) for sequence in normalized.values()}
if len(lengths) != 1:
    raise ValueError(f"MAFFT produced unequal alignment lengths for {snakemake.wildcards.gene}")

# Build the final per-gene record set. AMAS concatenation (rule
# amas_concatenate) glues every gene's FASTA together column-wise into one
# big supermatrix, matching records purely by sample header — it does not
# know or care which samples a given gene was actually aligned for. So each
# gene's FASTA must contain the exact same headers, in the exact same order,
# all the exact same length, or the supermatrix would be corrupted (columns
# would not line up, or AMAS would error on mismatched sample sets).
#
# `normalized` only has entries for samples where select_core_groups.py found
# exactly one copy of this gene (see that script's single-copy rule). Any
# analysis sample NOT in `normalized` simply never had this gene included in
# the input to MAFFT, so it has no aligned sequence to contribute here.
#
# Example: samples = [A, B, C], but this gene's family FASTA only had A and B
# (C was absent or had 2+ copies). After MAFFT, `normalized` = {A: "ACGT",
# B: "ACGA"} (alignment_length = 4). The loop below then emits:
#   A: "ACGT", B: "ACGA", C: "----"
# so sample C is still represented as an all-gap ("missing data") row for
# this gene, keeping row/column structure identical across every gene file.
alignment_length = lengths.pop()
records = {
    sample: normalized.get(sample, "-" * alignment_length)
    for sample in samples
}
write_fasta(snakemake.output.alignment, records)