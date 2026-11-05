def read_fasta(path):
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


max_gap = float(snakemake.params.gap_seq_max)
records = read_fasta(snakemake.input.alignment)

# Sequences are written ungapped so MAFFT realigns them from scratch.
kept = {}
removed = []
for name, sequence in records.items():
    gap_fraction = sequence.count("-") / len(sequence) if sequence else 1.0
    if gap_fraction <= max_gap:
        kept[name] = sequence.replace("-", "")
    else:
        removed.append((name, gap_fraction))

with open(snakemake.output.removed, "w") as out:
    out.write("gene\tsequence\tgap_fraction\n")
    for name, gap_fraction in removed:
        out.write(f"{snakemake.wildcards.gene}\t{name}\t{gap_fraction:.4f}\n")

if len(kept) < 2:
    raise ValueError(
        f"Fewer than 2 sequences with <= {max_gap:.0%} gaps remain for {snakemake.wildcards.gene}"
    )

with open(snakemake.output.fasta, "w") as out:
    for name, sequence in kept.items():
        out.write(f">{name}\n{sequence}\n")
