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


def count_variable_columns(sequences):
    """Count columns with at least two distinct states, ignoring gaps and N."""
    length = {len(sequence) for sequence in sequences}
    if len(length) != 1:
        raise ValueError("Sequences in the Panaroo alignment have unequal lengths")
    variable = 0
    for column in zip(*sequences):
        states = {base for base in (c.upper() for c in column) if base not in "-?N"}
        if len(states) > 1:
            variable += 1
    return variable


with open(snakemake.input.manifest) as handle:
    alignments = dict(line.rstrip("\n").split("\t") for line in handle if line.strip())

min_sites = int(snakemake.params.min_variable_sites)
variable, constant = [], []
for gene, path in alignments.items():
    records = read_fasta(path)
    n_sites = count_variable_columns(list(records.values())) if records else 0
    (variable if n_sites >= min_sites else constant).append(gene)

with open(snakemake.output.variable, "w") as out:
    out.write("".join(f"{gene}\n" for gene in variable))
with open(snakemake.output.constant, "w") as out:
    out.write("".join(f"{gene}\n" for gene in constant))
