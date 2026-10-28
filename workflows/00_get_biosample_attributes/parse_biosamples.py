import shutil
import xml.etree.ElementTree as ET

xml_files = snakemake.input
output_file = str(snakemake.output[0])
tmp_dir = str(snakemake.params.tmp_dir)

CORE_FIELDS = [
    "Biosample_Accession",
    "Bioproject_Accession",
    "Submission_Date",
    "Taxonomy_ID",
    "Taxonomy_Name",
    "Owner",
]

records = []
attribute_names = []
seen_attrs = set()

for xml_file in xml_files:
    biosample = ET.parse(xml_file).find(".//BioSample")
    if biosample is None:
        continue

    record = {
        "Biosample_Accession": biosample.get("accession", ""),
        "Submission_Date": biosample.get("submission_date", ""),
    }

    bioproject_acc = ""
    for link in biosample.findall("./Links/Link"):
        if link.get("target") == "bioproject":
            bioproject_acc = link.get("label", "")
            break
    record["Bioproject_Accession"] = bioproject_acc

    organism = biosample.find("./Description/Organism")
    record["Taxonomy_ID"] = organism.get("taxonomy_id", "") if organism is not None else ""
    record["Taxonomy_Name"] = organism.get("taxonomy_name", "") if organism is not None else ""

    owner_name = biosample.find("./Owner/Name")
    record["Owner"] = (owner_name.text or "") if owner_name is not None else ""

    for attr in biosample.findall("./Attributes/Attribute"):
        name = attr.get("attribute_name")
        if not name:
            continue
        if name not in seen_attrs:
            seen_attrs.add(name)
            attribute_names.append(name)
        record[name] = attr.text or ""

    records.append(record)

header = CORE_FIELDS + attribute_names

with open(output_file, "w") as out:
    out.write("\t".join(header) + "\n")
    for record in records:
        out.write("\t".join(record.get(col, "") for col in header) + "\n")

shutil.rmtree(tmp_dir)
