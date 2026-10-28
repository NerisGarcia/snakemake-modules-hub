#!/usr/bin/env python3
"""Summarise fastp JSON reports into a CSV table (one row per sample)."""
import argparse, csv, json, os, sys

COLS = ["sample", "reads_before", "reads_after", "gc_before", "gc_after",
        "passed_filter_reads", "low_quality_reads", "too_many_N_reads",
        "adapter_dimer_reads", "too_short_reads", "too_long_reads",
        "duplication_rate", "insert_size_peak",
        "adapter_trimmed_reads", "adapter_trimmed_bases"]


def row(path):
    d = json.load(open(path))
    s, f = d["summary"], d.get("filtering_result", {})
    a = d.get("adapter_cutting", {})
    sample = os.path.basename(path).removesuffix(".json").removesuffix(".fastp")
    return {
        "sample": sample,
        "reads_before": s["before_filtering"]["total_reads"],
        "reads_after": s["after_filtering"]["total_reads"],
        "gc_before": s["before_filtering"]["gc_content"],
        "gc_after": s["after_filtering"]["gc_content"],
        **{k: f.get(k, "") for k in COLS[5:11]},
        "duplication_rate": d.get("duplication", {}).get("rate", ""),
        "insert_size_peak": d.get("insert_size", {}).get("peak", ""),
        "adapter_trimmed_reads": a.get("adapter_trimmed_reads", ""),
        "adapter_trimmed_bases": a.get("adapter_trimmed_bases", ""),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("json", nargs="+", help="fastp JSON file(s)")
    p.add_argument("-o", "--output", help="output CSV (default: stdout)")
    p.add_argument("-d", "--delimiter", default=",", help="field delimiter (use '\\t' for TSV)")
    a = p.parse_args()
    delim = "\t" if a.delimiter == "\\t" else a.delimiter
    out = open(a.output, "w", newline="") if a.output else sys.stdout
    w = csv.DictWriter(out, fieldnames=COLS, delimiter=delim)
    w.writeheader()
    for j in a.json:
        w.writerow(row(j))


if __name__ == "__main__":
    main()
