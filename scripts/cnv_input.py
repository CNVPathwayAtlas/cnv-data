"""Read and parse the curated CNV input table (data/input/cnv_data.tsv).

The cnv-website build imports this module, so both repositories parse the input the same way.
"""
import re
from pathlib import Path

import pandas as pd

INPUT_FILE = Path(__file__).resolve().parent.parent / "data" / "input" / "cnv_data.tsv"

COLUMNS = ["locus", "extra", "pathway_name", "wikipathways_id", "orphacodes",
           "region_hg19", "region_hg38", "genes", "pubmed_id", "description"]

REGION_PATTERN = re.compile(r"^[cC]hr(\d+|[XY]):(\d+)-(\d+)$")
ORPHACODE_PATTERN = re.compile(r"^(\d+)(?:\(([^)]+)\))?$")


def load_input(path=INPUT_FILE):
    """Return the input table as a DataFrame of strings. Empty cells are empty strings."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return df.apply(lambda col: col.str.strip())


def split_list(value):
    """Split a list cell on ';' or ','. Return a list of non-empty items."""
    return [item.strip() for item in re.split(r"[;,]", value or "") if item.strip()]


def parse_orphacodes(value):
    """Return a list of (orphacode, cause) tuples. Cause is None when not given."""
    result = []
    for item in split_list(value):
        match = ORPHACODE_PATTERN.match(item)
        if match:
            result.append((match.group(1), match.group(2)))
    return result


def parse_region(value):
    """Return (chromosome, start, end) for 'ChrN:start-end', or (None, None, None)."""
    match = REGION_PATTERN.match(value or "")
    if not match:
        return None, None, None
    return match.group(1).upper(), int(match.group(2)), int(match.group(3))


def cnv_key(row):
    """Return the CNV identifier: locus, plus '-extra' when extra is given."""
    return f"{row['locus']}-{row['extra']}" if row["extra"] else row["locus"]
