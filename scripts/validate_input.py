"""Check the curated CNV input table. Exit with code 1 when an error is found.

Errors stop the data update. Warnings are printed and must be checked by a curator.
"""
import re
import sys

from cnv_input import COLUMNS, ORPHACODE_PATTERN, load_input, parse_region, split_list, cnv_key

# Region length ratio between the two builds that causes a warning
MAX_BUILD_LENGTH_RATIO = 1.25


def validate(df):
    """Return (errors, warnings) as two lists of strings."""
    errors, warnings = [], []

    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        return [f"Missing columns: {', '.join(missing)}"], []

    keys = df.apply(cnv_key, axis=1)
    for key in keys[keys.duplicated()].unique():
        errors.append(f"{key}: CNV identifier is not unique")

    for _, row in df.iterrows():
        key = cnv_key(row)
        if not row["locus"]:
            errors.append(f"{key}: locus is empty")

        if row["wikipathways_id"] and not re.fullmatch(r"WP\d+", row["wikipathways_id"]):
            errors.append(f"{key}: incorrect WikiPathways ID '{row['wikipathways_id']}'")

        for item in split_list(row["orphacodes"]):
            if not ORPHACODE_PATTERN.match(item):
                errors.append(f"{key}: incorrect orphacode item '{item}'")

        for pmid in split_list(row["pubmed_id"]):
            if not pmid.isdigit():
                errors.append(f"{key}: incorrect PubMed ID '{pmid}'")

        if not split_list(row["genes"]):
            warnings.append(f"{key}: gene list is empty")

        regions = {}
        for build in ("hg19", "hg38"):
            value = row[f"region_{build}"]
            chrom, start, end = parse_region(value)
            if chrom is None:
                errors.append(f"{key}: incorrect region_{build} '{value}'")
                continue
            if start >= end:
                errors.append(f"{key}: region_{build} start is not less than end")
                continue
            if start == 0:
                warnings.append(f"{key}: region_{build} start is 0 (coordinates start at 1)")
            regions[build] = (chrom, end - start)

        if len(regions) == 2:
            (chrom19, len19), (chrom38, len38) = regions["hg19"], regions["hg38"]
            if chrom19 != chrom38:
                errors.append(f"{key}: chromosome is different between hg19 and hg38")
            ratio = max(len19, len38) / min(len19, len38)
            if ratio > MAX_BUILD_LENGTH_RATIO:
                warnings.append(
                    f"{key}: region length hg19 {len19 / 1e6:.2f} Mb, hg38 {len38 / 1e6:.2f} Mb "
                    f"(ratio {ratio:.2f})"
                )

    return errors, warnings


def report(errors, warnings):
    for message in warnings:
        print(f"[WARN] {message}")
    for message in errors:
        print(f"[ERROR] {message}")


def main():
    errors, warnings = validate(load_input())
    report(errors, warnings)
    print(f"Validation: {len(errors)} errors, {len(warnings)} warnings")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
