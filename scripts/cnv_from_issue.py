"""Add or edit one CNV row in data/input/cnv_data.tsv from the body of a 'new-cnv' or 'edit-cnv' issue form.

Usage: ISSUE_BODY="..." python scripts/cnv_from_issue.py [report file]

A body with a 'CNV identifier' field edits that CNV: empty fields keep the current value, '-' makes a field empty.
A body without it adds a new CNV.
The script writes a Markdown report for the issue comment and exits with code 1 when the input has an error.
It writes the CNV identifier to $GITHUB_OUTPUT as 'key'.
"""
import os
import re
import sys

from cnv_input import COLUMNS, INPUT_FILE, cnv_key, load_input
from check_sources import check_row
from validate_input import validate

# Issue form label -> input column
LABELS = {
    "Locus": "locus",
    "Extra": "extra",
    "Pathway name": "pathway_name",
    "WikiPathways ID": "wikipathways_id",
    "ORPHAcodes": "orphacodes",
    "Region GRCh37 (hg19)": "region_hg19",
    "Region GRCh38 (hg38)": "region_hg38",
    "Genes": "genes",
    "PubMed IDs": "pubmed_id",
    "Description": "description",
}


def parse_body(body):
    """Return {label: value} from the '### Label' sections of an issue form body."""
    fields = {}
    for section in re.split(r"^### ", body.replace("\r\n", "\n"), flags=re.M)[1:]:
        label, _, value = section.partition("\n")
        value = value.strip()
        fields[label.strip()] = "" if value == "_No response_" else value
    return fields


def normalize(column, value):
    """Return the value in the format of the input table."""
    if column in ("genes", "pubmed_id"):
        return ";".join(re.split(r"[\s,;]+", value)).strip(";")
    if column == "orphacodes":
        items = (re.sub(r"\s*\(\s*", "(", item.strip()) for item in re.split(r"[\n;]+", value))
        return ";".join(item for item in items if item)
    if column.startswith("region_"):
        return re.sub(r"[\s,]", "", value)
    return " ".join(value.split())


def write_report(path, title, row, old, errors, warnings, notes=()):
    lines = [f"## {title}", "", "| Column | Value |", "|---|---|"]
    for c in COLUMNS:
        value = row[c] if old is None or old[c] == row[c] else f"~~{old[c]}~~ → **{row[c] or '(empty)'}**"
        lines.append(f"| `{c}` | {value.replace('|', '&#124;')} |")
    lines += [""] + [f"- :x: {m}" for m in errors] + [f"- :warning: {m}" for m in warnings]
    if notes:
        lines += ["", "Found in the sources:"] + [f"- {m}" for m in notes]
    lines += ["", "The check failed. Edit the issue to correct the errors." if errors
              else "The check passed. A pull request changes `data/input/cnv_data.tsv`."]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    report_path = sys.argv[1] if len(sys.argv) > 1 else "report.md"
    fields = parse_body(os.environ["ISSUE_BODY"])
    target = fields.get("CNV identifier", "").strip()
    df = load_input()
    keys = list(df.apply(cnv_key, axis=1))
    lines = INPUT_FILE.read_text(encoding="utf-8").splitlines()

    if target:
        if target not in keys:
            row = {c: "" for c in COLUMNS}
            write_report(report_path, f"Edit CNV `{target}`", row, None,
                         [f"{target}: CNV identifier is not in the input table"], [])
            sys.exit(1)
        index = keys.index(target)
        old = df.iloc[index].to_dict()
        row = dict(old)
        for label, column in LABELS.items():
            value = fields.get(label, "").strip()
            if value:
                row[column] = "" if value == "-" else normalize(column, value)
        lines[index + 1] = "\t".join(row[c] for c in COLUMNS)
        title = f"Edit CNV `{target}`"
    else:
        old, index = None, len(df)
        row = {column: normalize(column, fields.get(label, "")) for label, column in LABELS.items()}
        lines.append("\t".join(row[c] for c in COLUMNS))
        title = f"Add CNV `{cnv_key(row)}`"

    INPUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    key = cnv_key(row)
    df = load_input()
    errors, warnings = ([m for m in messages if m.startswith(f"{key}:")] for messages in validate(df))
    if len(df) <= index or df.iloc[index].to_dict() != row:
        errors.append(f"{key}: the row is not read back with the same values. Remove quote characters and try again.")
    source_errors, source_warnings, notes = check_row(row)
    errors += source_errors
    warnings += source_warnings
    write_report(report_path, title, row, old, errors, warnings, notes)

    if "GITHUB_OUTPUT" in os.environ:
        with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as f:
            f.write(f"key={key}\naction={'Edit' if target else 'Add'}\n")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
