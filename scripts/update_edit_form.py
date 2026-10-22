"""Write the CNV identifiers of data/input/cnv_data.tsv as the options of the 'Edit a CNV' issue form.
The workflow update-edit-form.yml runs this script when the input table changes.

Usage: python scripts/update_edit_form.py
"""
from pathlib import Path

from cnv_input import cnv_key, load_input

FORM = Path(__file__).resolve().parent.parent / ".github" / "ISSUE_TEMPLATE" / "edit-cnv.yml"
START, END = "        # CNV options start\n", "        # CNV options end\n"


def main():
    text = FORM.read_text(encoding="utf-8")
    head, rest = text.split(START)
    _, tail = rest.split(END)
    options = "".join(f'        - "{key}"\n' for key in load_input().apply(cnv_key, axis=1))
    FORM.write_text(head + START + options + END + tail, encoding="utf-8")


if __name__ == "__main__":
    main()
