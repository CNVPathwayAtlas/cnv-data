"""Write the release of each data source of the latest data version.

Usage:
  python scripts/source_notes.py            print the GitHub release notes (Markdown)
  python scripts/source_notes.py --zenodo   write the notes and the version to .zenodo.json

The notes are made from data/latest/sources.json and version.txt.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCES_FILE = ROOT / "data" / "latest" / "sources.json"
VERSION_FILE = ROOT / "version.txt"
ZENODO_FILE = ROOT / ".zenodo.json"

ORPHADATA_PRODUCTS = {"definitions": "product 1, disorders", "phenotypes": "product 4, HPO phenotypes",
                      "prevalence": "product 9, prevalence"}


def source_rows(sources):
    """Return (source, release, URL) rows for each source."""
    rows = []
    for name, label in ORPHADATA_PRODUCTS.items():
        info = sources.get("orphadata", {}).get(name, {})
        rows.append((f"Orphadata ({label})", f"{info.get('release_date', '')[:10]}, version {info.get('version', '')}",
                     info.get("url", "")))
    hgnc = sources.get("hgnc", {})
    rows.append(("HGNC", hgnc.get("last_modified", ""), hgnc.get("url", "")))
    hpo = sources.get("hpo", {})
    rows.append(("HPO", hpo.get("release", "not recorded"), hpo.get("url", "")))
    for wp, info in sorted(sources.get("wikipathways", {}).items()):
        rows.append((f"WikiPathways {wp}", info.get("revision", ""), info.get("url", "")))
    return rows


def markdown(version, sources):
    lines = [f"Data version {version}. Sources retrieved on {sources.get('retrieved', '')[:10]}.", "",
             "| Source | Release | File |", "|---|---|---|"]
    lines += [f"| {source} | {release} | {url} |" for source, release, url in source_rows(sources)]
    lines += ["", "The SHA-256 checksum of each file is in `data/latest/sources.json`."]
    return "\n".join(lines) + "\n"


def zenodo_notes(version, sources):
    items = "".join(f"<li>{source}: {release}</li>" for source, release, _ in source_rows(sources))
    return (f"<p>Data version {version}. Sources retrieved on {sources.get('retrieved', '')[:10]}.</p>"
            f"<ul>{items}</ul><p>The URL and SHA-256 checksum of each file are in data/latest/sources.json.</p>")


def write_zenodo():
    version = VERSION_FILE.read_text().strip()
    sources = json.loads(SOURCES_FILE.read_text())
    metadata = json.loads(ZENODO_FILE.read_text())
    metadata.update({"version": version, "notes": zenodo_notes(version, sources)})
    ZENODO_FILE.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")


def main():
    if "--zenodo" in sys.argv:
        write_zenodo()
    else:
        print(markdown(VERSION_FILE.read_text().strip(), json.loads(SOURCES_FILE.read_text())), end="")


if __name__ == "__main__":
    main()
