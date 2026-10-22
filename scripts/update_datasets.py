"""Download Orphadata, HGNC, WikiPathways and HPO data and build the processed CNV datasets.

The script does these steps:
1. Validate the curated input (data/input/cnv_data.tsv).
2. Read the release identifiers of each source and compare them with the last run.
   When no source has a new release, stop (use --force to continue).
3. Download and process the sources.
4. When the processed data is different from the latest data, write a new data version.
"""
import argparse
import datetime
import hashlib
import json
import re
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

from cnv_input import INPUT_FILE, load_input, parse_orphacodes, split_list, cnv_key
from source_notes import write_zenodo
from validate_input import validate, report

ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = ROOT / "data" / "processed"
LATEST_DIR = ROOT / "data" / "latest"
PATHWAYS_DIR = ROOT / "data" / "pathways"
VERSION_FILE = ROOT / "version.txt"

ORPHADATA_FILES = {
    "definitions": "https://www.orphadata.com/data/xml/en_product1.xml",
    "phenotypes": "https://www.orphadata.com/data/xml/en_product4.xml",
    "prevalence": "https://www.orphadata.com/data/xml/en_product9_prev.xml",
}
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
WIKIPATHWAYS_URL = "https://www.wikipathways.org/wikipathways-assets/pathways/{wp}/{wp}.{ext}"
HPO_RELEASE_URL = "https://api.github.com/repos/obophenotype/human-phenotype-ontology/releases/latest"
HPO_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/download/{release}/hp.json"
HPO_ROOT = "HP:0000118"  # Phenotypic abnormality: its direct children are the top-level categories
# Categories that are not the primary category of a term when the term has an other category
HPO_SECONDARY = {"HP:0025142", "HP:0002664"}  # Constitutional symptom, Neoplasm

TIMEOUT = 300

HGNC_COLUMNS = ["hgnc_id", "symbol", "name", "entrez_id", "ensembl_gene_id", "uniprot_ids",
                "omim_id", "locus_group", "prev_symbol", "alias_symbol"]

FREQUENCIES = {
    "Obligate (100%)": "Phenotypes_Obligate(100%)",
    "Very frequent (99-80%)": "Phenotypes_Very_frequent(99-80%)",
    "Frequent (79-30%)": "Phenotypes_Frequent(79-30%)",
    "Occasional (29-5%)": "Phenotypes_Occasional(29-5%)",
    "Very rare (<4-1%)": "Phenotypes_Very_rare(<4-1%)",
    "Excluded (0%)": "Phenotypes_Excluded(0%)",
}

# GPML node types that are genes or gene products
GENE_NODE_TYPES = {"GeneProduct", "Protein", "Rna"}

# GPML xref database name -> HGNC column used for the mapping
XREF_TO_HGNC_COLUMN = {
    "ensembl": "ensembl_gene_id",
    "entrez gene": "entrez_id",
    "ncbi gene": "entrez_id",
    "hgnc": "symbol",
    "hgnc accession": "hgnc_id",
    "uniprot-trembl": "uniprot_ids",
    "uniprot-swissprot": "uniprot_ids",
    "uniprot": "uniprot_ids",
}


def download(url, path):
    print(f"Downloading {url}")
    with requests.get(url, stream=True, timeout=TIMEOUT) as response:
        response.raise_for_status()
        with open(path, "wb") as f:
            for chunk in response.iter_content(chunk_size=1 << 20):
                f.write(chunk)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


# Release detection

def orphadata_release(url):
    """Return the JDBOR date and version from the start of an Orphadata XML file."""
    response = requests.get(url, headers={"Range": "bytes=0-1023"}, timeout=TIMEOUT)
    response.raise_for_status()
    match = re.search(r'<JDBOR date="([^"]+)" version="([^"]+)"', response.text)
    if not match:
        raise ValueError(f"No JDBOR release information in {url}")
    return {"release_date": match.group(1), "version": match.group(2)}


def hgnc_release():
    response = requests.head(HGNC_URL, timeout=TIMEOUT)
    response.raise_for_status()
    return {"last_modified": response.headers.get("Last-Modified", ""),
            "etag": response.headers.get("ETag", "").strip('"')}


def wikipathways_release(gpml_path):
    """Return the revision and last-modified attributes of a GPML file."""
    root = ET.parse(gpml_path).getroot()
    return {"revision": root.get("Version", ""), "last_modified": root.get("Last-Modified", "")}


def hpo_release():
    """Return the tag and date of the latest HPO release on GitHub."""
    response = requests.get(HPO_RELEASE_URL, timeout=TIMEOUT)
    response.raise_for_status()
    data = response.json()
    return {"release": data["tag_name"], "published": data["published_at"]}


def release_ids(sources):
    """Return only the release identifiers of a sources record, for comparison."""
    return {"orphadata": {k: v.get("release_date") for k, v in sources.get("orphadata", {}).items()},
            "hgnc": sources.get("hgnc", {}).get("etag"),
            "wikipathways": {k: v.get("revision") for k, v in sources.get("wikipathways", {}).items()},
            "hpo": sources.get("hpo", {}).get("release")}


def load_previous_sources():
    path = LATEST_DIR / "sources.json"
    return json.loads(path.read_text()) if path.exists() else {}


# Orphadata parsing

def parse_definitions(path, orphacodes):
    """Return {orphacode: {name, definition, omim}} from Orphadata product 1."""
    data = {}
    for disorder in ET.parse(path).getroot().iter("Disorder"):
        code = disorder.findtext("OrphaCode")
        if code not in orphacodes:
            continue
        definition = ""
        for section in disorder.findall("SummaryInformationList/SummaryInformation[@lang='en']"
                                        "/TextSectionList/TextSection"):
            if section.findtext("TextSectionType/Name[@lang='en']") == "Definition":
                definition = (section.findtext("Contents") or "").strip()
                break
        omim = [ref.findtext("Reference") for ref in disorder.findall("ExternalReferenceList/ExternalReference")
                if ref.findtext("Source") == "OMIM" and ref.findtext("Reference")]
        data[code] = {"name": disorder.findtext("Name[@lang='en']") or "",
                      "definition": definition, "omim": omim}
    return data


def parse_phenotypes(path, orphacodes):
    """Return {orphacode: {frequency column: ['term (HP:id)', ...]}} from Orphadata product 4."""
    phenotypes = defaultdict(lambda: defaultdict(list))
    for disorder in ET.parse(path).getroot().iter("Disorder"):
        code = disorder.findtext("OrphaCode")
        if code not in orphacodes:
            continue
        for assoc in disorder.findall("HPODisorderAssociationList/HPODisorderAssociation"):
            column = FREQUENCIES.get(assoc.findtext("HPOFrequency/Name[@lang='en']"))
            term = assoc.findtext("HPO/HPOTerm")
            hpo_id = assoc.findtext("HPO/HPOId")
            if column and term and hpo_id:
                phenotypes[code][column].append(f"{term} ({hpo_id})")
    return phenotypes


def parse_prevalence(path, orphacodes):
    """Return {orphacode: (class, [pmid, ...])} with the first worldwide prevalence class."""
    prevalence = {code: ("Unknown", []) for code in orphacodes}
    for disorder in ET.parse(path).getroot().iter("Disorder"):
        code = disorder.findtext("OrphaCode")
        if code not in orphacodes:
            continue
        for prev in disorder.findall("PrevalenceList/Prevalence"):
            geo = prev.findtext("PrevalenceGeographic/Name[@lang='en']")
            class_ = prev.findtext("PrevalenceClass/Name[@lang='en']")
            if geo == "Worldwide" and class_:
                pmids = re.findall(r"(\d+)\[PMID\]", prev.findtext("Source") or "")
                prevalence[code] = (class_, pmids)
                break
    return prevalence


def build_orphadata_table(tmp, orphacodes):
    definitions = parse_definitions(tmp / "definitions.xml", orphacodes)
    phenotypes = parse_phenotypes(tmp / "phenotypes.xml", orphacodes)
    prevalence = parse_prevalence(tmp / "prevalence.xml", orphacodes)

    rows = []
    for code in sorted(orphacodes, key=int):
        if code not in definitions:
            print(f"[WARN] ORPHA:{code} is not in the Orphadata release")
        info = definitions.get(code, {"name": "", "definition": "", "omim": []})
        prev_class, prev_pmids = prevalence[code]
        row = {"OrphaCode": code, "Name": info["name"], "Definition": info["definition"]}
        row.update({column: "; ".join(phenotypes[code][column]) for column in FREQUENCIES.values()})
        row.update({"Prevalence": prev_class, "Prevalence_pmid": ", ".join(prev_pmids),
                    "OMIM": "; ".join(info["omim"])})
        rows.append(row)
    return pd.DataFrame(rows)


# HPO categories

def load_hpo(path):
    """Return ({id: label}, {id: [parent ids]}, {old id: current id}) from an HPO obographs JSON file."""
    graph = json.loads(Path(path).read_text())["graphs"][0]
    short = lambda iri: iri.rsplit("/", 1)[-1].replace("_", ":")
    labels, parents, current = {}, defaultdict(list), {}
    for node in graph["nodes"]:
        if not node["id"].startswith("http://purl.obolibrary.org/obo/HP_"):
            continue
        hpo_id = short(node["id"])
        meta = node.get("meta", {})
        for prop in meta.get("basicPropertyValues", []):
            if prop["pred"].endswith("IAO_0100001"):  # replaced_by
                current[hpo_id] = prop["val"]
            elif prop["pred"].endswith("#hasAlternativeId"):
                current[prop["val"]] = hpo_id
        if not meta.get("deprecated"):
            labels[hpo_id] = node.get("lbl", "")
    for edge in graph["edges"]:
        if edge["pred"] == "is_a":
            parents[short(edge["sub"])].append(short(edge["obj"]))
    return labels, parents, current


def build_hpo_categories(orphadata, labels, parents, current, release):
    """Return the top-level HPO categories of the phenotypes in the Orphadata table.

    The categories of a term are its ancestors that are direct children of HP:0000118.
    The primary category is the nearest category. Constitutional symptom and Neoplasm are primary
    only when the term has no other category. A tie is resolved by the lowest HPO ID.
    """
    top = {c for c, ps in parents.items() if HPO_ROOT in ps and c in labels}
    hpo_ids = sorted({i for column in FREQUENCIES.values() for value in orphadata[column]
                      for i in re.findall(r"HP:\d{7}", value)})

    def distances(hpo_id):
        found, frontier, depth, seen = {}, [hpo_id], 0, set()
        while frontier:
            for term in frontier:
                if term in top and term not in found:
                    found[term] = depth
            seen.update(frontier)
            frontier = [p for t in frontier for p in parents.get(t, []) if p not in seen]
            depth += 1
        return found

    terms, not_mapped = {}, []
    for hpo_id in hpo_ids:
        term = current.get(hpo_id, hpo_id)
        found = distances(term) if term in labels else {}
        if not found:
            not_mapped.append(hpo_id)
            continue
        primary = min(found, key=lambda c: (c in HPO_SECONDARY, found[c], c))
        terms[hpo_id] = {"current_id": term, "label": labels[term], "primary": primary,
                         "categories": sorted(found)}
    if not_mapped:
        print(f"[WARN] {len(not_mapped)} HPO IDs are not in HPO {release} or have no category: {', '.join(not_mapped)}")

    primaries = sorted({t["primary"] for t in terms.values()})
    categories = {c: {"label": labels[c], "phenotypes": sorted(i for i, t in terms.items() if t["primary"] == c)}
                  for c in primaries}
    return {"version": 2, "hpo_release": release, "categories": categories,
            "phenotypeToCategory": {i: t["primary"] for i, t in terms.items()},
            "terms": terms, "not_mapped": not_mapped}


# HGNC mapping

def load_hgnc(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
    return df[HGNC_COLUMNS]


def resolve_symbols(df_input, hgnc):
    """Map each input gene symbol to an approved HGNC symbol.

    The order is: approved symbol, then a unique previous symbol, then a unique alias symbol.
    """
    approved = dict(zip(hgnc["symbol"], hgnc["hgnc_id"]))
    lookup = {"previous": defaultdict(set), "alias": defaultdict(set)}
    for kind, column in (("previous", "prev_symbol"), ("alias", "alias_symbol")):
        for symbols, approved_symbol in zip(hgnc[column], hgnc["symbol"]):
            for symbol in filter(None, symbols.split("|")):
                lookup[kind][symbol].add(approved_symbol)

    rows = []
    for _, cnv in df_input.iterrows():
        for symbol in split_list(cnv["genes"]):
            match_type, approved_symbol = "not_found", ""
            if symbol in approved:
                match_type, approved_symbol = "approved", symbol
            else:
                for kind in ("previous", "alias"):
                    candidates = lookup[kind][symbol]
                    if len(candidates) == 1:
                        match_type, approved_symbol = kind, next(iter(candidates))
                        break
                    if len(candidates) > 1:
                        match_type = f"{kind}_ambiguous"
                        break
            rows.append({"cnv": cnv_key(cnv), "input_symbol": symbol, "match_type": match_type,
                         "symbol": approved_symbol, "hgnc_id": approved.get(approved_symbol, "")})
    return pd.DataFrame(rows)


# WikiPathways

def build_xref_index(hgnc):
    """Return {hgnc column: {identifier: set(hgnc_id)}} for the xref columns used in GPML."""
    index = {}
    for column in set(XREF_TO_HGNC_COLUMN.values()):
        mapping = defaultdict(set)
        for value, hgnc_id in zip(hgnc[column], hgnc["hgnc_id"]):
            for identifier in filter(None, value.split("|")):
                mapping[identifier].add(hgnc_id)
        index[column] = mapping
    return index


def attr(element, *names):
    """Return the first attribute value found. GPML 2013a and GPML 2021 use different names."""
    for name in names:
        if element.get(name):
            return element.get(name)
    return ""


def parse_pathway_genes(wp_id, gpml_path, index, hgnc_by_id):
    """Return one row per gene node in the pathway, with the HGNC gene when the xref maps to one gene."""
    rows = []
    for node in ET.parse(gpml_path).getroot().iter():
        if not node.tag.endswith("DataNode") or attr(node, "Type", "type") not in GENE_NODE_TYPES:
            continue
        xref = next((child for child in node if child.tag.endswith("Xref")), None)
        database = attr(xref, "Database", "dataSource") if xref is not None else ""
        identifier = attr(xref, "ID", "identifier") if xref is not None else ""
        column = XREF_TO_HGNC_COLUMN.get(database.lower())
        if column == "hgnc_id" and identifier and not identifier.startswith("HGNC:"):
            identifier = f"HGNC:{identifier}"
        label = attr(node, "TextLabel", "textLabel")
        hits = index[column].get(identifier, set()) if column else set()
        hgnc_id, mapped_by = (next(iter(hits)), "xref") if len(hits) == 1 else ("", "")
        # When the xref does not map to one gene, use the label only when it is an approved symbol
        if not hgnc_id and len(index["symbol"].get(label, set())) == 1:
            hgnc_id, mapped_by = next(iter(index["symbol"][label])), "label"
        gene = hgnc_by_id.get(hgnc_id, {})
        rows.append({"wikipathways_id": wp_id, "label": label, "xref_database": database,
                     "xref_id": identifier, "hgnc_id": hgnc_id, "symbol": gene.get("symbol", ""),
                     "locus_group": gene.get("locus_group", ""), "mapped_by": mapped_by})
    return rows


# Output

def write_outputs(out_dir, tables, pathway_files, sources):
    for name, data in tables.items():
        if name.endswith(".json"):
            (out_dir / name).write_text(json.dumps(data, indent=2) + "\n")
            continue
        sep = "," if name.endswith(".csv") else "\t"
        data.to_csv(out_dir / name, sep=sep, index=False, lineterminator="\n")
    (out_dir / "pathways").mkdir()
    for path in pathway_files:
        shutil.copy(path, out_dir / "pathways" / path.name)
    (out_dir / "sources.json").write_text(json.dumps(sources, indent=2) + "\n")


def files_equal(new_dir, tables):
    """Return True when the new tables and pathway files are equal to the latest data."""
    for name in tables:
        stem, ext = name.rsplit(".", 1)
        latest = sorted(LATEST_DIR.glob(f"{stem}_*.{ext}"))
        if not latest or latest[-1].read_bytes() != (new_dir / name).read_bytes():
            return False
    new_pathways = {p.name: p.read_bytes() for p in (new_dir / "pathways").iterdir()}
    old_pathways = {p.name: p.read_bytes() for p in PATHWAYS_DIR.glob("*")}
    return new_pathways == old_pathways


def publish(new_dir, tables, today):
    """Write the new version to data/processed, data/latest, data/pathways and version.txt."""
    for directory in (PROCESSED_DIR, LATEST_DIR):
        directory.mkdir(parents=True, exist_ok=True)
        for path in directory.iterdir():
            if path.is_file() or path.is_symlink():
                path.unlink()

    for name in list(tables) + ["sources.json"]:
        stem, ext = name.rsplit(".", 1)
        target = PROCESSED_DIR / f"{stem}_{today}.{ext}"
        shutil.copy(new_dir / name, target)
        link = LATEST_DIR / ("sources.json" if name == "sources.json" else target.name)
        link.symlink_to(Path("..") / "processed" / target.name)

    shutil.rmtree(PATHWAYS_DIR, ignore_errors=True)
    shutil.copytree(new_dir / "pathways", PATHWAYS_DIR)
    VERSION_FILE.write_text(f"v_{today}")
    write_zenodo()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--force", action="store_true",
                        help="process the data and write a new version also when nothing is changed")
    args = parser.parse_args()

    df_input = load_input()
    errors, warnings = validate(df_input)
    report(errors, warnings)
    if errors:
        sys.exit("Input validation failed. Correct data/input/cnv_data.tsv.")

    orphacodes = {code for value in df_input["orphacodes"] for code, _ in parse_orphacodes(value)}
    wp_ids = sorted({wp for wp in df_input["wikipathways_id"] if wp})
    previous = load_previous_sources()
    input_changed = previous.get("input_sha256") != sha256(INPUT_FILE)

    with tempfile.TemporaryDirectory() as tmp_name:
        tmp = Path(tmp_name)

        releases = {"orphadata": {name: orphadata_release(url) for name, url in ORPHADATA_FILES.items()},
                    "hgnc": hgnc_release(), "wikipathways": {}, "hpo": hpo_release()}
        pathway_files = []
        for wp in wp_ids:
            for ext in ("gpml", "svg"):
                download(WIKIPATHWAYS_URL.format(wp=wp, ext=ext), tmp / f"{wp}.{ext}")
                pathway_files.append(tmp / f"{wp}.{ext}")
            releases["wikipathways"][wp] = wikipathways_release(tmp / f"{wp}.gpml")

        if release_ids(releases) == release_ids(previous) and not input_changed and not args.force:
            print("No new source release and no input change. Data is not changed.")
            return

        for name, url in ORPHADATA_FILES.items():
            download(url, tmp / f"{name}.xml")
            releases["orphadata"][name].update({"url": url, "sha256": sha256(tmp / f"{name}.xml")})
        download(HGNC_URL, tmp / "hgnc.txt")
        releases["hgnc"].update({"url": HGNC_URL, "sha256": sha256(tmp / "hgnc.txt")})
        for wp in wp_ids:
            releases["wikipathways"][wp]["url"] = WIKIPATHWAYS_URL.format(wp=wp, ext="gpml")
        hpo_url = HPO_URL.format(release=releases["hpo"]["release"])
        download(hpo_url, tmp / "hp.json")
        releases["hpo"].update({"url": hpo_url, "sha256": sha256(tmp / "hp.json")})

        orphadata = build_orphadata_table(tmp, orphacodes)
        hgnc = load_hgnc(tmp / "hgnc.txt")
        gene_mapping = resolve_symbols(df_input, hgnc)
        for _, row in gene_mapping[gene_mapping["match_type"] != "approved"].iterrows():
            print(f"[WARN] {row['cnv']}: gene '{row['input_symbol']}' -> {row['match_type']} {row['symbol']}")

        index = build_xref_index(hgnc)
        hgnc_by_id = hgnc.set_index("hgnc_id", drop=False).to_dict("index")
        pathway_genes = pd.DataFrame(
            [row for wp in wp_ids for row in parse_pathway_genes(wp, tmp / f"{wp}.gpml", index, hgnc_by_id)])

        used_ids = (set(gene_mapping["hgnc_id"]) | set(pathway_genes["hgnc_id"])) - {""}
        hgnc_subset = hgnc[hgnc["hgnc_id"].isin(used_ids)].sort_values(
            "hgnc_id", key=lambda ids: ids.str[5:].astype(int))

        hpo_categories = build_hpo_categories(orphadata, *load_hpo(tmp / "hp.json"), releases["hpo"]["release"])

        tables = {"orphadata_filtered.csv": orphadata, "hgnc_filtered.tsv": hgnc_subset,
                  "gene_mapping.tsv": gene_mapping, "pathway_genes.tsv": pathway_genes,
                  "hpo_categories.json": hpo_categories}

        new_dir = tmp / "out"
        new_dir.mkdir()
        sources = {"retrieved": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                   "input_sha256": sha256(INPUT_FILE), **releases}
        write_outputs(new_dir, tables, pathway_files, sources)

        if files_equal(new_dir, tables) and not args.force:
            print("Processed data is not changed. No new version.")
            return

        today = datetime.date.today().isoformat()
        publish(new_dir, tables, today)

    print(f"Update complete. Version: v_{today}")


if __name__ == "__main__":
    main()
