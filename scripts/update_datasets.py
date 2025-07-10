import os
import re
import requests
import xml.etree.ElementTree as ET
from collections import defaultdict
import pandas as pd
import tempfile
import datetime

# Directories
DATA_PROCESSED_DIR = os.path.join("data", "processed")
DATA_LATEST_DIR = os.path.join("data", "latest")
INPUT_DIR = os.path.join("data", "input")
VERSION_FILE = "version.txt"

# Orphadata URLs
ORPHADATA_FILES = {
    "definitions": "https://www.orphadata.com/data/xml/en_product1.xml",
    "phenotypes": "https://www.orphadata.com/data/xml/en_product4.xml",
    "prevalence": "https://www.orphadata.com/data/xml/en_product9_prev.xml",
}

# HGNC TSV URL
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

# Excel input file containing Orphacodes column
EXCEL_FILE = os.path.join(INPUT_DIR, "cnv_data.xlsx")

def ensure_dirs():
    os.makedirs(DATA_PROCESSED_DIR, exist_ok=True)
    os.makedirs(DATA_LATEST_DIR, exist_ok=True)
    os.makedirs(INPUT_DIR, exist_ok=True)

# Extract orphacodes
def load_orphacodes_from_excel():
    df = pd.read_excel(EXCEL_FILE, engine='openpyxl')
    orphacodes = set()
    pattern = re.compile(r"^(\d+)") 

    for entry in df['Orphacodes'].dropna():
        parts = [code.strip() for code in str(entry).split(",")]
        for part in parts:
            match = pattern.match(part)
            if match:
                orphacodes.add(match.group(1))
    return orphacodes

# Orphadata Parsing Functions
def parse_definitions(path, orphacodes):
    root = ET.parse(path).getroot()
    data = {}
    for disorder in root.findall(".//Disorder"):
        code = disorder.findtext("OrphaCode")
        if code in orphacodes:
            definition = ""
            for summary_info in disorder.findall("SummaryInformationList/SummaryInformation"):
                lang = summary_info.attrib.get("lang", "")
                if lang == "en":
                    for text_section in summary_info.findall("TextSectionList/TextSection"):
                        section_type = text_section.find("TextSectionType/Name[@lang='en']")
                        if section_type is not None and section_type.text == "Definition":
                            contents = text_section.findtext("Contents")
                            if contents:
                                definition = contents.strip()
                                break
                if definition:
                    break
            data[code] = definition
    return data

def parse_phenotypes(path, orphacodes, target_freq="Very frequent (99-80%)"):
    root = ET.parse(path).getroot()
    phenos = defaultdict(list)
    for disorder in root.findall(".//Disorder"):
        code = disorder.findtext("OrphaCode")
        if code not in orphacodes:
            continue
        assoc_list = disorder.find("HPODisorderAssociationList")
        if assoc_list is None:
            continue
        for assoc in assoc_list.findall("HPODisorderAssociation"):
            freq_name = assoc.findtext("HPOFrequency/Name[@lang='en']")
            if freq_name == target_freq:
                hpo_term = assoc.findtext("HPO/HPOTerm") or ""
                hpo_id = assoc.findtext("HPO/HPOId") or ""
                if hpo_term and hpo_id:
                    phenos[code].append(f"{hpo_term} ({hpo_id})")
    return phenos

def parse_prevalence(path, orphacodes):
    root = ET.parse(path).getroot()
    prev_data = defaultdict(list)
    pmid_pattern = re.compile(r'(\d+)\[PMID\]')
    for disorder in root.findall(".//Disorder"):
        code = disorder.findtext("OrphaCode")
        if code not in orphacodes:
            continue
        prev_list = disorder.find("PrevalenceList")
        if prev_list is None:
            continue
        for prev in prev_list.findall("Prevalence"):
            val = prev.findtext("ValMoy") or ""
            class_ = prev.findtext("PrevalenceClass/Name[@lang='en']") or ""
            source_text = prev.findtext("Source") or ""
            pmids = pmid_pattern.findall(source_text)
            prevalence_str = val
            if class_:
                prevalence_str += f" ({class_})"
            if pmids:
                prevalence_str += " (PMID:" + ",".join(pmids) + ")"
            prev_data[code].append(prevalence_str)
    return prev_data

def parse_omim(path, orphacodes):
    root = ET.parse(path).getroot()
    omim_map = defaultdict(list)
    for disorder in root.findall(".//Disorder"):
        code = disorder.findtext("OrphaCode")
        if code not in orphacodes:
            continue
        ext_refs = disorder.find("ExternalReferenceList")
        if ext_refs is None:
            continue
        for ext_ref in ext_refs.findall("ExternalReference"):
            source = ext_ref.findtext("Source")
            ref = ext_ref.findtext("Reference")
            if source == "OMIM" and ref:
                omim_map[code].append(ref)
    return omim_map

def save_combined_csv(defs, phenos, prevs, omims, orphacodes, output_path):
    rows = []
    for code in orphacodes:
        rows.append({
            "OrphaCode": code,
            "Definition": defs.get(code, ""),
            "Phenotypes": "; ".join(phenos.get(code, [])),
            "Prevalence": "; ".join(prevs.get(code, [])),
            "OMIM": "; ".join(omims.get(code, []))
        })
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)

# HGNC download and processing
def download_hgnc(tmpdir):
    print(f"Downloading HGNC file from {HGNC_URL}...")
    r = requests.get(HGNC_URL)
    r.raise_for_status()
    hgnc_path = os.path.join(tmpdir, "hgnc_complete_set.txt")
    with open(hgnc_path, "wb") as f:
        f.write(r.content)
    return hgnc_path

def filter_hgnc(hgnc_path):
    df = pd.read_csv(hgnc_path, sep="\t", dtype=str)
    columns_of_interest = [
        "symbol",
        "name",
        "entrez_id",
        "ensembl_gene_id",
        "uniprot_ids",
    ]
    return df[columns_of_interest].copy()

def save_hgnc(filtered_df, output_path):
    filtered_df.to_csv(output_path, sep="\t", index=False)

# Utility functions
def write_version_file(version):
    with open(VERSION_FILE, "w") as f:
        f.write(version)

def update_symlink(src_path, symlink_dir):
    os.makedirs(symlink_dir, exist_ok=True)

    # Make both paths absolute first
    src_path = os.path.abspath(src_path)
    symlink_dir = os.path.abspath(symlink_dir)
    symlink_path = os.path.join(symlink_dir, os.path.basename(src_path))

    # Compute relative path
    relative_src = os.path.relpath(src_path, symlink_dir)

    if os.path.islink(symlink_path) or os.path.exists(symlink_path):
        os.remove(symlink_path)

    os.symlink(relative_src, symlink_path)

def clean_latest_dir():
    for filename in os.listdir(DATA_LATEST_DIR):
        path = os.path.join(DATA_LATEST_DIR, filename)
        if os.path.islink(path) or os.path.isfile(path):
            os.remove(path)

# Main 
def main():
    ensure_dirs()
    orphacodes = load_orphacodes_from_excel()

    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Using temporary directory: {tmpdir}")

        # Download Orphadata files
        for name, url in ORPHADATA_FILES.items():
            out_path = os.path.join(tmpdir, f"{name}.xml")
            print(f"Downloading Orphadata {name}...")
            r = requests.get(url)
            r.raise_for_status()
            with open(out_path, "wb") as f:
                f.write(r.content)

        # Parse Orphadata files
        def_path = os.path.join(tmpdir, "definitions.xml")
        defs = parse_definitions(def_path, orphacodes)
        omims = parse_omim(def_path, orphacodes)
        phenos = parse_phenotypes(os.path.join(tmpdir, "phenotypes.xml"), orphacodes)
        prevs = parse_prevalence(os.path.join(tmpdir, "prevalence.xml"), orphacodes)

        # Download and filter HGNC
        hgnc_path = download_hgnc(tmpdir)
        filtered_hgnc = filter_hgnc(hgnc_path)

    today = datetime.datetime.now().strftime("%Y-%m-%d")

    # Save Orphadata processed CSV
    csv_path = os.path.join(DATA_PROCESSED_DIR, f"orphadata_filtered_{today}.csv")
    save_combined_csv(defs, phenos, prevs, omims, orphacodes, csv_path)

    # Save HGNC TSV
    hgnc_out_path = os.path.join(DATA_PROCESSED_DIR, f"hgnc_filtered_{today}.csv")
    save_hgnc(filtered_hgnc, hgnc_out_path)

    # Clean previous versions in the /latest folder
    clean_latest_dir()

    # Update latest symlinks
    update_symlink(csv_path, DATA_LATEST_DIR)
    update_symlink(hgnc_out_path, DATA_LATEST_DIR)

    # Save version
    write_version_file("v_" + today)

    print("Update complete!")
    print(f"Version: {today}")
    print(f"Processed files saved to {DATA_PROCESSED_DIR}")
    print(f"Latest symlinks updated in {DATA_LATEST_DIR}")

if __name__ == "__main__":
    main()
