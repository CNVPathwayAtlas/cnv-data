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
INPUT_DIR = "input"
VERSION_FILE = "version.txt"

# Orphadata URLs (XML files)
ORPHADATA_FILES = {
    "definitions": "https://www.orphadata.com/data/xml/en_product1.xml",
    "phenotypes": "https://www.orphadata.com/data/xml/en_product4.xml",
    "prevalence": "https://www.orphadata.com/data/xml/en_product9_prev.xml",
    # OMIM data is also in en_product1.xml (same as definitions)
    "omim": "https://www.orphadata.com/data/xml/en_product1.xml",
}

# HGNC complete dataset URL (TSV)
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

def ensure_dirs():
    os.makedirs(DATA_PROCESSED_DIR, exist_ok=True)
    os.makedirs(DATA_LATEST_DIR, exist_ok=True)
    os.makedirs(INPUT_DIR, exist_ok=True)

def load_orphacodes():
    path = os.path.join(INPUT_DIR, "orphacodes.txt")
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())

# --- Orphadata parsing functions ---

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
    # OMIM references are in the definitions XML file (en_product1.xml)
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

def save_combined_tsv(defs, phenos, prevs, omims, orphacodes, output_path):
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("OrphaCode\tDefinition\tPhenotypes\tPrevalence\tOMIM\n")
        for code in orphacodes:
            line = [
                code,
                defs.get(code, ""),
                "; ".join(phenos.get(code, [])),
                "; ".join(prevs.get(code, [])),
                "; ".join(omims.get(code, []))
            ]
            f.write("\t".join(line) + "\n")

def save_combined_excel(defs, phenos, prevs, omims, orphacodes, output_path):
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
    df.to_excel(output_path, index=False)

# --- HGNC download and processing ---

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
    filtered_df = df[columns_of_interest].copy()
    return filtered_df

def save_hgnc(filtered_df, output_path):
    filtered_df.to_csv(output_path, sep="\t", index=False)

# --- Utility functions ---

def write_version_file(version):
    with open(VERSION_FILE, "w") as f:
        f.write(version)

def update_symlink(src_path, symlink_dir):
    symlink_path = os.path.join(symlink_dir, os.path.basename(src_path))
    if os.path.islink(symlink_path) or os.path.exists(symlink_path):
        os.remove(symlink_path)
    os.symlink(os.path.abspath(src_path), symlink_path)

def main():
    ensure_dirs()
    orphacodes = load_orphacodes()

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
        defs = parse_definitions(os.path.join(tmpdir, "definitions.xml"), orphacodes)
        phenos = parse_phenotypes(os.path.join(tmpdir, "phenotypes.xml"), orphacodes)
        prevs = parse_prevalence(os.path.join(tmpdir, "prevalence.xml"), orphacodes)
        omims = parse_omim(os.path.join(tmpdir, "omim.xml"), orphacodes)

        # Download and filter HGNC
        hgnc_path = download_hgnc(tmpdir)
        filtered_hgnc = filter_hgnc(hgnc_path)

    today = datetime.datetime.now().strftime("%Y-%m-%d")

    # Save Orphadata processed files
    tsv_path = os.path.join(DATA_PROCESSED_DIR, f"orphadata_filtered_{today}.tsv")
    xlsx_path = os.path.join(DATA_PROCESSED_DIR, f"orphadata_filtered_{today}.xlsx")
    save_combined_tsv(defs, phenos, prevs, omims, orphacodes, tsv_path)
    save_combined_excel(defs, phenos, prevs, omims, orphacodes, xlsx_path)

    # Save HGNC filtered file
    hgnc_out_path = os.path.join(DATA_PROCESSED_DIR, f"hgnc_filtered_{today}.tsv")
    save_hgnc(filtered_hgnc, hgnc_out_path)

    # Update latest symlinks
    update_symlink(tsv_path, DATA_LATEST_DIR)
    update_symlink(xlsx_path, DATA_LATEST_DIR)
    update_symlink(hgnc_out_path, DATA_LATEST_DIR)

    # Save version file
    write_version_file("v_" + today)

    print("Update complete!")
    print(f"Version: {today}")
    print(f"Processed files saved to {DATA_PROCESSED_DIR}")
    print(f"Latest symlinks updated in {DATA_LATEST_DIR}")

if __name__ == "__main__":
    main()