"""Check one CNV row of the input table against the online sources.

The checks are:
- locus and regions: the chromosome of the locus and of the regions is the same, the region is in the
  chromosome, and the region overlaps the cytogenetic bands of the locus (UCSC cytoBand, hg19 and hg38)
- WikiPathways: the pathway ID exists, and the pathway name is the same as the WikiPathways title
- ORPHAcodes: each ORPHAcode is in Orphadata (product 1)
- genes: each symbol is an approved HGNC symbol, and the HGNC location is on the CNV chromosome
- PubMed IDs: each PubMed ID exists

check_row() returns (errors, warnings, notes) as three lists of strings.
"""
import gzip
import io
import re
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests

from cnv_input import cnv_key, parse_orphacodes, parse_region, split_list
from update_datasets import HGNC_URL, ORPHADATA_FILES, TIMEOUT, WIKIPATHWAYS_URL, download, resolve_symbols

CYTOBAND_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/{build}/database/cytoBand.txt.gz"
PUBMED_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
LOCUS_PATTERN = re.compile(r"^(\d+|X|Y)([pq][\d.]*)(?:-([pq]?[\d.]*))?$")
BAND_PATTERN = re.compile(r"^(\d+|X|Y)([pq][\d.]*)")


def load_cytobands(build):
    """Return the UCSC cytoBand table of a build: chromosome (no 'chr'), start, end, band."""
    response = requests.get(CYTOBAND_URL.format(build=build), timeout=TIMEOUT)
    response.raise_for_status()
    df = pd.read_csv(io.BytesIO(gzip.decompress(response.content)), sep="\t", header=None,
                     names=["chromosome", "start", "end", "band", "stain"], dtype={"chromosome": str})
    df["chromosome"] = df["chromosome"].str.replace("chr", "", regex=False)
    return df


def parse_locus(locus):
    """Return (chromosome, first band, last band) for e.g. '15q11.2-q13.1', or None."""
    match = LOCUS_PATTERN.match(locus)
    if not match:
        return None
    chrom, first, last = match.groups()
    if not last:
        last = first
    elif last[0] not in "pq":
        last = first[0] + last
    return chrom, first, last


def band_span(bands, chrom, first, last):
    """Return (start, end) of the bands from first to last on a chromosome, or None for an unknown band.

    A band also includes its sub-bands, e.g. q11.2 includes q11.21, q11.22 and q11.23.
    """
    rows = bands[bands["chromosome"] == chrom]
    spans = []
    for band in (first, last):
        match = rows[rows["band"].str.startswith(band)]
        if match.empty:
            return None
        spans.append((match["start"].min(), match["end"].max()))
    return min(s for s, _ in spans), max(e for _, e in spans)


def bands_at(bands, chrom, start, end):
    """Return the band names that the region start and end are in, e.g. '16p11.2-p12.1'."""
    rows = bands[bands["chromosome"] == chrom]
    names = [rows[(rows["start"] < pos) & (rows["end"] >= pos)]["band"] for pos in (start, end)]
    first, last = (n.iloc[0] if not n.empty else "?" for n in names)
    return chrom + first if first == last else f"{chrom}{first}-{last}"


def check_location(row, cytobands):
    errors, warnings = [], []
    key = cnv_key(row)
    locus = parse_locus(row["locus"])
    if locus is None:
        warnings.append(f"{key}: locus '{row['locus']}' is not in the format 15q11.2 or 15q11.2-q13.1; the locus is not checked")
    for build, column in (("hg19", "region_hg19"), ("hg38", "region_hg38")):
        chrom, start, end = parse_region(row[column])
        if chrom is None:
            continue
        bands = cytobands[build]
        length = bands[bands["chromosome"] == chrom]["end"].max()
        if pd.isna(length):
            errors.append(f"{key}: {column} chromosome {chrom} is not in {build}")
            continue
        if end > length:
            errors.append(f"{key}: {column} end {end:,} is after the end of chr{chrom} in {build} ({length:,})")
            continue
        if locus is None:
            continue
        locus_chrom, first, last = locus
        if chrom != locus_chrom:
            errors.append(f"{key}: {column} is on chr{chrom}, but the locus {row['locus']} is on chr{locus_chrom}")
            continue
        span = band_span(bands, chrom, first, last)
        if span is None:
            warnings.append(f"{key}: a band of the locus {row['locus']} is not in the {build} cytoBand table; {column} is not checked")
            continue
        in_region = bands_at(bands, chrom, start, end)
        if end < span[0] or start > span[1]:
            errors.append(f"{key}: {column} ({in_region}) does not overlap the locus {row['locus']} in {build}")
        elif start < span[0] or end > span[1]:
            warnings.append(f"{key}: {column} ({in_region}) is not fully in the locus {row['locus']} in {build}")
    return errors, warnings


def check_wikipathways(row):
    errors, warnings, notes = [], [], []
    key, wp, name = cnv_key(row), row["wikipathways_id"], row["pathway_name"]
    if not wp:
        warnings.append(f"{key}: no WikiPathways ID")
        return errors, warnings, notes
    response = requests.get(WIKIPATHWAYS_URL.format(wp=wp, ext="gpml"), timeout=TIMEOUT)
    if response.status_code == 404:
        errors.append(f"{key}: WikiPathways ID {wp} does not exist")
        return errors, warnings, notes
    response.raise_for_status()
    title = ET.fromstring(response.content).get("Name", "")
    notes.append(f"{wp}: {title}")
    if not name:
        warnings.append(f"{key}: pathway name is empty (WikiPathways title: '{title}')")
    elif name.strip().lower() != title.strip().lower():
        warnings.append(f"{key}: pathway name '{name}' is not the WikiPathways title '{title}'")
    return errors, warnings, notes


def check_orphacodes(row, tmp):
    errors, notes = [], []
    key = cnv_key(row)
    items = parse_orphacodes(row["orphacodes"])
    if not items:
        return errors, [f"{key}: no ORPHAcodes"], notes
    path = Path(tmp) / "product1.xml"
    download(ORPHADATA_FILES["definitions"], path)
    names = {}
    for _, element in ET.iterparse(path):
        if element.tag == "Disorder":
            names[element.findtext("OrphaCode")] = element.findtext("Name[@lang='en']") or ""
            element.clear()
    for code, cause in items:
        if code in names:
            notes.append(f"ORPHA:{code} ({cause or 'no cause'}): {names[code]}")
        else:
            errors.append(f"{key}: ORPHAcode {code} is not in Orphadata")
    return errors, [], notes


def check_genes(row, tmp):
    warnings = []
    key = cnv_key(row)
    if not split_list(row["genes"]):
        return warnings
    path = Path(tmp) / "hgnc.tsv"
    download(HGNC_URL, path)
    hgnc = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
    location = dict(zip(hgnc["symbol"], hgnc["location"]))
    chrom, _, _ = parse_region(row["region_hg19"])
    for _, gene in resolve_symbols(pd.DataFrame([row]), hgnc).iterrows():
        symbol = gene["input_symbol"]
        if gene["match_type"] == "not_found":
            warnings.append(f"{key}: gene {symbol} is not an HGNC symbol")
            continue
        if gene["match_type"] != "approved":
            warnings.append(f"{key}: gene {symbol} is a {gene['match_type'].replace('_', ' ')} symbol"
                            + (f" of {gene['symbol']}" if gene["symbol"] else ""))
        match = BAND_PATTERN.match(location.get(gene["symbol"], ""))
        if chrom and match and match.group(1) != chrom:
            warnings.append(f"{key}: gene {gene['symbol'] or symbol} is on {location[gene['symbol']]}, not on chr{chrom}")
    return warnings


def check_pubmed(row):
    key, ids = cnv_key(row), split_list(row["pubmed_id"])
    ids = [i for i in ids if i.isdigit()]
    if not ids:
        return []
    response = requests.get(PUBMED_URL, params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"}, timeout=TIMEOUT)
    response.raise_for_status()
    result = response.json()["result"]
    return [f"{key}: PubMed ID {i} does not exist" for i in ids if "error" in result.get(i, {"error": True})]


def check_row(row):
    """Return (errors, warnings, notes) for one input row. A source that cannot be read gives a warning."""
    errors, warnings, notes = [], [], []
    key = cnv_key(row)

    def run(name, check):
        try:
            result = check()
        except (requests.RequestException, ET.ParseError, ValueError) as e:
            warnings.append(f"{key}: {name} could not be checked ({e.__class__.__name__})")
            return
        for target, values in zip((errors, warnings, notes), result):
            target.extend(values)

    with tempfile.TemporaryDirectory() as tmp:
        run("the locus and the regions",
            lambda: check_location(row, {b: load_cytobands(b) for b in ("hg19", "hg38")}))
        run("the WikiPathways ID", lambda: check_wikipathways(row))
        run("the ORPHAcodes", lambda: check_orphacodes(row, tmp))
        run("the genes", lambda: ([], check_genes(row, tmp)))
        run("the PubMed IDs", lambda: (check_pubmed(row),))
    return errors, warnings, notes
