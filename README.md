# Copy Number Variation Data

**cnv-data** provides the data for the [CNVPathwayAtlas](https://github.com/CNVPathwayAtlas/cnv-website).

## Licence
The data in this repository is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/) (see [LICENSE](LICENSE)).
Orphadata content is © Orphanet, licensed under CC BY 4.0. HGNC data is from the HUGO Gene Nomenclature Committee. Pathways are from [WikiPathways](https://www.wikipathways.org/) (CC0). Phenotype categories are from the [Human Phenotype Ontology](https://hpo.jax.org/) (see its licence terms).

## Citing
For citation details of our work and the data resources we used, refer to the [Documentation](https://cnvpathwayatlas.github.io/cnv-website/documentation/#how-to-cite).

## Contents

| Path | Content |
|---|---|
| `data/input/cnv_data.tsv` | Curated CNV table (edit this file) |
| `data/processed/orphadata_filtered_<date>.csv` | Orphanet name, definition, phenotypes (HPO), prevalence and OMIM IDs per ORPHAcode |
| `data/processed/hgnc_filtered_<date>.tsv` | HGNC records of the CNV genes and the pathway genes |
| `data/processed/gene_mapping_<date>.tsv` | Input gene symbol → approved HGNC symbol (`approved`, `previous`, `alias`, `not_found`) |
| `data/processed/pathway_genes_<date>.tsv` | Gene nodes of each WikiPathways pathway, mapped to HGNC |
| `data/processed/hpo_categories_<date>.json` | Top-level HPO categories (children of HP:0000118) of each Orphadata phenotype, from the HPO release in `sources.json` |
| `data/processed/sources_<date>.json` | Release date, URL and SHA-256 of each source file |
| `data/latest/` | Links to the latest processed files |
| `data/pathways/` | GPML and SVG file of each WikiPathways pathway |
| `version.txt` | Data version (`v_<date>`) |

## Automatic update
The workflow [update-data.yml](.github/workflows/update-data.yml) runs on the first day of each month (05:00 UTC), and when `data/input/` or `scripts/` changes on `main`. It:
1. validates `data/input/cnv_data.tsv`,
2. reads the release date of each Orphadata product, the HGNC file version, the revision of each pathway and the latest HPO release,
3. stops when no source has a new release and the input is not changed,
4. otherwise downloads and processes the sources, and writes a new data version only when the processed data is changed.

A new version also writes the release of each source to `.zenodo.json` ([source_notes.py](scripts/source_notes.py)). A new `version.txt` creates a GitHub release with the same list of source releases and starts the cnv-website rebuild for that commit.
The workflows use the secret `CNV_RELEASE_AUTOMATION` (a token with write access to both repositories).

## Installation
```sh
conda create -n data_env python=3.11 -y
conda activate data_env
pip install -r requirements.txt
```

## Add or edit a CNV with an issue form
Open a [new issue](../../issues/new/choose) and select one of these forms:
- **Add a CNV**: complete one field for each column.
- **Edit a CNV**: select the CNV in the list (e.g. `16p11.2-distal`) and complete only the fields to edit. An empty field keeps the current value, `-` makes a field empty. A list field replaces the full current list.

The workflow [cnv-issue.yml](.github/workflows/cnv-issue.yml) runs [cnv_from_issue.py](scripts/cnv_from_issue.py). It:
1. adds or edits the row in `data/input/cnv_data.tsv` and checks it with `validate_input.py`,
   and with [check_sources.py](scripts/check_sources.py) against the online sources:
   - locus and regions: same chromosome, region in the chromosome, region overlaps the bands of the locus (UCSC cytoBand, hg19 and hg38)
   - WikiPathways ID exists and the pathway name is the WikiPathways title
   - each ORPHAcode is in Orphadata
   - each gene is an approved HGNC symbol on the CNV chromosome
   - each PubMed ID exists
2. writes the result as a comment on the issue (an edit shows the old and the new value),
3. when the check passes, opens a pull request from the branch `cnv/issue-<number>`.

The workflow [update-edit-form.yml](.github/workflows/update-edit-form.yml) runs [update_edit_form.py](scripts/update_edit_form.py) when `data/input/cnv_data.tsv` changes on `main`. It writes the CNV list of the Edit a CNV form.

Edit the issue to correct an error: the workflow runs again and updates the pull request. Merge the pull request to start the data update.
Anyone can submit a form. Only a maintainer with write access to the repository can merge the pull request, so every change needs a maintainer approval.
The repository must have the labels `new-cnv` and `edit-cnv`, and the setting **Settings → Actions → General → Allow GitHub Actions to create and approve pull requests** must be on.

## Add or change a CNV
Edit [data/input/cnv_data.tsv](data/input/cnv_data.tsv). Columns:

| Column | Content |
|---|---|
| `locus` | Cytogenetic band, e.g. `16p11.2` |
| `extra` | Sub-region or key gene, e.g. `distal` or `NPHP1` (optional). It is part of the CNV page name |
| `pathway_name` | WikiPathways pathway name |
| `wikipathways_id` | e.g. `WP4657` |
| `orphacodes` | `code(cause)` items separated by `;`, e.g. `261330(microdeletion);261337(microduplication)` |
| `region_hg19` | GRCh37 region, e.g. `Chr16:29592751-30190593` |
| `region_hg38` | GRCh38 region |
| `genes` | HGNC symbols separated by `;` |
| `pubmed_id` | PubMed IDs separated by `;` |
| `description` | Free text |

Check the file before a push:
```sh
python scripts/validate_input.py
```

## Update the datasets manually
```sh
python scripts/update_datasets.py          # only when a source or the input is changed
python scripts/update_datasets.py --force  # always write a new version
```
