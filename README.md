# Copy Number Variation Data 

**cnv-data** is a repository that provides data for the [CNVPathwayAtlas](https://github.com/CNVPathwayAtlas/cnv-website)

## Copyright information
This data is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

## Citing
For complete citation details of our work or the data resources we used, please refer to the [Documentation](https://cnvpathwayatlas.github.io/cnv-website/documentation/#how-to-cite).

## Main features of this repository

- Download the latest ORPHADATA files
- Download the latest HGNC gene information
- Filter and merge data into clean outputs

## Requirements
- [conda](https://docs.conda.io/en/latest/)

## Installation

### Setup conda environment
```sh
conda create -n data_env python=3.11 -y
conda activate data_env
```

### Install dependencies
```sh
pip install -r requirements.txt
```

## Add a new CNV (optional)
To add a new CNV or modify existing ones (e.g. add a disease identifier like an ORPHAcode, or a description), edit [cnv_data.xlsx](cnv-data/input/cnv_data.xlsx)

## Update the datasets
```sh
python scripts/update_datasets.py
```
