# Copy Number Variation Data 

**cnv-data** is a repository that provides data for the CNV booklet.

## Main Features

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

### Update the datasets
```sh
python scripts/update_datasets.py
```