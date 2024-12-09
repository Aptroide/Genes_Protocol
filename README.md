
## Description

This project processes genetic data related to cancer and depression studies. It uses the cBioPortal API to fetch mutation profiles and combines frequency data from multiple CSV files.

## Setup

1. Clone the repository.
2. Install the required dependencies:
    ```sh
    pip install -r requirements.txt
    ```

## Usage

### Configure Important file and Variables before run it
On the [main](main.py), you must change the config dictionary:

```py
config = {
    'file_path': csv File,
    'output_path': txt File,
    'num_cancer_studies': Integer,
    'num_genes_per_study': Integer or None
}
```
Where:
- `file_path`: 
    - Contains depression-related genes data
    - Input file with gene frequencies for depression
    - Required for comparing against cancer genes

- `output_path`: 
    - Path to file containing list of cancer study IDs
    - Each line has one study ID from cBioPortal
    - Required to know which studies process

- `num_cancer_studies`: 
    - Maximum number of cancer studies to process
    - Limits the total studies analyzed
    - Set (e.g., 1000) to process all available studies (depends on the studies available on cBioPortal)
    - Can be reduced for testing (e.g., 10)

- `num_genes_per_study`: 
    - Number of genes to analyze per cancer study
    - `None` means analyze all genes in study
    - Can be set to specific number to limit analysis

### Process Studies
The script processes studies listed in [study_ids.txt](study_ids.txt) and combines frequency data from CSV files in the `Frecuencias` folder.

### Combine CSV Files
The `load_and_combine_csv_files` function in [modules/proc_data.py](modules/proc_data.py) loads and combines CSV files with the suffix `_freqAltas.csv`.

### Clustering
We performe a K-mean clustering and a DBSCAN clustering.

## Output

- `comparacion_genes_depresion_Final.csv`: Combined and filtered frequency data.

- `times.json`: JSON file containing processing times.
