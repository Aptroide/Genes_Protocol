
## Description

This project processes genetic data related to cancer and depression studies. It uses the cBioPortal API to fetch mutation profiles and combines frequency data from multiple CSV files.

## Setup

1. Clone the repository.
2. Install the required dependencies:
    ```sh
    pip install -r requirements.txt
    ```

## Usage

### Initialize cBioPortal Client

The `initialize_cbioportal` function in [main.py](main.py) initializes and configures the cBioPortal client.

```py
def initialize_cbioportal():
    cbioportal = SwaggerClient.from_url(
        'https://www.cbioportal.org/api/v2/api-docs',
        config={
            "validate_requests": False,
            "validate_responses": False,
            "validate_swagger_spec": False,
            "use_models": False,
        }
    )
    return cbioportal
```

### Process Studies
The script processes studies listed in [study_ids.txt](study_ids.txt) and combines frequency data from CSV files in the `Frecuencias` folder.

The variables `num_cancer_studies` and `num_genes_per_study` control the number of studies and the number of genes, respectively.

```py
if __name__ == '__main__':
    cbioportal = initialize_cbioportal()
    file_path = 'depre_freqAltas.csv'
    output_path = 'study_ids.txt'

    with open(output_path, 'r') as file:
        study_ids_from_file = [line.strip() for line in file.readlines()]

    study_ids = study_ids_from_file
    os.makedirs('Frecuencias', exist_ok=True)
    eliminar_archivos_carpeta('Frecuencias')

    num_cancer_studies = 1000 # 1000 si quiero todos los estudios
    num_genes_per_study = None #None si quiero todos los genes

    study_ids = study_ids[:num_cancer_studies]

    for study_id in study_ids:
        try:
            data_extration.process_study(cbioportal, study_id, num_genes_per_study)
        except Exception as e:
            print(f"Error processing study '{study_id}': {e}")
            continue
```

### Combine CSV Files
The `load_and_combine_csv_files` function in [modules/proc_data.py](modules/proc_data.py) loads and combines CSV files with the suffix `_freqAltas.csv`.

```py
def load_and_combine_csv_files(folder_path, depre_file):
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('_freqAltas.csv')]
    master_df = pd.DataFrame()

    for csv_file in csv_files:
        study_id = csv_file.replace('_freqAltas.csv', '')
        df = pd.read_csv(os.path.join(folder_path, csv_file))
        df.set_index('Gene', inplace=True)
        df = df[['Frecuencia (%)']].rename(columns={'Frecuencia (%)': study_id})
        master_df = master_df.join(df, how='outer')
    
    return master_df
```

## Output

- `comparacion_genes_depresion_Final.csv`: Combined and filtered frequency data.

- `times.json`: JSON file containing processing times.
