import pandas as pd
import os
import numpy as np

def load_and_combine_csv_files(folder_path, depre_file):
    """
    Carga todos los archivos CSV con el sufijo '_freqAltas.csv' de la carpeta especificada,
    extrayendo la columna 'Frecuencia (%)' y combinándola en un solo DataFrame.

    Parámetros:
    folder_path (str): Ruta a la carpeta que contiene los archivos CSV.

    Retorna:
    pd.DataFrame: DataFrame combinado con las frecuencias de todos los archivos CSV.
    """
    # Obtener la lista de archivos CSV en la carpeta
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('_freqAltas.csv')]

    # Crear un DataFrame vacío para almacenar los datos
    master_df = pd.DataFrame()

    # Iterar sobre cada archivo CSV y agregar los datos al DataFrame maestro
    for csv_file in csv_files:
        # Obtener el study_id del nombre del archivo
        study_id = csv_file.replace('_freqAltas.csv', '')

        # Leer el archivo CSV
        df = pd.read_csv(os.path.join(folder_path, csv_file))

        # Establecer el índice del DataFrame como el nombre del gen
        df.set_index('Gene', inplace=True)

        # Seleccionar solo la columna de frecuencia y renombrarla con el study_id
        df = df[['Frecuencia (%)']].rename(columns={'Frecuencia (%)': study_id})

        # Unir el DataFrame actual con el DataFrame maestro
        master_df = master_df.join(df, how='outer')

                # # Guardar el DataFrame maestro en un archivo CSV
        output_filename = 'genes_cancer_depre_Final.csv'
        master_df.to_csv(output_filename)

    # original_df = pd.read_csv(path)
    original_df = master_df.copy()
    # Crear una nueva columna 'depre' en el DataFrame original
    original_df['depre'] = np.nan

    frecuencia_df = pd.read_csv(depre_file)

    # Iterar sobre el dataset de frecuencias
    for index, row in frecuencia_df.iterrows():
        gene = row['Gene']
        frecuencia = row['Frecuencia (%)']

        # Verificar si el gen existe en el dataset original
        if gene in original_df.index:
            # Si existe, agregar la frecuencia en la columna 'depre'
            original_df.at[gene, 'depre'] = frecuencia
        else:
            # Si no existe, crear un nuevo diccionario para la fila
            new_row = {col: np.nan for col in original_df.columns if col != 'depre'}
            new_row['depre'] = frecuencia

            # Crear un DataFrame con la nueva fila
            new_row_df = pd.DataFrame([new_row], index=[gene])

            # Concatenar el nuevo DataFrame al original
            original_df = pd.concat([original_df, new_row_df])

    return original_df