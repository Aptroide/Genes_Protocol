import os
from modules import data_extration, val_csv, proc_data, dim_red, k_means, dbscan
from bravado.client import SwaggerClient
import pandas as pd
import time
import json

def procesar_dataframe(data_filtered):
    # Crear una copia del DataFrame original
    df_filtrado = data_filtered.copy()
    # Reiniciar el índice para asegurarse de que empieza desde 0
    df_filtrado = df_filtrado.reset_index()
    # Renombrar la columna que contiene los genes a "Gene"
    df_filtrado = df_filtrado.rename(columns={"index": "Gene"})
    # Asegurarse de que el índice sea numérico nuevamente
    df_filtrado = df_filtrado.reset_index(drop=True)
    # Rellenar los valores NaN con 0.0
    df_filtrado.fillna(0.0, inplace=True)
    # Transponer el DataFrame
    df_transpuesto = df_filtrado.transpose()

    # Usar la primera fila como encabezado
    df_transpuesto.columns = df_transpuesto.iloc[0]  # Establece la primera fila como encabezado
    df_transpuesto = df_transpuesto[1:]  # Elimina la primera fila, ya que ahora es el encabezado

    # Restablecer el índice y renombrarlo a 'Cancer'
    df_transpuesto.reset_index(inplace=True)
    df_transpuesto = df_transpuesto.rename(columns={'index': 'Cancer'})

    # Eliminar el nombre del índice
    df_transpuesto.index.name = None

    # Convertir todas las columnas (excepto 'Cancer') a tipo float64
    for column in df_transpuesto.columns[1:]:
        df_transpuesto[column] = pd.to_numeric(df_transpuesto[column], errors='coerce')

    return df_transpuesto

def eliminar_archivos_carpeta(carpeta):
    for archivo in os.listdir(carpeta):
        archivo_path = os.path.join(carpeta, archivo)
        if os.path.isfile(archivo_path):
            os.remove(archivo_path)
        elif os.path.isdir(archivo_path):
            os.rmdir(archivo_path)  # Elimina la carpeta y su contenido

def initialize_cbioportal():
    """
    Inicializa y configura el cliente de cBioPortal.
    """
    cbioportal = SwaggerClient.from_url(
        'https://www.cbioportal.org/api/v2/api-docs',
        config={
            "validate_requests": False,
            "validate_responses": False,
            "validate_swagger_spec": False,
            "use_models": False,
        }
    )

    # Normalizar los nombres de los métodos
    for a in dir(cbioportal):
        if not a.startswith('_'):
            cbioportal.__setattr__(a.replace(' ', '_').lower(), cbioportal.__getattr__(a))

    return cbioportal

cbioportal = initialize_cbioportal()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    file_path = 'depre_freqAltas.csv'
    # Leer la lista de study_ids desde el archivo local
    output_path = 'study_ids.txt'

    with open(output_path, 'r') as file:
        study_ids_from_file = [line.strip() for line in file.readlines()]

    study_ids = study_ids_from_file
    # study_ids = fetch_all_studies()

    os.makedirs('Frecuencias', exist_ok=True)
    eliminar_archivos_carpeta('Frecuencias')

    num_cancer_studies = 1000 # 1000 si quiero todos los estudios
    num_genes_per_study = None #None si quiero todos los genes

    study_ids = study_ids[:num_cancer_studies]

    aux = 0
    print(f"Total de estudios a procesar: {len(study_ids)}")

    start_time = time.time()
    for study_id in study_ids:
        if aux == num_cancer_studies:
            break
        try:
            data_extration.process_study(cbioportal, study_id, num_genes_per_study)
            aux += 1
        except Exception as e:
            print(f"Se produjo un error al procesar el estudio '{study_id}': {e}")
            continue  # Continuar con el siguiente estudio en caso de error
    end_time = time.time()
    elapsed_time_seconds = end_time - start_time
    data_extration_time = elapsed_time_seconds 
    data_p = data_extration_time / 60
    print(f"Procesamiento de {aux} estudios completo en {data_p:.2f} minutos.")

    print("\nValidando formato de archivo...")
    

    if val_csv.validate_csv_file(file_path):
        # cBioportal Data extration
        print("Archivo Correcto.")
        depresion = pd.read_csv(file_path)
        genes_depresion = depresion['Gene'].tolist()
        print("\nCombinando archivos de frecuencias...")

        start_time = time.time()
        master_df = proc_data.load_and_combine_csv_files('/home/fernando/Documents/HPC/Final_Protocol/Genes_Protocol/Frecuencias', file_path)

        print(f"Total de estudios con al menos una mutación: {master_df.shape}")

    
        filtered_df = master_df.loc[genes_depresion]
        # Eliminar las columnas que contienen solo NaN
        filtered_df = filtered_df.dropna(axis=1, how='all')

        # Guardar el DataFrame filtrado en un archivo CSV
        filtered_output_filename = 'comparacion_genes_depresion_Final.csv'
        filtered_df.to_csv(filtered_output_filename)
        print(f"Datos filtrados tienen una dimension de '{filtered_df.shape}'.")
        end_time = time.time()
        data_proc_time = end_time - start_time
    

        # Clustering
        print("\nIniciando pre-procesameinto y cambio de dimensión...")
        df_transpuesto = procesar_dataframe(filtered_df)
        features = df_transpuesto.iloc[:, 1:]
        start_time = time.time()
        proj2d, kpca_features, tsne_features = dim_red.dim_reduction(df_transpuesto, features)
        print("Pre-procesamiento completo. Dimesiones hechas: umap, kpca y tsne")
        end_time = time.time()
        dim_red_time = end_time - start_time

        
        print("\nIniciando Clustering")
        print("\nUsando UMAP:")
        print("\nIniciando K-Means...")
        k_means.c_kmeans(proj2d)
        k_kmeans = input("Ingrese el numero de clusters: ")
        k_kmeans = int(k_kmeans)
        start_time = time.time()
        range_n_clusters = [k_kmeans]
        aux, aux_df, k_means_labels_umap = k_means.plot_kmeans_silhouette_plotly(proj2d, df_transpuesto, range_n_clusters)
        depre_clusters = k_means_labels_umap[k_means_labels_umap['Cancer'] == 'depre']
        print("Clusters for 'depre' in K-Means UMET:")
        a =depre_clusters['Cluster']
        print(a)
        end_time = time.time()
        k_means_umap_time = end_time - start_time

        start_time = time.time()
        print("\nIniciando DBSCAN...")
        r, dbscan_labels_umap = dbscan.optimize_dbscan(proj2d, df_transpuesto)
        print("Clusters for 'depre' in DBSCAN UMET:")
        print(dbscan_labels_umap[dbscan_labels_umap['Cancer'] == 'depre']['Cluster'])
        end_time = time.time()
        dbscan_umap_time = end_time - start_time

        print("\nUsando KPCA:")
        print("\nIniciando K-Means...")
        k_means.c_kmeans(kpca_features)
        k_kmeans = input("Ingrese el numero de clusters: ")
        k_kmeans = int(k_kmeans)
        range_n_clusters = [k_kmeans]
        start_time = time.time()
        aux, aux_df, k_means_labels_kpca = k_means.plot_kmeans_silhouette_plotly(kpca_features, df_transpuesto, range_n_clusters)
        depre_clusters = k_means_labels_kpca[k_means_labels_kpca['Cancer'] == 'depre']
        print("Clusters for 'depre' in K-Means KPCA:")
        a =depre_clusters['Cluster']
        print(a)
        end_time = time.time()
        k_means_kpca_time = end_time - start_time

        start_time = time.time()
        print("\nIniciando DBSCAN...")
        r, dbscan_labels_kpca = dbscan.optimize_dbscan(kpca_features, df_transpuesto)
        print("Clusters for 'depre' in DBSCAN KPCA:")
        print(dbscan_labels_kpca[dbscan_labels_kpca['Cancer'] == 'depre']['Cluster'])
        end_time = time.time()
        dbscan_kpca_time = end_time - start_time

        print("\nUsando t-SNE:")
        print("\nIniciando K-Means...")
        k_means.c_kmeans(tsne_features)
        k_kmeans = input("Ingrese el numero de clusters: ")
        k_kmeans = int(k_kmeans)
        range_n_clusters = [k_kmeans]
        start_time = time.time()
        aux, aux_df, k_means_labels_tsne = k_means.plot_kmeans_silhouette_plotly(tsne_features, df_transpuesto, range_n_clusters)
        depre_clusters = k_means_labels_tsne[k_means_labels_tsne['Cancer'] == 'depre']
        print("Clusters for 'depre' in K-Means t-SNE:")
        a =depre_clusters['Cluster']
        print(a)
        end_time = time.time()
        k_means_tsne_time = end_time - start_time

        start_time = time.time()
        print("\nIniciando DBSCAN...")
        r, dbscan_labels_t_sne = dbscan.optimize_dbscan(tsne_features, df_transpuesto)
        print("Clusters for 'depre' in DBSCAN t-SNE:")
        print(dbscan_labels_t_sne[dbscan_labels_t_sne['Cancer'] == 'depre']['Cluster'])  
        end_time = time.time()
        dbscan_tsne_time = end_time - start_time

    
        # Crear un diccionario con los tiempos
        times_dict = {
            "data_proc_time": data_proc_time,
            "data_extration_time": data_extration_time,
            "dim_red_time": dim_red_time,
            "k_means_umap_time": k_means_umap_time,
            "dbscan_umap_time": dbscan_umap_time,
            "k_means_kpca_time": k_means_kpca_time,
            "dbscan_kpca_time": dbscan_kpca_time,
            "k_means_tsne_time": k_means_tsne_time,
            "dbscan_tsne_time": dbscan_tsne_time
        }

        # Guardar el diccionario en un archivo JSON
        with open('times.json', 'w') as json_file:
            json.dump(times_dict, json_file, indent=4)

        print("Tiempos guardados en 'times.json'.")
    
    else:
        print("El archivo no cumple con las restricciones.")