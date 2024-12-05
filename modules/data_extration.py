from collections import defaultdict
from bravado.exception import HTTPNotFound, HTTPError
import pandas as pd

# def fetch_all_studies():
#     """
#     Obtiene todos los estudios disponibles utilizando la API de cBioPortal.

#     Returns:
#         list: Lista de identificadores de estudios.
#     """
#     try:
#         studies = cbioportal.Studies.getAllStudiesUsingGET().result()
#         study_ids = [study['studyId'] for study in studies]
#         print(f"Total number of studies: {len(study_ids)}")
#         return study_ids
#     except Exception as e:
#         print(f"Error al obtener la lista de estudios: {e}")
#         return []

def get_mutation_profile_id(cbioportal, study_id):
    """
    Obtiene el ID del perfil molecular de mutaciones para un estudio dado.

    Args:
        cbioportal: Cliente de cBioPortal ya configurado.
        study_id: ID del estudio.

    Returns:
        ID del perfil molecular de mutaciones o None si no se encuentra.
    """
    try:
        molecular_profiles = cbioportal.molecular_profiles.getAllMolecularProfilesInStudyUsingGET(
            studyId=study_id
        ).response().result

        # Filtrar perfiles de mutaciones
        mutation_profiles = [mp['molecularProfileId'] for mp in molecular_profiles if 'mutation' in mp['molecularProfileId'].lower()]

        if mutation_profiles:
            return mutation_profiles[0]  # Asumir que hay al menos un perfil de mutaciones
        else:
            print(f"No se encontró un perfil molecular de mutaciones para el estudio '{study_id}'.")
            return None
    except HTTPError as e:
        print(f"Error al obtener los perfiles moleculares para el estudio '{study_id}': {e}")
        return None

def get_sample_list_id(cbioportal, study_id):
    """
    Obtiene el ID de la lista de muestras "all" para un estudio dado.

    Args:
        cbioportal: Cliente de cBioPortal ya configurado.
        study_id: ID del estudio.

    Returns:
        ID de la lista de muestras "all" o None si no se encuentra.
    """
    try:
        sample_lists = cbioportal.sample_lists.getAllSampleListsInStudyUsingGET(
            studyId=study_id
        ).response().result

        # Filtrar listas de muestras con 'all' en el ID o nombre
        all_sample_lists = [sl['sampleListId'] for sl in sample_lists if 'all' in sl['sampleListId'].lower()]

        if all_sample_lists:
            return all_sample_lists[0]  # Asumir que hay al menos una lista de muestras "all"
        else:
            print(f"No se encontró una lista de muestras 'all' para el estudio '{study_id}'.")
            return None
    except HTTPError as e:
        print(f"Error al obtener las listas de muestras para el estudio '{study_id}': {e}")
        return None

def process_study(cbioportal, study_id, max_gen=None):
    """
    Procesa un estudio dado su study_id: obtiene mutaciones, calcula estadísticas y guarda un CSV.

    Args:
        cbioportal: Cliente de cBioPortal ya configurado.
        study_id: ID del estudio.
    """
    # print(f"\nProcesando el estudio: {study_id}")

    # Obtener el ID del perfil molecular de mutaciones
    molecular_profile_id = get_mutation_profile_id(cbioportal, study_id)
    if not molecular_profile_id:
        print(f"Saltando el estudio '{study_id}' debido a la falta de un perfil molecular de mutaciones.")
        return

    # Obtener el ID de la lista de muestras "all"
    sample_list_id = get_sample_list_id(cbioportal, study_id)
    if not sample_list_id:
        print(f"Saltando el estudio '{study_id}' debido a la falta de una lista de muestras 'all'.")
        return

    # Obtener todas las muestras del estudio
    try:
        samples_response = cbioportal.samples.getAllSamplesInStudyUsingGET(
            studyId=study_id
        ).response().result
        total_samples = len(samples_response)
        # print(f"Total de muestras perfiladas en el estudio '{study_id}': {total_samples}")
    except HTTPNotFound as e:
        print(f"Error: Estudio '{study_id}' no encontrado.")
        return
    except HTTPError as e:
        print(f"HTTP Error al obtener muestras para el estudio '{study_id}': {e}")
        return

    # Obtener todas las mutaciones en el perfil molecular y lista de muestras
    try:
        mutations_response = cbioportal.mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
            molecularProfileId=molecular_profile_id,
            sampleListId=sample_list_id,
            projection="DETAILED"
        ).response().result
        # print(f"Total de mutaciones obtenidas: {len(mutations_response)}")
    except HTTPNotFound as e:
        print(f"Error: Perfil molecular '{molecular_profile_id}' o lista de muestras '{sample_list_id}' no encontrada para el estudio '{study_id}'.")
        return
    except HTTPError as e:
        print(f"HTTP Error al obtener mutaciones para el estudio '{study_id}': {e}")
        return

    # Verificar si se obtuvieron mutaciones
    # if mutations_response:
    #     print("Hay mutaciones para procesar.")
    # else:
    #     print("No se obtuvieron mutaciones para este estudio.")
    #     return

    # Procesar las mutaciones para agregar datos por gen
    gene_mutation_count = defaultdict(int)
    gene_sample_set = defaultdict(set)

    for mut in mutations_response:
        gene_info = mut.get('gene', {})
        gene = gene_info.get('hugoGeneSymbol')  # Acceder al símbolo del gen correctamente
        sample_id = mut.get('sampleId')
        if gene and sample_id:
            gene_mutation_count[gene] += 1
            gene_sample_set[gene].add(sample_id)

    # Filtrar solo los genes que tienen mutaciones para optimizar
    genes_with_mutations = list(gene_mutation_count.keys())
    # print(f"Total de genes con al menos una mutación: {len(genes_with_mutations)}")

    if not genes_with_mutations:
        print(f"No hay genes con mutaciones para el estudio '{study_id}'.")
        return

    # Construir el DataFrame
    data = []
    for gene in genes_with_mutations:
        total_mutations = gene_mutation_count.get(gene, 0)
        samples_with_mutations = len(gene_sample_set.get(gene, set()))
        frecuencia = (samples_with_mutations / total_samples) * 100 if total_samples > 0 else 0
        data.append({
            'Gene': gene,
            '# Mut': total_mutations,
            '#': samples_with_mutations,
            'Frecuencia (%)': round(frecuencia, 2)
        })

    df = pd.DataFrame(data)

    if max_gen:
        # Ordenar el DataFrame por Frecuencia descendente
        df_sorted = df.sort_values(by="Frecuencia (%)", ascending=False).reset_index(drop=True)
        # Seleccionar los primeros 5000 genes
        df_final = df_sorted.head(max_gen)
        # Guardar el DataFrame en un archivo CSV
        output_filename = f"/home/fernando/Documents/HPC/Final_Protocol/Genes_Protocol/Frecuencias/{study_id}_freqAltas.csv"
        df_final.to_csv(output_filename, index=False)
        # print(f"Datos guardados en '{output_filename}'.")
    else:
        # Guardar el DataFrame en un archivo CSV
        output_filename = f"/home/fernando/Documents/HPC/Final_Protocol/Genes_Protocol/Frecuencias/{study_id}_freqAltas.csv"
        df.to_csv(output_filename, index=False)
        # print(f"Datos guardados en '{output_filename}'.")
