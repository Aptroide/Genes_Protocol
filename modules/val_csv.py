import os
import pandas as pd

def validate_csv_file(file_path):
    """
    Valida un archivo CSV asegurando que cumpla con ciertas restricciones:
    1. El nombre del archivo debe seguir un patrón específico.
    2. El archivo debe contener las columnas 'Gene' y 'Frecuencia (%)'.
    3. Los tipos de datos deben ser correctos.
    4. La columna 'Gene' debe tener valores únicos.
    5. La columna 'Frecuencia (%)' debe estar en el rango [0, 100] y no tener valores nulos.

    Parámetros:
    file_path (str): Ruta completa del archivo CSV a validar.

    Retorna:
    bool: True si el archivo cumple todas las restricciones, False en caso contrario.
    """
    # Verificar si el archivo existe
    if not os.path.exists(file_path):
        print(f"Error: El archivo {file_path} no existe.")
        return False

    # Verificar que el archivo tenga el sufijo adecuado
    if not file_path.endswith('_freqAltas.csv'):
        print(f"Error: El archivo debe tener el sufijo '_freqAltas.csv'.")
        return False

    # Leer el archivo CSV
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error al leer el archivo CSV: {e}")
        return False

    # Verificar que las columnas necesarias estén presentes
    required_columns = ['Gene', 'Frecuencia (%)']
    if not all(col in df.columns for col in required_columns):
        print(f"Error: Las columnas {required_columns} son necesarias.")
        return False

    # Verificar que la columna 'Gene' sea única
    if df['Gene'].duplicated().any():
        print("Error: La columna 'Gene' contiene valores duplicados.")
        return False

    # Verificar que 'Frecuencia (%)' sea numérica y esté en el rango [0, 100]
    if not pd.api.types.is_numeric_dtype(df['Frecuencia (%)']):
        print("Error: La columna 'Frecuencia (%)' debe ser numérica.")
        return False

    if not df['Frecuencia (%)'].between(0, 100).all():
        print("Error: Los valores de 'Frecuencia (%)' deben estar en el rango de 0 a 100.")
        return False

    # Verificar que no haya valores nulos en las columnas esenciales
    if df[required_columns].isnull().any().any():
        print("Error: No deben existir valores nulos en las columnas 'Gene' o 'Frecuencia (%).")
        return False

    return True