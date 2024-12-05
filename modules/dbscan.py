import itertools
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_samples, silhouette_score
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px

def get_scores_and_labels(combinations, X):
  scores =[]
  all_labels_list = []

  for i, (eps, num_samples) in enumerate(combinations):
    dbscan = DBSCAN(eps=eps, min_samples=num_samples)
    labels = dbscan.fit_predict(X)
    labels_set = set(labels)
    num_clusters = len(labels_set) - (1 if -1 in labels_set else 0)
    if num_clusters < 2 or num_clusters > 50:
      scores.append(-10)
      all_labels_list.append('bad')
      c = (eps, num_samples)
      # print(f"Combination {c} on iteration {i+1} of {len(combinations)} has {num_clusters} clusters. Moving on")
      continue
    scores.append(silhouette_score(X, labels))
    all_labels_list.append(labels)
    # print(f"Index: {i}, Score: {scores[-1]}, NumClusters: {num_clusters}")

  best_index = np.argmax(scores)
  best_parameters = combinations[best_index]
  best_labels = all_labels_list[best_index]
  best_score = scores[best_index]

  return {'best_epsilons': best_parameters[0],
          'best_min_samples': best_parameters[1],
          'best_score': best_score}


def optimize_dbscan(proj2d, df_transpuesto):
    # Definir los valores para eps y min_samples
    eps_values = np.linspace(0.1, 1.3, num=100)
    min_samples_values = np.arange(2, 20, step=3)

    # Generar todas las combinaciones posibles de (eps, min_samples)
    combinations = list(itertools.product(eps_values, min_samples_values))

    # Obtener los mejores valores de epsilon y min_samples usando get_scores_and_labels
    best_dict = get_scores_and_labels(combinations, proj2d)

    # Extraer los mejores valores de epsilon y min_samples
    eps_values = [best_dict['best_epsilons']]
    min_samples_values = [best_dict['best_min_samples']]

    # Obtener los resultados y las etiquetas de DBSCAN mediante plot_dbscan_silhouette_plotly
    results, dbscan_labels = plot_dbscan_silhouette_plotly(proj2d, df_transpuesto, eps_values, min_samples_values)

    return results, dbscan_labels

def plot_dbscan_silhouette_plotly(proj_2d, df_transpuesto, eps_values, min_samples_values):
    """
    Create silhouette analysis plots for DBSCAN clustering visualization using Plotly.
    This function returns the DataFrame corresponding to the best silhouette score.
    """

    results = {}
    best_silhouette_score = -1  # Start with a very low value
    best_df = None  # Placeholder for the DataFrame with the best silhouette score
    best_eps = None
    best_min_samples = None

    # Crear una copia del DataFrame para no modificar el original
    df_augmented = df_transpuesto.copy()
        # Inicializar el clusterizador

    for eps in eps_values:
        for min_samples in min_samples_values:
            # Inicializar el clusterizador
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            cluster_labels = dbscan.fit_predict(proj_2d)

            # Calcular las pun tuaciones de silueta
            if len(set(cluster_labels)) > 1:  # La puntuación de silueta no está definida para un solo clúster
                silhouette_avg = silhouette_score(proj_2d, cluster_labels)
                sample_silhouette_values = silhouette_samples(proj_2d, cluster_labels)
            else:
                silhouette_avg = -1
                sample_silhouette_values = np.zeros(len(proj_2d))

            # Almacenar los resultados
            results[(eps, min_samples)] = {
                'labels': cluster_labels,
                'silhouette_score': silhouette_avg
            }

            print(f"For eps = {eps}, min_samples = {min_samples}, The average silhouette_score is : {silhouette_avg:.3f}")

            # Verificar si el puntaje de silueta es el mejor
            if silhouette_avg > best_silhouette_score:
                best_silhouette_score = silhouette_avg
                best_eps = eps
                best_min_samples = min_samples

                # Actualizar el DataFrame con las etiquetas de los clústeres
                df_augmented['Cluster'] = cluster_labels

                # Reordenar las columnas para que 'Cluster' esté justo después de 'Cancer'
                cols = list(df_augmented.columns)
                cancer_idx = cols.index('Cancer')
                # Insertar 'Cluster' después de 'Cancer'
                cols.insert(cancer_idx + 1, cols.pop(cols.index('Cluster')))
                df_augmented = df_augmented[cols]

            # Crear subplots
            fig = make_subplots(rows=1, cols=2,
                                subplot_titles=("Silhouette Plot", "Clustered Data"),
                                column_widths=[0.4, 0.6])

            # Silhouette plot
            y_lower = 10
            unique_clusters = np.unique(cluster_labels)
            for i in unique_clusters:
                if i == -1:
                    continue  # Omitir puntos de ruido
                ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
                ith_cluster_silhouette_values.sort()

                size_cluster_i = ith_cluster_silhouette_values.shape[0]
                y_upper = y_lower + size_cluster_i

                # Añadir área rellena para los valores de silueta de cada clúster
                y_pts = np.arange(y_lower, y_upper)
                fig.add_trace(
                    go.Scatter(
                        x=ith_cluster_silhouette_values,
                        y=y_pts,
                        fill='tozerox',
                        mode='none',
                        name=f'Cluster {i}',
                        showlegend=False
                    ),
                    row=1, col=1
                )

                y_lower = y_upper + 10

            # Añadir línea vertical para la puntuación promedio de silueta
            fig.add_vline(x=silhouette_avg, line_dash="dash", line_color="red",
                          annotation_text=f"Average silhouette score: {silhouette_avg:.3f}",
                          row=1, col=1)

            # Scatter plot con clústeres
            scatter_df = pd.DataFrame({
                'Dim1': proj_2d[:, 0],
                'Dim2': proj_2d[:, 1],
                'Cancer Type': df_transpuesto['Cancer'],
                'Cluster': cluster_labels.astype(str)  # Convertir a string para manejar el color
            })

            # Definir colores para los clústeres, incluyendo ruido
            unique_clusters_str = scatter_df['Cluster'].unique()
            color_discrete_map = {str(cluster): px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)]
                                  for i, cluster in enumerate(unique_clusters_str)}


            fig_scatter = px.scatter(
                scatter_df,
                x='Dim1', y='Dim2',
                color='Cluster',
                color_discrete_map=color_discrete_map,
                hover_data={'Cancer Type': True, 'Cluster': True},
                labels={'Cluster': 'Cluster'}
            )

            # Añadir las trazas del scatter plot al subplot
            for trace in fig_scatter.data:
                fig.add_trace(trace, row=1, col=2)

            # Actualizar el layout
            fig.update_layout(
                title_text=f"Silhouette Analysis for DBSCAN Clustering (eps = {eps}, min_samples = {min_samples})",
                height=600,
                width=1200
            )

            # Actualizar los ejes
            fig.update_xaxes(title_text="Silhouette coefficient values", range=[-0.1, 1], row=1, col=1)
            fig.update_yaxes(showticklabels=False, title_text="Cluster label", row=1, col=1)
            fig.update_xaxes(title_text="First dimension", row=1, col=2)
            fig.update_yaxes(title_text="Second dimension", row=1, col=2)

            fig.show()

    # Después del loop, ahora df_augmented tiene las etiquetas de clúster correspondientes al mejor silueta promedio
    print(f"Best silhouette score found: {best_silhouette_score:.3f} for eps = {best_eps}, min_samples = {best_min_samples}")

    # Retornar los resultados y el DataFrame que corresponde a la mejor puntuación de silueta
    return results, df_augmented