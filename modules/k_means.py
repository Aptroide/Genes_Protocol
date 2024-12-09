from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_samples, silhouette_score
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
import pandas as pd

def c_kmeans(d_data, dim_red_name):
    """
    Plot the Elbow Method for KMeans clustering.

    Parameters:
    d_data : array-like
        The data to cluster
    dim_red_name : str
        The name of the dimension reduction method used
    """
    
    wcss=[]
    for i in range(1,6):
        kmeans=KMeans(n_clusters=i, init='k-means++',random_state=42)
        kmeans.fit(d_data)
        wcss.append(kmeans.inertia_)

    plt.plot(range(1,6),wcss)
    plt.title('The Elbow Method for '+dim_red_name)
    plt.xlabel('Number of Clusters')
    plt.ylabel('WCSS')
    plt.show()

def plot_kmeans_silhouette_plotly(proj_2d, df_transpuesto, range_n_clusters):
    """
    Create silhouette analysis plots for KMeans clustering visualization using Plotly
    and add a 'Grupo' column to the original DataFrame based on the best clustering.

    Parameters:
    proj_2d : array-like of shape (n_samples, 2)
        The 2D projected data to cluster
    df_transpuesto : pandas DataFrame
        The original dataframe containing the Cancer labels
    range_n_clusters : array-like
        The range of number of clusters to try

    Returns:
    dict: Dictionary with cluster labels and silhouette scores for each n_clusters value
    pandas.DataFrame: Original DataFrame with an added 'Grupo' column for the best clustering
    """
    results = {}
    summary = []  # Lista para almacenar los resultados

    for n_clusters in range_n_clusters:
        # Inicializar el clusterer
        clusterer = KMeans(n_clusters=n_clusters, random_state=10)
        cluster_labels = clusterer.fit_predict(proj_2d)

        # Calcular las puntuaciones de silueta
        silhouette_avg = silhouette_score(proj_2d, cluster_labels)
        sample_silhouette_values = silhouette_samples(proj_2d, cluster_labels)

        # Almacenar resultados
        results[n_clusters] = {
            'labels': cluster_labels,
            'silhouette_score': silhouette_avg
        }

        # Añadir al resumen
        summary.append({
            'n_clusters': n_clusters,
            'silhouette_avg': silhouette_avg
        })

        print(f"For n_clusters = {n_clusters}, The average silhouette_score is : {silhouette_avg:.3f}")

        # Crear subplots
        fig = make_subplots(rows=1, cols=2,
                            subplot_titles=("Silhouette Plot", "Clustered Data"),
                            column_widths=[0.4, 0.6])

        # Silhouette plot
        y_lower = 10
        for i in range(n_clusters):
            ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            # Añadir área rellenada para cada cluster
            y_pts = np.arange(y_lower, y_upper)
            fig.add_trace(
                go.Scatter(
                    x=ith_cluster_silhouette_values,
                    y=y_pts,
                    fill='tozerox',
                    name=f'Cluster {i}',
                    showlegend=False
                ),
                row=1, col=1
            )

            y_lower = y_upper + 10

        # Añadir línea vertical para la puntuación de silueta promedio
        fig.add_vline(x=silhouette_avg, line_dash="dash", line_color="red",
                     annotation_text=f"Average silhouette score: {silhouette_avg:.3f}",
                     row=1, col=1)

        # Scatter plot con clusters
        centers = clusterer.cluster_centers_

        # Crear scatter plot coloreado por clusters
        fig.add_trace(
            go.Scatter(
                x=proj_2d[:, 0],
                y=proj_2d[:, 1],
                mode='markers',
                marker=dict(
                    color=cluster_labels,
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title='Cluster')
                ),
                text=[f'Cluster {label}<br>Cancer: {cancer}'
                      for label, cancer in zip(cluster_labels, df_transpuesto.Cancer)],
                hoverinfo='text',
                name='Data points'
            ),
            row=1, col=2
        )

        # Añadir centros de clusters
        fig.add_trace(
            go.Scatter(
                x=centers[:, 0],
                y=centers[:, 1],
                mode='markers+text',
                marker=dict(
                    color='white',
                    size=15,
                    line=dict(color='black', width=2),
                    symbol='diamond'
                ),
                text=[f'Center {i}' for i in range(n_clusters)],
                name='Centroids',
                showlegend=True
            ),
            row=1, col=2
        )

        # Actualizar layout
        fig.update_layout(
            title_text=f"Silhouette Analysis for KMeans Clustering (n_clusters = {n_clusters})",
            height=600,
            width=1200,
            showlegend=True
        )

        # Actualizar ejes
        fig.update_xaxes(title_text="Silhouette coefficient values", range=[-0.1, 1], row=1, col=1)
        fig.update_yaxes(showticklabels=False, title_text="Cluster label", row=1, col=1)
        fig.update_xaxes(title_text="First dimension", row=1, col=2)
        fig.update_yaxes(title_text="Second dimension", row=1, col=2)

        fig.show()

    # Convertir el resumen a DataFrame
    results_df = pd.DataFrame(summary)

    # Identificar el n_clusters con la mejor puntuación de silueta
    best_n_clusters = results_df.loc[results_df['silhouette_avg'].idxmax(), 'n_clusters']
    best_silhouette = results_df['silhouette_avg'].max()
    print(f"\nEl mejor número de clusters es {best_n_clusters} con una puntuación de silueta de {best_silhouette:.3f}")

    # Obtener las etiquetas del mejor clustering
    best_labels = results[best_n_clusters]['labels']
    # Crear una copia del DataFrame original para no modificarlo en su lugar
    df_with_clusters = df_transpuesto.copy()

    # Agregar la columna 'Grupo' con las etiquetas del mejor clustering
    df_with_clusters['Cluster'] = best_labels


    # Reordenar las columnas para que 'Cluster' esté justo después de 'Cancer'
    cols = list(df_with_clusters.columns)
    cancer_idx = cols.index('Cancer')
    # Insertar 'Cluster' después de 'Cancer'
    cols.insert(cancer_idx + 1, cols.pop(cols.index('Cluster')))
    df_with_clusters = df_with_clusters[cols]


    return results, results_df, df_with_clusters