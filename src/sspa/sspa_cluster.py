import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


def sspa_cluster(mat, pathways, min_entity=2, projection=False):
    pathway_matrices = []
    pathway_ids = []
    for pathway, compounds in pathways.items():
        single_pathway_matrix = mat.drop(mat.columns.difference(compounds), axis=1)
        if single_pathway_matrix.shape[1] >= min_entity:
            pathway_matrices.append(single_pathway_matrix.values)
            pathway_ids.append(pathway)

    if projection:
        scores = []
        for m in pathway_matrices:
            kmeans = KMeans(n_clusters=2).fit(m)
            centroids1 = kmeans.cluster_centers_[0]
            centroids2 = kmeans.cluster_centers_[1]

            vec = centroids1 - centroids2
            unit_vec = vec / np.linalg.norm(vec)
            proj_data = unit_vec.dot(m.T)
            scores.append(proj_data)

    else:
        scores = []
        for m in pathway_matrices:
            kmeans = KMeans(n_clusters=2)
            new_data = kmeans.fit_transform(m)
            scores.append(new_data[:, 0])

    scores_df = pd.DataFrame(scores, columns=mat.index, index=pathways.keys())
    return scores_df