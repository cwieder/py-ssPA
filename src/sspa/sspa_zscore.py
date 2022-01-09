# Cite Lee et al. 2008
import pandas as pd
import numpy as np
import scipy.stats as stats

def sspa_zscore(mat, pathways):
    """
    Lee at al 2008 z-score method for single sample pathway analysis

    :param mat: pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
    Do not include metadata columns
    :param pathways: Dictionary of pathway identifiers (keys) and corresponding list of pathway entities (values).
    Entity identifiers must match those in the matrix columns

    :return: pandas DataFrame of pathway scores derived using the z-score method. Columns represent pathways and rows represnt samples.
    """

    mat_t = mat.T
    pathway_activities = []

    for k, v in pathways.items():
        pathway_mat = mat_t.iloc[mat_t.index.isin(v), :]
        pathway_mat = pathway_mat.to_numpy(dtype=float)
        zscores = stats.zscore(pathway_mat, axis=1)
        sum_zscore = np.sum(zscores, axis=0)
        pathway_act = sum_zscore / np.sqrt(pathway_mat.shape[0])
        pathway_activities.append(pathway_act)

    pathway_activities_df = pd.DataFrame(pathway_activities, columns=mat.index, index=pathways.keys())
    return pathway_activities_df

