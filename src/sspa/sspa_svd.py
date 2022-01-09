import pandas as pd
import numpy as np

def sspa_svd(mat, paths):

    """
    Tomfohr et al 2007 SVD/PLAGE method for single sample pathway analysis

    :param mat: pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
    Do not include metadata columns
    :param pathways: Dictionary of pathway identifiers (keys) and corresponding list of pathway entities (values).
    Entity identifiers must match those in the matrix columns

    :return: pandas DataFrame of pathway scores derived using the PLAGE method. Columns represent pathways and rows represnt samples.
    """

    # Transpose the matrix for SVD

    mat_t = mat.T
    pathway_activities = []

    # Create pathway matrices
    for k, v in paths.items():
        pathway_mat = mat_t.iloc[mat_t.index.isin(v), :]
        pathway_mat = pathway_mat.to_numpy(dtype=float)

        # s = singular values
        # u = left singular vector
        # v = right singular vector
        u, s, vh = np.linalg.svd(pathway_mat)

        pathway_activities.append(vh[0])

    pathway_activities_df = pd.DataFrame(pathway_activities, columns=mat.index, index=paths.keys())
    return pathway_activities_df
