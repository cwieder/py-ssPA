import pandas as pd
import numpy as np
import sspa.utils as utils

def sspa_svd(mat, pathway_df, min_entity=2):

    """
    Tomfohr et al 2004 SVD/PLAGE method for single sample pathway analysis

    Args:
        mat (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
        Do not include metadata columns
        pathways (pd.DataFrame): Dictionary of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        min_entity (int): minimum number of metabolites mapping to pathways for ssPA to be performed

    Returns:
        pandas DataFrame of pathway scores derived using the PLAGE method. Columns represent pathways and rows represent samples.
    """

    pathways = utils.pathwaydf_to_dict(pathway_df)

    pathway_activities = []

    # Create pathway matrices
    pathway_ids = []
    for pathway, compounds in pathways.items():
        single_pathway_matrix = mat.drop(mat.columns.difference(compounds), axis=1)
        if single_pathway_matrix.shape[1] >= min_entity:
            pathway_ids.append(pathway)
            pathway_mat = single_pathway_matrix.T.values

            # s = singular values
            # u = left singular vector
            # v = right singular vector
            u, s, vh = np.linalg.svd(pathway_mat)

            pathway_activities.append(vh[0])

    pathway_activities_df = pd.DataFrame(pathway_activities, columns=mat.index, index=pathway_ids).T
    return pathway_activities_df