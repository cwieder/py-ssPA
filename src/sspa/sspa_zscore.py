# Cite Lee et al. 2008
import pandas as pd
import numpy as np
import scipy.stats as stats
import sspa.utils as utils

def sspa_zscore(mat, pathway_df, min_entity=2):
    """
    Lee at al 2008 z-score method for single sample pathway analysis

    Args:
        mat (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
        Do not include metadata columns
        pathways (pd.DataFrame): Dictionary of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        min_entity (int): minimum number of metabolites mapping to pathways for ssPA to be performed


    Returns:
        pandas DataFrame of pathway scores derived using the z-score method. Columns represent pathways and rows represent samples.
    """

    pathways = utils.pathwaydf_to_dict(pathway_df)

    pathway_activities = []
    pathway_ids = []
    for pathway, compounds in pathways.items():
        single_pathway_matrix = mat.drop(mat.columns.difference(compounds), axis=1)
        if single_pathway_matrix.shape[1] >= min_entity:
            pathway_mat = single_pathway_matrix.T.values
            pathway_ids.append(pathway)

            zscores = stats.zscore(pathway_mat, axis=1)
            sum_zscore = np.sum(zscores, axis=0)
            pathway_act = sum_zscore / np.sqrt(pathway_mat.shape[0])
            pathway_activities.append(pathway_act)

    pathway_activities_df = pd.DataFrame(pathway_activities, columns=mat.index, index=pathway_ids).T
    return pathway_activities_df

