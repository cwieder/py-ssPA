import pandas as pd
from sklearn.decomposition import KernelPCA
import sspa.utils as utils

def sspa_kpca(mat, pathway_df, min_entity=2):

    """
    Kernel PCA method for single sample pathway analysis

    Args:
        mat (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
        Do not include metadata columns
        pathways (pd.DataFrame): Dictionary of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        min_entity (int): minimum number of metabolites mapping to pathways for ssPA to be performed

    Returns:
        pandas DataFrame of pathway scores derived using the kPCA method. Columns represent pathways and rows represent samples.
    """

    pathways = utils.pathwaydf_to_dict(pathway_df)

    pathway_matrices = []
    pathway_ids = []
    for pathway, compounds in pathways.items():
        single_pathway_matrix = mat.drop(mat.columns.difference(compounds), axis=1)
        if single_pathway_matrix.shape[1] >= min_entity:
            pathway_matrices.append(single_pathway_matrix.values)
            pathway_ids.append(pathway)

    scores = []
    for n, m in enumerate(pathway_matrices):
        kpca = KernelPCA(n_components=2, kernel="rbf")
        new_data = kpca.fit_transform(m)
        scores.append(new_data[:, 0])
    scores_df = pd.DataFrame(scores, columns=mat.index, index=pathway_ids).T
    return scores_df