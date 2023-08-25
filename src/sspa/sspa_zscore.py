# Cite Lee et al. 2008
import pandas as pd
import numpy as np
import scipy.stats as stats
import sspa.utils as utils
from sklearn.utils.validation import check_is_fitted
from sklearn.base import BaseEstimator

class sspa_zscore(BaseEstimator):
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

    def __init__(self, pathway_df, min_entity=2):
        self.pathway_df = pathway_df
        self.min_entity = min_entity
        self.pathways = utils.pathwaydf_to_dict(pathway_df)
        self.pathways_filt = {}
        self.fitted_models = []
        self.pathway_ids = []

    def fit(self, X, y=None):
        """
        Fit the model with X.
        
        Args:
            X (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
            Do not include metadata columns
            Returns: 
            self : object
        """
        self.X_ = X
        self.y_ = y

        self.is_fitted_ = True
        return self
    
    def transform(self, X, y=None):
        """
        Transform X.

        Args:
            X (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
            Do not include metadata columns
            Returns: 
            pandas DataFrame of pathway scores derived using the z-score method. Columns represent pathways and rows represent samples.
        """
        check_is_fitted(self, 'is_fitted_')

        scores = []
        pathway_ids = []

        for pathway, compounds in self.pathways.items():
            single_pathway_matrix = X.drop(X.columns.difference(compounds), axis=1)
            if single_pathway_matrix.shape[1] >= self.min_entity:
                pathway_ids.append(pathway)
                pathway_mat = single_pathway_matrix.T.values

                zscores = stats.zscore(pathway_mat, axis=1)
                sum_zscore = np.sum(zscores, axis=0)
                pathway_act = sum_zscore / np.sqrt(pathway_mat.shape[0])
                scores.append(pathway_act)

        pathway_activities_df = pd.DataFrame(scores, columns=X.index, index=pathway_ids).T
        return pathway_activities_df

    def fit_transform(self, X, y=None):
        """
        Fit the model with X and transform X.

        Args:
            X (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
            Do not include metadata columns
            Returns: 
            pandas DataFrame of pathway scores derived using the z-score method. Columns represent pathways and rows represent samples.
        """
        self.fit(X)
        return self.transform(X)

