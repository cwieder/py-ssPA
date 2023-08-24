import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import sspa.utils as utils
from sklearn.utils.validation import check_is_fitted
from sklearn.base import BaseEstimator

class sspa_ssClustPA(BaseEstimator):
    """
    K-means based clustering method for single sample pathway analysis

    Args:
        mat (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
        Do not include metadata columns
        pathways (pd.DataFrame): Dictionary of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        min_entity (int): minimum number of metabolites mapping to pathways for ssPA to be performed

    Returns:
        pandas DataFrame of pathway scores derived using the ssClustPA/(proj) method. Columns represent pathways and rows represent samples.
    """

    def __init__(self, pathway_df, min_entity=2, random_state=0):
        self.pathway_df = pathway_df
        self.min_entity = min_entity
        self.pathways = utils.pathwaydf_to_dict(pathway_df)
        self.pathways_filt = {}
        self.fitted_models = []
        self.pathway_ids = []
        self.random_state = random_state

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

        for pathway, compounds in self.pathways.items():
            single_pathway_matrix = X.drop(X.columns.difference(compounds), axis=1)
            if single_pathway_matrix.shape[1] >= self.min_entity:
                self.pathway_ids.append(pathway)
                kmeans = KMeans(n_clusters=2, random_state=self.random_state, n_init='auto')
                self.fitted_models.append(kmeans.fit(single_pathway_matrix.to_numpy()))

        self.pathways_filt = {k: v for k, v in self.pathways.items() if k in self.pathway_ids}
        self.is_fitted_ = True
        return self
    
    def transform(self, X, y=None):
        """
        Transform X.

        Args:
            X (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
            Do not include metadata columns
            Returns: 
            pandas DataFrame of pathway scores derived using the ssClustPA/(proj) method. Columns represent pathways and rows represent samples.
        """
        check_is_fitted(self, 'is_fitted_')

        scores = []
        for n, compounds in enumerate(self.pathways_filt.values()):
            single_pathway_matrix = X.drop(X.columns.difference(compounds), axis=1)
            centroids1 = self.fitted_models[n].cluster_centers_[0]
            centroids2 = self.fitted_models[n].cluster_centers_[1]
            vec = centroids1 - centroids2
            unit_vec = vec / np.linalg.norm(vec)
            proj_data = unit_vec.dot(single_pathway_matrix.T)
            scores.append(proj_data)

        scores_df = pd.DataFrame(scores, columns=X.index, index=self.pathway_ids).T
        return scores_df
    
    def fit_transform(self, X, y=None):
        """
        Fit the model with X and transform X.
        """
        self.fit(X)
        return self.transform(X)
    
    def fit_transform_(self, X, y=None):
        """
        Fit the model with X and transform X.

        Args:
            X (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
            Do not include metadata columns
            Returns: 
            pandas DataFrame of pathway scores derived using the ssClustPA/(proj) method. Columns represent pathways and rows represent samples.
        """

        self._X = X
        self._y = y
        scores = []
        for pathway, compounds in self.pathways.items():
            single_pathway_matrix = X.drop(X.columns.difference(compounds), axis=1)
            if single_pathway_matrix.shape[1] >= self.min_entity:
                self.pathway_ids.append(pathway)

                kmeans = KMeans(n_clusters=2, random_state=self.random_state, n_init='auto').fit(single_pathway_matrix)
                centroids1 = kmeans.cluster_centers_[0]
                centroids2 = kmeans.cluster_centers_[1]

                vec = centroids1 - centroids2
                unit_vec = vec / np.linalg.norm(vec)
                proj_data = unit_vec.dot(single_pathway_matrix.T)
                scores.append(proj_data)
        scores_df = pd.DataFrame(scores, columns=X.index, index=self.pathway_ids).T
        self.is_fitted_ = True
        return scores_df

