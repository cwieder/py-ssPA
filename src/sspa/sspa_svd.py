import pandas as pd
import numpy as np
import sspa.utils as utils
from sklearn.decomposition import PCA
from sklearn.utils.validation import check_is_fitted
from sklearn.base import BaseEstimator


class sspa_SVD(BaseEstimator):
    """
    Tomfohr et al 2005 PLAGE (SVD) method for single sample pathway analysis

    Args:
        pathway_df (pd.DataFrame): pandas DataFrame of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        min_entity (int): minimum number of metabolites mapping to pathways for ssPA to be performed

    """
    def __init__(self, pathway_df, min_entity=2, random_state=0):
        self.pathway_df = pathway_df
        self.min_entity = min_entity
        self.pathways = utils.pathwaydf_to_dict(pathway_df)
        self.pathways_filt = {}
        self.fitted_models = []
        self.pathway_ids = []
        self.random_state = random_state
        self.molecular_importance = {}

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
                pca = PCA(n_components=1, random_state=self.random_state)
                self.fitted_models.append(pca.fit(single_pathway_matrix.to_numpy()))

                # use loadings for PC1 molecular importances within the pathway
                loadings = pca.components_[0]
                self.molecular_importance[pathway] = pd.DataFrame(loadings, index=single_pathway_matrix.columns, columns=['PC1_Loadings'])

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
            self : object
        """
    
            # Check if fit has been called
        check_is_fitted(self, 'is_fitted_')

        # For each fitted model, transform the data
        scores = []
        for n, compounds in enumerate(self.pathways_filt.values()):
            single_pathway_matrix = X.drop(X.columns.difference(compounds), axis=1)
            new_data = self.fitted_models[n].transform(single_pathway_matrix.to_numpy())
            scores.append(new_data[:, 0])
        scores_df = pd.DataFrame(scores, columns=X.index, index=self.pathway_ids).T

        return scores_df
    
    def fit_transform(self, X, y=None):
            
            """
            Fit the model with X and transform X.
    
            Args:
                X (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
                Do not include metadata columns
                Returns: 
                self : object
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
            self : object
        """
        self.X_ = X
        self.y_ = y

        scores = []
        for pathway, compounds in self.pathways.items():
            single_pathway_matrix = X.drop(X.columns.difference(compounds), axis=1)
            if single_pathway_matrix.shape[1] >= self.min_entity:
                self.pathway_ids.append(pathway)
                pca = PCA(n_components=1, random_state=self.random_state)
                scores.append(pca.fit_transform(single_pathway_matrix)[:, 0])

                # use loadings for PC1 molecular importances within the pathway
                loadings = pca.components_[0]
                self.molecular_importance[pathway] = pd.DataFrame(loadings, index=single_pathway_matrix.columns, columns=['PC1_Loadings'])

        scores_df = pd.DataFrame(scores, columns=X.index, index=self.pathway_ids).T
        self.is_fitted_ = True
        return scores_df
