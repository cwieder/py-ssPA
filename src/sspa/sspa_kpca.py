import pandas as pd
from sklearn.decomposition import KernelPCA
import sspa.utils as utils
from sklearn.utils.validation import check_is_fitted

class sspa_kpca():
    """
    Kernel PCA method for single sample pathway analysis

    Args:
        pathway_df (pd.DataFrame): pandas DataFrame of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        min_entity (int): minimum number of metabolites mapping to pathways for ssPA to be performed

    """
    def __init__(self, pathway_df, min_entity=2):
        self.pathway_df = pathway_df
        self.min_entity = min_entity
        self.pathways = utils.pathwaydf_to_dict(pathway_df)
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

        self._X = X
        for pathway, compounds in self.pathways.items():
            single_pathway_matrix = X.drop(X.columns.difference(compounds), axis=1)
            if single_pathway_matrix.shape[1] >= self.min_entity:
                self.pathway_ids.append(pathway)
                kpca = KernelPCA(n_components=2, kernel="rbf")
                self.fitted_models.append(kpca.fit(single_pathway_matrix))
         
        return self

    def transform(self, X, y=None):
        """
        Transform X.

        Args:
            X (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
            Do not include metadata columns
            Returns: 
            pandas DataFrame of pathway scores derived using the kPCA method. Columns represent pathways and rows represent samples.
        """

        # Check if fit has been called
        check_is_fitted(self)

        # For each fitted model, transform the data
        scores = []
        for x in self.fitted_models:
            new_data = x.transform(X)
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
            pandas DataFrame of pathway scores derived using the kPCA method. Columns represent pathways and rows represent samples.
        """

        self.fit(X)
        return self.transform(X)
        
