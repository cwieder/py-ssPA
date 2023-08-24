import pandas as pd
import sspa.utils as utils
import gseapy
from sklearn.utils.validation import check_is_fitted
from sklearn.base import BaseEstimator

class sspa_ssGSEA(BaseEstimator):
    """
    Barbie et al ssGSEA method for single sample pathway analysis. 

    Uses the ssgsea function of the gseapy package (https://github.com/zqfang/GSEApy) as a backend. 

    All credit for ssGSEA code goes to developers of the GSEAPY python package (credit: 
    Zhuoqing Fang, Xinyuan Liu, Gary Peltz, GSEApy: 
    a comprehensive package for performing gene set enrichment analysis in Python,
    Bioinformatics, 2022;, btac757, https://doi.org/10.1093/bioinformatics/btac757)
    
    Args:
        pathway_df (pd.DataFrame): pandas DataFrame of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        min_entity (int): minimum number of metabolites mapping to pathways for ssPA to be performed


    Returns:
        pandas DataFrame of pathway scores derived using the ssGSEA method. Columns represent pathways and rows represent samples.
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
            pandas DataFrame of pathway scores derived using the ssGSEA method. Columns represent pathways and rows represent samples.
        """
        check_is_fitted(self, 'is_fitted_')

        ssgsea_res = gseapy.ssgsea(data=X.T,
                gene_sets=self.pathways,
                min_size=self.min_entity,
                outdir=None,
                sample_norm_method='rank', # choose 'custom' will only use the raw value of `data`
                no_plot=True)
        ssgsea_scores = ssgsea_res.res2d.pivot(index='Term', columns='Name', values='NES').T
        res_df = pd.DataFrame(ssgsea_scores, index=X.index)
        res_df = res_df.astype(float)
        return res_df
    
    def fit_transform(self, X, y=None):
        """
        Fit the model with X and transform X.

        Args:
            X (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
            Do not include metadata columns
            Returns: 
            pandas DataFrame of pathway scores derived using the ssGSEA method. Columns represent pathways and rows represent samples.
        """
        self.fit(X)
        return self.transform(X)
