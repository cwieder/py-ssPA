import numpy as np
import pandas as pd
import sspa.utils as utils
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# for rpy2
base = importr('base')

def sspa_fgsea(mat, metadata, pathway_df, min_entity=2):
    """
    Function for gene/metabolite set enrichment analysis using fGSEA method (Korotkevich et al.)
    Args:
        mat (pd.DataFrame): dataframe containing input metabolomics data
        metadata (pd.Series): series containing phenotype metadata e.g 'COVID', 'NON-COVID'
        pathway_df (pd.DataFrame): pathway dataframe containing compound identifiers
        min_entity (int): minimum number of metabolites mapping to a pathway for it to be tested
    Returns:
        DataFrame of GSEA results for each pathway, p-value, q-value, ES, NES, leading egde genes/metabolites
    """

    pathway_names = pathway_df["Pathway_name"].to_dict()
    pathways = utils.pathwaydf_to_dict(pathway_df)
    compounds_present = mat.columns.tolist()
    pathways = {k: v for k, v in pathways.items() if len([i for i in compounds_present if i in v]) >= min_entity}

    # Get rankings - SNR
    mat['Target'] = pd.factorize(metadata)[0]
    # Check user has only input two classes of samples 
    if len(set(mat['Target'].tolist())) > 2:
        raise ValueError('More than two metadata classes detected. Only two metadata classes are supported in GSEA.')

    class_a = mat.loc[mat["Target"] == 0]
    class_a.drop(['Target'], axis=1, inplace=True)
    class_b = mat.loc[mat["Target"] != 0]
    class_b.drop(['Target'], axis=1, inplace=True)

    a_mean = np.mean(class_a, axis=0)
    b_mean = np.mean(class_b, axis=0)
    a_std = np.std(class_a, axis=0)
    b_std = np.std(class_b, axis=0)

    means = a_mean.subtract(b_mean)
    stds = a_std + b_std

    snr = means / stds
    snr_dict = snr.to_dict()
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_list = ro.ListVector(pathways)  # pathways
        r_ranks = ro.ListVector(snr_dict) # ranks 

    ro.r('''
    unl <- function(v){
    test <- unlist(v)
    return(test)}
    ''')

    r_unl = ro.globalenv['unl']
    r_ranks = r_unl(r_ranks)

    fgsea_r = importr('fgsea')
    fgsea_res = fgsea_r.fgsea(pathways=r_list, stats=r_ranks)
    
    df = pd.DataFrame(fgsea_res)
    df = df.T
    df = df.apply(pd.to_numeric, errors="ignore")
    df.columns = ["ID", "P-value", "P-adjust", "log2err", "ES", "NES", "coverage", "leadingEdge"]
    df["Pathway_name"] = df["ID"].map(pathway_names)
    
    return df
