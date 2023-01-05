import numpy as np
import pandas as pd
import sspa.utils as utils
import gseapy

def sspa_gsea(mat, metadata, pathway_df, ranking_metric='signal_to_noise', min_entity=2):
    """Run GSEA using gseapy package by zqfang (https://github.com/zqfang/GSEApy)

    Args:
        mat (pd.DataFrame): dataframe containing input metabolomics data
        metadata (pd.Series): series containing phenotype metadata e.g 'COVID', 'NON-COVID'
        pathway_df (pd.DataFrame): GMT-like pathway dataframe containing compound identifiers
        ranking_metric (str): Ranking metric for molecules in GSEA. Default is signal-to-noise ratio. 
            Other options are 't_test' and see GSEApy package https://github.com/zqfang/GSEApy/blob/2b5419e14615b6fd19a575ff065256dc7099bbec/gseapy/gsea.py#L135 for more options. 
        min_entity (int, optional): minimum number of molecules mapping to pathways for GSEA to be performed. Defaults to 2.
    """
    
    pathway_names = pathway_df["Pathway_name"].to_dict()
    pathways = utils.pathwaydf_to_dict(pathway_df)
    compounds_present = mat.columns.tolist()
    pathways = {k: v for k, v in pathways.items() if len([i for i in compounds_present if i in v]) >= min_entity}

    gsea_res = gseapy.gsea(data=mat.T, 
                 gene_sets=pathways, 
                 cls=metadata,
                 min_size=min_entity,
                 permutation_type='phenotype',
                 permutation_num=1000, # reduce number to speed up test
                 outdir=None,  # do not write output to disk
                 method=ranking_metric)

    res_df = gsea_res.res2d
    res_df = res_df.rename(columns={'Term': 'Pathway_ID', 'NOM p-val': 'P-value', 'FDR q-val': 'P-adjust FDR', 'FWER p-val': 'P-adjust FWER',
     'Lead_genes': 'Leading_edge', 'Gene %': 'Entity %'})
    res_df = res_df.drop(['Name'], axis=1)
    name_col = res_df['Pathway_ID'].map(pathway_names)
    res_df.insert(1, 'Pathway_name', name_col)

    return res_df