import pandas as pd
import sspa.utils as utils
import gseapy
# import rpy2.robjects as ro
# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.conversion import localconverter

# # for rpy2
# base = importr('base')

def sspa_ssGSEA(mat, pathway_df, min_entity=2):

    """
    Barbie et al ssGSEA method for single sample pathway analysis. 
    Uses the ssgsea function of the gseapy package (https://github.com/zqfang/GSEApy) as a backend. 

    Args:
        mat (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
        Do not include metadata columns
        pathways (pd.DataFrame): Dictionary of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        min_entity (int): minimum number of metabolites mapping to pathways for ssPA to be performed


    Returns:
        pandas DataFrame of pathway scores derived using the ssGSEA method. Columns represent pathways and rows represent samples.
    """

    pathways = utils.pathwaydf_to_dict(pathway_df)
    compounds_present = mat.columns.tolist()
    pathways = {k: v for k, v in pathways.items() if len([i for i in compounds_present if i in v]) >= min_entity}

    ssgsea_res = gseapy.ssgsea(data=mat.T,
               gene_sets=pathways,
               min_size=min_entity,
               outdir=None,
               sample_norm_method='rank', # choose 'custom' will only use the raw value of `data`
               no_plot=True)

    ssgsea_scores = ssgsea_res.res2d.pivot(index='Term', columns='Name', values='NES').T
    res_df = pd.DataFrame(ssgsea_scores, index=mat.index, columns=pathways.keys())
    res_df = res_df.astype(float)
    return res_df

    # with localconverter(ro.default_converter + pandas2ri.converter):
    #     r_mat = ro.conversion.py2rpy(mat.T)
    # r_mat = base.as_matrix(r_mat)  # abundance matrix
    # row_vec = base.as_character(mat.columns.tolist())
    # r_mat.rownames = row_vec
    # r_list = ro.ListVector(pathways)  # pathways
    # gsva_r = importr('GSVA')
    # ssgsea_res = gsva_r.gsva(r_mat, r_list, method='ssgsea')
    # with localconverter(ro.default_converter + pandas2ri.converter):
    #     ssgsea_df = ro.conversion.rpy2py(ssgsea_res)
    # ssgsea_res_df = pd.DataFrame(ssgsea_df.T, columns=pathways.keys(), index=mat.index.tolist())

    # return ssgsea_res_df