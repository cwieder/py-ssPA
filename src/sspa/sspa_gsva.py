import pandas as pd
import sspa.utils as utils
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# for rpy2
base = importr('base')

def sspa_gsva(mat, pathway_df, min_entity=2):

    """
    Hanzelmann et al GSVA method for single sample pathway analysis. 
    This is an rpy2 wrapper script to run the R implementation of GSVA.

    :param mat: pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
    Do not include metadata columns
    :param pathways: Dictionary of pathway identifiers (keys) and corresponding list of pathway entities (values).
    Entity identifiers must match those in the matrix columns

    :return: pandas DataFrame of pathway scores derived using the GSVA method. Columns represent pathways and rows represnt samples.
    """

    pathways = utils.pathwaydf_to_dict(pathway_df)
    compounds_present = mat.columns.tolist()
    pathways = {k: v for k, v in pathways.items() if len([i for i in compounds_present if i in v]) >= min_entity}

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_mat = ro.conversion.py2rpy(mat.T)
    r_mat = base.as_matrix(r_mat)  # abundance matrix
    row_vec = base.as_character(mat.columns.tolist())
    r_mat.rownames = row_vec
    r_list = ro.ListVector(pathways)  # pathways
    gsva_r = importr('GSVA')
    gsva_res = gsva_r.gsva(r_mat, r_list)
    with localconverter(ro.default_converter + pandas2ri.converter):
        gsva_df = ro.conversion.rpy2py(gsva_res)
    gsva_res_df = pd.DataFrame(gsva_df.T, columns=pathways.keys(), index=mat.index.tolist())

    return gsva_res_df