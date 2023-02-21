import pandas as pd
import numpy as np
import pkg_resources
import scipy.stats as stats
import statsmodels.api as sm

def load_example_data(omicstype="metabolomics", processed=True):
    """
    Loads example datasets

    Args:
        omicstype (str): type of omics for example data. 
            Available options are "metabolomics" or "transcriptomics". 
            Metabolomics data are from Su et al 2020 https://doi.org/10.1016/j.cell.2020.10.037.
            Transcriptomics data - TO BE IMPLEMENTED
        processed (bool): Load processed (normalised, scaled) or raw data

    Returns:
        pre-processed omics data matrix consisting of m samples and n entities (metabolites/genes) in the form of a pandas DataFrame. 
        Contains one of more metadata columns at the end.
    """

    if omicstype == "metabolomics":
        if processed:
            stream = pkg_resources.resource_stream(__name__, 'example_data/Su_covid_metabolomics_processed.csv')
            f = pd.read_csv(stream, index_col=0, encoding='latin-1')
            return f
        else:
            stream = pkg_resources.resource_stream(__name__, 'example_data/Su_metab_data_raw.csv')
            f = pd.read_csv(stream, index_col=0, encoding='latin-1')
            return f


def t_tests(matrix, classes, multiple_correction_method, testtype="ttest"):
    """
    Performs two-sample independent t-tests

    Args:
        matrix (pd.DataFrame): processed sample-by-compound metabolomics dataframe
        classes (pd.Series): pandas series containing phenotype metadata (e.g. 'COVID', 'NON-COVID')
        multiple_correction_method (str): see https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for options
        testtype (str): Default is t-test, "mwu" also available to implement the Mann Whitney U test

    Returns:
        pd.DataFrame containing p-values and corrected p-values for each metabolite
    """
    matrix_copy = matrix.copy(deep=True)
    metabolites = matrix_copy.columns.tolist()
    matrix_copy['Target'] = pd.factorize(classes)[0]

    # Check user has only input two classes of samples 
    if len(set(matrix_copy['Target'].tolist())) > 2:
        raise ValueError('More than two metadata classes detected. Only two metadata classes are supported in ORA.')

    disease = matrix_copy[matrix_copy["Target"] == 0]
    disease = disease.drop(['Target'], axis=1)
    ctrl = matrix_copy[matrix_copy["Target"] != 0]
    ctrl = ctrl.drop(['Target'], axis=1)
    matrix_copy = matrix_copy.drop(['Target'], axis=1)

    # disease = matrix.loc[matrix["Target"] == 0]
    # disease.drop(['Target'], axis=1, inplace=True)
    # ctrl = matrix.loc[matrix["Target"] != 0]
    # ctrl.drop(['Target'], axis=1, inplace=True)
    # matrix = matrix.drop(['Target'], axis=1)
    
    if testtype == "mwu":
        pvalues = stats.mannwhitneyu(disease, ctrl, axis=0)[1]
    else:
        pvalues = stats.ttest_ind(disease, ctrl)[1]

    padj = sm.stats.multipletests(pvalues, 0.05, method=multiple_correction_method)
    results = pd.DataFrame(zip(metabolites, pvalues, padj[1]),
                           columns=["Entity", "P-value", "P-adjust"])
    return results


def pathwaydf_to_dict(df):
    """
    Converts pathway dataframe to dictionary, with pathway IDs as keys and metabolite lists as values
    Args:
        df (pd.DataFrame): Pandas DataFrame containing pathways 
    Returns: 
        python dict pathway representation
    """
    pathways_df = df.drop(["Pathway_name"], axis=1)
    pathway_dict = {}

    for pathway in pathways_df.index:
        pathway_compounds = list(set(pathways_df.loc[pathway, :].tolist()))
        pathway_compounds = [str(i) for i in pathway_compounds if str(i) not in ["None", np.nan, 'nan']]

        if len(pathway_compounds) > 1:
            pathway_dict[pathway] = pathway_compounds
    return pathway_dict