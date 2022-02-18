import pandas as pd
import pkg_resources
import scipy.stats as stats
import statsmodels.api as sm

def load_example_data(omicstype="metabolomics", processed=True):
    """
    Loads example datasets

    :param omicstype: type of omics for example data. Available options are "metabolomics" or "transcriptomics"
    Metabolomics data are from Su et al 2020 https://doi.org/10.1016/j.cell.2020.10.037

    :return: pre-processed omics data matrix consisting of m samples and n entities (metabolites/genes)
    in the form of a pandas DataFrame. Contains one of more metadata columns at the end.
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
    metabolites = matrix.columns.tolist()
    matrix['Target'] = pd.factorize(classes)[0]
    disease = matrix.loc[matrix["Target"] == 0]
    disease.drop(['Target'], axis=1, inplace=True)
    ctrl = matrix.loc[matrix["Target"] != 0]
    ctrl.drop(['Target'], axis=1, inplace=True)
    if testtype == "mwu":
        pvalues = stats.mannwhitneyu(disease, ctrl, axis=0)[1]
    else:
        pvalues = stats.ttest_ind(disease, ctrl)[1]

    padj = sm.stats.multipletests(pvalues, 0.05, method=multiple_correction_method)
    results = pd.DataFrame(zip(metabolites, pvalues, padj[1]),
                           columns=["Entity", "P-value", "P-adjust"])
    return results

def pathwaydf_to_dict(df):
    pathways_df = df.drop("Pathway_name", axis=1)
    pathway_dict = {}

    for pathway in df.index:
        pathway_compounds = list(set(pathways_df.loc[pathway, :].tolist()))
        pathway_compounds = [str(i) for i in pathway_compounds if str(i) != "None"]

        cpds = pathway_compounds[1:]
        if len(cpds) > 1:
            pathway_dict[pathway] = cpds
    return pathway_dict