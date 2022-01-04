import pandas as pd
import pkg_resources

def load_example_data(omicstype="metabolomics"):
    """
    Loads example datasets

    :param omicstype: type of omics for example data. Available options are "metabolomics" or "transcriptomics"
    Metabolomics data are from Su et al 2020 https://doi.org/10.1016/j.cell.2020.10.037

    :return: pre-processed omics data matrix consisting of m samples and n entities (metabolites/genes)
    in the form of a pandas DataFrame. Contains one of more metadata columns at the end.
    """

    if omicstype == "metabolomics":
        stream = pkg_resources.resource_stream(__name__, 'example_data/Su_covid_metabolomics_processed.csv')
        f = pd.read_csv(stream, index_col=0, encoding='latin-1')
        return f