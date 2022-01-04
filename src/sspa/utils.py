import pandas as pd
import pkg_resources

def load_example_data(omicstype="metabolomics"):
    """
    Loads example datasets
    :param omicstype: type of omics for example data. Available options are "metabolomics" or "transcriptomics"
    :return: pre-processed omics data matrix consisting of m samples and n entities (metabolites/genes)
    """

    if omicstype == "metabolomics":
        stream = pkg_resources.resource_stream(__name__, 'example_data/Su_covid_metabolomics_processed.csv')
        f = pd.read_csv(stream, index_col=0, encoding='latin-1')
        return f