import pandas as pd
import pkg_resources

def load_example_data(omicstype="metabolomics"):
    if omicstype == "metabolomics":
        stream = pkg_resources.resource_stream(__name__, 'example_data/Su_covid_metabolomics_processed.csv')
        f = pd.read_csv(stream, index_col=0, encoding='latin-1')
        return f