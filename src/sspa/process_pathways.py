import pandas as pd
import pkg_resources
import download_KEGG


class ProcessPathways:
    """
    Process pathway database files for use with sspa.

    :param infile : For Reactome, release 78 is available for all organisms and for KEGG release 98 (homo sapiens) is available.
     Otherwise, specify the full path to the pathway input file. For reactome pathways, download
     the ChEBI2Reactome_All_levels file from https://reactome.org/download-data and specify the full filepath to this.
    :param organism : the organism requried when using Reactome. Full name must be specified as in Reacome db e.g. ('Homo sapiens' or 'Danio rerio')

    :returns pathway DataFrame in GMT format. Each row represents a pathway. First column is ID, second is full name, rest of columns represent entities using DB specific identifiers. 
    """

    def __init__(self, infile=None):
        self.infile = infile

    def process_reactome(self, organism, download_latest=False):
        # Process CHEBI to reactome data

        if self.infile == None or self.infile == "R78":
            stream = pkg_resources.resource_stream(__name__, 'pathway_databases/ChEBI2Reactome_All_Levels_R78.txt')
            f = pd.read_csv(stream, sep="\t", header=None, encoding='latin-1')
        else:
            f = pd.read_csv(self.infile, sep="\t", header=None)
        f.columns = ['CHEBI', 'pathway_ID', 'link', 'pathway_name', 'evidence_code', 'species']
        f_filt = f[f.species == organism]
        name_dict = dict(zip(f_filt['pathway_ID'], f_filt['pathway_name']))

        groups = f_filt.groupby(['pathway_ID'])['CHEBI'].apply(list).to_dict()
        df = pd.DataFrame.from_dict(groups, orient='index', dtype="object")
        pathways_df = df.dropna(axis=0, how='all', subset=df.columns.tolist()[1:])
        pathways_df = df.dropna(axis=1, how='all')

        # Remove duplicated compounds
        mask = pathways_df.apply(pd.Series.duplicated, 1) & pathways_df.astype(bool)
        pathways_df = pathways_df.mask(mask, None)
        pathways_df["Pathway_name"] = pathways_df.index.map(name_dict)
        pathways_df.insert(0, 'Pathway_name', pathways_df.pop('Pathway_name'))

        return pathways_df

    def process_kegg(self, organism, download_latest=False):

        if download_latest:
            pathways_df = download_KEGG(organism)
            return pathways_df

        else:
            if self.infile == None or self.infile == "R98":
                stream = pkg_resources.resource_stream(__name__, 'pathway_databases/KEGG_human_pathways_compounds_R98.csv')
                pathways_df = pd.read_csv(stream, index_col=0, encoding='latin-1')
            else:
                pathways_df = pd.read_csv(self.infile, index_col=0)

            pathways_df = pathways_df.dropna(axis=0, how='all', subset=pathways_df.columns.tolist()[1:])
            pathways_df = pathways_df.dropna(axis=1, how='all')

            # Remove duplicated compounds
            mask = pathways_df.apply(pd.Series.duplicated, 1) & pathways_df.astype(bool)
            pathways_df = pathways_df.mask(mask, None)

            return pathways_df

    def process_gmt(self):
        pathways_df = pd.read_csv(self.infile, index_col=0)
        pathways_df = pathways_df.dropna(axis=0, how='all', subset=pathways_df.columns.tolist()[1:])
        pathways_df = pathways_df.dropna(axis=1, how='all')

        return pathways_df
