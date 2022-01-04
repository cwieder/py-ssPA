import pandas as pd
import pkg_resources


class ProcessPathways:
    """
    Process pathway database files for use with sspa.

    :param infile : For Reactome, release 78 is available for all organisms and for KEGG release 98 (homo sapiens) is available.
     Otherwise, specify the full path to the pathway input file. For reactome pathways, download
     the ChEBI2Reactome_All_levels file from https://reactome.org/download-data and specify the full filepath to this.
    :param organism : the organism requried when using Reactome. Full name must be specified as in Reacome db e.g. ('Homo sapiens' or 'Danio rerio')

    :returns pathway entity dictionary (composed of pathway IDs and the entities each pathway is composed of),
     and pathway name mapping dict (composed of pathway ID: pathway name)
    """

    def __init__(self, infile, organism):
        self.infile = infile
        self.organism = organism

    def process_reactome(self):
        # Process CHEBI to reactome data

        if self.infile == "R78":
            stream = pkg_resources.resource_stream(__name__, 'pathway_databases/ChEBI2Reactome_All_Levels_R78.txt')
            f = pd.read_csv(stream, sep="\t", header=None, encoding='latin-1')
        else:
            f = pd.read_csv(self.infile, sep="\t", header=None)
        f.columns = ['CHEBI', 'pathway_ID', 'link', 'pathway_name', 'evidence_code', 'species']
        f_filt = f[f.species == self.organism]
        name_dict = dict(zip(f_filt['pathway_ID'], f_filt['pathway_name']))

        groups = f_filt.groupby(['pathway_ID'])['CHEBI'].apply(list).to_dict()
        df = pd.DataFrame.from_dict(groups, orient='index', dtype="object")

        pathways_df = df.dropna(axis=0, how='all', subset=df.columns.tolist()[1:])
        pathways = pathways_df.index.tolist()
        pathway_dict = {}

        for pathway in pathways:
            pathway_compounds = list(set(pathways_df.loc[pathway, :].tolist()))
            pathway_compounds = [str(i) for i in pathway_compounds if str(i) != "None"]

            cpds = pathway_compounds[1:]
            if len(cpds) > 1:
                pathway_dict[pathway] = cpds
        return pathway_dict, name_dict

    def process_kegg(self):
        if self.infile == "R98":
            stream = pkg_resources.resource_stream(__name__, 'pathway_databases/KEGG_human_pathways_compounds_R98.csv')
            f = pd.read_csv(stream, index_col=0, encoding='latin-1')
        else:
            f = pd.read_csv(self.infile, index_col=0)
        name_dict = dict(zip(f.index, f['Pathway_name']))
        pathway_dict = {k: list(set(f.loc[k, '0':].tolist())) for k in list(name_dict.keys())}
        pathway_dict = {k: [i for i in v if pd.notnull(i)] for k, v in pathway_dict.items()}
        pathway_dict = {k: v for k, v in pathway_dict.items() if len(v) > 2}

        return pathway_dict, name_dict
        # remove dupes
