import pandas as pd

class ProcessPathways:
    def __init__(self, infile, organism):
        self.infile = infile
        self.organism = organism

    def process_reactome(self):
        # Process CHEBI to reactome data
        f = pd.read_csv("pathway_databases/" + self.infile, sep="\t", header=None)
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
        f = pd.read_csv("pathway_databases/" + self.infile, index_col=0)
        name_dict = dict(zip(f.index, f['Pathway_name']))
        pathway_dict = {k: list(set(f.loc[k, '0':].tolist())) for k in list(name_dict.keys())}
        pathway_dict = {k: [i for i in v if pd.notnull(i)] for k, v in pathway_dict.items()}
        pathway_dict = {k: v for k, v in pathway_dict.items() if len(v) > 2}

        return pathway_dict, name_dict
        # remove dupes