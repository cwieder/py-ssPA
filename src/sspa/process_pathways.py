import pandas as pd
import pkg_resources
import sspa.download_pathways 

def process_reactome(organism, infile=None, download_latest=False, filepath=None):
    '''
    Function to load Reactome pathways 
    Args:
        organism (str): Reactome organism name
        infile (str): default None, provide a Reactome pathway file to process into the GMT-style dataframe 
        download_latest (Bool): Downloads the latest version of Reactome metabolic pathways
        filepath (str): filepath to save pathway file to, default is None - save to variable
    Returns: 
        GMT-like pd.DataFrame containing Reactome pathways
    '''

    # Process CHEBI to reactome data

    if download_latest:
        pathways_df = sspa.download_pathways.download_reactome(organism, filepath)
        return pathways_df

    else:
        if infile == None or infile == "R78":
            stream = pkg_resources.resource_stream(__name__, 'pathway_databases/ChEBI2Reactome_All_Levels_R78.txt')
            f = pd.read_csv(stream, sep="\t", header=None, encoding='latin-1')
        else:
            f = pd.read_csv(infile, sep="\t", header=None)
        f.columns = ['CHEBI', 'pathway_ID', 'link', 'pathway_name', 'evidence_code', 'species']
        f_filt = f[f.species == organism]
        name_dict = dict(zip(f_filt['pathway_ID'], f_filt['pathway_name']))

        groups = f_filt.groupby(['pathway_ID'])['CHEBI'].apply(list).to_dict()
        groups = {k: list(set(v)) for k, v in groups.items()}
        df = pd.DataFrame.from_dict(groups, orient='index', dtype="object")
        pathways_df = df.dropna(axis=0, how='all', subset=df.columns.tolist()[1:])
        pathways_df = df.dropna(axis=1, how='all')

        pathways_df["Pathway_name"] = pathways_df.index.map(name_dict)
        pathways_df.insert(0, 'Pathway_name', pathways_df.pop('Pathway_name'))

        return pathways_df

def process_kegg(organism, infile=None, download_latest=False, filepath=None):
    '''
    Function to load KEGG pathways 
    Args:
        organism (str): KEGG organism code
        infile (str): default None, provide a KEGG pathway file to process into the GMT-style dataframe 
        download_latest (Bool): Downloads the latest version of KEGG metabolic pathways
        filepath (str): filepath to save pathway file to, default is None - save to variable
    Returns: 
        GMT-like pd.DataFrame containing KEGG pathways
    '''
    if download_latest:
        pathways_df = sspa.download_pathways.download_KEGG(organism, filepath)
        return pathways_df

    else:
        if infile == None or infile == "R98":
            stream = pkg_resources.resource_stream(__name__, 'pathway_databases/KEGG_human_pathways_compounds_R98.csv')
            pathways_df = pd.read_csv(stream, index_col=0, encoding='latin-1')
        else:
            pathways_df = pd.read_csv(infile, index_col=0)

        pathways_df = pathways_df.dropna(axis=0, how='all', subset=pathways_df.columns.tolist()[1:])
        pathways_df = pathways_df.dropna(axis=1, how='all')

        # Remove duplicated compounds
        mask = pathways_df.apply(pd.Series.duplicated, 1) & pathways_df.astype(bool)
        pathways_df = pathways_df.mask(mask, None)

        return pathways_df

def process_gmt(infile):
    '''
    Function to load pathways from a custom GMT-like file
    Args:
        infile (str): default None, provide a GMT pathway file to process into the GMT-style dataframe, file ending can be .csv or .gmt
    Returns: 
        GMT-like pd.DataFrame containing pathways
    '''
    if infile[-4:] == ".csv":
        pathways_df = pd.read_csv(infile, index_col=0)
    elif infile[-4:] == ".gmt":
        input_gmt = []
        with open(infile, "r") as f:
            for i in f:
                input_gmt.append(i.strip("\n").split("\t"))
        pathways_df = pd.DataFrame(input_gmt)
        pathways_df = pathways_df.rename({0:"Pathway_ID", 1:"Pathway_name"}, axis=1)
        pathways_df.index = pathways_df["Pathway_ID"]
        pathways_df = pathways_df.drop("Pathway_ID", axis=1)

    pathways_df = pathways_df.dropna(axis=0, how='all', subset=pathways_df.columns.tolist()[1:])
    pathways_df = pathways_df.dropna(axis=1, how='all')
    
    return pathways_df

