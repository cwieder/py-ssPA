import pandas as pd
from importlib import resources as resources
import sspa.download_pathways 

def process_reactome(organism, infile=None, download_latest=False, filepath=None, omics_type='metabolomics', identifiers=None, sep="\t"):
    '''
    Function to load Reactome pathways 
    Args:
        organism (str): Reactome organism name
        infile (str): default None, provide a Reactome pathway file to process into the GMT-style dataframe 
        download_latest (Bool): Downloads the latest version of Reactome metabolic pathways
        filepath (str): filepath to save pathway file to, default is None - save to variable
        omics_type(str): If using download_latest, specify type of omics pathways to download. Options are 'metabolomics', 'proteomics', 'transcriptomics', or 'multiomics'
        identifiers (list): list of identifiers to download for multi-omics pathways, default is None (download all). Options are 'chebi', 'uniprot', 'gene_symbol'
        sep (str): separator for the input file, default is '\t'
    Returns: 
        GMT-like pd.DataFrame containing Reactome pathways
    '''

    # Process CHEBI to reactome data

    if download_latest:
        pathways_df = sspa.download_pathways.download_reactome(organism, filepath, omics_type, identifiers)
        return pathways_df
    
    else:
        if omics_type != 'metabolomics':
            raise ValueError('Proteomics/multi-omics pathways only accessible when download_latest=True')
        if infile == None or infile == "R78":
            stream = resources.files(__name__).joinpath('pathway_databases/ChEBI2Reactome_All_Levels_R78.txt').open('rb')
            f = pd.read_csv(stream, sep=sep, header=None, encoding='latin-1')
            f.columns = ['CHEBI', 'pathway_ID', 'link', 'pathway_name', 'evidence_code', 'species']
            f_filt = f[f.species == organism]
            name_dict = dict(zip(f_filt['pathway_ID'], f_filt['pathway_name']))

            groups = f_filt.groupby(['pathway_ID'])['CHEBI'].apply(list).to_dict()
            groups = {k: list(set(v)) for k, v in groups.items()}
            df = pd.DataFrame.from_dict(groups, orient='index', dtype="object")
        else:
            df = pd.read_csv(infile, index_col=0,sep=sep, header=None)

        
        pathways_df = df.dropna(axis=0, how='all', subset=df.columns.tolist()[1:])
        pathways_df = df.dropna(axis=1, how='all')

        # Remove duplicated compounds
        mask = pathways_df.apply(pd.Series.duplicated, 1) & pathways_df.astype(bool)
        pathways_df = pathways_df.mask(mask, None)

        pathways_df = pathways_df.rename(columns={1: "Pathway_name"})

        return pathways_df

def process_kegg(organism, infile=None, download_latest=False, filepath=None, omics_type='metabolomics',sep="\t"):
    '''
    Function to load KEGG pathways 
    Args:
        organism (str): KEGG organism code
        infile (str): default None, provide a KEGG pathway file to process into the GMT-style dataframe 
        download_latest (Bool): Downloads the latest version of KEGG metabolic pathways
        filepath (str): filepath to save pathway file to, default is None - save to variable
        sep (str): separator for the input file, default is '\t'
    Returns: 
        GMT-like pd.DataFrame containing KEGG pathways
    '''
    if download_latest:
        pathways_df = sspa.download_pathways.download_KEGG(organism, filepath, omics_type)

        pathways_df = pathways_df.rename(columns={1: "Pathway_name"})
        return pathways_df

    else:
        if omics_type != 'metabolomics':
            raise ValueError('Proteomics/multi-omics pathways only accessible when download_latest=True')
        if infile == None or infile == "R98":
            stream = resources.files(__name__).joinpath('pathway_databases/KEGG_human_pathways_compounds_R98.csv').open('rb')
            pathways_df = pd.read_csv(stream, index_col=0, encoding='latin-1')
        else:
            pathways_df = pd.read_csv(infile, index_col=0, sep=sep, header=None)

        pathways_df = pathways_df.dropna(axis=0, how='all', subset=pathways_df.columns.tolist()[1:])
        pathways_df = pathways_df.dropna(axis=1, how='all')

        # Remove duplicated compounds
        mask = pathways_df.apply(pd.Series.duplicated, 1) & pathways_df.astype(bool)
        pathways_df = pathways_df.mask(mask, None)
        pathways_df = pathways_df.rename(columns={1: "Pathway_name"})

        return pathways_df
    

def process_pathbank(organism, infile=None, download_latest=False, filepath=None, omics_type='metabolomics', sep="\t"):
    '''
    Function to load PathBank pathways 
    Args:
        infile (str): default None, provide a PathBank pathway file to process into the GMT-style dataframe 
        download_latest (Bool): Downloads the latest version of PathBank metabolic pathways
        filepath (str): filepath to save pathway file to, default is None - save to variable
        omics_type(str): If using download_latest, specify type of omics pathways to download. Options are 'metabolomics', 'proteomics', or 'multiomics'
        sep (str): separator for the input file, default is '\t'
    Returns: 
        GMT-like pd.DataFrame containing PathBank pathways
    '''
    if download_latest:
        pathways_df = sspa.download_pathways.download_pathbank(organism, filepath, omics_type)
        return pathways_df

    else:
        if infile:
            pathways_df = pd.read_csv(infile, sep=sep,index_col=0)
        else: 
            print('Set download_latest=True to download latest version of PathBank pathways or provide a saved PathBank pathway .csv/.gmt file to load')


def process_gmt(infile, sep=","):
    '''
    Function to load pathways from a custom GMT-like file
    Args:
        infile (str): default None, provide a GMT pathway file to process into the GMT-style dataframe, file ending can be .csv or .gmt
        sep (str): separator for the input file, default is ','
    Returns: 
        GMT-like pd.DataFrame containing pathways
    '''
    if infile[-4:] == ".csv":
        pathways_df = pd.read_csv(infile, index_col=0, dtype='object', sep=sep)
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
    pathways_df = pathways_df.astype('object')
    return pathways_df

