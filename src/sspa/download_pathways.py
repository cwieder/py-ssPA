# Get latest KEGG metabolic pathways using rest API
# Similar to KEGGREST R in principle

import requests
import re
import pandas as pd
import warnings
import json
from tqdm import tqdm

def download_KEGG(organism, filepath=None):
    '''
    Function for KEGG pathway download
    Args:
        organism (str): KEGG 3 letter organism code
        filepath (str): filepath to save pathway file to, default is None - save to variable
    Returns: 
        GMT-like pd.DataFrame containing KEGG pathways
    '''
    print("Beginning KEGG download...")
    # get all pathways
    url = 'http://rest.kegg.jp/list/pathway/'+organism
    # change organism name
    data = requests.get(url)
    pathways = data.text
    pathways = pathways.split("\n")
    pathways = filter(None, pathways)
    pathway_dict = dict()
    
    for path in pathways:
        path = path.split("\t")
        name = path[1]
        pathid = re.search(r"(.*)", path[0]).group(1)
        pathway_dict[pathid] = name

    # get compounds for each pathway
    base_url = 'http://rest.kegg.jp/get/'

    pathway_ids = [*pathway_dict]
    pathway_names = list(pathway_dict.values())
    pathway_compound_mapping = dict()

    for index,i in enumerate(tqdm(pathway_ids)):
        complist = []
        current_url = base_url + "pathway:" +i
        # parse the pathway description page
        page = requests.get(current_url)
        lines = page.text.split("\n")

        try:
            cpds_start = [lines.index(i) for i in lines if i.startswith("COMPOUND")][0]
            reference_start = [lines.index(i) for i in lines if i.startswith("REFERENCE") or i.startswith("REL_PATHWAY")][0]
            cpds_lines = lines[cpds_start:reference_start]
            first_cpd = cpds_lines.pop(0).split()[1]
            complist.append(first_cpd)
            complist = complist + [i.split()[0] for i in cpds_lines]
            pathway_compound_mapping[i] = list(set(complist))
        except IndexError:
            pathway_compound_mapping[i] = []

    # get release details
    release_data = requests.get('http://rest.kegg.jp/info/kegg')
    version_no = release_data.text.split()[9][0:3]

    # create GMT style file
    df = pd.DataFrame.from_dict(pathway_compound_mapping, orient='index')
    df.insert(0, 'Pathway_name', pathway_names)

    if filepath:
        fpath = filepath + "/KEGG_" + organism + "_pathways_compounds_R" + str(version_no) + ".gmt"
        df.to_csv(fpath, sep="\t", header=False)
        print("KEGG DB file saved to " + fpath)
    print("Complete!")

    return df

def download_reactome(organism, filepath=None):
    '''
    Function for Reactome pathway download
    Args:
        organism (str): Reactome organism name
        filepath (str): filepath (str): filepath to save pathway file to, default is None - save to variable
    Returns: 
        GMT-like pd.DataFrame containing Reactome pathways
    '''
    print("Beginning Reactome download...")

     # get all pathways
    url = 'https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt'
    f = pd.read_csv(url, sep="\t", header=None)
    
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

    # get release details
    release_data = requests.get('https://reactome.org/download/current/reactome_stable_ids.txt')
    version_no = release_data.text.split()[6]

    if filepath:
        fpath = filepath + "/Reactome_" + "_".join(organism.split())+ "_pathways_compounds_R" + str(version_no) + ".gmt"
        pathways_df.to_csv(fpath, sep="\t", header=False)
        print("Reactome DB file saved to " + fpath)
    print("Complete!")

    return pathways_df

class MetExplorePaths:
    '''
    Class for downloading metexplore metabolic models in the form of pathways with mapped identifiers

    Attributes:
        model (str): identifier of genome scale metabolic model available on MetExplore
        id_type (str): identifier type for the model pathways
        filepath (str): filepath to save the pathway file to, default is None (save to variable)
        nMappedID (int): Number of metabolites mapping to the selected identifier type
        nMetab (int): Number of metabolites in the model
        pathways (pd.DataFrame): GMT-like format pathway DataFrame

    '''
    def __init__(self, model, id_type, filepath=None):
        self.model = model
        self.id_type = id_type
        self.filepath = filepath
        self.nMappedID = None
        self.nMetab = None
        self.pathways = None
        # downloads pathways on object instantiation
        self.download_metexplore()

    def download_metexplore(self):
        '''
        Function to download MetExplore pathways
        '''
        warnings.filterwarnings("ignore")
        metexploreURL = "https://metexplore.toulouse.inrae.fr/metexplore-api/"+str(self.model)+"/pathwaymetabolite/"+str(self.id_type)+"/"
        stats_nmapped_url = "https://metexplore.toulouse.inrae.fr/metexplore-api/stat/"+str(self.model)+"/"+str(self.id_type)+"/"
        stats_nmetab_url = "https://metexplore.toulouse.inrae.fr/metexplore-api/stat/"+str(self.model)+"/nbMetab/"
        
        stats_nmetab = requests.get(stats_nmetab_url, verify=False)
        stats_nmapped = requests.get(stats_nmapped_url, verify=False)
        data_api = requests.get(metexploreURL, verify=False)
        pathways_json = data_api.json()

        pathways_df = pd.DataFrame.from_dict(pathways_json, orient='columns')
        del(pathways_df["id"]) 
        
        pathways = pathways_df.merge(pathways_df.Metabolites.apply(pd.Series), right_index=True, left_index=True)
        del(pathways["Metabolites"])
        
        pathways = pathways.fillna("None")
        pathways.rename(columns={'name':"Pathway_name"}, inplace = True)
        
        pathways = pathways.fillna("None")
        pathways = pathways.set_index('dbIdentifier')
        pathways = pathways.drop(columns=['len'])
        pathways.rename(columns={'name':"Pathway_name"}, inplace = True)

        #if file path provided save gmt to drive
        if self.filepath:
            fpath = self.filepath + "/MetExplorePathways_" + str(self.model) + "_" + str(self.id_type) + ".gmt"
            pathways.to_csv(fpath, sep="\t", header=False)
            print("MetExplore metabolic network pathways file saved to " + fpath)

        self.pathways = pathways
        self.nMappedID = stats_nmapped.text.split("\n")[2]
        self.nMetab = stats_nmetab.text.split("\n")[2]
        
        print("Complete!")
        return pathways
