# Get latest KEGG metabolic pathways using rest API
# Similar to KEGGREST R in principle

import requests
import re
import pandas as pd
import warnings
from tqdm import tqdm
import zipfile
import requests
import io

def download_KEGG(organism, filepath=None, omics_type='metabolomics'):
    '''
    Function for KEGG pathway download
    Args:
        organism (str): KEGG 3 letter organism code
        filepath (str): filepath to save pathway file to, default is None - save to variable
        omics_type(str): type of omics pathways to download (metabolomics or multiomics)
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

    # get release details
    release_data = requests.get('http://rest.kegg.jp/info/kegg')
    version_no = release_data.text.split()[9][0:3]

    if omics_type == 'metabolomics':
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

        # remove empty pathway entries
        pathway_compound_mapping = {k: v for k, v in pathway_compound_mapping.items() if v}

        # create GMT style file
        df = pd.DataFrame.from_dict(pathway_compound_mapping, orient='index')
        # map pathway names onto first column
        df.insert(0, 'Pathway_name', df.index.map(pathway_dict.get))

        if filepath:
            fpath = filepath + "/KEGG_" + organism + "_pathways_compounds_R" + str(version_no) + ".gmt"
            df.to_csv(fpath, sep="\t", header=False)
            print("KEGG DB file saved to " + fpath)
        print("Complete!")

        return df
        

    if omics_type == 'multiomics':
        pathway_mapping = dict()

        for index,i in enumerate(tqdm(pathway_ids)):
            complist = []
            genelist = []
            current_url = base_url + "pathway:" +i
            # parse the pathway description page
            page = requests.get(current_url)
            lines = page.text.split("\n")

            try:
                genes_start = [lines.index(i) for i in lines if i.startswith("GENE")][0]
                cpds_start = [lines.index(i) for i in lines if i.startswith("COMPOUND")][0]
                reference_start = [lines.index(i) for i in lines if i.startswith("REFERENCE") or i.startswith("REL_PATHWAY")][0]
                genes_lines = lines[genes_start:cpds_start]
                cpds_lines = lines[cpds_start:reference_start]

                first_cpd = cpds_lines.pop(0).split()[1]
                complist.append(first_cpd)
                complist = complist + [i.split()[0] for i in cpds_lines]
                first_gene = genes_lines.pop(0).split()[1]
                genelist.append(first_gene)
                genelist = genelist + [i.split()[0] for i in genes_lines]
                pathway_mapping[i] = list(set(complist)) + list(set(genelist))
            except IndexError:
                pathway_mapping[i] = []

        # remove empty pathway entries
        pathway_mapping = {k: v for k, v in pathway_mapping.items() if v}

        # create GMT style file
        df = pd.DataFrame.from_dict(pathway_mapping, orient='index')
        # map pathway names onto first column
        df.insert(0, 'Pathway_name', df.index.map(pathway_dict.get))

        if filepath:
            fpath = filepath + "/KEGG_" + organism + "_pathways_multiomics_R" + str(version_no) + ".gmt"
            df.to_csv(fpath, sep="\t", header=False)
            print("KEGG DB file saved to " + fpath)
        print("Complete!")

        return df

def download_reactome(organism, filepath=None, omics_type='metabolomics', identifiers=None):
    '''
    Function for Reactome pathway download
    Args:
        organism (str): Reactome organism name
        filepath (str): filepath (str): filepath to save pathway file to, default is None - save to variable
        omics_type(str): type of omics pathways to download. 
        Options are 'metabolomics' (ChEBI identifiers), 'proteomics' (UniProt identifiers), 'transcriptomics' (Gene Symbol), or 'multiomics' (ChEBI, UniProt and Gene Symbol identifiers)
        identifiers (list): list of identifiers to download for multi-omics pathways, default is None (download all). Options are 'chebi', 'uniprot', 'gene_symbol'
    Returns: 
        GMT-like pd.DataFrame containing Reactome pathways
    '''
    print("Beginning Reactome download...")

    # get release details
    release_data = requests.get('https://reactome.org/download/current/reactome_stable_ids.txt')
    version_no = release_data.text.split()[6]

    # get all pathways
    if omics_type == 'metabolomics':
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

        if filepath:
            fpath = filepath + "/Reactome_" + "_".join(organism.split())+ "_pathways_ChEBI_R" + str(version_no) + ".gmt"
            pathways_df.to_csv(fpath, sep="\t", header=False)
            print("Reactome DB file saved to " + fpath)
        
        print("Complete!")
        return pathways_df

    if omics_type == 'proteomics':
        url = 'https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
        f = pd.read_csv(url, sep="\t", header=None)
        
        f.columns = ['UniProt', 'pathway_ID', 'link', 'pathway_name', 'evidence_code', 'species']
        f_filt = f[f.species == organism]
        name_dict = dict(zip(f_filt['pathway_ID'], f_filt['pathway_name']))
        groups = f_filt.groupby(['pathway_ID'])['UniProt'].apply(list).to_dict()
        groups = {k: list(set(v)) for k, v in groups.items()}
        df = pd.DataFrame.from_dict(groups, orient='index', dtype="object")
        pathways_df = df.dropna(axis=0, how='all', subset=df.columns.tolist()[1:])
        pathways_df = df.dropna(axis=1, how='all')
        pathways_df["Pathway_name"] = pathways_df.index.map(name_dict)
        pathways_df.insert(0, 'Pathway_name', pathways_df.pop('Pathway_name'))

        if filepath:
            fpath = filepath + "/Reactome_" + "_".join(organism.split())+ "_pathways_UniProt_R" + str(version_no) + ".gmt"
            pathways_df.to_csv(fpath, sep="\t", header=False)
            print("Reactome DB file saved to " + fpath)
        
        print("Complete!")
        return pathways_df
    
    if omics_type == 'transcriptomics':
        url_gmt = 'https://reactome.org/download/current/ReactomePathways.gmt.zip'
        resp_gmt = requests.get(url_gmt)

        input_gmt = []
        with zipfile.ZipFile(io.BytesIO(resp_gmt.content)) as zip_ref:
            gmt_filename = [f for f in zip_ref.namelist() if f.endswith('.gmt')][0]

            # Extract the GMT file contents
            with zip_ref.open(gmt_filename) as gmt_file:
                for i in gmt_file:
                    i = i.decode('utf-8').strip()
                    input_gmt.append(i.strip("\n").split("\t"))
        pathways_df = pd.DataFrame(input_gmt)
        pathways_df = pathways_df.rename({0:"Pathway_name", 1:"Pathway_ID"}, axis=1)
        pathways_df.index = pathways_df["Pathway_ID"]
        pathways_df = pathways_df.drop("Pathway_ID", axis=1)
        pathways_df = pathways_df.dropna(axis=0, how='all', subset=pathways_df.columns.tolist()[1:])
        pathways_df = pathways_df.dropna(axis=1, how='all')

        if filepath:
            fpath = filepath + "/Reactome_" + "_".join(organism.split())+ "_pathways_GeneSymbol_R" + str(version_no) + ".gmt"
            pathways_df.to_csv(fpath, sep="\t", header=False)
            print("Reactome DB file saved to " + fpath)
        
        print("Complete!")
        return pathways_df

    if omics_type == 'multiomics':
        
        if organism != 'Homo sapiens' and ('gene_symbol' in identifiers):
            print('WARNING: Reactome does not provide gene_symbols for a non-Human organism, use UniProt instead')
            
        name_df = pd.read_csv('https://reactome.org/download/current/ReactomePathways.txt', sep="\t", header=None)
        url_prot = 'https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
        url_metab = 'https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt'
        url_gmt = 'https://reactome.org/download/current/ReactomePathways.gmt.zip'
        
        # read the chebi and uniprot files
        pathway_dicts = {'chebi': None, 'uniprot': None, 'gene_symbol': None}
        for k, url in {'uniprot': url_prot, 'chebi': url_metab}.items():
            f = pd.read_csv(url, sep="\t", header=None)
            f.columns = ['molecule_ID', 'pathway_ID', 'link', 'pathway_name', 'evidence_code', 'species']
            f_filt = f[f.species == organism]
            groups = f_filt.groupby(['pathway_ID'])['molecule_ID'].apply(list).to_dict()
            groups = {k: list(set(v)) for k, v in groups.items()}
            pathway_dicts[k] = groups

        # read the GMT with gene symbols
        url_gmt = 'https://reactome.org/download/current/ReactomePathways.gmt.zip'
        resp_gmt = requests.get(url_gmt)

        input_gmt = []
        gene_dict = {}
        with zipfile.ZipFile(io.BytesIO(resp_gmt.content)) as zip_ref:
            gmt_filename = [f for f in zip_ref.namelist() if f.endswith('.gmt')][0]

            # Extract the GMT file contents
            with zip_ref.open(gmt_filename) as gmt_file:
                for i in gmt_file:
                    i = i.decode('utf-8').strip()
                    line = i.strip("\n").split("\t")
                    gene_dict[line[1]] = line[2:]
        pathway_dicts['gene_symbol'] = gene_dict

        # combine all dicts by key
        # combine the required dictionaries based on the identifiers
        if identifiers:
            required_dicts = [pathway_dicts[k] for k in identifiers]
        else:
            required_dicts = list(pathway_dicts.values())

        combined_dict = {}
        for k in set().union(*required_dicts):
            combined_dict[k] = list(d.get(k) for d in required_dicts)
        # remove none values
        combined_dict = {k: [x for x in v if x is not None] for k, v in combined_dict.items()}
        # flatten the list of lists
        combined_dict = {k: [item for sublist in v for item in sublist] for k, v in combined_dict.items()}

        # create the dataframe
        reactome_mo = pd.DataFrame.from_dict(combined_dict, orient='index', dtype="object")

        # add the pathway names
        reactome_mo["Pathway_name"] = reactome_mo.index.map(dict(zip(name_df[0], name_df[1])))
        reactome_mo.insert(0, 'Pathway_name', reactome_mo.pop('Pathway_name'))

        if filepath:
            fpath = filepath + "/Reactome_" + "_".join(organism.split())+ "_pathways_multiomics_R" + str(version_no) + ".gmt"
            reactome_mo.to_csv(fpath, sep="\t", header=False)
            print("Reactome DB file saved to " + fpath)

        print("Complete!")
        return reactome_mo


# class MetExplorePaths:
#     '''
#     Class for downloading metexplore metabolic models in the form of pathways with mapped identifiers (DEPRECATED)

#     Attributes:
#         model (str): identifier of genome scale metabolic model available on MetExplore
#         id_type (str): identifier type for the model pathways
#         filepath (str): filepath to save the pathway file to, default is None (save to variable)
#         nMappedID (int): Number of metabolites mapping to the selected identifier type
#         nMetab (int): Number of metabolites in the model
#         pathways (pd.DataFrame): GMT-like format pathway DataFrame

#     '''
#     def __init__(self, model, id_type, filepath=None):
#         self.model = model
#         self.id_type = id_type
#         self.filepath = filepath
#         self.nMappedID = None
#         self.nMetab = None
#         self.pathways = None
#         # downloads pathways on object instantiation
#         self.download_metexplore()

#     def download_metexplore(self):
#         '''
#         Function to download MetExplore pathways
#         '''
#         warnings.filterwarnings("ignore")
#         metexploreURL = "https://metexplore.toulouse.inrae.fr/metexplore-api/"+str(self.model)+"/pathwaymetabolite/"+str(self.id_type)+"/"
#         stats_nmapped_url = "https://metexplore.toulouse.inrae.fr/metexplore-api/stat/"+str(self.model)+"/"+str(self.id_type)+"/"
#         stats_nmetab_url = "https://metexplore.toulouse.inrae.fr/metexplore-api/stat/"+str(self.model)+"/nbMetab/"
        
#         stats_nmetab = requests.get(stats_nmetab_url, verify=False)
#         stats_nmapped = requests.get(stats_nmapped_url, verify=False)
#         data_api = requests.get(metexploreURL, verify=False)
#         pathways_json = data_api.json()

#         pathways_df = pd.DataFrame.from_dict(pathways_json, orient='columns')
#         del(pathways_df["id"]) 
        
#         pathways = pathways_df.merge(pathways_df.Metabolites.apply(pd.Series), right_index=True, left_index=True)
#         del(pathways["Metabolites"])
        
#         pathways = pathways.fillna("None")
#         pathways.rename(columns={'name':"Pathway_name"}, inplace = True)
        
#         pathways = pathways.fillna("None")
#         pathways = pathways.set_index('dbIdentifier')
#         pathways = pathways.drop(columns=['len'])
#         pathways.rename(columns={'name':"Pathway_name"}, inplace = True)

#         #if file path provided save gmt to drive
#         if self.filepath:
#             fpath = self.filepath + "/MetExplorePathways_" + str(self.model) + "_" + str(self.id_type) + ".gmt"
#             pathways.to_csv(fpath, sep="\t", header=False)
#             print("MetExplore metabolic network pathways file saved to " + fpath)

#         self.pathways = pathways
#         self.nMappedID = stats_nmapped.text.split("\n")[2]
#         self.nMetab = stats_nmetab.text.split("\n")[2]
        
#         print("Complete!")
#         return pathways


def download_pathbank(organism, filepath=None, omicstype='metabolomics'):
    '''
    Function for PathBank pathway download
    Args:
        organism (str): PathBank organism name
        filepath (str): filepath to save pathway file to, default is None - save to variable
        omics_type(str): type of omics pathways to download. 
        Options are 'metabolomics' (ChEBI identifiers), 'proteomics' (UniProt identifiers), or 'multiomics' (ChEBI and UniProt identifiers)
    '''
    organisms = ['Homo sapiens', 'Escherichia coli', 'Mus musculus', 'Arabidopsis thaliana',
    'Saccharomyces cerevisiae', 'Bos taurus', 'Caenorhabditis elegans',
    'Rattus norvegicus', 'Drosophila melanogaster', 'Pseudomonas aeruginosa']
    if organism not in organisms:
        raise ValueError('Organism must be one of '+ ", ".join(organisms))

    version_no = None
    pathway_names = pd.read_csv('https://pathbank.org/downloads/pathbank_all_pathways.csv.zip', compression='zip', sep=',', header=0)
    name_dict = dict(zip(pathway_names['SMPDB ID'], pathway_names['Name']))


    if omicstype == 'metabolomics':
        metabolites_url = 'https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip'
        chebi_pathways = pd.read_csv(metabolites_url, compression='zip', sep=',', header=0, dtype=str)
        chebi_pathways = chebi_pathways[chebi_pathways['Species'] == organism]

        # reformat to gmt style such that each row contains chebi ID per pathway and each column is a chebi id
        chebi_pathways = chebi_pathways.groupby(['PathBank ID', 'Pathway Name'])['ChEBI ID'].apply(list).reset_index()
        chebi_pathways_gmt = pd.DataFrame(chebi_pathways['ChEBI ID'].values.tolist(), index=chebi_pathways['PathBank ID'])
        chebi_pathways_gmt['Pathway_name'] = chebi_pathways_gmt.index.map(name_dict)
        chebi_pathways_gmt.insert(0, 'Pathway_name', chebi_pathways_gmt.pop('Pathway_name'))

        if filepath:
            fpath = filepath + "/Pathbank_" + "_".join(organism.split())+ "_pathways_ChEBI" + ".gmt"
            chebi_pathways_gmt.to_csv(fpath, sep="\t", header=False)
            print("Pathbank DB file saved to " + fpath)

        print("Complete!")
        return chebi_pathways_gmt
    
    if omicstype == 'proteomics':
        proteins_url = 'https://pathbank.org/downloads/pathbank_all_proteins.csv.zip'
        uniprot_pathways = pd.read_csv(proteins_url, compression='zip', sep=',', header=0, dtype=str)
        uniprot_pathways = uniprot_pathways[uniprot_pathways['Species'] == organism]

        # reformat to gmt style such that each row contains uniprot ID per pathway and each column is a uniprot id
        uniprot_pathways = uniprot_pathways.groupby(['PathBank ID', 'Pathway Name'])['Uniprot ID'].apply(list).reset_index()
        uniprot_pathways_gmt = pd.DataFrame(uniprot_pathways['Uniprot ID'].values.tolist(), index=uniprot_pathways['PathBank ID'])
        uniprot_pathways_gmt['Pathway_name'] = uniprot_pathways_gmt.index.map(name_dict)
        uniprot_pathways_gmt.insert(0, 'Pathway_name', uniprot_pathways_gmt.pop('Pathway_name'))

        if filepath:
            fpath = filepath + "/Pathbank_" + "_".join(organism.split())+ "_pathways_UniProt" + ".gmt"
            uniprot_pathways_gmt.to_csv(fpath, sep="\t", header=False)
            print("Pathbank DB file saved to " + fpath)

        print("Complete!")
        return uniprot_pathways_gmt
    
    if omicstype == 'multiomics':
        metabolites_url = 'https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip'
        chebi_pathways = pd.read_csv(metabolites_url, compression='zip', sep=',', header=0, dtype=str)
        chebi_pathways = chebi_pathways[chebi_pathways['Species'] == organism]

        # reformat to gmt style such that each row contains chebi ID per pathway and each column is a chebi id
        chebi_pathways = chebi_pathways.groupby(['PathBank ID', 'Pathway Name'])['ChEBI ID'].apply(list).reset_index()
        chebi_pathways_gmt = pd.DataFrame(chebi_pathways['ChEBI ID'].values.tolist(), index=chebi_pathways['PathBank ID'])

        proteins_url = 'https://pathbank.org/downloads/pathbank_all_proteins.csv.zip'
        uniprot_pathways = pd.read_csv(proteins_url, compression='zip', sep=',', header=0, dtype=str)
        uniprot_pathways = uniprot_pathways[uniprot_pathways['Species'] == organism]

        # reformat to gmt style such that each row contains uniprot ID per pathway and each column is a uniprot id
        uniprot_pathways = uniprot_pathways.groupby(['PathBank ID', 'Pathway Name'])['Uniprot ID'].apply(list).reset_index()
        uniprot_pathways_gmt = pd.DataFrame(uniprot_pathways['Uniprot ID'].values.tolist(), index=uniprot_pathways['PathBank ID'])

        multiomics_pathways_gmt = pd.concat([chebi_pathways_gmt, uniprot_pathways_gmt], axis=1)
        multiomics_pathways_gmt['Pathway_name'] = multiomics_pathways_gmt.index.map(name_dict)
        multiomics_pathways_gmt.insert(0, 'Pathway_name', multiomics_pathways_gmt.pop('Pathway_name'))

        if filepath:
            fpath = filepath + "/Pathbank_" + "_".join(organism.split())+ "_pathways_multiomics" + ".gmt"
            multiomics_pathways_gmt.to_csv(fpath, sep="\t", header=False)
            print("Pathbank DB file saved to " + fpath)

        print("Complete!")
  
        return multiomics_pathways_gmt
    