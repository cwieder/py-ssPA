# Get latest KEGG metabolic pathways using rest API
# Similar to KEGGREST R in principle

import requests
import re
import pandas as pd

def download_KEGG(organism, filepath=None):
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
        pathid = re.search(r"path:(.*)", path[0]).group(1)
        pathway_dict[pathid] = name

    # get compounds for each pathway
    base_url = 'http://rest.kegg.jp/get/'

    pathway_ids = [*pathway_dict]
    pathway_names = list(pathway_dict.values())
    pathway_compound_mapping = dict()

    for i in pathway_ids:
        complist = []
        current_url = base_url + "pathway:" + i

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
