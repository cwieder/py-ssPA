# Get latest KEGG metabolic pathways using rest API
# Similar to KEGGREST R in principle

import requests
import re
import pandas as pd

def download_KEGG(organism, export=False, filepath=None):
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
        print(i)
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

    if export and filepath:
        df.to_csv(filepath + "KEGG_" + organism + "_pathways_compounds_" + str(version_no) + ".csv")
    print("Complete!")

    return df
