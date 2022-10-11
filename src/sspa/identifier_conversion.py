import requests
import pandas as pd
import numpy as np

def identifier_conversion(input_type, compound_list):
    """
    Use Metaboanalyst API for identifier conversion
    Args:
        input_type (str): identifier type present in input data - any of ('name', 'hmdb', 'pubchem', 'chebi', 'metlin', 'kegg')
        compound_list (list): list of identifiers in the data
    Returns:
        (pd.DataFrame) Dataframe containing identifier matches 
    """
    url = "http://api.xialab.ca/mapcompounds"

    resps = []
    for i in range(0, len(compound_list), 100):
        cpds = compound_list[i:i+100]
        input_cpds = {"queryList": cpds, "inputType": input_type}
        headers = {
            # 'Content-Type': "application/json",
            'cache-control': "no-cache",
            }

        response = requests.request("POST", url, json=input_cpds, headers=headers)
        resp_dict = response.json()
        resps.append(resp_dict)

    res_dict = {k: [d[k] for d in resps] for k in resps[0]}
    res_dict = {k: [item for sublist in v for item in sublist] for k, v in res_dict.items()}

    resp_df = pd.DataFrame(res_dict)
    return resp_df

def map_identifiers(query_df, output_id_type, matrix):
    """
    Map desired identifiers to input data
    Args:
        query_df (pd.DataFrame): DataFrame obtained using the identifier_conversion function containing ID mappings
        output_id_type (str): Any of ('Match', 'HMDB', 'PubChem', 'ChEBI', 'KEGG', 'METLIN','SMILES')
        matrix (pd.DataFrame): sample-by-compound metabolomics data matrix
    Returns:
        Sample-by-compound metabolomics data matrix with mapped identifiers, any compounds without a matching ID will be dropped
    """


    # output id type can be any of ['HMDB', 'PubChem', 'ChEBI', 'KEGG', 'METLIN','SMILES']
    cpd_mapping_dict = dict(zip(query_df["Query"].tolist(), query_df[output_id_type].tolist()))
    cpd_mapping_dict = {k: v for k, v in cpd_mapping_dict.items() if v not in ["NA", np.nan, None, "None"]}
    renamed_mat = matrix.drop([i for i in matrix.columns if i not in cpd_mapping_dict.keys()], axis=1)
    renamed_mat = renamed_mat.rename(cpd_mapping_dict, axis=1)
    return renamed_mat