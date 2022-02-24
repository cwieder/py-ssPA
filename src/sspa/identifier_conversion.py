import requests
import pandas as pd

def identifier_conversion(input_type, compound_list):
    url = "http://api.xialab.ca/mapcompounds"

    cpd_list = ";".join(compound_list)+";"

    input_cpds = {"queryList": cpd_list, "inputType": input_type}
    headers = {
        'Content-Type': "application/json",
        'cache-control': "no-cache",
        }

    response = requests.request("POST", url, json=input_cpds, headers=headers)

    resp_dict = response.json()
    resp_df = pd.DataFrame(resp_dict)
    return resp_df

def map_identifiers(query_df, output_id_type, matrix):
    # output id type can be any of ['HMDB', 'PubChem', 'ChEBI', 'KEGG', 'METLIN','SMILES']
    cpd_mapping_dict = dict(zip(query_df["Query"].tolist(), query_df[output_id_type].tolist()))
    cpd_mapping_dict = {k: v for k, v in cpd_mapping_dict.items() if pd.notnull(v)}
    renamed_mat = matrix.drop([i for i in matrix.columns if i not in cpd_mapping_dict.keys()], axis=1)
    renamed_mat = renamed_mat.rename(cpd_mapping_dict, axis=1)
    return renamed_mat