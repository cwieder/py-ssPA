import requests
import pandas as pd
import numpy as np

def identifier_conversion(input_type, compound_list):
    """
    Use Metaboanalyst API for identifier conversion
    Args:
    input_type: identifier type present in input data
    compound_list: list of identifiers in the data
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
    """Fetches rows from a Smalltable.

    Retrieves rows pertaining to the given keys from the Table instance
    represented by table_handle.  String keys will be UTF-8 encoded.

    Args:
        table_handle: An open smalltable.Table instance.
        keys: A sequence of strings representing the key of each table
          row to fetch.  String keys will be UTF-8 encoded.
        require_all_keys: If True only rows with values set for all keys will be
          returned.

    Returns:
        A dict mapping keys to the corresponding table row data
        fetched. Each row is represented as a tuple of strings. For
        example:

        {b'Serak': ('Rigel VII', 'Preparer'),
         b'Zim': ('Irk', 'Invader'),
         b'Lrrr': ('Omicron Persei 8', 'Emperor')}

        Returned keys are always bytes.  If a key from the keys argument is
        missing from the dictionary, then that row was not found in the
        table (and require_all_keys must have been False).

    Raises:
        IOError: An error occurred accessing the smalltable.
    """



    # output id type can be any of ['HMDB', 'PubChem', 'ChEBI', 'KEGG', 'METLIN','SMILES']
    cpd_mapping_dict = dict(zip(query_df["Query"].tolist(), query_df[output_id_type].tolist()))
    cpd_mapping_dict = {k: v for k, v in cpd_mapping_dict.items() if v not in ["NA", np.nan, None, "None"]}
    renamed_mat = matrix.drop([i for i in matrix.columns if i not in cpd_mapping_dict.keys()], axis=1)
    renamed_mat = renamed_mat.rename(cpd_mapping_dict, axis=1)
    return renamed_mat