#


### identifier_conversion
[source](https://github.com/cwieder/py-ssPA\blob\master\src/sspa/identifier_conversion.py\#L5)
```python
.identifier_conversion(
   input_type, compound_list
)
```

---
Use Metaboanalyst API for identifier conversion

**Args**

* **input_type** (str) : identifier type present in input data - any of ('name', 'hmdb', 'pubchem', 'chebi', 'metlin', 'kegg')
* **compound_list** (list) : list of identifiers in the data


**Returns**

(pd.DataFrame) Dataframe containing identifier matches 

----


### map_identifiers
[source](https://github.com/cwieder/py-ssPA\blob\master\src/sspa/identifier_conversion.py\#L35)
```python
.map_identifiers(
   query_df, output_id_type, matrix
)
```

---
Map desired identifiers to input data

**Args**

* **query_df** (pd.DataFrame) : DataFrame obtained using the identifier_conversion function containing ID mappings
* **output_id_type** (str) : Any of ('Match', 'HMDB', 'PubChem', 'ChEBI', 'KEGG', 'METLIN','SMILES')
* **matrix** (pd.DataFrame) : sample-by-compound metabolomics data matrix


**Returns**

Sample-by-compound metabolomics data matrix with mapped identifiers, any compounds without a matching ID will be dropped
