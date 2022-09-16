#


## MetExplorePaths
[source](https://github.com/cwieder/py-ssPA/blob/master/src/sspa/download_pathways.py/#L118)
```python 
MetExplorePaths(
   model, id_type, filepath = None
)
```


---
Class for downloading metexplore metabolic models in the form of pathways with mapped identifiers


**Attributes**

* **model** (str) : identifier of genome scale metabolic model available on MetExplore
* **id_type** (str) : identifier type for the model pathways
* **filepath** (str) : filepath to save the pathway file to, default is None (save to variable)
* **nMappedID** (int) : Number of metabolites mapping to the selected identifier type
* **nMetab** (int) : Number of metabolites in the model
* **pathways** (pd.DataFrame) : GMT-like format pathway DataFrame



**Methods:**


### .download_metexplore
[source](https://github.com/cwieder/py-ssPA/blob/master/src/sspa/download_pathways.py/#L141)
```python
.download_metexplore()
```

---
Function to download MetExplore pathways

----


### download_KEGG
[source](https://github.com/cwieder/py-ssPA/blob/master/src/sspa/download_pathways.py/#L10)
```python
.download_KEGG(
   organism, filepath = None
)
```

---
Function for KEGG pathway download

**Args**

* **organism** (str) : KEGG 3 letter organism code
* **filepath** (str) : filepath to save pathway file to, default is None - save to variable


**Returns**

GMT-like pd.DataFrame containing KEGG pathways

----


### download_reactome
[source](https://github.com/cwieder/py-ssPA/blob/master/src/sspa/download_pathways.py/#L77)
```python
.download_reactome(
   organism, filepath = None
)
```

---
Function for Reactome pathway download

**Args**

* **organism** (str) : Reactome organism name
* **filepath** (str): filepath (str) : filepath to save pathway file to, default is None - save to variable


**Returns**

GMT-like pd.DataFrame containing Reactome pathways
