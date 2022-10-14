#


### process_reactome
[source](https://github.com/cwieder/py-ssPA\blob\master\src/sspa/process_pathways.py\#L5)
```python
.process_reactome(
   organism, infile = None, download_latest = False, filepath = None
)
```

---
Function to load Reactome pathways 

**Args**

* **organism** (str) : Reactome organism name
* **infile** (str) : default None, provide a Reactome pathway file to process into the GMT-style dataframe 
* **download_latest** (Bool) : Downloads the latest version of Reactome metabolic pathways
* **filepath** (str) : filepath to save pathway file to, default is None - save to variable


**Returns**

GMT-like pd.DataFrame containing Reactome pathways

----


### process_kegg
[source](https://github.com/cwieder/py-ssPA\blob\master\src/sspa/process_pathways.py\#L44)
```python
.process_kegg(
   organism, infile = None, download_latest = False, filepath = None
)
```

---
Function to load KEGG pathways 

**Args**

* **organism** (str) : KEGG organism code
* **infile** (str) : default None, provide a KEGG pathway file to process into the GMT-style dataframe 
* **download_latest** (Bool) : Downloads the latest version of KEGG metabolic pathways
* **filepath** (str) : filepath to save pathway file to, default is None - save to variable


**Returns**

GMT-like pd.DataFrame containing KEGG pathways

----


### process_gmt
[source](https://github.com/cwieder/py-ssPA\blob\master\src/sspa/process_pathways.py\#L75)
```python
.process_gmt(
   infile
)
```

---
Function to load pathways from a custom GMT-like file

**Args**

* **infile** (str) : default None, provide a GMT pathway file to process into the GMT-style dataframe, file ending can be .csv or .gmt


**Returns**

GMT-like pd.DataFrame containing pathways
