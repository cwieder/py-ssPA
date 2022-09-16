#


### load_example_data
[source](https://github.com/cwieder/py-ssPA/blob/master/src/sspa/utils.py/#L7)
```python
.load_example_data(
   omicstype = 'metabolomics', processed = True
)
```

---
Loads example datasets


**Args**

* **omicstype** (str) : type of omics for example data. 
    Available options are "metabolomics" or "transcriptomics". 
    Metabolomics data are from Su et al 2020 https://doi.org/10.1016/j.cell.2020.10.037.
    Transcriptomics data - TO BE IMPLEMENTED
* **processed** (bool) : Load processed (normalised, scaled) or raw data


**Returns**

pre-processed omics data matrix consisting of m samples and n entities (metabolites/genes) in the form of a pandas DataFrame. 
Contains one of more metadata columns at the end.

----


### t_tests
[source](https://github.com/cwieder/py-ssPA/blob/master/src/sspa/utils.py/#L34)
```python
.t_tests(
   matrix, classes, multiple_correction_method, testtype = 'ttest'
)
```

---
Performs two-sample independent t-tests


**Args**

* **matrix** (pd.DataFrame) : processed sample-by-compound metabolomics dataframe
* **classes** (pd.Series) : pandas series containing phenotype metadata (e.g. 'COVID', 'NON-COVID')
* **multiple_correction_method** (str) : see https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for options
* **testtype** (str) : Default is t-test, "mwu" also available to implement the Mann Whitney U test


**Returns**

pd.DataFrame containing p-values and corrected p-values for each metabolite

----


### pathwaydf_to_dict
[source](https://github.com/cwieder/py-ssPA/blob/master/src/sspa/utils.py/#L66)
```python
.pathwaydf_to_dict(
   df
)
```

---
Converts pathway dataframe to dictionary, with pathway IDs as keys and metabolite lists as values

**Args**

* **df** (pd.DataFrame) : Pandas DataFrame containing pathways 


**Returns**

python dict pathway representation
