# sspa
![sspa_logo](sspa_logo.png)


## Single sample pathway analysis tools for omics data
**This package is still in development stages and not yet available on pip.**

Install via pip

```
pip install sspa
```

### Input data

Omics matrix in the form of a Pandas DataFrame. DataFrame must contain rows representing samples and columns representing compound/gene identifiers. Identifiers must match those in the pathway database desired (i.e. ChEBI compounds for Reactome metabolite pathways, or KEGG compound IDs for KEGG metaboltie pathways.)

Data matrix must be scaled prior to use with sspa (each feature must have mean = 0 and SD = 1)


### Basic usage

```
import sspa
```

We will import and process the metabolite pathways from the Reactome database into a python dictionary. We can specify any of the Reactome organism names.
This returns a dictionary of pathways and compounds and a pathway name mapping dictionary

```
reactome_pathways, reactome_pathway_names = sspa.ProcessPathways("R78", "Homo sapiens").process_reactome()
```

Load some example metabolomics data in the form of a pandas DataFrame:

```
example_data = sspa.load_example_data(omicstype="metabolomics")
```

Generate pathway scores using kPCA method

```
kpca_scores = sspa.sspa_kpca(example_data.iloc[:, :-2], reactome_pathways)
```

### Available single sample pathway analysis methods:
- kPCA
- ssClustPA and ssClustPA(proj)
- z-score (Lee et al. 2008)
- SVD (PLAGE) (Tomfohr et al. 2007)

```
# kPCA
sspa.sspa_kpca(mat, pathways_dictionary)

# ssClustPA
sspa.sspa_cluster(mat, pathways_dictionary)

# ssClustPA(proj)
sspa.sspa_cluster(mat, pathways_dictionary, projection=True)

# z-score
sspa.sspa_zscore(mat, pathways_dictionary)

# SVD (PLAGE)
sspa.sspa_svd(mat, pathways_dictionary)
```

### Additional available conventional pathway analysis methods
Over representation analysis (ORA)

cutoff_thresh: threshold on the FDR adjusted P-value to select differential metabolites

```
# ORA
sspa.sspa_ora(mat, metadata_column, pathways_dictionary, cutoff_thresh=0.05):
```

## License
GNU GPL 3.0

## Citation
If you are using this tool in your work, please cite: Wieder et al 2022 (manuscript in preparation).
If you are using the methods SVD (PLAGE), z-score, or ORA, please cite the original publications alongside this tool.