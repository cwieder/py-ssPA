# sspa
![sspa_logo](sspa_logo.png)


## Single sample pathway analysis tools for omics data

Install via pip

```
pip install sspa
```

Basic usage

```
import sspa
```

We will import the metabolite pathways from the Reactome database. We can specify any of the Reactome organism names.
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

### Supported single sample pathway analysis methods:
- kPCA
- ssClustPA and ssClustPA(proj)

```
# kPCA
sspa.sspa_kpca(mat, pathways_dictionary)

# ssClustPA
sspa.sspa_cluster(mat, pathways_dictionary)

# ssClustPA(porj)
sspa.sspa_cluster(mat, pathways_dictionary, projection=True)

```

## License
GNU GPL 3.0

## Citation
If you are using this tool in your work, please cite: Wieder et al 2021 (manuscript in preparation)