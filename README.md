# sspa
![sspa_logo](sspa_logo.png)

[![DOI](https://zenodo.org/badge/442446643.svg)](https://zenodo.org/badge/latestdoi/442446643)
[![](https://readthedocs.org/projects/pip/badge/?version=latest&style=flat)](https://readthedocs.org/projects/py-ssPA/badge/)

## Single sample pathway analysis tools for omics data

Full walkthrough notebook available on Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1X_ZBjXYxBVv43LcArRv2aB2HJkUp8A4x?usp=sharing)

Documentation is available on our [Read the Docs page](https://cwieder.github.io/py-ssPA/)

## Quickstart
```
pip install sspa
```
Load Reactome pathways
```
reactome_pathways  = sspa.process_reactome(organism="Homo sapiens")
```

Load some example metabolomics data in the form of a pandas DataFrame:

```
covid_data_processed = sspa.load_example_data(omicstype="metabolomics", processed=True)
```

Generate pathway scores using kPCA method

```
kpca_scores = sspa.sspa_kpca(covid_data_processed, reactome_pathways)
```

## Loading pathways 
```
# Pre-loaded pathways
# Reactome v78
reactome_pathways  = sspa.process_reactome(organism="Homo sapiens")

# KEGG v98
kegg_human_pathways  = sspa.process_kegg(organism="hsa")
```

Load a custom GMT file (extension .gmt or .csv)
```
custom_pathways = sspa.process_gmt("wikipathways-20220310-gmt-Homo_sapiens.gmt")
```

Download latest version of pathways
```
# download KEGG latest
kegg_mouse_latest = sspa.process_kegg("mmu", download_latest=True, filepath=".")

# download Reactome latest
reactome_mouse_latest = sspa.process_reactome("Mus musculus", download_latest=True, filepath=".")
```

## Identifier harmonization 
```
# download the conversion table
compound_names = processed_data.columns.tolist()
conversion_table = sspa.identifier_conversion(input_type="name", compound_list=compound_names)

# map the identifiers to your dataset
processed_data_mapped = sspa.map_identifiers(conversion_table, output_id_type="ChEBI", matrix=processed_data)
```

## Conventional pathway analysis
ORA
```
ora = sspa.sspa_ora(processed_data_mapped, covid_data["Group"], reactome_pathways, 0.05, custom_background=None)

# perform ORA 
ora_res = ora.over_representation_analysis()

# get t-test results
ora.ttest_res

# obtain list of differential molecules input to ORA
ora.DA_molecules
```

GSEA
```
sspa.sspa_fgsea(processed_data_mapped, covid_data['Group'], reactome_pathways)
```

## Single sample pathway analysis methods
```
# ssclustPA
ssclustpa_proj_res = sspa.sspa_cluster(processed_data_mapped, reactome_pathways)

# kPCA
kpca_scores = sspa.sspa_kpca(processed_data_mapped, reactome_pathways)

# z-score
zscore_res = sspa.sspa_zscore(processed_data_mapped, reactome_pathways)

# SVD (PLAGE)
svd_res = sspa.sspa_svd(processed_data_mapped, reactome_pathways)

# GSVA
gsva_res = sspa.sspa_gsva(processed_data_mapped, reactome_pathways)
```


## License
GNU GPL 3.0

## Citation
[![DOI](https://zenodo.org/badge/442446643.svg)](https://zenodo.org/badge/latestdoi/442446643)

If you are using this tool in your work, please cite: <div class="csl-entry">Wieder, C., Lai, R. P. J., &#38; Ebbels, T. (2022). Single sample pathway analysis in metabolomics : performance evaluation and application. <i>BioRxiv</i>, 2022.04.11.487976. https://doi.org/10.1101/2022.04.11.487976</div>

If you are using the methods SVD (PLAGE), z-score, or ORA, please cite the original publications alongside this tool.