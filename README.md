# sspa
![sspa_logo](sspa_logo.png)

[![PyPI version](https://badge.fury.io/py/sspa.svg)](https://badge.fury.io/py/sspa)
[![DOI](https://zenodo.org/badge/442446643.svg)](https://zenodo.org/badge/latestdoi/442446643)
[![ssPA docs](https://github.com/cwieder/py-sspa/actions/workflows/sspa-docs.yml/badge.svg)](https://cwieder.github.io/py-ssPA/)
[![Downloads](https://pepy.tech/badge/sspa)](https://pepy.tech/project/sspa)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Single sample pathway analysis toolkit
sspa provides a Python interface for metabolomics pathway analysis. In addition to conventional methods over-representation analysis (ORA) and gene/metabolite set enrichment analysis (GSEA), it also provides a wide range of single-sample pathway analysis (ssPA) methods. 

## Features
- Over-representation analysis
- Metabolite set enrichment analysis (based on GSEA)
- Single-sample pathway analysis
- Compound identifier conversion
- Pathway database download (KEGG, Reactome, and MetExplore metabolic networks)

Although this package is designed to provide a user-friendly interface for metabolomics pathway analysis, the methods are also applicable to other datatypes such as normalised RNA-seq data. 

## Documentation and tutorials
This README provides a quickstart guide to the package and its functions. For new users we **highly recommend following our full walkthrough notebook tutorial** available on Google Colab which provides a step-by-step guide to using the package.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1PdJueJkpdkpplE2ieNoO8CKcEmtXF1-x?usp=sharing)

Click the link above and save a copy of the Colab notebook to your Google Drive. Alternatively, you can download the notebook from the Colab tutorial as an '.ipynb' file and run it locally using Jupyter Notebook or Jupyter Lab.

Documentation is available on our [Read the Docs page](https://cwieder.github.io/py-ssPA/). This includes a function API reference. 

## Quickstart
```
pip install sspa
```
Load Reactome pathways
```python
reactome_pathways  = sspa.process_reactome(organism="Homo sapiens")
```

Load some example metabolomics data in the form of a pandas DataFrame:

```python
covid_data_processed = sspa.load_example_data(omicstype="metabolomics", processed=True)
```

Generate pathway scores using kPCA method

```python
kpca_scores = sspa.sspa_kpca(reactome_pathways, min_entity=2).fit_transform(covid_data_processed.iloc[:, :-2])
```

## Loading example data
Note we provide processed and non-processed versins of the COVID example metabolomics dataset ([Su et al. 2020, Cell](https://data.mendeley.com/datasets/tzydswhhb5/5)). The processed version (set `processed=True`) already has ChEBI identifiers as column names, whereas the non-processed version has metabolite names. 

```python
covid_data = sspa.load_example_data(omicstype="metabolomics", processed=False)
```

Here we demonstrate some simple pre-processing for this dataset in order to enable conventional and ssPA pathway analysis:

```python
# Keep only metabolites (exclude metadata columns)
covid_values = covid_data.iloc[:, :-2]

# Remove metabolites with too many NA values
data_filt = covid_values.loc[:, covid_values.isin([' ', np.nan, 0]).mean() < 0.5]

# Impute using the median
imputed_mat = data_filt.fillna(data_filt.median())

# Log transform the data
log2_mat = np.log2(imputed_mat)

# Standardise the data (metabolite values) using z-score (mean=0, std=1) by subtracting the mean and dividing by the standard deviation
processed_data = (log2_mat - log2_mat.mean(axis=0)) / log2_mat.std(axis=0)
```

## Loading pathways 
```python
# Pre-loaded pathways
# Reactome v78
reactome_pathways  = sspa.process_reactome(organism="Homo sapiens")

# KEGG v98
kegg_human_pathways  = sspa.process_kegg(organism="hsa")
```

Load a custom GMT file (extension .gmt or .csv)
```python
custom_pathways = sspa.process_gmt("wikipathways-20220310-gmt-Homo_sapiens.gmt")
```

Download latest version of pathways
```python
# download KEGG latest
kegg_mouse_latest = sspa.process_kegg("mmu", download_latest=True, filepath=".")

# download Reactome latest
reactome_mouse_latest = sspa.process_reactome("Mus musculus", download_latest=True, filepath=".")
```

## Identifier harmonization 
```python
# download the conversion table
compound_names = processed_data.columns.tolist()
conversion_table = sspa.identifier_conversion(input_type="name", compound_list=compound_names)

# map the identifiers to your dataset
processed_data_mapped = sspa.map_identifiers(conversion_table, output_id_type="ChEBI", matrix=processed_data)
```

## Conventional pathway analysis
Over-representation analysis (ORA)
```python
ora = sspa.sspa_ora(processed_data_mapped, covid_data["Group"], reactome_pathways, 0.05, DA_testtype='ttest', custom_background=None)

# perform ORA 
ora_res = ora.over_representation_analysis()

# get t-test results
ora.ttest_res

# obtain list of differential molecules input to ORA
ora.DA_test_res
```

Gene Set Enrichment Analysis (GSEA), applicable to any type of omics data

```python
sspa.sspa_gsea(processed_data_mapped, covid_data['Group'], reactome_pathways)
```

## Single sample pathway analysis methods
All ssPA methods now have a `fit()`, `transform()` and `fit_transform()` method for compatibility with SciKitLearn. This allows integration of ssPA transformation with various machine learning functions in SKLearn such as `Pipeline` and `GridSearchCV`. Specifically for `sspa.sspa_ssClustPA`, `sspa.sspa_SVD`, and `sspa.sspa_KPCA` methods the model can be fit on the training data and the test data is transformed using the fitted model.

```python
# ssclustPA
ssclustpa_res = sspa.sspa_ssClustPA(reactome_pathways, min_entity=2).fit_transform(processed_data_mapped)

# kPCA 
kpca_scores = sspa.sspa_kpca(reactome_pathways, min_entity=2).fit_transform(processed_data_mapped)

# z-score (Lee et al. 2008)
zscore_res = sspa.sspa_zscore(reactome_pathways, min_entity=2).fit_transform(processed_data_mapped)

# SVD (PLAGE, Tomfohr et al. 2005)
svd_res = sspa.sspa_svd(reactome_pathways, min_entity=2).fit_transform(processed_data_mapped)

# ssGSEA (Barbie et al. 2009)
ssgsea_res = sspa.sspa_ssGSEA(reactome_pathways, min_entity=2).fit_transform(processed_data_mapped)
```


## License
GNU GPL 3.0

## Citing us
[![DOI](https://zenodo.org/badge/442446643.svg)](https://zenodo.org/badge/latestdoi/442446643)

If you found this package useful, please consider citing us:

ssPA package
```
@article{Wieder22a,
   author = {Cecilia Wieder and Nathalie Poupin and Clément Frainay and Florence Vinson and Juliette Cooke and Rachel PJ Lai and Jacob G Bundy and Fabien Jourdan and Timothy MD Ebbels},
   doi = {10.5281/ZENODO.6959120},
   month = {8},
   title = {cwieder/py-ssPA: v1.0.4},
   url = {https://zenodo.org/record/6959120},
   year = {2022},
}
```


Single-sample pathway analysis in metabolomics
```
@article{Wieder2022,
   author = {Cecilia Wieder and Rachel P J Lai and Timothy M D Ebbels},
   doi = {10.1186/s12859-022-05005-1},
   issn = {1471-2105},
   issue = {1},
   journal = {BMC Bioinformatics},
   pages = {481},
   title = {Single sample pathway analysis in metabolomics: performance evaluation and application},
   volume = {23},
   url = {https://doi.org/10.1186/s12859-022-05005-1},
   year = {2022},
}

```

## Contributing
Read our [contributor's guide](https://github.com/cwieder/py-ssPA/blob/main/CONTRIBUTING.md) to get started

### Contributors
We are grateful for our contributors who help develop and maintain py-ssPA:
- Maëlick Brochut [@mbrochut](https://github.com/mbrochut)

## News and updates
<details>
<summary>Read more</summary>

### [v1.0.0] - 25/08/23
- Add compatability with SciKitLearn by implementing `fit()`, `transform()` and `fit_transform()` methods for all ssPA methods. This allows integration of ssPA transformation with various machine learning functions in SKLearn such as `Pipeline` and `GridSearchCV`. Specifically for `sspa.sspa_ssClustPA`, `sspa.sspa_SVD`, and `sspa.sspa_KPCA` methods the model can be fit on the training data and the test data is transformed using the fitted model. 
- Fixed ID conversion bug in `sspa.map_identifiers()` due to MetaboAnalyst API URL change

### [v0.2.4] - 04/07/23
Enable the download of multi-omics (ChEBI and UniProt) Reactome pathways for multi-omics integration purposes. Enable `omics_type='multiomics'` to download:
```
reactome_mouse_latest_mo = sspa.process_reactome("Mus musculus", download_latest=True, filepath=".", omics_type='multiomics')
```

### [v0.2.3] - 23/06/23
- @mbrochut Bug fix in KEGG pathway downloader 
- @mbrochut Add tqdm progress bar for long KEGG downloads

### [v0.2.1] - 05/01/23
- Removal of rpy2 dependency for improved compatibility across systems
- Use [GSEApy](https://github.com/zqfang/GSEApy) as backend for GSEA and ssGSEA 
- Minor syntax changes
   - `ora.ttest_res` is now `ora.DA_test_res` (as we can implement t-test or MWU tests)
   - `sspa_fgsea()` is now `sspa_gsea()` and uses gseapy as the backend rather than R fgsea
   - `sspa_gsva()` is temporarily deprecated due to the need for the rpy2 compatability - use the [GSVA R package](https://bioconductor.org/packages/release/bioc/html/GSVA.html)

</details>
<!-- - Allow download of gene/protein pathways from KEGG and Reactome -->