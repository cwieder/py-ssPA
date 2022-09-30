# sspa
![sspa_logo](sspa_logo.png)

[![DOI](https://zenodo.org/badge/442446643.svg)](https://zenodo.org/badge/latestdoi/442446643)
[![](https://readthedocs.org/projects/pip/badge/?version=latest&style=flat)](https://readthedocs.org/projects/py-ssPA/badge/)

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
   author = {Cecilia Wieder and Rachel P J Lai and Timothy Ebbels},
   doi = {10.1101/2022.04.11.487976},
   journal = {bioRxiv},
   month = {4},
   pages = {2022.04.11.487976},
   publisher = {Cold Spring Harbor Laboratory},
   title = {Single sample pathway analysis in metabolomics : performance evaluation and application},
   url = {https://www.biorxiv.org/content/10.1101/2022.04.11.487976v1 https://www.biorxiv.org/content/10.1101/2022.04.11.487976v1.abstract},
   year = {2022},
}
```