# Single sample pathway analysis toolkit
sspa provides a Python interface for metabolomics pathway analysis. In addition to conventional methods over-representation analysis (ORA) and gene/metabolite set enrichment analysis (GSEA), it also provides a wide range of single-sample pathway analysis (ssPA) methods. 

## Features
- Over-representation analysis
- Metabolite set enrichment analysis (based on GSEA)
- Single-sample pathway analysis
- Compound identifier conversion
- Pathway database download (KEGG, Reactome, and MetExplore metabolic networks)

``` mermaid 
graph LR 
  A[Download pathways ] --> C;
  B[Input data] --> G[Identifier conversion];
  G -->C[Pathway analysis];
  C --> D[ORA];
  C --> E[GSEA];
  C --> F[ssPA];
```

!!! note

    Although this package is designed to provide a user-friendly interface for metabolomics pathway analysis, the methods are also applicable to other datatypes such as normalised RNA-seq data. Gene/protein pathway collections can be input as .GMT files (see tutorials).



## Installation
```
pip install sspa
```

## Citing us
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