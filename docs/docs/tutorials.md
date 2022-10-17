
Please see our interactive walkthrough tutorial for a detailed overview on the package functionalities:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1X_ZBjXYxBVv43LcArRv2aB2HJkUp8A4x?usp=sharing)

A local version of this tutorial is available to download from [GitHub](https://github.com/cwieder/py-ssPA/blob/main/sspa_walkthrough_local.ipynb).

## Quickstart quide

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

## Input data
The input data should consist of an $m×n$ pandas DataFrame of values corresponding to the abundance of all annotated metabolites, either as absolute (concentrations) or relative (peak intensities) quantifications. Rows represent $m$ samples and columns represent $n$ annotated compounds. 

Sample metadata is necessary for some types of pathway analysis, particularly the conventional methods over-representation analysis (ORA) and gene set enrichment analysis (GSEA). Sample metadata should contain details of the clinical outcome or phenotype used for comparison and should be represented as a column either as part of the metabolite abundance DataFrame, or as a separate column/Series, as long as the sample identifiers can be matched between data structures. 

??? example

    Below is a condensed example of the input data format:

    | sample\_id | spermidine | 1-methylnicotinamide | 12,13-DiHOME | alpha-ketoglutarate | kynurenate | Group |
    | ---------- | ---------- | -------------------- | ------------ | ------------------- | ---------- | ----- |
    | 1004596    | \-0.756979 | 0.552163             | \-0.317382   | 0.726321            | \-0.608606 | A     |
    | 1008097    | 0.079818   | \-0.839393           | 0.49128      | \-1.867786          | \-0.044496 | A     |
    | 1008631    | 0.978372   | \-1.281277           | \-0.199487   | 0.355229            | 0.014784   | B     |
    | 1012545    | \-0.93754  | \-0.242391           | 1.63653      | 2.080704            | \-0.31561  | B     |
    | 1022407    | \-0.652496 | \-0.110733           | 0.814461     | \-0.886903          | 0.409608   | B     |

    The index are sample identifiers and the column names are metabolite names/identifiers. Metdata columns can also be added at the end of the DataFrame.

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

!!! note

    Downloading the lastest version of KEGG pathways can take up to ten minutes. 

## Identifier harmonization 
```
# download the conversion table
compound_names = processed_data.columns.tolist()
conversion_table = sspa.identifier_conversion(input_type="name", compound_list=compound_names)

# map the identifiers to your dataset
processed_data_mapped = sspa.map_identifiers(conversion_table, output_id_type="ChEBI", matrix=processed_data)
```

!!! warning

    This step requries active internet connection and access to the MetaboAnalyst API. If you are experiencing errors please check you are able to query the [API](https://www.metaboanalyst.ca/docs/APIs.xhtml).

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
ssclustpa_res = sspa.sspa_ssClustPA(processed_data_mapped, reactome_pathways)

# kPCA
kpca_scores = sspa.sspa_kpca(processed_data_mapped, reactome_pathways)

# z-score
zscore_res = sspa.sspa_zscore(processed_data_mapped, reactome_pathways)

# SVD (PLAGE)
svd_res = sspa.sspa_svd(processed_data_mapped, reactome_pathways)

# GSVA
gsva_res = sspa.sspa_gsva(processed_data_mapped, reactome_pathways)
```