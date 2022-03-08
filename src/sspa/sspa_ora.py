import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import sspa.utils as utils

def sspa_ora(mat, metadata_column, pathways, cutoff_thresh=0.05, custom_bgset=None):
    t_test_res = utils.t_tests(mat, metadata_column, "fdr_bh")
    DA_compounds = t_test_res[t_test_res["P-adjust"] <= cutoff_thresh]["Entity"].tolist()
    if not custom_bgset:
        bg_set = mat.columns.to_list()
        ora_res = over_representation_analysis(DA_compounds, bg_set, pathways)
    else:
        ora_res = over_representation_analysis(DA_compounds, custom_bgset, pathways)
    ora_res = ora_res.rename(columns={"Pathway_ID": "ID"})
    return ora_res

def over_representation_analysis(DA_list, background_list, pathways_df):
    """
    Function for over representation analysis using Fisher exact test (right tailed)
    :param DA_list: List of differentially abundant metabolite IDENTIFIERS
    :param background_list: background list of IDENTIFIERS
    :param pathways_df: pathway dataframe containing compound identifiers
    :return: DataFrame of ORA results for each pathway, p-value, q-value, hits ratio
    """

    pathway_names = pathways_df["Pathway_name"].to_dict()
    pathway_dict = utils.pathwaydf_to_dict(pathways_df)

    # Remove pathways not present in the dataset
    compounds_present = DA_list
    pathways_present = {k: v for k, v in pathway_dict.items() if len([i for i in compounds_present if i in v]) >= 1}

    pathways_with_compounds = []
    pvalues = []
    pathway_ratio = []
    pathway_count = 0
    pathway_coverage = []

    for pathway in pathways_present:
        # perform ORA for each pathway
        pathway_compounds = pathway_dict[pathway]
        pathway_compounds = [i for i in pathway_compounds if str(i) != "nan"]
        if not pathway_compounds or len(pathway_compounds) < 2:
            # ignore pathway if contains no compounds or has less than 3 compounds
            continue
        else:

            DA_in_pathway = len(set(DA_list) & set(pathway_compounds))
            # k: compounds in DA list AND pathway
            DA_not_in_pathway = len(np.setdiff1d(DA_list, pathway_compounds))
            # K: compounds in DA list not in pathway
            compound_in_pathway_not_DA = len(set(pathway_compounds) & set(np.setdiff1d(background_list, DA_list)))
            # not DEM compounds present in pathway
            compound_not_in_pathway_not_DA = len(
                np.setdiff1d(np.setdiff1d(background_list, DA_list), pathway_compounds))
            # compounds in background list not present in pathway
            if DA_in_pathway == 0 or (compound_in_pathway_not_DA + DA_in_pathway) < 2:
                # ignore pathway if there are no DEM compounds in that pathway
                continue
            else:
                pathway_count += 1
                # Create 2 by 2 contingency table
                pathway_ratio.append(str(DA_in_pathway) + "/" + str(compound_in_pathway_not_DA + DA_in_pathway))
                pathway_coverage.append(
                    str(compound_in_pathway_not_DA + DA_in_pathway) + "/" + str(len(pathway_compounds)))
                pathways_with_compounds.append(pathway)
                contingency_table = np.array([[DA_in_pathway, compound_in_pathway_not_DA],
                                              [DA_not_in_pathway, compound_not_in_pathway_not_DA]])
                # Run right tailed Fisher's exact test
                oddsratio, pvalue = stats.fisher_exact(contingency_table, alternative="greater")
                pvalues.append(pvalue)
    try:
        padj = sm.stats.multipletests(pvalues, 0.05, method="fdr_bh")
        results = pd.DataFrame(
            zip(pathways_with_compounds, pathway_ratio, pathway_coverage, pvalues,
                padj[1]),
            columns=["Pathway_ID",  "Hits", "Coverage", "P-value", "P-adjust"])
        results["Pathway_name"] = results["Pathway_ID"].map(pathway_names)
        results.insert(1, 'Pathway_name', results.pop('Pathway_name'))

    except ZeroDivisionError:
        padj = [1] * len(pvalues)
        results = pd.DataFrame(zip(pathways_with_compounds, pathway_ratio, pvalues, padj),
                               columns=["Pathway_ID", "Hits", "Coverage", "P-value", "P-adjust"])
        results["Pathway_name"] = results["Pathway_ID"].map(pathway_names)
        results.insert(1, 'Pathway_name', results.pop('Pathway_name'))
        
    return results