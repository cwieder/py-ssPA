import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import sspa.utils as utils

"""
Module docstring
"""

class sspa_ora:
    """
    Class for overrepresentation analysis 

    Attributes:
        mat (pd.DataFrame): dataframe containing input metabolomics data
        metadata (pd.Series): series containing phenotype metadata e.g 'COVID', 'NON-COVID'
        pathways (pd.DataFrame): pathway dataframe containing compound identifiers
        DA_cutoff (float): pFDR cutoff for selecting differential metabolites e.g. 0.05 or 0.01
        DA_testtype (str): Test type for selecing differential metabolites. Can either be 'ttest' (default) for the independent students T-test or 'mwu' for the Mann Whitney U test, both as implemented in SciPy. 
        custom_background (list): background list of identifiers, default is to use annotated compounds in input data (i.e. mat.columns)
    """
    def __init__(self, mat, metadata, pathways, DA_cutoff, DA_testtype='ttest', custom_background=None):
        self.data = mat
        self.metadata = metadata
        self.pathways = pathways
        self.threshold = DA_cutoff
        self.testtype = DA_testtype
        self.background_set = custom_background if custom_background is not None else mat.columns.to_list()

        self.DA_test_res = utils.t_tests(self.data.copy(deep=True), self.metadata, "fdr_bh", testtype=self.testtype)
        self.DA_molecules = self.DA_test_res[self.DA_test_res["P-adjust"] <= self.threshold]["Entity"].tolist()
        self.results = []

    def over_representation_analysis(self):

        """
        Function for over representation analysis using Fisher exact test (right tailed)
        Returns:
            DataFrame of ORA results for each pathway, p-value, q-value, hits ratio
        """

        pathway_names = self.pathways["Pathway_name"].to_dict()
        pathway_dict = utils.pathwaydf_to_dict(self.pathways)

        # Remove pathways not present in the dataset
        compounds_present = self.DA_molecules
        pathways_present = {k: v for k, v in pathway_dict.items() if len([i for i in compounds_present if i in v]) >= 1}

        pathways_with_compounds = []
        pvalues = []
        pathway_ratio = []
        pathway_coverage = []

        for pathway in pathways_present:
            # perform ORA for each pathway
            pathway_compounds = pathway_dict[pathway]
            pathway_compounds = [i for i in pathway_compounds if str(i) != "nan"]
            if not pathway_compounds or len(pathway_compounds) < 2:
                # ignore pathway if contains no compounds or has less than 3 compounds
                continue
            else:
                DA_in_pathway = len(set(self.DA_molecules) & set(pathway_compounds))
                # k: compounds in DA list AND pathway
                DA_not_in_pathway = len(np.setdiff1d(self.DA_molecules, pathway_compounds))
                # K: compounds in DA list not in pathway
                compound_in_pathway_not_DA = len(set(pathway_compounds) & set(np.setdiff1d(self.background_set, self.DA_molecules)))
                # not DEM compounds present in pathway
                compound_not_in_pathway_not_DA = len(
                    np.setdiff1d(np.setdiff1d(self.background_set, self.DA_molecules), pathway_compounds))
                # compounds in background list not present in pathway
                if DA_in_pathway == 0 or (compound_in_pathway_not_DA + DA_in_pathway) < 2:
                    # ignore pathway if there are no DEM compounds in that pathway
                    continue
                else:
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
                columns=["ID",  "Hits", "Coverage", "P-value", "P-adjust"])
            results["Pathway_name"] = results["ID"].map(pathway_names)
            results.insert(1, 'Pathway_name', results.pop('Pathway_name'))

        except ZeroDivisionError:
            padj = [1] * len(pvalues)
            results = pd.DataFrame(zip(pathways_with_compounds, pathway_ratio, pvalues, padj),
                                columns=["ID", "Hits", "Coverage", "P-value", "P-adjust"])
            results["Pathway_name"] = results["ID"].map(pathway_names)
            results.insert(1, 'Pathway_name', results.pop('Pathway_name'))
            
        self.results = results
        return results