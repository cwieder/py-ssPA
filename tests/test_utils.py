from sspa.utils import pathwaydf_to_dict, t_tests, load_example_data
import pandas as pd
from pandas.testing import assert_frame_equal
from io import StringIO

class TestUtils():
    dummy_pathway_df_data = """,Pathway_name,0,1,2,3,4\nR-HSA-1059683,Interleukin-6 signaling,30616,456216,,,\nR-HSA-109581,Apoptosis,61120,4705,456216,28494,36080\nR-HSA-109582,Hemostasis,15366,91144,15377,15378,15379\nR-HSA-109606,Intrinsic Pathway for Apoptosis,456216,28494,36080,15377,43474"""
    dummy_pathway_df = pd.read_csv(StringIO(dummy_pathway_df_data), index_col=0, dtype='object', sep=",")
    dummy_metab_data = """sample_id,1372,16610,72665,27823,30915\n
    1004596,-0.8224695827246687,0.3701686486032686,0.5372458799876,0.3622428848480213,-0.5783614704105667\n
    1008097,0.1492912149182915,-0.8026382628405486,-1.220562888082628,-0.774416961826044,-0.1236628685191437\n
    1008631,1.012770830569165,-1.0332911743536763,0.1456895140374768,0.5240497759384892,-0.0735033496314952\n
    1012545,-0.9795544012788546,-0.3958535807477214,2.449097066915658,0.8667518532258393,-0.3540404921571244\n
    1022407,-0.7398166507127327,-0.2887098484131596,-0.7825521081279915,-0.5261277379794408,0.2880482185138752
    """
    dummy_metab = pd.read_csv(StringIO(dummy_metab_data), index_col=0)
    dummy_classes = ['CTRL', 'CASE', 'CASE', 'CTRL', 'CTRL']
    
    def test_pathwaydf_to_dict_base(self):
        expected = {'R-HSA-1059683': ['30616', '456216'],
    'R-HSA-109581': ['36080', '28494', '61120', '4705', '456216'],
    'R-HSA-109582': ['15366', '15377', '91144', '15379', '15378'],
    'R-HSA-109606': ['36080', '43474', '28494', '15377', '456216']}
        actual = pathwaydf_to_dict(self.dummy_pathway_df)
        print(actual)
        assert expected.keys() == actual.keys()
        assert [set(i) for i in expected.values()] == [set(i) for i in actual.values()]

    def test_ttests(self):
        expected_ttest_data = """,Entity,P-value,P-adjust\n
    0,1372,0.02354850918927161,0.11774254594635805\n
    1,16610,0.08512755262079287,0.21281888155198217\n
    2,72665,0.40434649622894775,0.6739108270482462\n
    3,27823,0.6494831059231535,0.752080138459922\n
    4,30915,0.752080138459922,0.752080138459922
    """
        expected = pd.read_csv(StringIO(expected_ttest_data), index_col=0, dtype={'Entity':str, 'P-value':float, 'P-adjust':float})
        actual = t_tests(self.dummy_metab, self.dummy_classes, 'fdr_bh', testtype='ttest')
        assert_frame_equal(actual, expected)

    # def test_loadexampledata(self):
    #     expected_met_raw = pd.read_csv('../src/sspa/example_data/Su_metab_data_raw.csv', index_col=0)
    #     actual_met_raw = load_example_data(omicstype='metabolomics', processed=False)
    #     assert_frame_equal(expected_met_raw, actual_met_raw)