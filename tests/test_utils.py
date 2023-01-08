from sspa.utils import pathwaydf_to_dict, t_tests, load_example_data
import pandas as pd
from pandas.testing import assert_frame_equal

class TestUtils():
    dummy_pathway_df_data = '''
        ,Pathway_name,0,1,2,3,4
        R-HSA-1059683,Interleukin-6 signaling,30616,456216,,,
        R-HSA-109581,Apoptosis,61120,4705,456216,28494,36080
        R-HSA-109582,Hemostasis,15366,91144,15377,15378,15379
        R-HSA-109606,Intrinsic Pathway for Apoptosis,456216,28494,36080,15377,43474
    '''
    dummy_pathway_df = pd.read_csv(pd.compat.StringIO(dummy_pathway_df_data), index_col=0, dtype='object')
    # dummy_metab = pd.read_csv('test_data/mini_metabolomics.csv', index_col=0)
    # dummy_classes = ['CTRL', 'CASE', 'CASE', 'CTRL', 'CTRL']

    def test_pathwaydf_to_dict_base(self):
        expected = {'R-HSA-1059683': ['30616', '456216'],
    'R-HSA-109581': ['36080', '28494', '61120', '4705', '456216'],
    'R-HSA-109582': ['15366', '15377', '91144', '15379', '15378'],
    'R-HSA-109606': ['36080', '43474', '28494', '15377', '456216']}
        
        actual = pathwaydf_to_dict(self.dummy_pathway_df)
        assert expected.keys() == actual.keys()
        assert [set(i) for i in expected.values()] == [set(i) for i in actual.values()]

    # def test_ttests(self):
    #     expected = pd.read_csv('test_data/ttest_expected.csv', index_col=0, dtype={'Entity':str, 'P-value':float, 'P-adjust':float})
    #     actual = t_tests(self.dummy_metab, self.dummy_classes, 'fdr_bh', testtype='ttest')
    #     assert_frame_equal(actual, expected)

    # def test_loadexampledata(self):
    #     expected_met_raw = pd.read_csv('../src/sspa/example_data/Su_metab_data_raw.csv', index_col=0)
    #     actual_met_raw = load_example_data(omicstype='metabolomics', processed=False)
    #     assert_frame_equal(expected_met_raw, actual_met_raw)