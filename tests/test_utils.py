from sspa.utils import pathwaydf_to_dict
from sspa.utils import t_tests
import pandas as pd

def test_pathwaydf_to_dict_base():
    dummy_pathway_df = pd.read_csv('test_data/dummy_pathway_df.csv', index_col=0, dtype='object')
    expected = {'R-HSA-1059683': ['30616', '456216'],
 'R-HSA-109581': ['36080', '28494', '61120', '4705', '456216'],
 'R-HSA-109582': ['15366', '15377', '91144', '15379', '15378'],
 'R-HSA-109606': ['36080', '43474', '28494', '15377', '456216']}
    actual = pathwaydf_to_dict(dummy_pathway_df)
    assert expected.keys() == actual.keys()
    assert [set(i) for i in expected.values()] == [set(i) for i in actual.values()]

def test_ttests():
    pass