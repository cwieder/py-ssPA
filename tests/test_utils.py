from sspa.utils import pathwaydf_to_dict

def test_pathwaydf_to_dict_base():
    dummy_pathway_df = pd.DataFrame({
        "" #...
    })
    expected = [] #...
    actual = pathwaydf_to_dict(dummy_pathway_df)
    assert expected == actual
