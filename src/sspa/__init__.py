from pkg_resources import get_distribution
__version__ = get_distribution('sspa').version

from .process_pathways import process_reactome, process_kegg, process_gmt
from .sspa_cluster import sspa_ssClustPA
from .sspa_kpca import sspa_kpca
from .sspa_zscore import sspa_zscore
from .sspa_svd import sspa_svd
from .utils import load_example_data, t_tests
from .sspa_ora import sspa_ora
from .sspa_gsea import sspa_gsea
from .sspa_ssGSEA import sspa_ssGSEA
from .download_pathways import download_KEGG, download_reactome, MetExplorePaths
from .identifier_conversion import identifier_conversion, map_identifiers