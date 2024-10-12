"""scCyclone Analysis in Python."""

from __future__ import annotations



import scanpy.plotting as pl
import scanpy.tools as tl
import scanpy.preprocessing as pp
import scanpy.get as get
from scanpy import datasets, experimental, external, get, logging, metrics, queries
from scanpy.neighbors import Neighbors
from scanpy.readwrite import read, read_10x_h5, read_10x_mtx, read_visium, write
from scanpy._settings import Verbosity, settings



from anndata import (
    AnnData,
    concat,
    read_csv,
    read_excel,
    read_h5ad,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
    read_umi_tools,
)


from .read  import generate_Iso_adata, generate_PSI_adata, generate_Gene_adata, generate_IF_adata
from . import tools as ctl
from . import get as cget
from . import plotting as cpl




__all__ = [
    "AnnData",
    "concat",
    "read_csv",
    "read_excel",
    "read_h5ad",
    "read_hdf",
    "read_loom",
    "read_mtx",
    "read_text",
    "read_umi_tools",
    "read",
    "read_10x_h5",
    "read_10x_mtx",
    "read_visium",
    "write",
    "datasets",
    "experimental",
    "external",
    "get",
    "logging",
    "metrics",
    "queries",
    "pl",
    "pp",
    "tl",
    "Verbosity",
    "settings",
    "Neighbors",
    "set_figure_params",
    "generate_Iso_adata",
    "generate_PSI_adata",
    "generate_PSI_adata",
    "generate_Gene_adata",
    "generate_IF_adata",
    "cget"
    "cpl",
    "ctl",
    "cget",
]





name = "scCyclone"
__version__ = "1.0.1"
# scCyclone_logo="""
#      _______.  ______   ______ ____    ____  ______  __        ______   .__   __.  _______ 
#     /       | /      | /      |\   \  /   / /      ||  |      /  __  \  |  \ |  | |   ____|
#    |   (----`|  ,----'|  ,----' \   \/   / |  ,----'|  |     |  |  |  | |   \|  | |  |__   
#     \   \    |  |     |  |       \_    _/  |  |     |  |     |  |  |  | |  . `  | |   __|  
# .----)   |   |  `----.|  `----.    |  |    |  `----.|  `----.|  `--'  | |  |\   | |  |____ 
# |_______/     \______| \______|    |__|     \______||_______| \______/  |__| \__| |_______|
                                                                                           
                                                                                                                                                       
# """
# print(scCyclone_logo)
# print(f'Version: {__version__}, Author: Dawn')

