"""scCyclone Analysis in Python."""

from __future__ import annotations



from .read  import generate_Iso_adata, generate_PSI_adata, generate_Gene_adata, generate_IF_adata
from . import tools as tl
from . import get as get
from . import plotting as pl



__all__ = [
    "generate_Iso_adata",
    "generate_PSI_adata",
    "generate_PSI_adata",
    "generate_Gene_adata",
    "generate_IF_adata",
    "get"
    "pl",
    "tl",
    "get",
]




name = "scCyclone"
__version__ = "1.0.3"


scCyclone_logo="""
     _______.  ______   ______ ____    ____  ______  __        ______   .__   __.  _______ 
    /       | /      | /      |\   \  /   / /      ||  |      /  __  \  |  \ |  | |   ____|
   |   (----`|  ,----'|  ,----' \   \/   / |  ,----'|  |     |  |  |  | |   \|  | |  |__   
    \   \    |  |     |  |       \_    _/  |  |     |  |     |  |  |  | |  . `  | |   __|  
.----)   |   |  `----.|  `----.    |  |    |  `----.|  `----.|  `--'  | |  |\   | |  |____ 
|_______/     \______| \______|    |__|     \______||_______| \______/  |__| \__| |_______|
                                                                                           
                                                                                                                                                       
"""
print(scCyclone_logo)
print(f'Version: {__version__}, Author: Dawn')

