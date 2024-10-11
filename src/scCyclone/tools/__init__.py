from __future__ import annotations


from ._iso_annotation import add_sqanti3, add_cpc2, add_pfam, add_deepLoc2, add_custom, add_gtf
from ._rank_ifs_groups import rank_ifs_groups
from ._rank_psis_groups import rank_psis_groups
from ._rank_swiths_groups import rank_switchs_groups, rank_switch_consequences_groups


__all__ = [
    "add_sqanti3",
    "add_gtf",
    "add_cpc2",
    "add_pfam",
    "add_deepLoc2",
    "add_custom",
    "rank_ifs_groups",
    "rank_psis_groups",
    "rank_switchs_groups",
    "rank_switch_consequences_groups",
]
