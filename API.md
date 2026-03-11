# scCyclone API（总览）

> 先给出一个可读的 Markdown 版 API 清单，覆盖当前仓库中主要可调用接口。

## 1) 顶层入口（`scCyclone`）

- `generate_Iso_adata(data_path, mode, matrix_prefix="OUT.discovered_transcript_grouped_counts", chunk_size=10000)`
- `generate_PSI_adata(adata, event_info_path)`
- `generate_Gene_adata(adata, var_name="gene_name")`
- `generate_IF_adata(adata, var_name="gene_name", bulk=False, obs_name=None)`

别名子包：
- `scCyclone.get`
- `scCyclone.tl`（tools）
- `scCyclone.pl`（plotting）

---

## 2) 数据构建模块（`scCyclone.read`）

### `generate_Iso_adata(...)`
从 isoform 计数文件构建 `AnnData`。
- 支持 `mode="csv"`（读取 TSV/CSV）
- 支持 `mode="mtx"`（读取 matrix/barcodes/features）

### `generate_PSI_adata(...)`
根据事件信息计算 PSI 并生成新的 `AnnData`。

### `generate_Gene_adata(...)`
将转录本层表达聚合到基因层。

### `generate_IF_adata(...)`
生成 isoform fraction（IF）矩阵，支持单细胞或 bulk 分组。

---

## 3) 结果提取模块（`scCyclone.get` / `scCyclone.get.get`）

- `rank_ifs_groups_df(adata, group=None, key="rank_ifs_groups", pval_cutoff=0.05, min_dif=0, max_dif=1, dpr_cutoff=None, tpr_cutoff=None, rpr_cutoff=None, compare_abs=False)`
- `rank_switchs_groups_df(adata, key="rank_if_switchs_groups")`
- `rank_switch_consequences_groups_df(adata, key="rank_if_switch_consequences_groups")`
- `rank_psis_groups_df(adata, group=None, key="rank_psis_groups", pval_cutoff=None, min_dpsi=0, max_dpsi=1, dpr_cutoff=None, tpr_cutoff=None, rpr_cutoff=None, compare_abs=False)`
- `psis_rmaps_df(event_list, type, gtf_file)`
- `psis_modal_df(adata, groupby, event_list=None, valid_cells=5, groups=None, pkl="PSI_random_forest_model4.pkl")`

---

## 4) 分析工具模块（`scCyclone.tools`）

### 4.1 注释相关
- `add_sqanti3(adata, sqanti3_result_path)`
- `add_gtf(adata, sqanti3_gtf_path)`
- `add_cpc2(adata, cpc2_result_path)`
- `add_pfam(adata, pfam_result_path)`
- `add_deepLoc2(adata, deeploc2_result_path)`
- `add_custom(adata, custom_result_path, feature_colname)`

### 4.2 差异与统计
- `rank_ifs_groups(adata, groupby, groups=None, reference=None, key_added="rank_ifs_groups", valid_cells=5, n_bins=1000, random_seed=0, var_name="gene_name")`
- `rank_psis_groups(adata, groupby, groups=None, reference=None, key_added="rank_psis_groups", percent=0.05, valid_cells=5, n_bins=1000, random_seed=0, var_name="gene_id", two_tailed=False, method="permutation")`
- `rank_if_switchs_groups(adata, group, key="rank_ifs_groups", key_added="rank_if_switchs_groups", pval_cutoff=0.05, dpr_cutoff=0.1, tpr_cutoff=0.2, rpr_cutoff=0.2, abs_min_dif=0.1, abs_max_dif=1)`
- `rank_if_switch_consequences_groups(adata, var_name_list=["ORF_length", "CDS_length", "NMD_status", "coding"], key="rank_if_switchs_groups", key_added="rank_if_switch_consequences_groups")`

### 4.3 评分
- `isoform_entropy_score(adata, groupby, gene_list=None, groups=None, var_name="gene_name")`
- `splice_score(adata_iso, event_dict, score_column="score")`

---

## 5) 可视化模块（`scCyclone.plotting`）

- `dotplot_switch_consequences(adata, key="rank_if_switch_consequences_groups", colcmap="YlOrRd", vmin=0, vmax=1, dot_title=None, colorbar_title="fraction", cluster_order=None, **kwargs)`
- `transcript_structure(adata, gene_name, var_name="gene_name", output_path="./")`

---

## 6) 内部函数（高级/开发者）

以下函数位于下划线模块中，一般不建议作为稳定公共 API 依赖，但当前也可在源码中查看：

- `scCyclone.tools._rank_ifs_groups`: `_bonferroni_safe`, `_two_sided_pvalue`, `_compute_dif`
- `scCyclone.tools._rank_psis_groups`: `_dict_to_rec`, `_filter_event`, `_obs_stats_median_mean`, `_perm_pval_py`, `_pval_approx_normal`, `_compute_dpsi_observed`
- `scCyclone.tools._rank_if_switchs_groups`: `_convert_table`
- `scCyclone.tools._iso_entropy`: `_compute_score`
- `scCyclone.plotting._rank_ifs`: `_fig_show_save_or_axes`, `_get_values_to_plot`, `_rank_ifs_groups_plot`

---

## 7) 说明

- 本文档是“先给一个可读总览”的 Markdown 版本，便于快速查函数名和参数。
- 若需要，我下一步可以把每个函数补成「参数解释 + 返回值 + 最小示例」的详细版（按模块拆分成多文件）。
