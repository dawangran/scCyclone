
# scCyclone: Isoform-Resolved Single-Cell Long-Read Analysis in Python

[![PyPI](https://img.shields.io/pypi/v/scCyclone?logo=pypi)](https://pypi.org/project/scCyclone/)
[![Docs](https://readthedocs.org/projects/sccyclone/badge/?version=latest)](https://sccyclone.readthedocs.io/en/latest/)
![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

---

## 🔬 Overview

**scCyclone** is a comprehensive Python toolkit for **single-cell long-read transcriptomics** (ONT & PacBio).  
It provides isoform-aware quantification, transcript structure annotation, differential transcript usage (DTU), splicing event analysis (DSE), and RBP motif enrichment with rMAPS-compatible output.

Built on top of **Scanpy / AnnData**, scCyclone introduces a standardized multi-layer representation of:

- **Isoform × Cell** matrices  
- **Gene × Cell** matrices  
- **Isoform Fraction (IF)** matrices  
- **PSI (Percent Spliced In)** event-level matrices  

---

## ✨ Key Features

### **1. Multi-layer AnnData Generation**
Convert long-read quantification into structured AnnData matrices:

- `generate_Iso_adata` — isoform count matrix  
- `generate_Gene_adata` — aggregated gene expression  
- `generate_IF_adata` — isoform fraction per gene per cell  
- `generate_PSI_adata` — event-level PSI using SUPPA2-style IOE files  

---

### **2. SQANTI3-based Isoform Annotation**
Automatically incorporates:

- Structural categories (FSM, ISM, NIC, NNC…)  
- Exon/intron structure  
- TSS/TES information  
- CDS / UTR lengths  
- Splice site features  

---

### **3. Differential Transcript Usage (DTU)**
`tl.rank_ifs_groups` provides robust DTU testing using:

- Two-sided permutation tests  
- NA-aware denominator logic  
- Gene-wise normalization  
- Bonferroni & FDR correction  

---

### **4. Differential Splicing Event Analysis (DSE)**
Supports SUPPA2 event types:

- A3 / A5  
- SE  
- RI  
- MX  
- AF / AL  

Output: PSI matrix + statistical testing (API growing).

---

### **5. RBP Analysis & rMAPS Output**
Generate rMAPS-compatible files for motif-centric RBP enrichment analysis.

---

## 📦 Installation


```bash
git clone https://github.com/dawangran/scCyclone
cd scCyclone
pip install -e .
```

Requires Python ≥ 3.8.

---

## 🚀 Quick Start

### **1. Build Isoform AnnData**

```python
import scCyclone as scc

iso_file = "transcript_model_grouped_counts.tsv"
adata_iso = scc.generate_Iso_adata(iso_file)
```

---

### **2. Add SQANTI3 annotation & gene names**

```python
import pandas as pd
import scanpy as sc

scc.tl.add_sqanti3(adata_iso, "sqanti_classification.txt")

gene_map = pd.read_csv("gene_id_translate.txt", sep=" ", header=None)
gene_map = gene_map.rename(columns={0: "associated_gene", 1: "gene_name"})
adata_iso.var["gene_name"] = pd.merge(
    adata_iso.var, gene_map, on="associated_gene", how="left"
)["gene_name"]

adata_iso = adata_iso[:, adata_iso.var["gene_name"].notna()]
```

---

### **3. Build Gene / IF / PSI Layers**

```python
adata_gene = scc.generate_Gene_adata(adata_iso, var_name="gene_name")
adata_IF   = scc.generate_IF_adata(adata_iso, var_name="gene_name")
adata_PSI  = scc.generate_PSI_adata(adata_iso, "events.ioe")
```

---

### **4. Perform DTU (Differential Transcript Usage)**

```python
adata_iso.obs["condition"] = pd.Categorical(adata_iso.obs["condition"])

scc.tl.rank_ifs_groups(
    adata_iso,
    groupby="condition",
)
```

Results appear under:

```python
adata_iso.uns["rank_ifs_groups"]
```

---

## 📁 Input Requirements

| Input Type | Description |
|-----------|-------------|
| Isoform counts | IsoQuant `*grouped_counts.tsv` (recommended) |
| Annotation | SQANTI3 extended classification |
| Event definitions | SUPPA2 `.ioe` files |
| Gene mapping | Two-column table: associated_gene → gene_name |
| Metadata | Any per-cell metadata (condition/cluster/time) |

---

## 🧠 API Overview

### **Matrix Construction**

| Function | Description |
|----------|-------------|
| `generate_Iso_adata` | Build isoform × cell AnnData |
| `generate_Gene_adata` | Aggregate isoforms → gene |
| `generate_IF_adata` | Compute isoform fractions |
| `generate_PSI_adata` | Build event-level PSI matrix |

---

### **Tools (tl module)**

| Function | Description |
|----------|-------------|
| `tl.add_sqanti3` | Add SQANTI3 annotation |
| `tl.rank_ifs_groups` | Differential isoform usage |
| *(more coming)* | DSE testing, RBP analysis |

---

## 📊 Visualization

```python
import scanpy as sc

top_iso = adata_iso.uns["rank_ifs_groups"]["names"]["WT"][:20]

sc.pl.matrixplot(
    adata_IF,
    var_names=top_iso,
    groupby="condition",
    standard_scale="var"
)
```



---

## 🌏 简体中文简介

scCyclone 是一个面向 **单细胞长读长测序** 的完整分析工具包，支持：

- 构建 *转录本 / 基因 / IF / PSI* 多层级 AnnData  
- 单细胞层面的 **差异转录本使用（DTU）**  
- **可变剪接事件** PSI 计算与显著性分析  
- 结合 SQANTI3 的结构注释  
- 输出 rMAPS 兼容文件做 RBP motif 富集  

适合从“读表达量”升级到“读结构 + 调控机制”的研究。

---

## 📄 License

MIT License.
