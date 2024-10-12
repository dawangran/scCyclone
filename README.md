[![PyPI](https://img.shields.io/pypi/v/scCyclone?logo=PyPI)](https://pypi.org/project/scCyclone)
[![Documentation Status](https://readthedocs.org/projects/sccyclone/badge/?version=latest)](https://sccyclone.readthedocs.io/en/latest/?badge=latest)



# scCyclone: Single-Cell Cyclone Analysis in Python


scCyclone is a comprehensive Python package designed to analyze single-cell full-length transcriptome data.


## Features


* **Personalized Matrix Generation**: This feature allows users to create custom matrices tailored to their specific dataset requirements, providing flexibility in handling various types of single-cell sequencing data.

* **Differential Transcript Usage (DTU) Analysis**: This tool helps identify variations in transcript usage across different cell populations, enabling researchers to understand how gene expression patterns differ under different conditions.

* **Functional and Structural Analysis of Differential Transcripts**: This feature goes beyond simple identification of differential transcripts by providing a deeper analysis of their functional and structural implications. This can help researchers understand the biological significance of the observed transcript usage differences.

* **Differential Splicing Event (DSE) Analysis**: This feature is designed to detect differences in splicing patterns within single-cell data, which is crucial for identifying splice variants that may be associated with specific cellular processes or conditions.

* **Modal Analysis of Splicing Events**: This feature involves examining the distribution of splicing events to identify common patterns or central tendencies. This can provide insights into the most frequent or significant splicing patterns within the dataset.

* **RNA Binding Protein (RBP) Analysis with rMAPS Output**: This feature allows for the analysis of RNA binding protein (RBP) binding patterns and provides output files compatible with rMAPS（http://rmaps.cecsresearch.org/）, a popular tool for differential splicing analysis. This can help researchers understand the role of RBPs in regulating gene expression and splicing events.



## Getting Started

### Prerequisites

Before you begin, ensure you have met the following requirements:
- Python 3.8 or higher
- A compatible operating system (e.g., Windows, macOS, Linux)

* [scanpy](https://pythonhosted.org/gffutils/installation.html)
* [dotplot](https://github.com/kn-bibs/dotplot)


### Installation

To install the software, follow these steps:

```
pip install scCyclone
```

### Documentation
Read the documentation : https://sccyclone.readthedocs.io

## License
scCyclone is released under the MIT License.

## Contact
For any questions or feedback, please reach out to us at 605547565@qq.com.