.. scCyclone documentation master file, created by
   sphinx-quickstart on Sat Oct 12 15:08:57 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

scCyclone documentation
=======================

Add your content using ``reStructuredText`` syntax. See the
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
documentation for details.

scCyclone is a comprehensive Python package designed to analyze single-cell full-length transcriptome data.


* **Personalized Matrix Generation**: This feature allows users to create custom matrices tailored to their specific dataset requirements, providing flexibility in handling various types of single-cell sequencing data.

* **Differential Transcript Usage (DTU) Analysis**: This tool helps identify variations in transcript usage across different cell populations, enabling researchers to understand how gene expression patterns differ under different conditions.

* **Functional and Structural Analysis of Differential Transcripts**: This feature goes beyond simple identification of differential transcripts by providing a deeper analysis of their functional and structural implications. This can help researchers understand the biological significance of the observed transcript usage differences.

* **Differential Splicing Event (DSE) Analysis**: This feature is designed to detect differences in splicing patterns within single-cell data, which is crucial for identifying splice variants that may be associated with specific cellular processes or conditions.

* **Modal Analysis of Splicing Events**: This feature involves examining the distribution of splicing events to identify common patterns or central tendencies. This can provide insights into the most frequent or significant splicing patterns within the dataset.

* **RNA Binding Protein (RBP) Analysis with rMAPS Output**: This feature allows for the analysis of RNA binding protein (RBP) binding patterns and provides output files compatible with rMAPS（http://rmaps.cecsresearch.org/）, a popular tool for differential splicing analysis. This can help researchers understand the role of RBPs in regulating gene expression and splicing events.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   notebooks/tutorial