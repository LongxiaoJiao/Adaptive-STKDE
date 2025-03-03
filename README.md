# Adaptive-STKDE
Performing STKDE with adaptive bandwidths (ASTKDE) in R.
- **Introduction**

&nbsp;&nbsp;&nbsp;&nbsp;This repository corresponds to the paper “*Urban Road Collapse Disaster Hotspot Analysis in Zhengzhou: A Novel Adaptive Bandwidth Spatiotemporal Kernel Density Estimation Method*”.  

&nbsp;&nbsp;&nbsp;&nbsp;It contains the implementation of ASTKDE and the code for all experiments presented in the paper. Its goal is to illustrate the advantages of ASTKDE over the traditional fixed-bandwidth STKDE (TSTKDE), particularly for analyzing spatiotemporally inhomogeneous data.

&nbsp;&nbsp;&nbsp;&nbsp;Computation is extremely fast thanks to R’s broadcasting mechanism. In particular, parallel processing has also been implemented, allowing multi-threaded STKDE calculations for further performance improvements.

&nbsp;&nbsp;&nbsp;&nbsp;TSTKDE has already seen widespread use in areas such as criminology, traffic accident analysis, and epidemiological studies. However, ASTKDE has received far less attention. Our findings show that ASTKDE outperforms TSTKDE when working with spatiotemporally inhomogeneous data. By open-sourcing the code, we hope to foster broader adoption and application of ASTKDE.

2025-03-03
