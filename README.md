# Adaptive-STKDE
Performing STKDE with adaptive bandwidths (ASTKDE) in R.
- <h1><strong>Introduction</strong></h1>

&nbsp;&nbsp;&nbsp;&nbsp;This repository corresponds to the paper “*Urban Road Collapse Disaster Hotspot Analysis in Zhengzhou: A Novel Adaptive Bandwidth Spatiotemporal Kernel Density Estimation Method*”.  

&nbsp;&nbsp;&nbsp;&nbsp;It includes the implementation of ASTKDE and TSTKDE, along with the code for all experiments presented in the paper. Results demonstrate that ASTKDE outperforms TSTKDE, particularly in the analysis of spatiotemporally inhomogeneous data.

&nbsp;&nbsp;&nbsp;&nbsp;Computation is extremely fast thanks to R’s broadcasting mechanism. In particular, parallel processing has also been implemented, allowing multi-threaded STKDE calculations for further performance improvements.

&nbsp;&nbsp;&nbsp;&nbsp;TSTKDE has already seen widespread use in areas such as criminology, traffic accident analysis, and epidemiological studies. However, ASTKDE has received far less attention. By open-sourcing the code, we hope to foster broader adoption and application of ASTKDE.

&nbsp;&nbsp;&nbsp;&nbsp;The **"Code"** folder contains R scripts, divided into multiple source files. 

&nbsp;&nbsp;&nbsp;&nbsp;The **"testDATA"** folder includes the study area (a closed counterclockwise polygon) and a set of observation points. Note that the data provided is **NOT** the actual dataset used in the paper; it is for demonstration purposes only.

- <h1><strong>Functions</strong></h1>

**Data.R:** Provides essential variables that must be generated in R's Global Environment.

**STKDE.R:** Contains functions for generating ASTKDE & TSTKDE and performing other key computations.

**Synthetic.R:** Generates synthetic PDFs from the generic formula in the paper and calculates both global and local ISE.

**Prediction.R:** Performs prediction iterations on future unseen data.

**MonteCarlo.R:** Performs hotspot identification using Monte Carlo simulations.

**Visualization.R:** Plots contour lines, isosurfaces, and hotspot trajectories in a 3D view.

**Export.R:** Exports the density cube for further rendering.

2025-03-03
