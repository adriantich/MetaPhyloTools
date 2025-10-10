# MetaPhyloTools

<img src="logo.png" align="right" width="120"/>

<!-- badges: start -->
<!-- badges: end -->

## Overview

MetaPhyloTools is an R package that provides a collection of tools to process and analyse phylogeographic patterns from metabarcoding data, thus analysing hundreds of Operational Taxonomic Units at the same time. This package provides functions to obtain Molecular Operational Taxonomic Units (MOTUs) and the Exact Sequence Variants (ESV) that cluster into theese MOTUs from fasta files (in development). In addition, the package includes functions for data normalization, phylogenetic analysis, and visualization of results.

## Philosophy

The package is designed to work with standard data table formats (data.frames) as input and output publication-ready ggplot2 visualizations, making it easy to incorporate into existing R workflows and customize plots according to specific research needs.


## Installation

You can install the development version of MetaPhyloTools from GitHub:

```r
# install.packages("devtools")
devtools::install_github("adriantich/MetaPhyloTools")
```

## Main Functions

- `edenetwork_percolation()`: Calculate percolation thresholds for ecological networks
- `pairwise_djost()`: Calculate pairwise Jost's D differentiation indices
- `haplo_ggplot()`: Create haplotype network visualizations
- `haploNet_data()` & `haploNet_plot()`: Process and plot haplotype networks
- `rarefy_within_motu()`: Rarefy data within MOTUs
- `print_network()` & `print_network_ggplot()`: Network printing utilities

### Developing
- `declu_seq()`: Denoising + Clustering of Sequences
- `entropy_calculator()`: Computes the entropy values of each position in the codon
- `d_calibrator()`: rapid assessment of d curves for swarm clustering
- `fasta2tab()`: converts fasta files into table format using vsearch

## Quick Example

```r
library(MetaPhyloTools)

# Calculate percolation threshold
dist_mat <- as.matrix(dist(iris[1:10, 1:4]))
result <- edenetwork_percolation(dist_mat)
```

## License

GPL-3

## Author

Adrià Antich (a.antich@ceab.csic.es)