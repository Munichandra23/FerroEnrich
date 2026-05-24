# FerroEnrich

FerroEnrich is a web-based R Shiny application for integrated analysis of ferroptosis and senescence programs in liver transcriptomic datasets.

## Overview

FerroEnrich enables users to upload RNA-seq count matrices and metadata, perform differential expression analysis, identify ferroptosis-prone and ferroptosis-resistant genes, calculate a Ferroptosis Index Value (FIV), perform liver-focused ferroptosis and senescence module enrichment, visualize gene interaction networks, and evaluate ferroptosis-senescence cross-talk.

## Main Features

- Upload RNA-seq count matrix and sample metadata
- Run DESeq2-based differential expression analysis
- Identify ferroptosis-prone and ferroptosis-resistant genes
- Calculate Ferroptosis Index Value
- Perform liver ferroptosis and senescence module enrichment
- Visualize gene interaction networks
- Analyze ferroptosis-senescence cross-talk
- Download publication-quality plots and results

## Input Files

### Count Matrix

Rows should represent genes and columns should represent samples.

### Metadata

Rows should represent samples and columns should contain sample-level information such as condition or group.

## Required R Packages

FerroEnrich requires the following major R packages:

- shiny
- DESeq2
- ggplot2
- dplyr
- readr
- tidyr
- stringr
- clusterProfiler
- enrichplot
- AnnotationDbi
- STRINGdb
- visNetwork
- igraph
- DT
- ragg
- Cairo

## How to Run

Open R or RStudio and run:

```r
shiny::runApp()