2020.04

# RNA-Seq data analysis with R
Plateforme Genom'IC, Institut Cochin, Paris, France

This script is used on the Genom'IC plateform to analyse RNA-Seq data from normalization to supervised analysis. It is the most simple version of the script, which is modified and adapted for each project and each design.

Data used are obtained from `RSEM` software, and the script uses several packages for differential analysis and figures creation. The script is based on [DESEeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [ggplot2](https://rdrr.io/cran/ggplot2/)

## Packages
```
library(biomaRt)
library(tximport)
library(readr)
library(DESeq2)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(factoextra)
```
