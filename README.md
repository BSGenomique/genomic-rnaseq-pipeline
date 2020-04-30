2020.04

# RNA-Seq data analysis with R
Plateforme Genom'IC, Institut Cochin, Paris, France

This script is used on the Genom'IC plateform to analyse RNA-Seq data. It is the most simple version of the script, which is modified and adapted for each project. We chose to start with a PROJECT.conf and PROJECT.contrast files for the adaptability this kinf of files offer. The projects treated on the plateform are very specific, each with its own number of samples and own design, and this imply a very general script with the possibility to change each step and adapt each methods.

Data used are obtained from `RSEM` software, and the script uses several packages for differential analysis and figures creation. The script is based on [DESEeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [ggplot2](https://rdrr.io/cran/ggplot2/). 

## Pipeline big lines

We start with the RSEM results files, and import them in R with the [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) package. [DESEeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) is able to handle these data, so we use this package for the whole analysis, from normalization to differential expression.

We also use the [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) package, since we choose to align our data on the [Ensembl](https://www.ensembl.org/index.html) references. This choice has been made for the easy accession Ensembl provide to get fasta, gtf and APIs. This part of the script is often modified to add informations on the several matrices, informations needed by the project team, so specific to this very project.

At last the [ggplot2](https://rdrr.io/cran/ggplot2/) package is used to generate figures as PCA, Volcano & MA plots.

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
