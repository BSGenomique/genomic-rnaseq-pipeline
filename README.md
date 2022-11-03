2021.09

# RNA-Seq data analysis with R
Plateforme Genom'IC, Institut Cochin, Paris, France

This script is used on the Genom'IC plateform to analyse RNA-Seq data. It is the most simple version of the script, which is modified and adapted for each project. We chose to start with a PROJECT.conf and PROJECT.contrast files for the adaptability this kinf of files offer. The projects treated on the plateform are very specific, each with its own number of samples and own design, and this imply a very general script with the possibility to change each step and adapt each methods.

Data used are obtained from `RSEM` software, and the script uses several packages for differential analysis and figures creation. The script is based on [DESEeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [ggplot2](https://rdrr.io/cran/ggplot2/). 

## Main Pipeline

### Data & File configuration

To get ready, you have to decompress the Renv_Archive.zip in the working directory, then load it in R with the restore function
```{r eval=FALSE}
renv::restore()
```
Once the environnement is loaded, with the librairies, please copy your RSEM files in your working directory, then edit the two PROJECT.conf and PROJECT.contrast as follow, conserving the header.

PROJECT.conf: add one column by Condition you want to add in your design

File | Name | Condition
------------- | ------------- | -------------
RSEM Filename  | Sample name  | Condition, treatment, time

PROJECT.contrast: add one line for each comparison you want to make

Condition name  | Group to compare  | Control group
------------- | ------------- | -------------

You will also modify the following line to fit to your reference:
```{r eval=FALSE}
ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version = 101)
```
During the data importation, there a filter to reduce the data size (and the results matrix) and to improve the computation time. This threshold filter out genes where there are less than 3 samples with normalized counts greater than or equal to 10.

### Normalization

The normalization matrix is extracted from the `dds` object, and the `biomart` package is used to add information for each sample (gene names, description, and GO terms by default). If you need more informations, you can get the different attributes available with the `listAttributes(ensembl)` function, then modify the following line: 
rtyèrty
```
  IDsWithNamesDesc <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol','description','go_id','name_1006'), filters = 'ensembl_gene_id', values = ensemblIDs, mart = ensembl,useCache = FALSE)
```

### Unsupervised analysis

To make the clustering and the PCA, we need to transform the data with the `rld` or `vst` functions, then custom functions will allow you to create the corresponding figures. You can easily modify the source code of these functions in the GenomIC_RNASeq-pipeline_functions.r file, with the `ggplot2` package.

### Supervised analysis

A for loop is parsing the PROJECT.contrast file, to extract the good results matrix and to generate corresponding figures. Then, once the loop is over, do not forget to modify the report.rmd file to generate you project report.

## Main packages

RSEM    # Li, B., Dewey, C.N. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 12, 323 (2011). https://doi.org/10.1186/1471-2105-12-323

DESeq2  # Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8

ggplot2 # Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4


