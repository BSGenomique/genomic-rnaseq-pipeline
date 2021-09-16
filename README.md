2021.09

# RNA-Seq data analysis with R
Plateforme Genom'IC, Institut Cochin, Paris, France

This script is used on the Genom'IC plateform to analyse RNA-Seq data. It is the most simple version of the script, which is modified and adapted for each project. We chose to start with a PROJECT.conf and PROJECT.contrast files for the adaptability this kinf of files offer. The projects treated on the plateform are very specific, each with its own number of samples and own design, and this imply a very general script with the possibility to change each step and adapt each methods.

Data used are obtained from `RSEM` software, and the script uses several packages for differential analysis and figures creation. The script is based on [DESEeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [ggplot2](https://rdrr.io/cran/ggplot2/). 

## Main Pipeline

To get ready, you have to decompress the Renv_Archive.zip in the working directory, then load it in R with the restore function
```{r eval=FALSE}
renv::restore()
```
Once the environnement is loaded, with the librairies, please copy your RSEM files in your working directory, then edit the two PROJECT.conf and PROJECT.contrast as follow:

Table Header  | Second Header
------------- | -------------
Table Cell  | Cell 2
Cell 3  | Cell 4 
