source("GenomIC_RNASeq-pipeline_functions.r")


confFile <- "PROJET.conf"
contrastFile <- "PROJET.contrastes"
projectName <- "NGS21-104"

ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version = 101)
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 101)


# Step 1: Preparing the configuration table with conditions and factors to take in account

  currentDirectory <- getwd()
  configuration <- read.table(confFile,header=TRUE,sep="\t",colClasses = c("character","character","factor")) # Change the colClasses parameter to fit with the project design
  files <- file.path(currentDirectory, configuration$File)
  names(files) <- configuration$Name
  
# Step 2: Import the data

  txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
  txi.rsem$length <- ifelse(txi.rsem$length==0,1,txi.rsem$length) # Setting the length to 1, to avoid error for the dds object

# Step 3: Creating the DDS object

  dds <- DESeqDataSetFromTximport(txi.rsem, configuration, ~ Preparation + Condition)
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 3 # filter out genes where there are less than 3 samples with normalized counts greater than or equal to 10.
  dds <- dds[idx,]
  dds <- DESeq(dds) # no need to set the Reference group, as the script use the contrast parameter later
  
# Step 4: Save the normalized matrix
  countNorm<-counts(dds, normalized=TRUE)
  ensemblIDs <- row.names(countNorm)
  
  #Human
  # IDsWithNamesDesc <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol','description','go_id','name_1006'), filters = 'ensembl_gene_id', values = ensemblIDs, mart = ensembl,useCache = FALSE)
  #Mouse
  IDsWithNamesDesc <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol','description','go_id','name_1006'), filters = 'ensembl_gene_id', values = ensemblIDs, mart = ensembl,useCache = FALSE)
  
  
  rownames(IDsWithNamesDesc) <- make.names(IDsWithNamesDesc$ensembl_gene_id, unique=TRUE)
  dataMerged <- merge(as.data.frame(countNorm),IDsWithNamesDesc,by="row.names",all.x=TRUE)
  nameNorm <- paste(projectName,"_deseq2_NormalizedMatrix.tsv",sep="")
  write.table(dataMerged,file=nameNorm,sep='\t',row.names=F)
  
  
# Step 7: Unsupervised analysis
  
  make_nice_boxplot(dds,configuration,projectName)
  
  rld<- rlog(dds,blind=TRUE)
  
  make_nice_clusters(rld,configuration,projectName)

  make_nice_pca(dds,rld,projectName)
  
  # Warning! Before using this function, please modify it with the right Effect name, in the functions file
  make_PCA_noBatchEffect(dds,rld,projectName)
  
# Step 8: Supervised analysis

  contrasteList <- read.table(contrastFile,header=FALSE,sep="\t",colClasses = c("character","character","character"))
  nComp <- nrow(contrasteList)

# Step 8.5: with many factors, a proposed solution is to merge two of them ->

  # Many factors
  # dds$Condition <- factor(paste0(dds$Traitement, dds$Temps))
  # design(dds) <- ~ Condition
  # dds <- DESeq(dds)
  # resultsNames(dds)
  # 
  # configuration$Condition <- factor(paste0(configuration$Traitement, configuration$Temps))


  for (z in 1:nComp){
    
    ctrst <- unname(unlist(contrasteList[z,]))
    res <- results(dds, contrast=ctrst)

    name <- paste(projectName,"deseq2_results_contrast",ctrst[2],'vs',ctrst[3],sep="_")
    
    resLFC <- lfcShrink(dds, contrast=ctrst, res=res)
    nameMA <- paste(name,"MAplot.png",sep="_")
    png(filename=nameMA,width=7 ,height=7, units="in",res = 600 ) 
    plotMA(resLFC,alpha = 0.05, ylim=c(-10,10),main=paste(ctrst[2],'vs',ctrst[3],sep=" "))
    dev.off()
    
    res <- merge(as.data.frame(res),IDsWithNamesDesc,by="row.names",all.x=TRUE)
    
    row.names(res)<-res$Row.names
    res<-res[,2:10]
    
    fullData <- merge(as.data.frame(res),countNorm,by="row.names",all.x=TRUE)
    row.names(fullData)<-fullData$Row.names
    
    fullData<-fullData[order(fullData$padj),]
    write.table(as.data.frame(fullData),file=paste(name,".tsv",sep=""),sep='\t',row.names=F)
    
    make_nice_volcanoPlot(fullData,name,ctrst)
    make_nice_diffPlot(dds,fullData,configuration,name,ctrst) 
    
  }
  
  
render("report.Rmd", output_file=paste0(projectName,".html"))

######
# RSEM    # Li, B., Dewey, C.N. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 12, 323 (2011). https://doi.org/10.1186/1471-2105-12-323
# DESeq2  # Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8
# ggplot2 # Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4

######
#    Copyright (C) 2021  Benjamin Saintpierre - Plateforme GENOM'IC - Institut Cochin
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
