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

confFile <- "PROJET.conf"
contrastFile <- "PROJET.contrastes"
projectName <- "PROJET_reference"

make_nice_boxplot <- function(dds,configuration,projectName) {
  
  countNorm<-counts(dds, normalized=TRUE)
  pseudoNormCount <- as.data.frame(log2(countNorm+1),row.names=row.names(countNorm))
  pseudoNormCount$Ids <- row.names(pseudoNormCount)
  datNorm <- melt(pseudoNormCount, id.vars = "Ids", variable.name = "Samples", value.name = "count")
  colnames(datNorm)<-c("Ids", "Samples", "count")
  datNorm <- merge(datNorm, configuration, by.x = "Samples", by.y = "Name")
  
  countRaw<-counts(dds, normalized=FALSE)
  pseudoRawCount <- as.data.frame(log2(countRaw+1),row.names=row.names(countRaw))
  pseudoRawCount$Ids <- row.names(pseudoRawCount)
  datRaw <- melt(pseudoRawCount, id.vars = "Ids", variable.name = "Samples", value.name = "count")
  colnames(datRaw)<-c("Ids", "Samples", "count")
  datRaw <- merge(datRaw, configuration, by.x = "Samples", by.y = "Name")
  
  gDens <- ggplot(datNorm, aes(x = count, colour = Samples, fill = Samples)) + geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) + theme(legend.position = "top") + xlab(expression(log[2](count + 1)))
  ggsave(filename=paste(projectName,"_DensityNormCount.jpg",sep=""), plot=gDens)
  
  gNorm <- ggplot(data = datNorm, aes(x = Samples, y = count, fill = Condition)) + geom_boxplot() + xlab("Samples") + ylab("log2(Count+1)") + ggtitle("Normalized log2(counts) per sample") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename="boxplotNormCount.jpg", plot=gNorm)
  
  gRaw <- ggplot(data = datRaw, aes(x = Samples, y = count, fill = Condition)) + geom_boxplot() + xlab("Samples") + ylab("log2(Count+1)") + ggtitle("Raw log2(counts) per sample") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename="boxplotRawCount.jpg", plot=gRaw)
  
}

make_nice_clusters <- function(rld,configuration,projectName) {
  distsRL <- dist(t(assay(rld)))
  matDist <- as.matrix(distsRL)
  hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
  Conditions <- data.frame(configuration$Condition,row.names=configuration$Name)
  colnames(Conditions) <- c("condition")
  nameClustering <- paste(projectName,"deseq2_Unsupervised_clustering_euclidean-complete.png",sep="_")
  png(filename=nameClustering,width=7 ,height=7, units="in",res = 600 ) 
  pheatmap(matDist, 
           col=hmcol, 
           annotation_col = Conditions,
           show_rownames=F,
           show_colnames=T)
  dev.off()
  nameClustering <- paste(projectName,"deseq2_Unsupervised_clustering_euclidean-complete_no-sampleNames.png",sep="_")
  png(filename=nameClustering,width=7 ,height=7, units="in",res = 600 ) 
  pheatmap(matDist, 
           col=hmcol, 
           annotation_col = Conditions,
           show_rownames=F,
           show_colnames=F)
  dev.off()
  
  png(filename="clustering.png",width=7 ,height=7, units="in",res = 600 ) 
  pheatmap(matDist, 
           col=hmcol, 
           annotation_col = Conditions,
           show_rownames=F,
           show_colnames=T)
  dev.off()
}

make_nice_pca <- function(dds,rld,projectName) {
  
  ntop <- 500
  
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t( assay(rld)[select, ] )
  pc <- prcomp(mat)
  eig <- (pc$sdev)^2
  variance <- eig*100/sum(eig)
  
  PCAdata<-as.data.frame(pc$x[,1:3])
  PCAdata$condition <- dds$Condition
  
  nameACP1 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC2_top500varGenes.jpg",sep="")
  pca1 <- ggplot(PCAdata,aes(x=PC1,y=PC2,col=condition,label=rownames(PCAdata))) + geom_point(aes(shape=condition, color=condition), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC2: ",round(variance[2],1),"% variance"))
  ggsave(filename=nameACP1, plot=pca1)
  ggsave(filename="PCA1.jpg", plot=pca1)
  
  nameACP2 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC2vsPC3_top500varGenes.jpg",sep="")
  pca2 <- ggplot(PCAdata,aes(x=PC2,y=PC3,col=condition,label=rownames(PCAdata))) + geom_point(aes(shape=condition, color=condition), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC2: ",round(variance[2],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP2, plot=pca2)
  ggsave(filename="PCA2.jpg", plot=pca2)
  
  nameACP3 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC3_top500varGenes.jpg",sep="")
  pca3 <- ggplot(PCAdata,aes(x=PC1,y=PC3,col=condition,label=rownames(PCAdata))) + geom_point(aes(shape=condition, color=condition), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP3, plot=pca3)
  
  nameACP_contrib <- paste(projectName,"_deseq2_Unsupervised_PCA_contribution.png",sep="")
  contrib <- fviz_eig(pc)
  ggsave(filename=nameACP_contrib, plot=contrib)
  
}

make_nice_pca_3factors <- function(dds,rld,projectName) {
  
  ntop <- 500
  
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t( assay(rld)[select, ] )
  pc <- prcomp(mat)
  eig <- (pc$sdev)^2
  variance <- eig*100/sum(eig)
  
  # Condition 1
  
  PCAdata<-as.data.frame(pc$x[,1:3])
  PCAdata$Patient <- dds$Patient
  
  nameACP1 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC2_top500varGenes_Patient.jpg",sep="")
  pca1 <- ggplot(PCAdata,aes(x=PC1,y=PC2,col=Patient,label=rownames(PCAdata))) + geom_point(aes(shape=Patient, color=Patient), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC2: ",round(variance[2],1),"% variance"))
  ggsave(filename=nameACP1, plot=pca1)
  nameACP2 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC2vsPC3_top500varGenes_Patient.jpg",sep="")
  pca2 <- ggplot(PCAdata,aes(x=PC2,y=PC3,col=Patient,label=rownames(PCAdata))) + geom_point(aes(shape=Patient, color=Patient), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC2: ",round(variance[2],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP2, plot=pca2)
  nameACP3 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC3_top500varGenes_Patient.jpg",sep="")
  pca3 <- ggplot(PCAdata,aes(x=PC1,y=PC3,col=Patient,label=rownames(PCAdata))) + geom_point(aes(shape=Patient, color=Patient), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP3, plot=pca3)
  
  # Condition 2
  
  PCAdata<-as.data.frame(pc$x[,1:3])
  PCAdata$Traitement <- dds$Traitement
  
  nameACP1 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC2_top500varGenes_Traitement.jpg",sep="")
  pca1 <- ggplot(PCAdata,aes(x=PC1,y=PC2,col=Traitement,label=rownames(PCAdata))) + geom_point(aes(shape=Traitement, color=Traitement), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC2: ",round(variance[2],1),"% variance"))
  ggsave(filename=nameACP1, plot=pca1)
  nameACP2 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC2vsPC3_top500varGenes_Traitement.jpg",sep="")
  pca2 <- ggplot(PCAdata,aes(x=PC2,y=PC3,col=Traitement,label=rownames(PCAdata))) + geom_point(aes(shape=Traitement, color=Traitement), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC2: ",round(variance[2],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP2, plot=pca2)
  nameACP3 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC3_top500varGenes_Traitement.jpg",sep="")
  pca3 <- ggplot(PCAdata,aes(x=PC1,y=PC3,col=Traitement,label=rownames(PCAdata))) + geom_point(aes(shape=Traitement, color=Traitement), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP3, plot=pca3)
  
  # Condition 3
  
  PCAdata<-as.data.frame(pc$x[,1:3])
  PCAdata$Temps <- dds$Temps
  
  nameACP1 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC2_top500varGenes_Temps.jpg",sep="")
  pca1 <- ggplot(PCAdata,aes(x=PC1,y=PC2,col=Temps,label=rownames(PCAdata))) + geom_point(aes(shape=Temps, color=Temps), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC2: ",round(variance[2],1),"% variance"))
  ggsave(filename=nameACP1, plot=pca1)
  nameACP2 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC2vsPC3_top500varGenes_Temps.jpg",sep="")
  pca2 <- ggplot(PCAdata,aes(x=PC2,y=PC3,col=Temps,label=rownames(PCAdata))) + geom_point(aes(shape=Temps, color=Temps), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC2: ",round(variance[2],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP2, plot=pca2)
  nameACP3 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC3_top500varGenes_Temps.jpg",sep="")
  pca3 <- ggplot(PCAdata,aes(x=PC1,y=PC3,col=Temps,label=rownames(PCAdata))) + geom_point(aes(shape=Temps, color=Temps), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP3, plot=pca3)
  
}

make_nice_volcanoPlot <- function(fullData,name,ctrst) {
  
  
  newdata<-fullData[,c(1,3,7,9)]
  
  volcanoTitle <- paste(ctrst[2],'vs',ctrst[3],sep=" ")
  
  colnames(newdata)<-c("Id","log2FoldChange","padj","GeneSymbol")
  
  # add a column of NAs
  newdata$DEGs <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  newdata$DEGs[newdata$log2FoldChange > log2(1.5)] <- "UP"
  newdata$DEGs[newdata$log2FoldChange > log2(1.5) & newdata$padj < 0.05] <- "UP & DE"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  newdata$DEGs[newdata$log2FoldChange < -log2(1.5)] <- "DOWN"
  newdata$DEGs[newdata$log2FoldChange < -log2(1.5) & newdata$padj < 0.05] <- "DOWN & DE"

  nameVolcanoPlot<-paste(name,"_volcanoPlot.png",sep="")
  
  volcanoNoLabel<-ggplot(newdata) +
    aes(x = log2FoldChange, y = -log10(padj), colour = DEGs) +
    geom_point(shape = "circle", size = 1.5) +
    xlab("log2(Fold-Change") +
    ylab("-log10(padj)") +
    labs(title = volcanoTitle) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="grey50",linetype=3) + 
    geom_hline(yintercept=-log10(0.05), col="grey50",linetype=3) +
    scale_color_manual(values=c("#CCCCFF","blue3", "snow3", "#FFCCCC", "firebrick2"))
  
  ggsave(filename=nameVolcanoPlot, plot=volcanoNoLabel)
  
  top10DEG<-newdata[1:10,]
  
  nameVolcanoPlot<-paste(name,"_volcanoPlot_Top10.png",sep="")
  
  volcanoTop10<-ggplot(newdata) +
    aes(x = log2FoldChange, y = -log10(padj), colour = DEGs) +
    geom_point(shape = "circle", size = 1.5) +
    xlab("log2(Fold-Change") +
    ylab("-log10(padj)") +
    labs(title = volcanoTitle ,subtitle = "Top 10 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="grey50",linetype=3) + 
    geom_hline(yintercept=-log10(0.05), col="grey50",linetype=3) +
    scale_color_manual(values=c("#CCCCFF","blue3", "snow3", "#FFCCCC", "firebrick2"))+
    geom_label_repel(data=top10DEG, aes(label=GeneSymbol), color = 'black',show.legend=FALSE)
  ggsave(filename=nameVolcanoPlot, plot=volcanoTop10)
  
  
}

make_nice_diffPlot <- function(dds,fullData,configuration,name,ctrst) {
  
  diffTitle <- paste(ctrst[2],'vs',ctrst[3],sep=" ")
  
  countNorm<-counts(dds, normalized=TRUE)
  pseudoNormCount <- as.data.frame(log2(countNorm+1),row.names=row.names(countNorm))
  pseudoNormCount$Ids <- row.names(pseudoNormCount)
  
  ListTop_10 <- fullData[1:10,1]
  
  Top10DEGs <- subset(countNorm, row.names(countNorm) %in% ListTop_10)
  Top10DEGs <- melt(Top10DEGs)
  colnames(Top10DEGs) <- c("EnsemblID","Sample","NormalizedCount")
  
  top_10 <- merge(Top10DEGs, configuration, by.x = "Sample",by.y="Name")
  
  top_10 <- subset(top_10, Condition %in% ctrst)
  
  diffplot <- ggplot(top_10) +
    geom_point(aes(x = EnsemblID, y = NormalizedCount, color = Condition, shape=Condition)) +
    scale_y_log10() +
    xlab("Genes") +
    ylab("Normalized Counts") +
    labs(title = diffTitle ,subtitle = "Top 10 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  nameDiffPlot<-paste(name,"_diffPlot_Top10.png",sep="")
  ggsave(filename=nameDiffPlot, plot=diffplot)
}




# Step 1: Preparing the configuration table with conditions and factors to take in account

  currentDirectory <- getwd()
  configuration <- read.table(confFile,header=TRUE,sep="\t",colClasses = c("character","character","factor")) # Change the colClasses parameter to fit with the project design
  files <- file.path(currentDirectory, configuration$File)
  names(files) <- configuration$Name
  
  
  ensembl <- useEnsembl(biomart = "ensembl", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "uswest")

# Step 2: Import the data

  txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
  txi.rsem$length <- ifelse(txi.rsem$length==0,1,txi.rsem$length) # Setting the length to 1, to avoid error for the dds object

# Step 3: Creating the DDS object

  dds <- DESeqDataSetFromTximport(txi.rsem, configuration, ~ Condition)
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 3 # filter out genes where there are less than 3 samples with normalized counts greater than or equal to 10.
  dds <- dds[idx,]
  dds <- DESeq(dds) # no need to set the Reference group, as the script use the contrast parameter later
  
# Step 4: Save the normalized matrix
  countNorm<-counts(dds, normalized=TRUE)
  ensemblIDs <- row.names(countNorm)
  
  #Human
  IDsWithNamesDesc <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol','description'), filters = 'ensembl_gene_id', values = ensemblIDs, mart = ensembl,useCache = FALSE)
  rownames(IDsWithNamesDesc) <- make.names(IDsWithNamesDesc$ensembl_gene_id, unique=TRUE)
  dataMerged <- merge(as.data.frame(countNorm),IDsWithNamesDesc,by="row.names",all.x=TRUE)
  nameNorm <- paste(projectName,"_deseq2_NormalizedMatrix.tsv",sep="")
  write.table(dataMerged,file=nameNorm,sep='\t',row.names=F)
  
  
# Step 7: Unsupervised analysis
  
  make_nice_boxplot(dds,configuration,projectName)
  
  rld<- rlog(dds,blind=TRUE)
  
  make_nice_clusters(rld,configuration,projectName)

  make_nice_pca(dds,rld,projectName)
  
  make_nice_pca_3factors(dds,rld,projectName)
  
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
  
library(rmarkdown)
  
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
