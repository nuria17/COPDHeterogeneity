## CODE FOR COMPUTING A DENDROGRAM WHERE THE DISTANCES BETWEEN THE NODES ARE BASED ON THE FREQUENCY THAT THEY ARE IN THE SAME COMMUNITY
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat


### R LIBRARIES ###
library(tidyverse)
library(dplyr)
library(xlsx)
library(openxlsx)
library(gtools)
library(dendextend)
library(gplots)
library(RColorBrewer)
library(purrr)
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                       ##               Computes the distance matrix based on the frequency nodes are together when the resolution is tuned  ##
                       ##                                           and plots the dendrogram                                                 ##
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##


### FUNCTIONS ###
ClusteringParsing <- function(file, phenoall){
    phenoall$Clust=NULL
    tmp= read.table(file)
    clust=grep('ClusterID:',unlist(tmp))
    llista_clust_all=list()
    for (i in 1:(length(clust))){ if (i < (length(clust))) {llista_clust_all=append(llista_clust_all,list(tmp[(clust[i]+1):(clust[i+1]-1),]))} else {llista_clust_all=append(llista_clust_all,list(tmp[(clust[i]+1):dim(tmp)[[1]],]))}}
    for (e in 1:length(llista_clust_all)) {phenoall$Clust[phenoall$ID %in% llista_clust_all[[e]]]= paste0('Cluster', e)}
    phenoall$Clust= as.factor(phenoall$Clust)
    return(phenoall)
}
# Create ind x ind matrix to represent the frequency that the nodes are in the same community.
phenoall= read_excel('Data/Phenoall.xlsx')
files_reso= list.files('Data/MoltTi', pattern='p.*[^csv]$')
Matrix_sum= matrix(0, ncol= 164, nrow= 164)
for (i in 1:length(files_reso)){
    phenoall1= ClusteringParsing(file= paste0('/home/nuria/Desktop/Reviewers_comments/ResolutionTuning/',files_reso[[i]]), phenoall= phenoall)
    Matrix= matrix(0,ncol= length(phenoall1$ID), nrow= length(phenoall1$ID))
    colnames(Matrix)= phenoall1$ID
    rownames(Matrix)= phenoall1$ID
    for (b in 1:ncol(Matrix)){
        Matrix[,b]= ifelse(phenoall1$ID %in% phenoall1$ID[phenoall1$Clust == phenoall1$Clust[phenoall1$ID == colnames(Matrix)[[b]]]], 1, 0)
    }
    Matrix_sum= Matrix_sum + Matrix       
}

hc=hclust(dist(Matrix_sum), method= 'ward.D2')

rr = factor(phenoall1$Clust, levels= c(paste0('Cluster', 1:11)))
brew = brewer.pal(length(levels(rr)), "Paired")
cols_rr = brew2[rr]

png("Data/MolTi/Dendrogram.png", width = 1000, height = 1500)
dend= as.dendrogram(hc)
cols_rr <- cols_rr[order.dendrogram(dend)] 
dend = dend %>% 
        set("labels_colors", cols_rr) %>% 
        set("branches_lwd", 1) %>% set("labels_cex",1) %>%
        plot(horiz= TRUE, axes=FALSE)
legend("topleft", legend = levels(rr), fill = brew, cex = 2)
dev.off()

