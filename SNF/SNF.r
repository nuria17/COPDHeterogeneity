
## CODE FOR COMPUTING THE SIMILARITY NETWORK FUSION APPROACH AND THE ENRICHMENT ANALYSIS FOR mRNA, miRNA, METHYLATION
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat

#load the packages
library(tidyverse)
library(dplyr)
library(xlsx)
library(ggplot2)
library(usethis)
library(devtools)
library(WGCNA)
library(openxlsx)
library(dendextend)
library(gplots)
library(RColorBrewer)
library(SNFtool)
library(Cairo)
library(purrr)
library(ChAMP)
library(GSEABase)
library(GSVAdata)
library(GSVA)
library(fgsea)
library(genefilter)
library(ROSE)
library(methylGSA)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(readxl)
library(dendextend)
library(Spectrum)

#generate required directories
dir.created("Data/SNF")

map <- purrr::map
select <- dplyr::select
rename <- dplyr::rename
mutate <- dplyr::mutate
filter <- dplyr::filter

                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                       ##                       Compute the Similarity Network Fusion algorithm, retrieve the spectral                       ##
                       ##                       clustering groups and perform enrichment analysis between the groups                         ##
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

# Before running this script, DEanalysis.r must be run first in order to have the mRNA, miRNA and Methylation 
# genes/probes x individuals matrices.
# Function to compute SNF from genes/probes x matrices, dataframe with clinical data and string to indicate file suffix.
# Returns a list with the dataframe (including the spectral clustering groups), the affinity matrix, a matrix with the 
# concordance of NMI between the group partitions and the NMI between the partition and the clinical labels.
SimilarityNetworkFusion <- function(matrixA, matrixB, matrixC, dataframe, output){
    # Standard normalization to have mean 0 and SD 1
    Data1= standardNormalization(t(matrixA))
    Data2= standardNormalization(t(matrixB))
    Data3= standardNormalization(t(matrixC))
    # Calculate distance matrices(here we calculate Euclidean Distance)
    Dist1= (dist2(as.matrix(Data1), as.matrix(Data1)))^(1/2)
    Dist2= (dist2(as.matrix(Data2), as.matrix(Data2)))^(1/2)
    Dist3= (dist2(as.matrix(Data3), as.matrix(Data3)))^(1/2)
    # Construct similarity graphs
    W1= affinityMatrix(Dist1, K=15, sigma = 0.5)
    W2= affinityMatrix(Dist2, K=15, sigma= 0.5)
    if(!is.na(Dist3)){W3= affinityMatrix(Dist3, K=15, sigma= 0.5)} else {W3 = NULL}
    if(is.null(W3)){W_list= list(W1, W2)} else {W_list= list(W1,W2,W3)}
    # Fuse all the graphs,then the overall matrix can be computed by similarity network fusion(SNF).
    W= SNF(W_list, 20, 20)
    # Save the matrix that will be used to 
    W %>% write.table(paste0('Data/SNF/Matrix_', output))
    # Estimate number of clusters using eigengap
    eigengap_W= estimate_k(W,maxk = 10, showplots = T)
    # Cluster labels for each data point by spectral clustering.
    labels = spectralClustering(W, eigengap_W)
    dataframe$SpectClust= paste0('SNF', labels)
    dataframe %>% dplyr::select(ID,SpectClust,group5) %>% write.table(paste0('Data/SNF/Table_', output))
    # Return a matrix containing the NMIs between cluster assignments made with spectral clustering on all matrices provided
    cond=concordanceNetworkNMI(W_list, C=3)
    # Compute Normalized Mutual Information (NMI)
    NMI_clinical_labels= calNMI(dataframe$SpectClust, dataframe$group5)
    
    return(list(dataframe, W, cond, NMI_clinical_labels))
}

# SNF with mRNA and miRNA datasets
mrna_inter_micro= mrna_matrix[,colnames(mrna_matrix) %in% intersect(phenoall_mrna$ID, phenoall_micro$ID)]
micro_inter_mrna= micro_matrix[,colnames(micro_matrix) %in% intersect(phenoall_mrna$ID, phenoall_micro$ID)]
phenoall_micro_mrna= phenoall[phenoall$ID %in% intersect(phenoall_mrna$ID, phenoall_micro$ID),]
fun= SimilarityNetworkFusion(matrixA = mrna_inter_micro, matrixB = micro_inter_mrna, dataframe = phenoall_micro_mrna, output= 'micro_mrna')
phenoall_micro_mrna= fun[[1]]

# SNF with methylationf and miRNA datasets
methy_inter_micro= methy_matrix[,colnames(methy_matrix) %in% intersect(phenoall_methy$ID, phenoall_micro$ID)]
micro_inter_methy= micro_matrix[,colnames(micro_matrix) %in% intersect(phenoall_methy$ID, phenoall_micro$ID)]
phenoall_methy_micro= phenoall[phenoall$ID %in% intersect(phenoall_methy$ID, phenoall_micro$ID),]
fun2= SimilarityNetworkFusion(matrixA = methy_inter_micro, matrixB = micro_inter_methy, dataframe = phenoall_methy_micro, output = 'methy_micro')
phenoall_methy_micro= fun2[[1]]

# SNF with mRNA and methylation datasets
mrna_inter_methy= mrna_matrix[,colnames(mrna_matrix) %in% intersect(phenoall_mrna$ID, phenoall_methy$ID)]
methy_inter_mrna= methy_matrix[,colnames(methy_matrix) %in% intersect(phenoall_mrna$ID, phenoall_methy$ID)]
phenoall_mrna_methy= phenoall[phenoall$ID %in% intersect(phenoall_mrna$ID, phenoall_methy$ID),]
fun3= SimilarityNetworkFusion(mrna_inter_methy, methy_inter_mrna, phenoall_mrna_methy, 'mrna_methy')
phenoall_mrna_methy=fun3[[1]]

# SNF with mRNA, miRNA and methylation datasets
mrna_inter= mrna_matrix[,colnames(mrna_matrix) %in% Reduce(intersect, list(phenoall_mrna$ID, phenoall_methy$ID, phenoall_micro$ID))]
methy_inter= methy_matrix[,colnames(methy_matrix) %in% Reduce(intersect, list(phenoall_mrna$ID, phenoall_methy$ID, phenoall_micro$ID))]
micro_inter= micro_matrix[,colnames(micro_matrix) %in% Reduce(intersect, list(phenoall_mrna$ID, phenoall_methy$ID, phenoall_micro$ID))]
phenoall_inter= phenoall[phenoall$ID %in% iReduce(intersect, list(phenoall_mrna$ID, phenoall_methy$ID, phenoall_micro$ID)),]
fun4= SimilarityNetworkFusion(mrna_inter_methy, methy_inter_mrna, phenoall_mrna_methy, 'mrna_methy')
phenoall_inter=fun4[[1]]

# Plot to show the NMI between networks

colnames(fun[[3]])= c('Fused Network', 'mRNA Network', 'miRNA Network', 'Clinical groups')
rownames(fun[[3]])= c('Fused Network', 'mRNA Network', 'miRNA Network', 'Clinical groups')
png('Data/SNF/mRNA_micro_NMI.png')
corrplot(fun[[3]],  type="lower",is.corr = F,method='number', tl.cex = 1, col=colorRampPalette(c("white","pink", "red"))(10) )
dev.off()

colnames(fun2[[3]])= c('Fused Network', 'Methylation Network', 'miRNA Network', 'Clinical groups')
rownames(fun2[[3]])= c('Fused Network', 'Methylation Network', 'miRNA Network', 'Clinical groups')
png('Data/SNF/methy_miRNA_NMI.png')
corrplot(fun2[[3]],  type="lower",is.corr = F,method='number', tl.cex = 1, col=colorRampPalette(c("white","pink", "red"))(10) )
dev.off()

colnames(fun3[[3]])= c('Fused Network', 'mRNA Network', 'Methylation Network', 'Clinical groups')
rownames(fun3[[3]])= c('Fused Network', 'mRNA Network', 'Methylation Network', 'Clinical groups')
png('Data/SNF/mRNA_methy_NMI.png')
corrplot(fun3[[3]],  type="lower",is.corr = F,method='number', tl.cex = 1, col=colorRampPalette(c("white","pink", "red"))(10) )
dev.off()

colnames(fun4[[3]])= c('Fused Network', 'mRNA Network', 'miRNA Network','Methylation Network', 'Clinical groups')
rownames(fun4[[3]])= c('Fused Network', 'mRNA Network','miRNA Network', 'Methylation Network', 'Clinical groups')
png('Data/SNF/mRNA_miRNA_methy_NMI.png')
corrplot(fun4[[3]],  type="lower",is.corr = F,method='number', tl.cex = 1, col=colorRampPalette(c("white","pink", "red"))(10) )
dev.off()

# Clinical data with the patients that had the three omics profiled
phenoall=read.xlsx('Data/phenoall.xlsx')
phenoall_inter= phenoall[!is.na(phenoall$Metilacio_teixit) & !is.na(phenoall$small_teixit) & !is.na(phenoall$mRNA_teixit),]
# Run DEanalysis.r script with the resulting SNF groups on the three layers to have the 
# table with the DE genes and p-values.
# Pathways used in the enrichment analysis came from 'The Molecular Signature Database' (MSigDB)
# that can be loaded from the ChAMP library. The same pathway database was used for mRNA, miRNA and methylation. 
data(PathwayList)
# Enrichment analysis for mRNA
# GSEA
GSEA<- function(tt){
gene_list = (tt$logFC) * (-log(tt$adj.P.Val))
names(gene_list) = rownames(tt)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

fgRes <- fgsea::fgsea(pathways = PathwayList, 
                      stats = gene_list,
                      minSize=15,
                      maxSize=600,
                      nperm=1000) %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < !!0.05)

top=fgRes[order(fgRes$padj, -abs(fgRes$NES)),]
}
# From the topTable returned in the DE analysis, GSEA function returns the table with the enriched pathways.
gsea_mrna_SNF1SNF2=GSEA(diff_mrna_SNF1SNF2)
gsea_mrna_SNF1SNF3=GSEA(diff_mrna_SNF1SNF3)
gsea_mrna_SNF2SNF3=GSEA(diff_mrna_SNF2SNF3)

# Enrichment analysis for miRNA
# Fisher's extact test mapping miRNA probes to genes through
# miRTarBase_MTI database (discarding the weak miRNA-gene interactions)
mitar= read.xlsx('RawData/miRNA/miRTarBase_MTI.xlsx')
mitar=mitar[mitar$`Species.(miRNA)` == 'Homo sapiens',]
mitar=mitar[mitar$Support.Type != "Functional MTI (Weak)",]
# Map miRNA probes of the topTable from DE analysis to gene symbol
llista_micro_snf_genes=lapply(c(list(diff_micro_SNF1SNF2$table), list(diff_micro_SNF1SNF3$table), list(diff_micro_SNF2SNF3$table)), function(x) x %>% rename(miRNA= genes) %>% left_join(mitar, by= 'miRNA'))
llista_micro_snf_genes_uniq= purrr::map(llista_micro_snf_genes, function(x) dplyr::distinct(.data=x, x$miRNA, x$Target.Gene, .keep_all= T))

llista_micro_pathwaylist_clust=list()
  for (q in 1:length(llista_micro_snf_genes)){
    Genes_meta=llista_micro_snf_genes_uniq[[q]]$Target.Gene[llista_micro_snf_genes_uniq[[q]]$FDR < 0.05] %>% unique()
       if(length(Genes_meta) > 0){
      llista_micro= purrr::map(PathwayList, function(x) { if(dim(table(llista_micro_snf_genes_uniq[[q]]$Target.Gene[llista_micro_snf_genes_uniq[[q]]$FDR < 0.05] %in% x)) > 1) {fisher.test(table(ifelse(llista_micro_snf_genes_uniq[[q]]$FDR < 0.05, 'Si', 'No'), ifelse(llista_micro_snf_genes_uniq[[q]]$Target.Gene %in% x, 'Si-P', 'No-P')), alternative = 'greater')$p.value}})
    }
    else {
      llista_micro= 0
    }
    llista_micro_pathwaylist_clust= append(llista_micro_pathwaylist_clust, list(llista_micro))
}
# Adjust p-value of the Fisher's exact test
llista_micro_pathwaylist_clust_adj= purrr::map(llista_micro_pathwaylist_clust, function(x) p.adjust(x %>% unlist(), method = 'BH',n= length(PathwayList))) 

# Enrichment analysis for Methylation
# ebGSEA method from ChAMP library
# This method doesn't need to have the topTable with the DE probes, just the normalized beta values.
# Run the methylation pipeline from DEanalysis.r script to normalize and adjust for the batch effect. 
myCombat_inter_filt=varFilter(myCombat_inter, var.func = IQR, var.cutoff = 0.25, filterByQuantile = TRUE)
fgRes_SNF1_SNF2=champ.ebGSEA(beta=myCombat_inter_filt[, colnames(myCombat_inter_filt) %in% phenoall_inter$ID[phenoall_inter$SNF %in% c('SNF1', 'SNF2')]], pheno=phenoall_inter$SNF[phenoall_inter$SNF %in% c('SNF1', 'SNF2')], minN=5, adjPval=1, arraytype="EPIC")
fgRes_SNF1_SNF3=champ.ebGSEA(beta=myCombat_inter_filt[, colnames(myCombat_inter_filt) %in% phenoall_inter$ID[phenoall_inter$SNF %in% c('SNF1', 'SNF3')]], pheno=phenoall_inter$SNF[phenoall_inter$SNF %in% c('SNF1', 'SNF3')], minN=5, adjPval=1, arraytype="EPIC")
fgRes_SNF2_SNF3=champ.ebGSEA(beta=myCombat_inter_filt[, colnames(myCombat_inter_filt) %in% phenoall_inter$ID[phenoall_inter$SNF %in% c('SNF2', 'SNF3')]], pheno=phenoall_inter$SNF[phenoall_inter$SNF %in% c('SNF2', 'SNF3')], minN=5, adjPval=1, arraytype="EPIC")

# Clinical profiling
# The function returns a matrix with the same structure as the clinical tables (Table S1-S8).
# As input, it needs the dataframe with the clinical data. 
ClinicalData_SNF <- function(phenoall){
    phenoall$Paquetes[is.na(phenoall$Paquetes) & phenoall$SpectClust == 'no.fumador']= 0
    Matrix= matrix(ncol=7, nrow= 11)
    rownames(Matrix)= c('Age,years', 'Males, %', 'BMI, kg/m', 'Packs-year', 'Years since quitting smoking', 'FEV1/FVC, %', 'FEV1 % ref. Pre-BD','FEV1 % ref. Post-BD', 'DLCO % ref.', 'CT- Emphysema(Y/N)', 'Inhaled Corticosteroid treatment')
    colnames(Matrix)= c('N', 'SNF1', 'N', 'SNF2', 'N', 'SNF3', 'p-value')
    a=0
    c=0
    llista= list('SNF1', 'SNF2', 'SNF3')
    phenoall$SpectClust= phenoall$clustering
    phenoall$SpectClust= phenoall$clustering
    for(a in 1:length(llista)){
    d= a + c
    print(a)
    Matrix[1,d]= length(phenoall[phenoall$SpectClust == llista[a],]$Edad_all)
    Matrix[1, d+1]= paste0(as.character(round(mean(phenoall[phenoall$SpectClust == llista[a],]$Edad_all),2)), '±', as.character(sd(phenoall[phenoall$SpectClust == llista[a],]$Edad_all)))  
    Matrix[2, d]= length(phenoall[phenoall$SpectClust == llista[a],]$Genero_all[phenoall[phenoall$SpectClust == llista[a],]$Genero_all == 'HOMBRE' & !is.na(phenoall[phenoall$SpectClust == llista[a],]$Genero_all)])
    Matrix[2, d+1]= length(phenoall[phenoall$SpectClust == llista[a],]$Genero_all[phenoall[phenoall$SpectClust == llista[a],]$Genero_all == 'HOMBRE' & !is.na(phenoall[phenoall$SpectClust == llista[a],]$Genero_all)])/length(phenoall[phenoall$SpectClust == llista[a],]$ID)
    Matrix[3,d]= length(phenoall[phenoall$SpectClust == llista[a],]$BMI[!is.na(phenoall[phenoall$SpectClust == llista[a],]$BMI)])
    Matrix[3,d+1]= paste0(as.character(round(mean(phenoall[phenoall$SpectClust == llista[a],]$BMI[!is.na(phenoall[phenoall$SpectClust == llista[a],]$BMI)]),2)), '±', as.character(sd(phenoall[phenoall$SpectClust == llista[a],]$BMI[!is.na(phenoall[phenoall$SpectClust == llista[a],]$BMI)])))
    Matrix[4,d]= length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$Paquetes),]$Paquetes)
    Matrix[4, d+1]= paste0(as.character(round(mean(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$Paquetes),]$Paquetes),2)), '±', as.character(sd(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$Paquetes),]$Paquetes)))
    Matrix[5,d]= length(phenoall[phenoall$SpectClust == llista[a],]$Anos_Sin_Fumar[!is.na(phenoall[phenoall$SpectClust == llista[a],]$Anos_Sin_Fumar)])
    Matrix[5, d+1]= paste0(mean(phenoall[phenoall$SpectClust == llista[a],]$Anos_Sin_Fumar[!is.na(phenoall[phenoall$SpectClust == llista[a],]$Anos_Sin_Fumar)]), '±', sd(phenoall[phenoall$SpectClust == llista[a],]$Anos_Sin_Fumar[!is.na(phenoall[phenoall$SpectClust == llista[a],]$Anos_Sin_Fumar)]))
    Matrix[6, d]= length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$FEV1_FVC),]$FEV1_FVC)
    Matrix[6, d +1]= paste0(as.character(round(mean(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$FEV1_FVC),]$FEV1_FVC), 2)),'±',as.character(sd(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$FEV1_FVC),]$FEV1_FVC)))
    Matrix[7, d]= length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$pcFEV1_Pre_BD),]$pcFEV1_Pre_BD)
    Matrix[7, d+1]= paste0(as.character(round(mean(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$pcFEV1_Pre_BD),]$pcFEV1_Pre_BD),2)),'±',as.character(sd(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$pcFEV1_Pre_BD),]$pcFEV1_Pre_BD)))
    Matrix[8, d]= length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$pcFEV1_Post_BD),]$pcFEV1_Post_BD)
    Matrix[8, d+1]= paste0(as.character(mean(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$pcFEV1_Post_BD),]$pcFEV1_Post_BD)),'±',as.character(sd(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$pcFEV1_Post_BD),]$pcFEV1_Post_BD)))
    Matrix[9, d]= length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$DLCO_pc),]$DLCO_pc)
    Matrix[9, d+1]= paste0(as.character(mean(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$DLCO_pc),]$DLCO_pc)), '±', as.character(sd(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$DLCO_pc),]$DLCO_pc)))
    Matrix[10, d]= length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$TAC_Enfisema_Revisado),]$TAC_Enfisema_Revisado)
    Matrix[10, d+1]= length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$TAC_Enfisema_Revisado) & phenoall$TAC_Enfisema_Revisado == 'SI',]$TAC_Enfisema_Revisado)/length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$TAC_Enfisema_Revisado),]$TAC_Enfisema_Revisado)
    Matrix[11, d]= length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$CorticoidesInhalados),]$CorticoidesInhalados)
    Matrix[11, d+1]= length(phenoall[phenoall$SpectClust == llista[a] & !is.na(phenoall$CorticoidesInhalados) & phenoall$CorticoidesInhalados == 'SI',]$CorticoidesInhalados)/ length(phenoall[ phenoall$SpectClust == llista[a] & !is.na(phenoall$CorticoidesInhalados),]$CorticoidesInhalados)
    c= c +1
    }   

    #ANOVA
    Matrix[1,7]= unlist(summary(aov(phenoall$Edad_all  ~ phenoall$SpectClust))[[1]][5])[1]
    Matrix[3,7]= unlist(summary(aov(phenoall$BMI ~ phenoall$SpectClust))[[1]][5])[1]
    Matrix[4,7]= unlist(summary(aov(phenoall$Paquetes ~ phenoall$SpectClust))[[1]][5])[1]
    Matrix[5,7]= unlist(summary(aov(phenoall$Anos_Sin_Fumar ~ phenoall$SpectClust))[[1]][5])[1]
    Matrix[6,7]= unlist(summary(aov(phenoall$FEV1_FVC ~ phenoall$SpectClust))[[1]][5])[1]
    Matrix[7,7]= unlist(summary(aov(phenoall$pcFEV1_Pre_BD ~ phenoall$SpectClust))[[1]][5])[1]
    Matrix[8,7]= unlist(summary(aov(phenoall$pcFEV1_Post_BD ~ phenoall$SpectClust))[[1]][5])[1]
    Matrix[9,7]= unlist(summary(aov(phenoall$DLCO_pc ~ phenoall$SpectClust))[[1]][5])[1]
    #CHI-SQUARE
    Matrix[2,7]= prop.trend.test(x= table(phenoall$SpectClust, phenoall$Genero_all)[,1], n= lapply(list(1,2,3), function(x) sum(table(phenoall$SpectClust,phenoall$Genero_all)[x,])) %>% unlist())[[3]]
    Matrix[10,7]= prop.trend.test(x= table(phenoall$SpectClust, phenoall$TAC_Enfisema_Revisado)[,1], n= lapply(list(1,2,3), function(x) sum(table(phenoall$SpectClust,phenoall$TAC_Enfisema_Revisado)[x,])) %>% unlist())[[3]] #no tenim aquesta info pels no fumadors i el test dona error si els incloc
    Matrix[11,7]= prop.trend.test(x= table(phenoall$SpectClust, phenoall$CorticoidesInhalados)[,1], n= lapply(list(1,2,3), function(x) sum(table(phenoall$SpectClust, phenoall$CorticoidesInhalados)[x,])) %>% unlist())[[3]]
return(Matrix)
}  

# Save the clinical table
Matrix_all_SNF= ClinicalData_SNF(phenoall= phenoall_inter)
Matrix_all_SNF %>% write.csv('Data/SNF/Clinical_table_all_SNF.csv')

Matrix_micro_methy= ClinicalData_SNF(phenoall = phenoall_methy_micro)
Matrix_micro_methy %>% write.csv('Data/SNF/Clinical_table_mrna_SNF.csv')
Matrix_mrna_micro= ClinicalData_SNF(phenoall= phenoall_micro_mrna)
Matrix_mrna_micro %>% write.csv('Data/SNF/Clinical_table_micro_SNF.csv')
Matrix_mrna_methy= ClinicalData_SNF(phenoall = phenoall_mrna_methy)
Matrix_mrna_methy %>% write.csv('Data/SNF/Clinical_table_methy_SNF.csv')

Matrix_mrna= ClinicalData_SNF(phenoall = phenoall_mrna)
Matrix_mrna %>% write.csv('Data/SNF/CombSNF/Clinical_table_mrna_SNF.csv' )
Matrix_micro= ClinicalData_SNF(phenoall = phenoall_micro)
Matrix_micro %>% write.csv('Data/SNF/Clinical_table_micro_SNF.csv' )
Matrix_methy= ClinicalData_SNF(phenoall= phenoall_methy)
Matrix_methy %>% write.csv('Data/SNF/Clinical_table_methy_SNF.csv' )
