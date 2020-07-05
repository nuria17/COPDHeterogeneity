
## CODE FOR COMPUTING THE DE, ENRICHMENT AND CLINICAL ANALYSIS BETWEEN ALL THE MOLTI COMMUNITIES. 
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat

#load the packages
library(oligo)
  library(oligoData)
  library(tidyverse)
  library(simpleaffy)
  library(xlsx)
  library(ggplot2)
  library(limma)
  library(ggfortify)
  library(biomaRt)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(hgu219.db)
  library(GO.db)
  library(usethis)
  library(devtools)
  library(WGCNA)
  library(edgeR)
  library(openxlsx)
  library(sva)
  library(gtools)
  library(VennDiagram)
  library(dendextend)
  library(gplots)
  library(RColorBrewer)
  library(SNFtool)
  library(mixOmics)
  library(Cairo)
  library(purrr)
  library(igraph)
  library(GOstats)
  library(org.Hs.eg.db)
  library(DOSE)
  library(pathview)
  library(clusterProfiler)
  library(AnnotationHub)
  library(ensembldb)
  library(splines)
  library(devtools)
  library(ChAMP)
  library(GSEABase)
  library(GSVAdata)
  library(GSVA)
  library(ROSE)
  library(pROC)
  library(igraph)
  library(GeneAnswers)
  library(methylGSA)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(goeveg)
  library(disparityfilter)
  library(RNentropy)
  library(ggrepel)
  library(Category)
  library(SPIA)
  library(jaccard)

#generate required directories
dir.created("Data/MolTi")

map <- purrr::map
select <- dplyr::select
rename <- dplyr::rename
mutate <- dplyr::mutate
filter <- dplyr::filter
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                       ##                      Perform the DE, enrichment and clinical analysis for mRNA, miRNA and methylation              ##
                       ##                       for all the comparisons between the communities.                                             ##
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

# Before running this script, DEanalysis.r must be run to perform the normalization and the batch effect adjustment for each layer.
# The normalized matrices and the clinical data need to be loaded for the further analysis.

## mRNA ##
# Function to get in a dataframe the ID and the cluster group they belong in the communities detected by MolTi.
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
# Filtering out the least variables probes.
mrna= varFilter(data_rma_collapse, var.func = IQR, var.cutoff = 0.25, filterByQuantile = T)
# Pathways from the Molecular Signature Database
data(PathwayList)
# Function to perform the DE analysis and GSEA for all the possible comparisons between the communities detected by MolTi.
GSEA_mRNA <- function(phenoall_mrna){
  DE_mrna=list()
  Genes=list()
  Path= list()
  phenoall_mrna$Clust= phenoall_mrna$Clust %>% as.factor()
  for(i in 1:length(levels(phenoall_mrna$Clust))){
  fact= phenoall_mrna$Clust %>% as.factor()
  design <- model.matrix(~0 + fact)
  #Estimating surrogate variables using SVA algorithm
  sv <- sva(mrna,design)
  design <- cbind(design, sv$sv)
  colnames(design)[12:length(colnames(design))] = paste0('SV_', 1:sv$n.sv) 
  fit= lmFit(mrna, design)
  arg= c(unique((lapply(phenoall_mrna$Clust[phenoall_mrna$Clust != levels(phenoall_mrna$Clust)[[i]]], function(x) if (x != levels(phenoall_mrna$Clust)[[i]]){paste0('fact',levels(phenoall_mrna$Clust)[[i]],'-fact',x)}))), list(levels=design))
  contrast.matrix <- do.call(makeContrasts, arg)
  data.fit.con = contrasts.fit(fit,contrast.matrix)
  data.fit.eb <- eBayes(data.fit.con)
  tt= purrr::map(c(1:(length(levels(phenoall_mrna$Clust))-1)), function(x) topTable(data.fit.eb, coef=x, n=Inf))
  names(tt)= colnames(topTable(data.fit.eb))[c(1:(length(levels(phenoall_mrna$Clust))-1))]
  pathway= purrr::map(tt, function(x){gene_list= (x$logFC * (-log10(x$adj.P.Val)));names(gene_list) = rownames(x);gene_list = sort(gene_list, decreasing = T);
  gene_list = gene_list[!duplicated(names(gene_list))]; fgRes <- fgsea::fgsea(pathways = GO,stats = gene_list,minSize=15,maxSize=600,nperm=1000) %>% as.data.frame() %>% dplyr::filter(padj < !!0.05); fgRes[order(fgRes$padj, -abs(fgRes$NES)),]})
  Path= append(Path, list(pathway))
  names(Path)[[i]]= levels(phenoall_mrna$Clust)[[i]]
  
  }
  return(Path)
}
phenoall_mrna_molti= ClusteringParsing(file = 'Data/Molti/Community_file', phenoall = phenoall_mrna)
de_mrna_mult_path= GSEA_mRNA(phenoall_mrna = phenoall_mrna_molti)

## Methylation ##
# Function to perform the DE analysis with the methylation dataset for all the possible 
# comparisions between the communities detected by MolTi software.
# loading data to map CpG probes with genes.
data("probe.features.epic")
myCombat_molti_filt=varFilter(myCombat, var.func = IQR, var.cutoff = 0.25, filterByQuantile = TRUE)
DE_methylation <- function(phenoall_methy) {
    phenoall_methy$Clust= as.factor(phenoall_methy$Clust)
    gometh=list()
    gometh_kegg= list()
    RRA2= list()
    Table_methy= list()
    fgRes_methy= list()
    for (c in 1:length(levels(phenoall_methy$Clust))){
        fact= phenoall_methy$Clust %>% as.factor()
        design <- model.matrix(~0+ fact)
        fit <- lmFit(myCombat_molti_filt, design)
        arg= c(unique((lapply(phenoall_methy$Clust[phenoall_methy$Clust != levels(phenoall_methy$Clust)[[i]]], function(x) if (x != levels(phenoall_methy$Clust)[[i]]){paste0('fact',levels(phenoall_methy$Clust)[[i]],'-fact',x)}))), list(levels=design))
        contrast.matrix <- do.call(makeContrasts, arg)
        data.fit.con = contrasts.fit(fit,contrast.matrix)
        data.fit.eb <- eBayes(data.fit.con)
        tt= purrr::map(c(1:10), function(x) topTable(data.fit.eb, coef=x, n=Inf))
        names(tt)= colnames(topTable(data.fit.eb))[c(1:10)]
        Table_methy= append(Table_methy, list(tt))          
    }
    return(Table_methy)
}
phenoall_methy_molti= ClusteringParsing(file = 'Data/Molti/Community_file', phenoall = phenoall_methy)
DE_methy= DE_methylation(phenoall_methy = phenoall_methy_molti)
# Pathway enrichment analysis for methylation using ebGSEA
ebGSEA <- function(phenoall_methy) {
 phenoall_methy$Clust= as.factor(phenoall_methy$Clust)
 fgRes_methy= list()
 d=0
 llista= list()
 for (i in 1:length(levels(phenoall_methy$Clust))){
    for (c in 1:length(levels(phenoall_methy$Clust))){
        if (levels(phenoall_methy$Clust)[[i]] != levels(phenoall_methy$Clust)[[c]] & (!paste0(levels(phenoall_methy$Clust)[[c]],'-',levels(phenoall_methy$Clust)[[i]]) %in% llista)){
            d= d + 1
            fact= phenoall_methy$Clust[phenoall_methy$Clust %in% c(levels(phenoall_methy$Clust)[[i]], levels(phenoall_methy$Clust)[[c]])]
            fgRes=champ.ebGSEA(beta=myCombat_molti_filt[,colnames(myCombat_molti_filt) %in% phenoall_methy$ID[phenoall_methy$Clust %in% c(levels(phenoall_methy$Clust)[[i]], levels(phenoall_methy$Clust)[[c]])]], pheno=as.character(fact), minN=5, adjPval=0.05, arraytype="EPIC")
            fgRes_methy= append(fgRes_methy, list(fgRes))
            names(fgRes_methy)[[d]]= paste0(levels(phenoall_methy$Clust)[[c]], '-', levels(phenoall_methy$Clust)[[i]])
            print(paste0(levels(phenoall_methy$Clust)[[c]], '-', levels(phenoall_methy$Clust)[[i]]))	
            llista= append(llista, c(paste0(levels(phenoall_methy$Clust)[[i]],'-', levels(phenoall_methy$Clust)[[c]])))  
        }
   
    } 
 }
return(fgRes_methy)
}
DE_methy_myCombat025_all_mult_gene= ebGSEA(phenoall_methy = phenoall_methy_molti)
DE_methy_myCombat025_all_mult_gene=purrr::map(DE_methy_myCombat025_all_mult, function(y) {purrr::map(y, function(x) x[x$adj.P.Val < 0.05,])})
DE_methy_myCombat025_all_mult_gene= purrr::map(DE_methy_myCombat025_all_mult_gene, function(y) {purrr::map(y, function(x) x %>% rownames_to_column(var='probe') %>% left_join(probe.features %>% rownames_to_column(var= 'probe'), by= 'probe'))})

## miRNA ##
# Function to perform the DE analysis with the miRNA dataset for all the comparisons between the communities detected with MolTi software.
# Database to perform the miRNA probe mapping to gene symbol.
mitar= read.xlsx('Data/miRNA/miRTarBase_MTI.xlsx')
mitar=mitar[mitar$`Species.(miRNA)` == 'Homo sapiens',]
mitar=mitar[mitar$Support.Type != "Functional MTI (Weak)",]

DE_mirna <- function(phenoall_micro){
    dge <- DGEList(counts= micro_counts %>% dplyr::select(-micro) %>% as.matrix(), group= phenoall_micro$Clust, genes= micro_counts$micro)
    dge <- calcNormFactors(dge, method = 'TMM')
    dge$logCPM <- cpm(dge, log= TRUE, prior.count = 0.25, normalized.lib.sizes = TRUE)
    cpmcutoff <- round(10/min(dge$sample$lib.size/1e+06), digits = 1)
    nsamplescutoff <- min(table(dge$samples$group))
    mask <- rowSums(cpm(dge,  log= TRUE, prior.count = 0.25, normalized.lib.sizes = TRUE) > cpmcutoff) >= nsamplescutoff
    dge <- dge[mask,]
    dge$logCPM <- dge$logCPM[mask,]
    phenoall_micro$Clust= factor(phenoall_micro$Clust)
    rownames(dge$logCPM)= dge$genes$genes
    top_micro=list()
    DE_micro=list()
    for (i in 1:length(levels(phenoall_micro$Clust))){
        fact= phenoall_micro$Clust
        design= model.matrix(~ 0 + fact)
        sva <- sva(dge$logCPM, design)
        design <- cbind(design,sva$sv)
        colnames(design)[11:length(colnames(design))]= paste0('SV_', 1:sva$n.sv)
        fit <- lmFit(dge$logCPM, design)
        arg= c(unique((lapply(phenoall_micro$Clust[phenoall_micro$Clust != levels(phenoall_micro$Clust)[[i]]], function(x) if (x != levels(phenoall_micro$Clust)[[i]]){paste0('fact',levels(phenoall_micro$Clust)[[i]],'-fact',x)}))), list(levels=design))
        contrast.matrix <- do.call(makeContrasts, arg)
        data.fit.con = contrasts.fit(fit,contrast.matrix)
        data.fit.eb <- eBayes(data.fit.con)
        DEresults = decideTests(data.fit.eb,method='global',adjust.method="BH",p.value=0.05,lfc=1)
        tt= purrr::map(c(1:(length(levels(phenoall_micro$Clust))-1)), function(x) topTable(data.fit.eb, coef=x, n=Inf))
        names(tt)= colnames(topTable(data.fit.eb))[c(1:(length(levels(phenoall_micro$Clust))-1))]
        tt2= purrr::map(tt, function(x) x %>% rownames_to_column(var='miRNA') %>% left_join(mitar %>% dplyr::select(miRNA, Target.Gene, `Target.Gene.(Entrez.ID)`), by= 'miRNA')) 
        DE_micro= append(DE_micro, list(tt2))
        names(DE_micro)[[i]]= levels(phenoall_micro$Clust)[[i]]
    }
    return(DE_micro)
}
phenoall_micro_molti= ClusteringParsing(file = phenoall_micro)
micro_all_mult= DE_mirna(phenoall_micro = phenoall_micro_molti)
micro_all_mult=purrr::map(micro_all_mult, function(y) {purrr::map(y, function(x) dplyr::distinct(.data=x, x$miRNA, x$Target.Gene, .keep_all= T))})
# Fisher's exact test with the genes that have an FDR < 0.05 and the pathways.
llista_micro_pathwaylist_clust_mult=list()
for (q in 1:length(levels(phenoall_micro_molti$Clust))){
    llista_micro_GO_mult=purrr::map(micro_all_mult[[q]], function(y){ if (length(y$Target.Gene[y$adj.P.Val < 0.05]) > 0) {purrr::map(PathwayList, function(x) { if(dim(table(y$Target.Gene[y$adj.P.Val < 0.05] %in% x)) > 1) {fisher.test(table(ifelse(y$adj.P.Val < 0.05, 'Si', 'No'), ifelse(y$Target.Gene %in% x, 'Si-P', 'No-P')), alternative = 'greater')$p.value}})}})
    llista_micro_pathwaylist_clust_mult= append(llista_micro_pathwaylist_clust_mult, list(llista_micro_GO_mult))
    names(llista_micro_pathwaylist_clust_mult)[[q]]=levels(phenoall$Clust)[[q]]
}
llista_micro_pathwaylist_clust_mult_adj= purrr::map(llista_micro_pathwaylist_clust_mult, function(y) {purrr::map(y,function(x) p.adjust(x %>% unlist(),method = 'BH', n = length(PathwayList)))})

# Function to keep the pathways for every community that were deregulated in more than one comparison
Paths_more_than_one_comp <- function(omics_pathways, input_type){
omic_path_more=list()
omic_path_temp=list()
for (i in 1:length(omics_pathways)){
  for (a in 1:length(omics_pathways[[i]])){
    all= omics_pathways[[i]][names(omics_pathways[[i]]) != names(omics_pathways[[i]])[[a]]]
    all_path= ifelse(input_type == 'mRNA', purrr::map(all, function(x) x$pathway) %>% unlist(), ifelse(input_type == 'miRNA',purrr::map(all, function(x) names(x)) %>% unlist(),purrr::map(all, function(x) x$`Rank(P)` %>% rownames()) %>% unlist()))
    temp= ifelse(input_type == 'mRNA', omics_pathways[[i]][[a]]$pval)[(omics_pathways[[i]][[a]]$pathway) %in% all_path & abs(omics_pathways[[i]][[a]]$NES) > 2],
    ifelse(input_type == 'miRNA',omics_pathways[[i]][[a]][names(omics_pathways)[[i]][[a]] %in% all_path],(omics_pathways[[i]][[a]]$`Rank(P)`)[rownames(omics_pathways[[i]][[a]]$`Rank(P)`) %in% all_path,]))
    
    omic_path_temp= append(omic_path_temp, list(temp))
    names(omic_path_temp)[[a]]= names(omic_path_full[[i]])[[a]]
  }
  omic_path_more= append(omic_path_more, list(omic_path_temp))
  names(omic_path_more)[[i]]= names(omic_path_full)[[i]]
  omic_path_temp=list()
}
return(omic_path_more)
}
methy_path_more= Paths_more_than_one_comp(DE_methy_myCombat025_all_mult_gene, input_type = 'Methy')

mrna_path_more=Paths_more_than_one_comp(de_mrna_mult_path, input_type = 'mRNA')

micro_path_more=Paths_more_than_one_comp(llista_micro_pathwaylist_clust_mult_adj, input_type = 'miRNA')

# Functions to obtain the pathways that were solely found in a particular community
Uniq_mrna<- function(mrna_path_more_NES_2){
  mrna_path_more_NES_2_uniq=list()
  for (i in 1:length(mrna_path_more_NES_2)){
    All= purrr::map(mrna_path_more_NES_2[names(mrna_path_more_NES_2) != names(mrna_path_more_NES_2)[[i]]], function(y) purrr::map(y, function(x) x)) %>% unlist() %>% unique()
    mrna_path_more_NES_2_uniq= append(mrna_path_more_NES_2_uniq, list(purrr::map(mrna_path_more_NES_2[[i]], function(x) x[!x$pathway %in% All,])))
  }
  return(mrna_path_more_NES_2_uniq)
}
mrna_path_more_uniq= Uniq_mrna(mrna_path_more)

Uniq_methy<- function(methy_path_more_clustmethy){
    methy_path_more_clustmethy_uniq=list()
  for (i in 1:length(methy_path_more_clustmethy)){
    All= purrr::map(methy_path_more_clustmethy[names(methy_path_more_clustmethy) != names(methy_path_more_clustmethy)[[i]]], function(y) purrr::map(y, function(x) rownames(x))) %>% unlist() %>% unique()
    methy_path_more_clustmethy_uniq= append(methy_path_more_clustmethy_uniq, list(purrr::map(methy_path_more_clustmethy[[i]], function(x) as.data.frame(x)[!rownames(x) %in% All,])))
  }
  return(methy_path_more_clustmethy_uniq)
}
methy_path_more_uniq= Uniq_methy(methy_path_more_clustmethy = methy_path_more)

Uniq_micro<- function(micro_path_more){
    micro_path_more_uniq=list()
  for (i in 1:length(micro_path_more)){
    All= purrr::map(micro_path_more[names(micro_path_more) != names(micro_path_more)[[i]]], function(y) purrr::map(y, function(x) names(x))) %>% unlist() %>% unique()
    micro_path_more_uniq= append(micro_path_more_uniq, list(purrr::map(micro_path_more[[i]], function(x) x[!names(x) %in% All])))
  }
  return(micro_path_more_uniq)
}
micro_path_more_uniq= Uniq_micro(micro_path_more)

# Function to compute Jaccard index between the communities
# JACCARD INDEX
jaccard_mat <- function(omics_path, phenoall){
    jaccard_matrix=matrix(nrow=length(omics_path), ncol= length(omics_path))
    for (i in 1:length(omics_path)){
        for (j in 1:length(omics_path)){
            if (i != j){
                setA= ifelse(names(PathwayList) %in% (omics_path[[i]][(names(omics_path[[i]]) != paste0('fact',names(omics_path)[[i]],'.fact', names(omics_path)[[j]])) & (names(omics_path[[i]]) != paste0('fact',names(omics_path)[[j]],'.fact', names(omics_path)[[i]]))] %>% unlist() %>% unique), 1, 0)
                setB= ifelse(names(PathwayList) %in% (omics_path[[j]][(names(omics_path[[i]]) != paste0('fact',names(omics_path)[[i]],'.fact', names(omics_path)[[j]])) & (names(omics_path[[i]]) != paste0('fact',names(omics_path)[[j]],'.fact', names(omics_path)[[i]]))] %>% unlist() %>% unique), 1, 0)
                jaccard_matrix[i,j]= jaccard(setA,setB)
            } else {
                jaccard_matrix[i,j]=1
            }
        }
    }
    rownames(jaccard_matrix)= levels(phenoall$Clust)
    colnames(jaccard_matrix)=levels(phenoall$Clust)
    return(jaccard_matrix)
}

jaccard_mrna_out= jaccard_mat(omics_path = mrna_path_more, phenoall = phenoall_mrna_molti)
jaccard_mrna_out=jaccard_mrna_out[c(1,4,5,6,7,8,9,10,11,2,3),c(1,4,5,6,7,8,9,10,11,2,3)]

jaccard_methy_out=jaccard_mat(methy_path_more, phenoall_methy_molti)
jaccard_methy_out=jaccard_methy_out[c(1,4,5,6,7,8,9,10,11,2,3),c(1,4,5,6,7,8,9,10,11,2,3)]

jaccard_micro_out=jaccard_mat(micro_path_more, phenoall_micro_molti)
jaccard_micro_out=jaccard_micro_out[c(1,4,5,6,7,8,9,10,11,2,3),c(1,4,5,6,7,8,9,10,11,2,3)]

# Plotting the corrplot with the Jaccard Index
pdf('Data/MolTi/Jaccard',height = 20, width = 15)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
corrplot(jaccard_mrna_out,  type="lower",is.corr = F,method='shade', tl.cex = 2.5 )
dev.off()

## Heatmap of the adjusted p-values ##
# Function to select the pathways the would be represented in the heatmap
Llista_pathways <- function(mrna_path_more_NES_2_names_, mrna_path_more_NES_2_uniq){
llista=list()
for (i in 1:length(mrna_path_more_NES_2_names_)){
  llista= c(llista, as.list(unique(unlist(mrna_path_more_NES_2_uniq[[i]]))[c(1:5)]))
  for (j in 1:length(mrna_path_more_NES_2_names_)){
    if (i != j){
      bla= c(unlist(mrna_path_more_NES_2_names_[[i]])[unlist(mrna_path_more_NES_2_names_[[i]]) %in% unlist(mrna_path_more_NES_2_names_[[j]])][c(1:5)])
      llista=c(llista,as.list(bla))
    }
  }
}
llista2= unlist(llista) %>% unique()
return(llista2)
}
llista_methy_path=Llista_pathways(methy_path_more, methy_path_more_uniq)
llista_micro_path=Llista_pathways(micro_path_more, micro_path_more_uniq)
llista_mrna_path= Llista_pathways(mrna_path_more, mrna_path_more_uniq)

# Function that creates a matrix with the pathways, communities and the corresponding adjusted p-values.
Clustering_matrix<- function(omics_path_more,omic_path_morepval,matrix){
  for(i in 1:length(omics_path_more)){
    if (length(unlist(omics_path_morepval[[i]])[unlist(omics_path_more[[i]]) %in% rownames(matrix)]) != 0){
      matrix[rownames(matrix)[which(rownames(matrix) %in% unlist(omics_path_more[[i]]))],i]= unlist(omics_path_morepval[[i]])[unlist(omics_path_more[[i]]) %in% rownames(matrix) & !duplicated(unlist(omics_path_more[[i]]))]
    } else {
      matrix[,i]= 1
    }
  }
  return(matrix)
}

Clust_mrna_all=matrix(ncol=11,nrow=62)
colnames(Clust_mrna_all)= c(sort(paste0('C', 1:11)))
rownames(Clust_mrna_all)= llista_mrna_path[!is.na(llista_mrna_path)]

Clust_methy_all=matrix(ncol=11, nrow=85)
colnames(Clust_methy_all)= c(sort(paste0('C', 1:11)))
rownames(Clust_methy_all)= llista_methy_path[!is.na(llista_methy_path)]

Clust_micro_all=matrix(ncol=11, nrow= 38)
colnames(Clust_micro_all)= c(sort(paste0('C', 1:11)))
rownames(Clust_micro_all)= llista_micro_path[!is.na(llista_micro_path)]

heatmap_mrna=Clustering_matrix(mrna_path_more,mrna_path_more_pvals, Clust_mrna_all)
heatmap_mrna[is.na(heatmap_mrna)]=1
heatmap_mrna=heatmap_mrna[,c(1,4,5,6,7,8,9,10,11,2,3)]
heatmap_mrna=heatmap_mrna[-c(6,7,16),] # delete repeated pathways

heatmap_methy= Clustering_matrix(methy_path_more, methy_path_more_pvals, Clust_methy_all)
heatmap_methy[is.na(heatmap_methy)]=1
heatmap_methy= heatmap_methy[-c(3,10,11,12,13,14,15,17,18,23,35,33,32,30),]
heatmap_methy= heatmap_methy[,c(1,4,5,6,7,8,9,10,11,2,3)]

heatmap_micro= Clustering_matrix(micro_path_more, micro_path_more_pval,Clust_micro_all)
heatmap_micro[is.na(heatmap_micro)]=1
heatmap_micro= heatmap_micro[-c(18,17),]
heatmap_micro= heatmap_micro[,c(1,4,5,6,7,8,9,10,11,2,3)]

## CLINICAL DATA PROFILING OF THE COMMUNITIES FOR ALL THE POSSIBLE COMBINATIONS BETWEEN THE LAYERS ##

llista_combinations=list()
llista_comb_mrna=list()
llista_comb_methy= list()
llista_comb_micro=list()
llista_comb_p=list()
phenoall_comb=list()
phenoall$Clust=NULL
for (a in 1:length(combinations)){
    llista_tab= list()
    llista_w= list()
    llista_mrna=list()
    llista_methy=list()
    llista_micro=list()
    llista_p=list()
    d=0
    c=0
    f=0
    phenoall= ClusteringParsing(file= paste0('Data/MolTi/', combinations[[a]]), phenoall = phenoall)
    print(combinations[[a]])
    llista_categ= list(Cardio=phenoall$Farmacos_Uso_Cardiovascular, Gender=phenoall$Genero_all.x, Emphysema=phenoall$CategoriaEnfisema, Cortis=phenoall$CortisInhalados_all, GOLD=phenoall$GOLD)
    llista_cont= list(FEV1=phenoall$pcFEV1_Pre_BD,FVC=phenoall$pc_FVC_Pre_BD ,Packs=phenoall$Paquetes, Age=phenoall$Edad_all.x, FEV1_FVC=phenoall$FEV1_FVC)
    for (i in 1:length(levels(phenoall$Clust))){
        # Fisher
        llista_tab= append(llista_tab, list(purrr::map(llista_categ,function(x) {fisher.test(table(phenoall$Clust == levels(phenoall$Clust)[[i]],x), alternative = 'two.sided')$p.value})))
        names(llista_tab)[[i]]= levels(phenoall$Clust)[[i]]
        # Wilcoxon
        llista_w= append(llista_w, list(purrr::map(llista_cont,function(x) {if((length(as.numeric(x)[phenoall$Clust == levels(phenoall$Clust)[[i]] & !is.na(x)])) > 0){wilcox.test(as.numeric(x) ~ phenoall$Clust == levels(phenoall$Clust)[[i]])$p.value} else {1}})))
        names(llista_w)[[i]]= levels(phenoall$Clust)[[i]]
        # Fisher with other layers
        phenoall_mrna_comb= phenoall_mrna[phenoall_mrna$ID %in% phenoall$ID,]
        phenoall_comb_mrna= phenoall[phenoall$ID %in% phenoall_mrna$ID,]    
        for(b in 1:length(levels(phenoall_mrna_comb$Clust))){
            d= d + 1
            llista_mrna= append(llista_mrna, fisher.test(table(phenoall_mrna_comb$Clust == levels(phenoall_mrna_comb$Clust)[[b]], phenoall_comb_mrna$Clust == levels(phenoall_comb_mrna$Clust)[[i]]), alternative = 'two.sided')$p.value)
            names(llista_mrna)[[d]]= paste0(levels(phenoall_mrna_comb$Clust)[[b]],'-', levels(phenoall_comb_mrna$Clust)[[i]])
        }
        phenoall_methy_comb= phenoall_methy[phenoall_methy$ID %in% phenoall$ID,]
        phenoall_comb_methy= phenoall[phenoall$ID %in% phenoall_methy$ID,]
        for(b in 1:length(levels(phenoall_methy_comb$Clust))){
            c= c + 1
            llista_methy= append(llista_methy, fisher.test(table(phenoall_methy_comb$Clust == levels(phenoall_methy_comb$Clust)[[b]], phenoall_comb_methy$Clust == levels(phenoall_comb_methy$Clust)[[i]]), alternative = 'two.sided')$p.value)
            names(llista_methy)[[c]]= paste0(levels(phenoall_methy_comb$Clust)[[b]],'-', levels(phenoall_comb_methy$Clust)[[i]])
        }
        phenoall_micro_comb= phenoall_micro[phenoall_micro$ID %in% phenoall$ID,]
        phenoall_comb_micro= phenoall[phenoall$ID %in% phenoall_micro$ID,]
        for(b in 1:length(levels(phenoall_micro_comb$Clust))){
            f= f + 1
            if(dim(table(phenoall_comb_micro$Clust == levels(phenoall_comb_micro$Clust)[[i]])) > 1){
            llista_micro= append(llista_micro, fisher.test(table(phenoall_micro_comb$Clust == levels(phenoall_micro_comb$Clust)[[b]], phenoall_comb_micro$Clust == levels(phenoall_comb_micro$Clust)[[i]]), alternative = 'two.sided')$p.value)
            } else{
                llista_micro= append(llista_micro, 1)
            }
            names(llista_micro)[[f]]= paste0(levels(phenoall_micro_comb$Clust)[[b]],'-', levels(phenoall_comb_micro$Clust)[[i]])
        }
        # Percentages of the clinical groups
        p= table(phenoall$Clust == levels(phenoall$Clust)[[i]], phenoall$group5)
        llista_p= append(llista_p, list(list((p[2]/(p[2] + p[4] + p[6]))*100,(p[4]/(p[2] + p[4] + p[6]))*100, (p[6]/(p[2] + p[4] + p[6]))*100)))

    }
    # Clinical variables
    test=purrr::map2(llista_tab, llista_w, function(x,y) c(x,y))
    llista_combinations=append(llista_combinations,list(test))
    names(llista_combinations)[[a]]= combinations[[a]]
    # Overlap with mRNA communities
    llista_comb_mrna= append(llista_comb_mrna, list(llista_mrna))
    names(llista_comb_mrna)[[a]]= combinations[[a]]
    llista_comb_methy= append(llista_comb_methy, list(llista_methy))
    names(llista_comb_methy)[[a]]= combinations[[a]]
    llista_comb_micro= append(llista_comb_micro, list(llista_micro))
    names(llista_comb_micro)[[a]]= combinations[[a]]
    # Percentages of the clinical groups
    tab=llista_p[[1]] %>% unlist()
    for (i in 2:length(llista_p)){tab= rbind(tab, llista_p[[i]] %>% unlist())}
    rownames(tab)= levels(phenoall$Clust)
    colnames(tab)= c('GOLD12', 'GOLD34', 'NON-SMOKERS')
    llista_comb_p= append(llista_comb_p, list(tab))
    names(llista_comb_p)[[a]]= combinations[[a]]
    phenoall_comb= append(phenoall_comb, list(phenoall))
    names(phenoall_comb)[[a]]= combinations[[a]]
    phenoall$Clust=NULL  
}
llista_combinations
llista_comb_p
llista_comb_mrna
llista_comb_methy
llista_comb_micro



