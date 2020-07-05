## CODE FOR CALCULATING THE PROBABILITY OF FINDING A COMMUNITY IN THE RANDOM NETWORKS THAT AT LEAST ENRICHED ONE CLINICAL FEATURE
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat


### R LIBRARIES ###

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
library(readxl)
library(dendextend)
library(Spectrum)
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                       ##                                         Performs the Fisher's exact test and Wilcoxon signed-rank test             ##
                       ##                                           for all the random communites.                                           ##
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##


# First, it is needed to load the clinical data
Random_communities_clinical_data <- function(random_communities_files, phenoall){
  counter= 0
  total_comm=0
  llista_probabilities=list()
  phenoall$Clust=NULL
  for (i in 1:length(random_communities_files)){
    tmp= read.table(paste0('Data/MolTi/Random_graphs/', random_communities_files[[i]]))
    print(random_communities_files[[i]])
    clust=grep('ClusterID:',unlist(tmp))
    llista_clust_all=list()
    for (i in 1:(length(clust))){ if (i < (length(clust))) {llista_clust_all=append(llista_clust_all,list(tmp[(clust[i]+1):(clust[i+1]-1),]))} else {llista_clust_all=append(llista_clust_all,list(tmp[(clust[i]+1):dim(tmp)[[1]],]))}}
    for (e in 1:length(llista_clust_all)) {phenoall$Clust[phenoall$ID %in% llista_clust_all[[e]]]= paste0('Cluster', e)}
    phenoall$Clust= as.factor(phenoall$Clust)
    llista_tab= list()
    llista_w= list()
    llista_categ= list(phenoall$Farmacos_Uso_Cardiovascular, phenoall$Genero_all.x, phenoall$CategoriaEnfisema, phenoall$CortisInhalados_all, phenoall$GOLD)
    llista_cont= list(phenoall$pcFEV1_Pre_BD,phenoall$pc_FVC_Pre_BD ,phenoall$Paquetes, phenoall$Edad_all.x, phenoall$FEV1_FVC)
      for (i in 1:length(levels(phenoall$Clust))){
        # Fisher
        llista_tab= append(llista_tab, list(purrr::map(llista_categ,function(x) {fisher.test(table(phenoall$Clust == levels(phenoall$Clust)[[i]],x), alternative = 'two.sided')$p.value})))
        names(llista_tab)[[i]]= levels(phenoall$Clust)[[i]]
        # Wilcoxon
        llista_w= append(llista_w, list(purrr::map(llista_cont,function(x) {if((length(as.numeric(x)[phenoall$Clust == levels(phenoall$Clust)[[i]] & !is.na(x)])) > 0){wilcox.test(as.numeric(x) ~ phenoall$Clust == levels(phenoall$Clust)[[i]])$p.value} else {1}})))
        names(llista_w)[[i]]= levels(phenoall$Clust)[[i]]
    }
    sig_categ= purrr::map(llista_tab, function(x) if(length(unlist(x)) > 0){x[x < 0.05][1]})
    sig_cont= purrr::map(llista_w, function(x){x[x < 0.05][1]})
    comm_sig_names= unique(c(names(unlist(sig_categ)), names(unlist(sig_cont))))
    comm_sig= length(unique(c(names(unlist(sig_categ)), names(unlist(sig_cont)))))
    llista_probabilities= append(llista_probabilities, list(comm_sig/length(levels(phenoall$Clust))))
    total_comm= total_comm + length(llista_tab)
    counter= counter + comm_sig
    phenoall$Clust=NULL
  }
  return(list(list(counter), list(total_comm), list(llista_probabilities)))
}
### CALLING FUNCTIONS ###
phenoall= read_excel('Data/Phenoall.xlsx')
random_communities_files= list.files('Data/MolTi/Random_graphs/', pattern= 'Random_clusters_all_([0-9])$')
Ran_comm_clinical=Random_communities_clinical_data(random_communities_files = random_communities_files, phenoall = phenoall)
