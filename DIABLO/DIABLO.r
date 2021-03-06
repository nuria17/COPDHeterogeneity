## CODE FOR COMPUTING THE DIABLO MODEL WITH THE mRNA, miRNA AND METHYLATION DATA BLOCKS. 
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat

#load the packages
library(mixOmics)
library(dplyr)
library(purrr)
library(tidyverse)
library(readxl)
library(xlsx)

#generate required directories
dir.created("Data/DIABLO")

map <- purrr::map
select <- dplyr::select
rename <- dplyr::rename
mutate <- dplyr::mutate
filter <- dplyr::filter

                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                       ##                       Compute the DIABLO model from MixOmics package, using mRNA, miRNA and methylation            ##
                       ##                       data blocks.                                                                                 ##
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

# Before running this script, DEanalysis.r must be run first in order to have the mRNA, miRNA and Methylation 
# genes/probes x individuals matrices.
# As input, it needs the matrices, the dataframe with the clinical groups and a string indicating the suffix of the output.
DIABLO_model <- function(matrixA, matrixB, matrixC=NULL, phenoall, output){
    Y= phenoall$group5 %>% as.factor() # outcome
    #  Non sparse analysis with block.splsda to examine the correlation between
    #  the different blocks via the modelled components and gettting the max correlation
    #  to determine the connection between them
    if (is.null(matrixC)){
        data= list(datatypeA= t(matrixA), datatypeB= t(matrixB))           
        param= block.splsda(data, Y)
        connection= cor(param$variates$datatypeA[,1], param$variates$datatypeB[,1])
    } else {
        data= list(datatypeA= t(matrixA), datatypeB= t(matrixB), datatypeC= t(matrixC))        
        param= block.splsda(data, Y)
        connection= max(cor(param$variates$datatypeA[,1], param$variates$datatypeB[,1]), cor(param$variates$datatypeA[,1], param$variates$datatypeC[,1]),cor(param$variates$datatypeB[,1], param$variates$datatypeC[,1]))
    }
    # The matrix design determines which blocks should be connected to 
    # maximize the correlation or covariance between components.      
    design = matrix(connection, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
    diag(design) = 0
    # Fit DIABLO model without variable selection and choose the number of components for the final model.
    sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, 
                           design = design)
    set.seed(123)
    perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
    png(paste0('Perf.diablo_', output))
    plot(perf.diablo)
    dev.off()
    ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
    # Choose the optimal number of variables to select in each data set using the tune function, for a grid of keepX values
    test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
    tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 5,
                              cpus = 10, dist = "centroids.dist")
    list.keepX = tune.TCGA$choice.keepX # number of features to select on each component and dataset
    # Final model : Projection to Latent Structures models (PLS) with sparse Discriminant Analysis
    sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)
    # Sample plots 
    # Plot DIABLO diagnostic plot to check whether the correlation
    # between components from each data set has been maximized as specified in the design matrix.
    png(paste0('PlotDIABLO_', output))
    plotDiablo(sgccda.res, ncomp = 1)
    dev.off()
    # This plot projects each sample into the space spanned by the components of each block. 
    png(paste0('PlotIndiv_', output))
    plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
    dev.off()
    # Arrow plot
    png(paste0('PlotArrow_', output))
    plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
    dev.off()
    # Clustered image map of multi-omics molecular signature
    png(paste0('cimDIABLO_', output))
    cimDiablo(sgccda.res)
    dev.off()
    # Performance of the model using 10-fold cross-validation repeated 10 times 
    set.seed(123) 
    perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, 
                   dist = 'centroids.dist')
    # Performance based on the weighted vote, where weight is defined according to the correlation
    # between the latent component associated to a particular data set and the outcome. 
    Performance=perf.diablo$WeightedVote.error.rate
    
    return(Performance)
}

phenoall= read_excel('Data/phenoall.xlsx')
phenoall_mrna= phenoall[!is.na(phenoall$mRNA_teixit),]
phenoall_micro=phenoall[!is.na(phenoall$small_teixit),]
phenoall_methy= phenoall[!is.na(phenoall$Metilacio_teixit),]

# DIABLO model with mRNA and miRNA matrices
mrna_inter_micro= mrna_matrix[,colnames(mrna_matrix) %in% intersect(phenoall_mrna$ID, phenoall_micro$ID)]
micro_inter_mrna= micro_matrix[,colnames(micro_matrix) %in% intersect(phenoall_mrna$ID, phenoall_micro$ID)]
phenoall_micro_mrna= phenoall[phenoall$ID %in% intersect(phenoall_mrna$ID, phenoall_micro$ID),]
Perf_mrna_micro= DIABLO_model(mrna_inter_micro, micro_inter_mrna,phenoall=phenoall_micro_mrna, output = 'mrna_micro')

# DIABLO model with methylation and miRNA matrices
methy_inter_micro= methy_matrix[,colnames(methy_matrix) %in% intersect(phenoall_methy$ID, phenoall_micro$ID)]
micro_inter_methy= micro_matrix[,colnames(micro_matrix) %in% intersect(phenoall_methy$ID, phenoall_micro$ID)]
phenoall_methy_micro= phenoall[phenoall$ID %in% intersect(phenoall_methy$ID, phenoall_micro$ID),]
Perf_methy_micro= DIABLO_model(methy_inter_micro, micro_inter_methy,phenoall_methy_micro, output = 'methy_micro' )

# DIABLO model with mRNA and methylation matrices
mrna_inter_methy= mrna_matrix[,colnames(mrna_matrix) %in% intersect(phenoall_mrna$ID, phenoall_methy$ID)]
methy_inter_mrna= methy_matrix[,colnames(methy_matrix) %in% intersect(phenoall_mrna$ID, phenoall_methy$ID)]
phenoall_mrna_methy= phenoall[phenoall$ID %in% intersect(phenoall_mrna$ID, phenoall_methy$ID),]
Perf_mrna_methy= DIABLO_model(mrna_inter_methy, methy_inter_mrna, phenoall_mrna_methy, output = 'mrna_methy')

# DIABLO model with mRNA, miRNA and methylation matrices
mrna_inter= mrna_matrix[,colnames(mrna_matrix) %in% Reduce(intersect, list(phenoall_mrna$ID, phenoall_methy$ID, phenoall_micro$ID))]
methy_inter= methy_matrix[,colnames(methy_matrix) %in% Reduce(intersect, list(phenoall_mrna$ID, phenoall_methy$ID, phenoall_micro$ID))]
micro_inter= micro_matrix[,colnames(micro_matrix) %in% Reduce(intersect, list(phenoall_mrna$ID, phenoall_methy$ID, phenoall_micro$ID))]
phenoall_inter= phenoall[phenoall$ID %in% Reduce(intersect, list(phenoall_mrna$ID, phenoall_methy$ID, phenoall_micro$ID),]
Perf_inter= DIABLO_model(mrna_inter, micro_inter, methy_inter, phenoall= phenoall_inter, output = 'inter')

