# COPDHeterogeneity

This repository contains the scripts used to perform the analysis described in the master thesis "New insights into molecular heterogeneity of Chronic Obstructive Pulmonary Disease", supervised by Rosa Faner, Rosalba Lepora and Jon Sanchez from IDIBAPS and Barcelona Supercomputing Center. The omics and clinical data that were used in this study are not currently available in the Gene Expression Omnibus (GEO) repository since the article has not been published in a peer-review journal yet. 

# Data

Documents showing how the data is organised and the different files used in the study are found in the Data/RawData directory. 

# Scripts to be run

1. DEAnalysis.r\
This is the first script that needs to be run. It performs a the differential expression (DE) analysis for mRNA (microarray), miRNA (transcriptomics) and Methylation (microarray) between the clinical groups (GOLD 1-2, GOLD 3-4 and NON-SMOKERS). From these analysis, we retrieve the genes/probes x individuals matrices that will be used to build the multi-omics models. Moreover, several diagnostic plots such as PCA, hierarchical clustering, heatmap and volcano plot are plotted to evaluate the data. 

2. DIABLO\
DIABLO.r builds a Data Integration Analysis for Biomarker discovery using Latent variable approaches for ‘Omics studies (DIABLO) model for all the possible omics combinations. It also returns the plots showed in the main paper and Supplementary Material 1 (individual plot, arrow plot, clustered image map).

3. SNF\
SNF.r builds a Similarity Network Fusion (SNF) model, performs the enrichment analysis using GSEA, ebGSEA and Fisher's exact test for mRNA, methylome and miRNA, respectively. It also returns the clinical tables, where the available clinical variables were evaluated to see possible differences between the SNF groups. NcolFormat.py transforms a individual x individual correlation/affinity matrix to a an unweighted graph in the ncol format through a cutoff given by the user to decide whether two nodes should be connected or not. Once it is in the suitable format, the file can be given to Cytoscape to represent the networks.

4. MOLTI\
DistanceClosure.py returns a two column file, where the pair of nodes that are connected through a backbone edge, are represented. GenerateRandomGraphs.py creates 10000 random graphs preserving the degree distribution of the original network and runs the MolTi software on each graph to recover the random communities. ProbabilityCalculation.py computes the probability of finding a community with the same size, the nodes and pair of nodes as the given original network. RandomNetworksClinicalData.r computes the probability of finding clinically relevant communities (clusters with at least one enriched feature) in the random networks. EnrichmentMolTi.r performs the differential expression and enrichment analysis for all the possible comparisons between the eleven communities retrieved with MolTi. It also plots the heatmap of the adjusted p-values of the pathways and the corr plots of the Jaccard Index. PairNodesTogether.py computes the number of pairs that are together in the same community when the resolution parameter is tuned. DendogramResolution.r generates a dendrogram plot, where the distances between the nodes are based on the frequency that nodes tend to be in the same community.  
