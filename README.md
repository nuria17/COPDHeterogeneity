# COPDHeterogeneity

This repository contains the scripts used to perform the analysis described in the master thesis "New insights into molecular heterogeneity of Chronic Obstructive Pulmonary Disease", supervised by Rosa Faner, Rosalba Lepora and Jon Sanchez from IDIBAPS and Barcelona Supercomputing Center. The omics and clinical data that were used in this study are not currently available in the Gene Expression Omnibus (GEO) repository since the article has not been published in a peer-review journal yet. 

# Data

Documents showing how the data is organised and the different files used in the study are found in the Data/RawData directory. 

# Scripts to be run

1. DEAnalysis.r \\
This is the first script that needs to be run. It performs a the differential expression (DE) analysis for mRNA (microarray), miRNA (transcriptomics) and Methylation (microarray) between the clinical groups (GOLD 1-2, GOLD 3-4 and NON-SMOKERS). From these analysis, we retrieve the genes/probes x individuals matrices that will be used to build the multi-omics models. Moreover, several diagnostic plots such as PCA, hierarchical clustering, heatmap and volcano plot are plotted to evaluate the data. 

2. DIABLO \\
DIABLO.r builds a Data Integration Analysis for Biomarker discovery using Latent variable approaches for ‘Omics studies (DIABLO) model for all the possible omics combinations. It also returns the plots showed in the main paper and Supplementary Material 1 (individual plot, arrow plot, clustered image map).

3. SNF \\
SNF.r builds a Similarity Network Fusion (SNF) model, performs the enrichment analysis using GSEA, ebGSEA and Fisher's exact test for mRNA, methylome and miRNA, respectively. It also returns the clinical tables, where the available clinical variables were evaluated to see possible differences between the SNF groups. 
