
## CODE FOR THE DIFFERENTIAL EXPRESSION ANALYSIS BETWEEN GOLD1-2, GOLD3-4 AND NON-SMOKERS INDIVIDUALS IN mRNA, miRNA AND METHYLOME
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat

#load the packages
library(tidyverse)
library(dplyr)
library(oligo)
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
library(purrr)
library(org.Hs.eg.db)
library(AnnotationHub)
library(ChAMP)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#generate required directories
dir.created("Data/mRNA")
dir.created("Data/miRNA")
dir.created("Data/Methylome")

map <- purrr::map
select <- dplyr::select
rename <- dplyr::rename
mutate <- dplyr::mutate
filter <- dplyr::filter

                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                       ##                       Generate genes/probes x individuals matrices for mRNA, miRNA and methylation                 ##
                       ##                       from the differential expression analysis between G1-2, G3-4 and non-smokers.                ##
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

### mRNA Data ###
#load cel files and clinical/array data
phenoall=read.xlsx('Data/phenoall.xlsx') # clinical data
celfiles =list.celfiles('Data/RawData/mRNA') # microarray 2018
celfiles80 <- list.celfiles('Data/RawData/mRNA/wetransfer-9a9283') # microarray 2019
celfiles16 <- list.celfiles('Data/RawData/mRNA/wetransfer-132b1d') # microarray 2019
files= phenoall %>% filter(!is.na(mRNA_teixit)) %>% select(ID, mRNATissue)
files= files[order(files$ID),]
# adding the origin of the files to evaluate a possible batch effect in further analysis. 
phenoall= phenoall %>% mutate(ORIGIN= ifelse(phenoall$mRNATissue %in% files$mRNATissue, '2018', ifelse(phenoall$ID %in% c(map(celfiles80, ~sub('_.*.CEL$', '', .x)), map(celfiles16, ~sub('_.*.CEL$', '', .x))), '2019', NA)))
phenoall$mRNA_teixit[!is.na(phenoall$ORIGIN)] <- 'SI'
# defining the phenotypes
phenoall= phenoall %>% dplyr::mutate(group5 = 'no.fumador',
group5 = ifelse(GOLD %in% c('1','2') & Group == 'EPOC EXFUMADOR','gold12',group5), 
group5 = ifelse(GOLD %in% c('1','2') & Group == 'EPOC FUMADOR','gold12_S',group5), 
group5 = ifelse(GOLD %in% c('3','4') & Group == 'EPOC FUMADOR','gold34_S',group5), 
group5 = ifelse(GOLD %in% c('3','4') & Group == 'EPOC EXFUMADOR','gold34',group5), 
group5 = ifelse(GOLD %in% c('NO') & Group == 'NO FUMADOR','no.fumador',group5), 
group5 = ifelse(GOLD %in% c('NO') & Group == 'FUMADOR SANO','fumador.sano',group5),
group5 = ifelse(GOLD %in% c('NO') & Group == 'EXFUMADOR SANO','ex.fumador.sano',group5))
phenoall= phenoall %>% filter(group5 %in% c('gold12', 'gold34', 'no.fumador'))
phenoall= phenoall[order(phenoall$ID),]
phenoall_mrna= phenoall %>% filter(!is.na(mRNA_teixit)) # 154 individuals
phenoall_micro= phenoall %>% filter(!is.na(small_teixit)) # 87 individuals
phenoall_methy= phenoall %>% filter(!is.na(Metilacio_teixit)) # 150 individuals
# read celfiles
celfiles_all= c(celfiles[celfiles %in% phenoall_mrna$mRNATissue], celfiles80[map(celfiles80,~sub('_.*.CEL$', '', .x)) %in% phenoall_mrna$ID], celfiles16[map(celfiles16,~sub('_.*.CEL$', '', .x)) %in% phenoall_mrna$ID])
affybatch<-read.affybatch(paste0('Data/RawData/mRNA', celfiles_all), phenoData = AnnotatedDataFrame(phenoall_mrna))
# plot the distribution of log base 2 intensities 
pmexp= exprs(affybatch)
sampleNames = vector(); logs = vector(); group= vector()
for (i in 1:dim(affybatch@assayData$exprs)[2])
{
sampleNames = c(sampleNames,rep(affybatch@phenoData@data[i,1],dim(pmexp)[1]))
group= c(group, rep(affybatch@phenoData$group5[i], dim(pmexp)[1]))
logs = c(logs,log2(pmexp[,i]))
}
logData = data.frame(logInt=logs,sampleName=sampleNames, Group= as.character(group))
ggplot2::ggplot(logData, aes(x=logData$logInt, color= logData$Group)) + geom_density()
# MA plots for each sample before normalising
for (i in 1:154){ oligo::MAplot(affybatch, which=i) }
# Normalization of microarray data using Robust Multiarray Average (RMA)
affybatch.r= affy::rma(affybatch)
for (i in 1:10){oligo::MAplot(affybatch.r, which=i,)}
# Get ENTREZ ID from probes names
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
dict_entrez_symbol= getBM(attributes=c('entrezgene_id',"external_gene_name", "hgnc_symbol"), values='*', mart= ensembl) %>% 
rename(entrez_id = entrezgene_id, symbol = external_gene_name)
# Get the probe identifiers that are mapped to an ENTREZ Gene ID
mapped_probes <- mappedkeys(hgu219ENTREZID)
xx <- as.list(hgu219ENTREZID[mapped_probes])
probe_entrez = data_frame(
probe_id = xx %>% unlist %>% names,
entrez_id = xx %>% unlist %>% as.numeric) %>% 
left_join(dict_entrez_symbol) %>% filter(!is.na(symbol)) %>% 
distinct(entrez_id, .keep_all=TRUE) %>% 
distinct(symbol, .keep_all=TRUE) %>% 
dplyr::select(-hgnc_symbol)
eqs = data_frame(probe_id = row.names(affybatch.r@assayData$exprs)) %>% left_join(probe_entrez)
# Collapse the rows of expression matrix by forming an average for each group of rows specified by symbol (grouping variable).
data_rm= affybatch.r@assayData$exprs
collapse.object=collapseRows(datET=as.matrix(data_rm), rowGroup=eqs$symbol, rowID=eqs$probe_id)
data_rma_collapse = collapse.object$datETcollapsed
# Remove genes exhibiting little variation across samples
collapse_filt = varFilter(as.matrix(data_rma_collapse), var.func=IQR, var.cutoff=0.25, filterByQuantile=TRUE)
# PCA
pca <- prcomp(t(collapse_filt),scale.=TRUE)
g <- ggbiplot::ggbiplot(pca, choices = 1:2,  obs.scale = 1, var.scale = 1, 
labels = phenoall_mrna$ID,
labels.size = 3, 
size = 10, 
groups = as.factor(phenoall_mrna$ORIGIN)) 
ellipse = FALSE, circle = TRUE, var.axes = F) + ggtitle('Batch:') 
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top') + scale_color_manual(values = c("red","coral4","chartreuse4","blue2"))
# Hierarchical clustering of samples
d <- as.dist(1 - cor(collapse_filt, method = "spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(phenoall_mrna$ORIGIN))
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- colnames(collapse_filt)
outcome <- phenoall_mrna$ID
names(outcome) <- colnames(collapse_filt)
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) {
  if (is.leaf(x)) {
    attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")]))  
    attr(x, "label") <- as.vector(labels[attr(x, "label")]) 
  }
  x
}, batch, outcome)

plot(sampleDendrogram, main = "Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(phenoall_mrna$ORIGIN))), fill=sort(unique(batch)))
# DE analysis
fact = phenoall_mrna %>% pull(group5) %>% as.factor
design <- model.matrix(~0 + fact)
# Estimating surrogate variables using SVA algorithm
sv <- sva(collapse_filt,design)
design <- cbind(design, sv$sv)
colnames(design)[4:length(colnames(design))] = paste0('SV_', 1:sv$n.sv)
arg= c('factno.fumador-factgold34','factno.fumador-factgold12', 'factgold12-factgold34', list(levels= design))
contrast.matrix <- do.call(makeContrasts, arg)
fit <- lmFit(collapse_filt, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
diff_no34= topTable(fit, n= 1000000, adjust.method = 'BH', coef = 'factno.fumador-factgold34') %>% 
  as.data.frame %>% 
  rownames_to_column('symbol') %>%
  mutate(pair='no.fumador_to_gold34') %>% 
  mutate(fc=round(logratio2foldchange(logFC),2)) %>% 
   mutate_at( vars(logFC, AveExpr, t, B) , funs(round(.,4))) %>%
   dplyr::rename(FDR=adj.P.Val) %>% dplyr::select(pair, everything())
diff_mrna_no12 = topTable(fit, n=1000000, adjust="BH", coef = "factno.fumador-factgold12") %>% 
  as.data.frame %>% 
  rownames_to_column('symbol') %>%
  mutate(pair='no.fumador_to_gold12') %>% 
  mutate(fc=round(logratio2foldchange(logFC),2)) %>% 
   mutate_at( vars(logFC, AveExpr, t, B) , funs(round(.,4))) %>%
   dplyr::rename(FDR=adj.P.Val) %>% dplyr::select(pair, everything())
diff_mrna_1234 = topTable(
    fit, n=1000000, adjust="BH", coef = "factgold12-factgold34") %>% 
  as.data.frame %>% 
  rownames_to_column('symbol') %>%
  mutate(pair='factgold34_to_gold12') %>% 
  mutate(fc=round(logratio2foldchange(logFC),2)) %>% 
   mutate_at( vars(logFC, AveExpr, t, B) , funs(round(.,4))) %>%
   dplyr::rename(FDR=adj.P.Val) %>% dplyr::select(pair, everything())
# Contrast GOLD1-2 and GOLD3-4 with FEV1 values
collapse_filt_cont= collapse_filt[,!colnames(collapse_filt) %in% phenoall$ID[is.na(phenoall$FEV1_Post_BD_L)]]
collapse_filt_cont= collapse_filt_cont[,colnames(collapse_filt_cont) %in% phenoall$ID[phenoall$group5 %in% c('gold12', 'gold34')]]
fev1_num= phenoall$FEV1_Post_BD_L[phenoall$group5 %in% c('gold12', 'gold34') & !is.na(phenoall$FEV1_Post_BD_L)] %>% as.numeric()
design2 <- model.matrix(~fev1_num)
fit3 <- lmFit(collapse_filt_cont, design2)
fit3 <- eBayes(fit3)
res3 <- decideTests(fit3)
diff_mrna_1234_cont = topTable(
    fit3, n=1000000, adjust="BH") %>% 
  as.data.frame %>% 
  rownames_to_column('symbol') %>%
  mutate(pair='factgold34_to_gold12') %>% 
  mutate(fc=round(logratio2foldchange(logFC),2)) %>% 
  mutate_at( vars(logFC, AveExpr, t, B) , funs(round(.,4))) %>%
   dplyr::rename(FDR=adj.P.Val) %>% dplyr::select(pair, everything())
# Volcano plots
diff_no34 = diff_no34 %>% mutate(Sig= ifelse(FDR < 0.05, 'Yes', 'No'))
g = ggplot(data=diff_no34, aes(x=logFC, y=-log10(FDR), col= Sig)) + 
  geom_point( size=0.1)+
theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = diff_no34[diff_no34$FDR < 0.05 & abs(diff_no34$logFC) > 1,], col = "black", na.rm = TRUE, hjust = 1, label= diff_no34$symbol[diff_no34$FDR < 0.05 & abs(diff_no34$logFC) > 1]) +
xlab("log2 fold change") + ylab("-log10 FDR") + guides(fill=TRUE)

diff_mrna_no12 = diff_mrna_no12 %>% mutate(Sig= ifelse(FDR < 0.05, 'Yes', 'No'))
g = ggplot(data= diff_mrna_no12, aes(x=logFC, y=-log10(FDR), col= Sig)) + 
  geom_point( size=0.1)+
theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = diff_mrna_no12[diff_mrna_no12$FDR < 0.05 & abs(diff_mrna_no12$logFC) > 1,], col = "black", na.rm = TRUE, hjust = 1, label= diff_mrna_no12$symbol[diff_mrna_no12$FDR < 0.05 & abs(diff_mrna_no12$logFC) > 1]) +
xlab("log2 fold change") + ylab("-log10 FDR")

  diff_mrna_1234_cont = diff_mrna_1234_cont %>% mutate(Sig= ifelse(diff_mrna_1234_cont$FDR < 0.05, 'Yes', 'No'))
  g = ggplot(data= diff_mrna_1234_cont, aes(x=logFC, y=-log10(adj.P.Val), col= Sig)) + 
    geom_point( size=0.1)+
  theme(legend.position = "none") +
    ggrepel::geom_text_repel(data = diff_mrna_1234_cont[diff_mrna_1234_cont$FDR < 0.05 & abs(diff_mrna_1234_cont$logFC) > 1,], col = "black", na.rm = TRUE, hjust = 1, label= diff_mrna_1234_cont$symbol[diff_mrna_1234_cont$FDR < 0.05 & abs(diff_mrna_1234_cont$logFC) > 1]) +
  xlab("log2 fold change") + ylab("-log10 FDR")
# Selecting the top genes in each comparison
diff_no34_sel= diff_no34 %>% filter(FDR < 0.05 & abs(logFC) > 0.4)
diff_no12_sel= diff_mrna_no12 %>% filter(FDR < 0.05 & abs(logFC) > 0.4)
diff_1234_sel= diff_mrna_1234_cont %>% filter(FDR < 0.05 & abs(logFC) > 0.4)
# Matrix with the DE genes that will be used to build SNF and DIABLO models
mrna_matrix= collapse_filt[rownames(collapse_filt) %in% unique(c(diff_no34_sel$symbol,diff_no12_sel$symbol, diff_1234_sel$symbol)),]
# Matrix with most variable genes (coefficient of variation) that will be used to build MolTi model.
llista= numeric(17788)
names(llista)= rownames(data_rma_collapse)
for (l in 1:dim(data_rma_collapse)[1]){
  llista[l]= cv(data_rma_collapse[l,])
}
hist(as.numeric(llista), breaks= 'FD')
sel_var= llista[llista > quantile(as.numeric(llista), probs= 0.75)]
molti_mrna= data_rma_collapse[rownames(data_rma_collapse) %in% names(sel_var),]
Cor_mrna=cor(molti_mrna, method= 'pearson')
Cor_mrna %>% write.table('Data/mRNA/Matrix_cor_mrna.txt')
# Clustering and heatmap representation of DE genes and samples
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
mrna_matrix_scale = mrna_matrix %>% t %>% as_tibble %>% mutate_all(.funs = funs(range01(rank(., na.last = F)))) %>% t
colnames(mrna_matrix)= phenoall_mrna$ID

dend_r <- mrna_matrix_scale %>% dist(method = "minkowski") %>% hclust(method = "ward.D2" ) %>% as.dendrogram %>% ladderize %>%
    color_branches(k=3) #genes
dend_c <- t(mrna_matrix_scale) %>% dist(method = "minkowski") %>% hclust(method = "ward.D2") %>% as.dendrogram %>% ladderize %>%
    color_branches(k=3) #samples

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)

rr = as.factor(phenoall_mrna$group5)
brew = brewer.pal(length(levels(rr)), "Set1")
cols_rr = brew[rr]

cr = 0.3
ww = 11
marg = 12
pdf(
'Data/mRNA/clustering_mRNA',
 width = 10, height = 10)
gplots::heatmap.2(as.matrix(mrna_pca_matrix_scale), 
          main = "heatmap mRNA",
          srtCol = 35,
          Rowv = dend_r,
          Colv = dend_c,
          trace="none", hline = NA, tracecol = "darkgrey",         
          margins =c(4,marg),      
          key.xlab = "",
          denscol = "grey",
         density.info = "density", cexRow = cr,
          col = mypalette, cexCol = 0.01, 
ColSideColors = cols_rr)
legend("topright", legend=levels(rr), fill=brew, bg = 'white', title="", cex=1)
dev.off(); gc()
### miRNA Data ###
# Trimming and alignment of the miRNA transcriptomic data was already performed, so I loaded the data.
load("Data/miRNA/micro_data.RData")
micro_counts = cat_comps5_N$as_mat_micro_full[[1]]
micro_counts= micro_counts[,order(colnames(micro_counts))]
fact_micro_= phenoall_micro$group5 %>% as.factor()
dge <- DGEList(counts= micro_counts %>% dplyr::select(-micro) %>% as.matrix(), group= phenoall_micro$group5, genes= micro_counts$micro)
dge <- calcNormFactors(dge, method = 'TMM') # Normalization
dge$logCPM <- cpm(dge, log= TRUE, prior.count = 0.25, normalized.lib.sizes = TRUE)
# Filtering lowly-expressed genes with a minimum CPM cutoff value of expression (10 reads per million in the sample with lowest depth)
# Keep only genes that have this minimum level of expresseion in at least as many samples as the smallest group of comparison
cpmcutoff <- round(10/min(dge$sample$lib.size/1e+06), digits = 1)
nsamplescutoff <- min(table(dge$samples$group))
mask <- rowSums(cpm(dge,  log= TRUE, prior.count = 0.25, normalized.lib.sizes = TRUE) > cpmcutoff) >= nsamplescutoff
dge <- dge[mask,]
dge$logCPM <- dge$logCPM[mask,]
design= model.matrix(~ 0 + fact_micro_)
sva <- sva(cpm, design)
design <- cbind(design,sva$sv)
colnames(design)[7:length(colnames(design))]= paste0('SV_', 1:sva$n.sv)
fit <- glmQLFit(dge, design)
arg= c('fact_micro_gold34-fact_micro_no.fumador','fact_micro_gold12-fact_micro_no.fumador', 'fact_micro_gold34-fact_micro_gold12', list(levels= design))
# Estimate Common,Trended and Tagwise Dispersion for Negative Binomial 
dge <- estimateGLMCommonDisp(dge, design)
dge<- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
contrast.matrix <- do.call(makeContrasts, arg)
contrast.matrix= contrast.matrix %>% as_tibble()
# Quasi-likelihood (QL) F-test
glf_SNF1SNF2 <-glmQLFTest(fit, contrast = contrast.matrix$`fact_micro_SNF1-fact_micro_SNF2`)
glf_SNF1SNF3 <-glmQLFTest(fit, contrast = contrast.matrix$`fact_micro_SNF1-fact_micro_SNF3`)
glf_SNF2SNF3 <-glmQLFTest(fit, contrast = contrast.matrix$`fact_micro_SNF2-fact_micro_SNF3`)
diff_micro_1234= topTags(glf1234, n= 1000000)
diff_micro_no34= topTags(glfno34, n= 1000000)
diff_micro_no12= topTags(glfno12, n= 1000000)
# Volcano plot
diff_micro_no34$table= diff_micro_no34$table %>% mutate(Sig= ifelse(FDR < 0.05, 'Yes', 'No'))
g = ggplot(data=diff_micro_no34$table, aes(x=logFC, y=-log10(FDR), col= Sig)) + 
  geom_point( size=0.1)+
theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = diff_micro_no34$table[diff_micro_no34$table$FDR < 0.05 & abs(diff_micro_no34$table$logFC) > 1,], col = "black", na.rm = TRUE, hjust = 1, label= diff_micro_no34$table$genes[diff_micro_no34$table$FDR < 0.05 & abs(diff_micro_no34$table$logFC) > 1]) + xlab("log2 fold change") + ylab("-log10 FDR") + guides(fill=TRUE)

diff_methy_no12_genes= diff_methy_no12_genes %>% mutate(Sig= ifelse(FDR < 0.05, 'Yes', 'No'))
g = ggplot(data=diff_methy_no12_genes, aes(x=logFC, y=-log10(FDR), col= Sig)) + 
  geom_point( size=0.1)+
theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = diff_methy_no12_genes[diff_methy_no12_genes$FDR < 0.05 & abs(diff_methy_no12_genes$logFC) > 1,], col = "black", na.rm = TRUE, hjust = 1, label= diff_methy_no12_genes$genes[diff_methy_no12_genes$FDR < 0.05 & abs(diff_methy_no12_genes$logFC) > 1]) + xlab("log2 fold change") + ylab("-log10 FDR") + guides(fill=TRUE)

diff_methy_1234= diff_methy_1234 %>% mutate(Sig= ifelse(FDR < 0.05, 'Yes', 'No'))
g = ggplot(data=diff_methy_1234, aes(x=logFC, y=-log10(FDR), col= Sig)) + 
  geom_point( size=0.1)+
theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = diff_methy_1234[diff_methy_1234$FDR < 0.05 & abs(diff_methy_1234$logFC) > 1,], col = "black", na.rm = TRUE, hjust = 1, label= diff_methy_1234$genes[diff_methy_1234$FDR < 0.05 & abs(diff_methy_1234$logFC) > 1]) + xlab("log2 fold change") + ylab("-log10 FDR") + guides(fill=TRUE)
# Hierarchical clustering
cpm_batch= cpm[,colnames(cpm) %in% phenoall_micro$ID]
d <- as.dist(1 - cor(cpm_batch, method = "spearman"))
sampleClustering <- hclust(d)
batch <- as.factor(phenoall_micro$group5)
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- colnames(cpm_batch)
outcome <- phenoall_micro$group5
names(outcome) <- colnames(cpm_batch)
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) {
  if (is.leaf(x)) {
    attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")]))  
    attr(x, "label") <- as.vector(labels[attr(x, "label")]) 
  }
  x
}, batch, outcome)
plot(sampleDendrogram, main = "Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(phenoall$CortisInhalados_all))), fill=sort(unique(batch)))
# PCA
cpm_batch= cpm[,colnames(cpm) %in% phenoall_micro$ID]
cpm_batch= cpm_batch[rownames(cpm_batch) %in% symbols_micro,]
pca <- prcomp(t(cpm_batch),scale.=TRUE)
g <- ggbiplot::ggbiplot(pca, choices = 1:2,  obs.scale = 1, var.scale = 1, 
labels = phenoall_micro$ID,
labels.size = 3, 
size = 10, 
                groups = as.factor(phenoall_micro$group5), 
ellipse = FALSE, circle = TRUE, var.axes = F) + ggtitle('microRNA:') 
g <- g + scale_color_discrete(name = '') 
g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top') + scale_color_manual(values = c("red","coral4","chartreuse4","blue2"))
# Selecting the top probes in each comparison
diff_micro_no34_sel= diff_micro_no34$table %>% filter(FDR < 0.05)
diff_micro_no12_sel= diff_micro_no12$table %>% filter(FDR < 0.05)
diff_micro_1234_sel= diff_micro_1234$table %>% filter(FDR < 0.05)
x1= unique(diff_micro_no34_sel$genes)
x2= unique(diff_micro_no12_sel$genes) 
x3= unique(diff_micro_1234_sel$genes)
symbols_micro= unique(c(x1,x2,x3))
# Matrix of DE probes that will be used to build the SNF and DIABLO model
micro_matrix= micro_matrix[rownames(micro_matrix) %in% symbols_micro,]
# Matrix of the most variables probes (coefficient of variation) that will used for MolTi
llista2= list()
for (w in 1:dim(cpm)[1]){
  llista2[w]= cv(cpm[w,])
  names(llista2)[w]= rownames(cpm)[w]
}
hist(as.numeric(llista2), breaks= 'FD')
sel_micro= llista2[llista2 > quantile(as.numeric(llista2), probs= 0.75)]
molti_micro= cpm[rownames(cpm) %in% names(sel_micro),]
Cor_micro= cor(molti_micro, method= 'pearson') 
Cor_micro %>% write.table('Data/miRNA/Matrix_cor_miRNA.txt')
# Clustering and heatmap representation of DE probes and samples
range01 <- function(x){(x-min(x))/(max(x)-min(x))} 
micro_matrix_scale = micro_matrix %>% t %>% as_tibble %>% mutate_all(.funs = funs(range01(rank(., na.last = F)))) %>% t
colnames(micro_matrix_scale)= phenoall_micro$ID

dend_r <- micro_matrix_scale %>% dist(method = "minkowski") %>% hclust(method = "ward.D2" ) %>% as.dendrogram %>% ladderize %>%
    color_branches(k=3) #genes
dend_c <- t(micro_matrix_scale) %>% dist(method = "minkowski") %>% hclust(method = "ward.D2") %>% as.dendrogram %>% ladderize %>%
    color_branches(k=3) #samples

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)

rr = as.factor(phenoall_micro$group5)
brew = brewer.pal(length(levels(rr)), "Set1")
cols_rr = brew[rr]

cr = 0.3
ww = 11
marg = 12
pdf(
'Data/miRNA/clustering_microRNA',
 width = 10, height = 10)
gplots::heatmap.2(as.matrix(micro_matrix_scale), 
          main = "heatmap microRNA",
          srtCol = 35,
          Rowv = dend_r,
          Colv = dend_c,
          trace="none", hline = NA, tracecol = "darkgrey",         
          margins =c(4,marg),      
          key.xlab = "",
          denscol = "grey",
         density.info = "density", cexRow = cr,
          col = mypalette, cexCol = 0.01, 
ColSideColors = cols_rr )
legend("topright", legend=levels(rr), fill=brew, bg = 'white', title="", cex=1)
dev.off(); gc()
### METHYLATION ###
# Loading and filtering data
# champ.load function filters probes with detection p-value < 0.01, probes with <3 beads in at least 5% of samples per probe,
# all SNP-related probes, multi-hit probes, probes located in chromosome X and Y. 
myLoad <- champ.load('RawData/Methylation', arraytype = 'EPIC')
champ.QC()  # Draws MDS plot and density plot of raw data.
CpG.GUI(CpG=rownames(myLoad$beta),arraytype="EPIC") # Check CpG distribution on chromosomes.
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",method = 'PBC') # Normalization
champ.SVD(beta=myNorm_inter,pd=myLoad$pd) # Singular value decomposition method (SVD) to assess number and significance of possible covariates.
# Hierarchical clustering to check for possible batch effect
d <- as.dist(1 - cor(myLoad$beta, method = "spearman"))
batch <- as.integer(factor(myLoad$pd$Array)) 
sampleClustering <- hclust(d)
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- myLoad$pd$Sample_Name
outcome <- as.character(myLoad$pd$Sample_Group)
names(outcome) <-  myLoad$pd$Sample_Name
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) {
  if (is.leaf(x)) {
    attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")]))
    attr(x, "label") <- as.vector(labels[attr(x, "label")])
      }
  x
}, batch, outcome)

plot(sampleDendrogram, main= 'Clustering')

plotMDS(myLoad$beta,labels= myLoad$pd$Sample_Name,col= as.integer(factor(myLoad$pd$Sample_Plate)), dim.plot = c(7,5))
#In the clustering, it can be seen clearly that the array plate has an important batch effect.
myCombat <-champ.runCombat(beta= myNorm, pd= myLoad$pd, batchname = c('Array'))
champ.SVD(beta=myCombat, pd=myLoad$pd)
myCombat_filt= varFilter(myCombat, var.func = IQR, var.cutoff = 0.25, filterByQuantile = TRUE)
fact= phenoall_methy %>% pull('group5') %>% as.factor()
design <- model.matrix(~0+ fact)
arg= c('factgold12-factgold34', 'factno.fumador-factgold12', 'factno.fumador-factgold34', list(levels=design))
contrast.matrix <- do.call(makeContrasts, arg)
fit <- lmFit(myCombat_filt, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
diff_methy_no34_2= topTable(fit,n=100000, adjust.method = 'BH', coef= 'factno.fumador-factgold34')
diff_methy_no12_2= topTable(fit,n= 100000, adjust.method = 'BH', coef= 'factno.fumador-factgold12')
diff_methy_1234_2= topTable(fit,n= 100000, adjust.method = 'BH', coef= 'factgold12-factgold34')
# Load from ChAMP package the dataframe with the mapping between CpG probes and genes
data("probe.features.epic")
# Mapping CpG probes to the genes
diff_methy_no34_genes= diff_methy_no34 %>% rownames_to_column(var= 'probeID') %>% left_join(probe.features %>% rownames_to_column(var= 'probeID') %>% select(probeID, gene), by= 'probeID')
diff_methy_no34_genes$gene <- droplevels(diff_methy_no34_genes$gene, exclude="")
diff_methy_no34_genes= diff_methy_no34_genes %>% filter(!is.na(gene)) %>% filter(!duplicated(gene))

diff_methy_no12_genes= diff_methy_no12 %>% rownames_to_column(var= 'probeID') %>% left_join(probe.features %>% rownames_to_column(var= 'probeID') %>% select(probeID, gene), by= 'probeID')
diff_methy_no12_genes$gene <- droplevels(diff_methy_no12_genes$gene, exclude= "")
diff_methy_no12_genes= diff_methy_no12_genes %>% filter(!is.na(gene)) %>% filter(!duplicated(gene))

diff_methy_1234_genes= diff_methy_1234 %>% rownames_to_column(var= 'probeID') %>% left_join(probe.features %>% rownames_to_column(var= 'probeID') %>% select(probeID, gene), by= 'probeID')
diff_methy_1234_genes$gene= droplevels(diff_methy_1234_genes$gene, exclude= "")
diff_methy_1234_genes= diff_methy_1234_genes %>% filter(!is.na(gene)) %>% filter(!duplicated(gene))
# Volcano plot
diff_methy_no34_genes= diff_methy_no34_genes %>% mutate(Sig= ifelse(adj.P.Val < 0.05, 'Yes', 'No'))
g = ggplot(data=diff_methy_no34_genes, aes(x=logFC, y=-log10(adj.P.Val), col= Sig)) + 
  geom_point( size=0.1)+
theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = diff_methy_no34_genes[diff_methy_no34_genes$adj.P.Val < 0.05 & abs(diff_methy_no34_genes$logFC) > 1,], col = "black", na.rm = TRUE, hjust = 1, label= diff_methy_no34_genes$gene[diff_methy_no34_genes$adj.P.Val < 0.05 & abs(diff_methy_no34_genes$logFC) > 1]) + xlab("log2 fold change") + ylab("-log10 adj.P.Val") + guides(fill=TRUE)

diff_methy_no12_genes= diff_methy_no12 %>% mutate(Sig= ifelse(adj.P.Val < 0.05, 'Yes', 'No'))
g = ggplot(data=diff_methy_no12, aes(x=logFC, y=-log10(adj.P.Val), col= Sig)) + 
  geom_point( size=0.1)+
theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = diff_methy_no12[diff_methy_no12$adj.P.Val < 0.05 & abs(diff_methy_no12$logFC) > 1,], col = "black", na.rm = TRUE, hjust = 1, label= diff_methy_no12_genes$gene[diff_methy_no12_genes$adj.P.Val < 0.05 & abs(diff_methy_no12_genes$logFC) > 1]) + xlab("log2 fold change") + ylab("-log10 adj.P.Val") + guides(fill=TRUE)


diff_methy_1234= diff_methy_1234 %>% mutate(Sig= ifelse(adj.P.Val < 0.05, 'Yes', 'No'))
g = ggplot(data=diff_methy_1234, aes(x=logFC, y=-log10(adj.P.Val), col= Sig)) + 
  geom_point( size=0.1)+
theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = diff_methy_1234[diff_methy_1234$adj.P.Val < 0.05 & abs(diff_methy_1234$logFC) > 1,], col = "black", na.rm = TRUE, hjust = 1, label= diff_methy_1234$genes[diff_methy_1234$adj.P.Val < 0.05 & abs(diff_methy_1234$logFC) > 1]) + xlab("log2 fold change") + ylab("-log10 FDR") + guides(fill=TRUE)

# Selecting top genes in each comparison
diff_methy_no34_venn = diff_methy_no34_genes %>% filter(adj.P.Val < 0.05)
diff_methy_no12_venn = diff_methy_no12_genes %>% filter(adj.P.Val < 0.05)
diff_methy_1234_venn = diff_methy_1234_genes %>% filter(adj.P.Val < 0.05)
x1= unique(diff_methy_no34_venn$probe)
x2= unique(diff_methy_no12_venn$probe)
x3= unique(diff_methy_1234_venn$probe)
symbols_met= unique(c(x1,x2,x3))
# Matrix of DE genes that will be used to build the SNF and DIABLO models
methy_matrix= myCombat_filt[rownames(myCombat_filt) %in% symbols_met,]
# Matrix with most variable genes (coefficient of variation) that will be used to build MolTi model.
llista3 <- numeric(23125)
names(llista3)= rownames(myCombat)
for (k in 1:dim(myCombat)[1]){
  llista3[k] <- cv(myCombat_molti[k,])
}
hist(llista3, xlim= c(0,0.4), breaks= 'FD')
sel_methy= llista3[llista3 > quantile(llista3, probs= 0.75)]
molti_methy= myCombat[rownames(myCombat) %in% names(sel_methy),]
Cor_methy=cor(molti_methy, method= 'pearson') 
Cor_methy %>% write.table('Data/Methylation/Matrix_cor_methy.txt')
# Clustering and heatmap of the DE probes and samples
range01 <- function(x){(x-min(x))/(max(x)-min(x))} #re-scales a numeric vector so it ranges from 1 to 10
methy_matrix_scale = methy_matrix %>% t %>% as_tibble %>% mutate_all(.funs = funs(range01(rank(., na.last = F)))) %>% t
      
dend_r <- methy_matrix_scale %>% dist(method = "minkowski") %>% hclust(method = "ward.D2" ) %>% as.dendrogram %>% ladderize %>%
    color_branches(k=3) #genes
dend_c <- t(methy_matrix_scale) %>% dist(method = "minkowski") %>% hclust(method = "ward.D2") %>% as.dendrogram %>% ladderize %>% color_branches(k=3) #samples

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)

rr = as.factor(phenoall_methy$group5)
brew = brewer.pal(length(levels(rr)), "Set1")
cols_rr = brew[rr]

cr = 0.3
ww = 11
marg = 12
pdf(
'Data/Methylation/clustering_methylation',
 width = 10, height = 10)
gplots::heatmap.2(as.matrix(methy_pca_matrix_scale), 
          main = "heatmap methylation",
          srtCol = 35,
          Rowv = dend_r,
          Colv = dend_c,
          trace="none", hline = NA, tracecol = "darkgrey",         
          margins =c(4,marg),      
          key.xlab = "",
          denscol = "grey",
         density.info = "density", cexRow = cr,
          col = mypalette, cexCol = 0.01, 
ColSideColors = cols_rr )
legend("topright", legend=levels(rr), fill=brew, bg = 'white', title="", cex=1)
dev.off(); gc()






