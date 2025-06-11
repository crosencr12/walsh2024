# Created a txt and csv files to point to of all htseq-counts data from all sequencing samples
## Make sure htseq-count data is in DESeq2 folder for processing using this structure

# For a new R directory/project:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("RColorBrewer")
install.packages("pheatmap")
BiocManager::install("DESeq2")
BiocManager::install('ComplexHeatmap')
BiocManager::install('EnhancedVolcano')
## Heat map resource using: https://link.springer.com/protocol/10.1007/978-1-0716-1084-8_17
install.packages("remotes")
remotes::install_github('YuLab-SMU/ggtree') ## needed for clusterprofiler
BiocManager::install("clusterProfiler")
BiocManager::install('org.Hs.eg.db')
# BiocManager::install('circlize')
# BiocManager::install('digest')
# BiocManager::install("biomaRt")


## Going to need these packages
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library("DESeq2")
library("ComplexHeatmap")
library('clusterProfiler')
library("circlize")
library('org.Hs.eg.db')
library('EnhancedVolcano')
# library(biomaRt)
# library(reshape2)
# library("digest")
# library("cluster")
# library(dplyr)
# library(dbplyr)
# library("PoiClaClu")

# BiocManager::install("ggtreeExtra")
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(enrichplot)
# library(GOSemSim)
# library(DOSE)
# library(ggtreeExtra)


## For lengthy printing
options(max.print=100000)

## htseq_counts files (from TT/RNAseq pipeline) should be in this folder for processing
parentDir = '~/htseq_counts'
# Any output files should go here
outputDir = '~/'

# 74 samples total (including 1xKD TT/RNAseq samples)
# metaData <- read.csv('htseqcounts_file_DESEQ2-4_metadata.csv', header = TRUE, sep = ",")
metaData <- read.csv('v2_htseqcounts_file_DESEQ2-4_metadata.csv', header = TRUE, sep = ",")

sampleNames <- metaData[,1]

fileNames <- c()
for (i in 1:length(sampleNames))
{
  fileNames <- append(fileNames, paste(sampleNames[i], "_gene_counts.txt", sep=""))	
}

sampleConditions <- metaData[,2]

sampleTable <- data.frame(sampleName = sampleNames, fileName = fileNames, condition = sampleConditions)

# sampleTable_24control <- sampleTable[1:10,]
sampleTable_RNA <- sampleTable[11:26,]
sampleTable_TR <- sampleTable[27:49,]
sampleTable_TT <- sampleTable[50:72,]
# sampleTable_TR.2xKD <- sampleTable[34:49,]
# sampleTable_TT.2xKD <- sampleTable[57:72,]
# sampleTable_TT.INF <- sampleTable[c(52:56,59:64,67:72),]
# sampleTable_TTALL <- sampleTable[c(3:10,50:72),]
## Just processing 72hr / 96hr TT on own
# sampleTable_TT72 <- sampleTable[65:72,]
# sampleTable_TT96 <- sampleTable[50:64,]

### Create formatted DESeq2 dataset
# dds_24control <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_24control, directory = parentDir, design= ~ condition)
dds_RNA <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_RNA, directory = parentDir, design= ~ condition)
dds_TR.all <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_TR, directory = parentDir, design= ~ condition)
dds_TT.all <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_TT, directory = parentDir, design= ~ condition)


### I need to strip the version extension (after the decimal) to generate gene names / easy annotation
# rownames(dds_24control) <- gsub("\\..*","",rownames(dds_24control))
rownames(dds_RNA) <- gsub("\\..*","",rownames(dds_RNA))
rownames(dds_TT.all) <- gsub("\\..*","",rownames(dds_TT.all))
rownames(dds_TR.all) <- gsub("\\..*","",rownames(dds_TR.all))


### Retain a copy with viral reads
# ddsVirus_24control <- dds_24control
ddsVirus_RNA <- dds_RNA
ddsVirus_TT <- dds_TT.all
ddsVirus_TR <- dds_TR.all


### Remove viral genes and only keep human genes for default dds
# dds_24control <- dds_24control[ grepl("ENSG", rownames(dds_24control)) ]
dds_RNA <- dds_RNA[ grepl("ENSG", rownames(dds_RNA)) ]
dds_TT.all <- dds_TT.all[ grepl("ENSG", rownames(dds_TT.all)) ]
dds_TR.all <- dds_TR.all[ grepl("ENSG", rownames(dds_TR.all)) ]


### Remove rows with less than sum of 1 in row
#### 11/16/23 -- changed threshold to 4 (maybe clean up data)
# This is more strict than removing lines with no counts, 
# less than previous (Remove rows with less than 10 counts in at less than 2 samples)
# dds_24control <- dds_24control[ rowSums(counts(dds_24control)) > 4, ]
dds_RNA <- dds_RNA[ rowSums(counts(dds_RNA)) > 4, ]
dds_TT.all <- dds_TT.all[ rowSums(counts(dds_TT.all)) > 4, ]
dds_TR.all <- dds_TR.all[ rowSums(counts(dds_TR.all)) > 4, ]


##----------DE analysis--------#
## Datasets = dds_24control , dds_RNA , dds_TT.2xKD , dds_TR.2xKD 
# Generate DESeq dataset object for each data set

# dds_24control <- DESeq(dds_24control)
dds_RNA <- DESeq(dds_RNA)
dds_TT.all <- DESeq(dds_TT.all)
dds_TR.all <- DESeq(dds_TR.all)

# ddsVirus_24control <- DESeq(ddsVirus_24control)
ddsVirus_RNA <- DESeq(ddsVirus_RNA)
ddsVirus_TT <- DESeq(ddsVirus_TT)
ddsVirus_TR <- DESeq(ddsVirus_TR)

##----------Normalized expression analysis--------##
## rlog transformation docs:
## https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/rlog
## https://compbiocore.github.io/deseq-workshop-1/assets/deseq_workshop_1.html
# rld_24control <- rlog(dds_24control, blind=FALSE)
rld_RNA <- rlog(dds_RNA, blind=FALSE)
rld_TT.all <- rlog(dds_TT.all, blind=FALSE)
rld_TR.all <- rlog(dds_TR.all, blind=FALSE)

##----------DE analysis--------#
## Set up contrasts
# 96hr
contrast1.TT <- c("condition", "TsiCinf96", "TsiCmock96")
# 72hr
contrast4.TT <- c("condition", "T72siCinf", "T72siCmock")
# 96hr
contrast1.R <- c("condition", "RsiCinf96", "RsiCmock96")
# 120hr
contrast4.R <- c("condition", "RsiCinf120", "RsiCmock120")
# 96hr
contrast1.TR <- c("condition", "TRsiCinf96", "TRsiCmock96")
# 72hr
contrast4.TR <- c("condition", "TR72siCinf", "TR72siCmock")

# This is the significance cutoff used for optimizing the independent filtering (by default it is set to 0.1). 
# If the adjusted p-value cutoff (FDR) will be a value other than 0.1 (for our final list of significant genes), alpha should be set to that value.
### Set thresholds
alpha_threshold <- 0.1 ## v10 script default = 0.1
padj.cutoff <- 0.05

# 96hr TT
res_C_mock96 <- results(dds_TT.all, contrast=contrast1.TT, alpha = alpha_threshold)
# 72hr TT
res_C_mock72 <- results(dds_TT.all, contrast=contrast4.TT, alpha = alpha_threshold)
# 96hr RNA
res_C_mock96.R <- results(dds_RNA, contrast=contrast1.R, alpha = alpha_threshold)
# 120hr RNA
res_C_mock120.R <- results(dds_RNA, contrast=contrast4.R, alpha = alpha_threshold)
# 96hr TR
res_C_mock96.TR <- results(dds_TR.all, contrast=contrast1.TR, alpha = alpha_threshold)
# 72hr TR
res_C_mock72.TR <- results(dds_TR.all, contrast=contrast4.TR, alpha = alpha_threshold)

################################################################################
##---------- Threshold analysis--------#
################################################################################


### 96hr TT
Sigres_C_mock96 <- subset(res_C_mock96, padj < padj.cutoff)
Sigres_UP_C_mock96 <- subset(Sigres_C_mock96, log2FoldChange > 0)
Sigres_DOWN_C_mock96 <- subset(Sigres_C_mock96, log2FoldChange < 0)
## top 1000 genes based on FC or basemean
topgenesdown96_FC <- head(Sigres_DOWN_C_mock96[order(Sigres_DOWN_C_mock96$log2FoldChange),], 1000)
topgenesup96_FC <- head(Sigres_UP_C_mock96[order(Sigres_UP_C_mock96$log2FoldChange, decreasing = TRUE),], 1000)
# topgenesdown96_baseMean <- head(Sigres_DOWN_C_mock96[order(Sigres_DOWN_C_mock96$baseMean, decreasing = TRUE),], 1000)
# topgenesup96_baseMean <- head(Sigres_UP_C_mock96[order(Sigres_UP_C_mock96$baseMean),], 1000)

### 72hr TT
Sigres_C_mock72 <- subset(res_C_mock72, padj < padj.cutoff)
Sigres_UP_C_mock72 <- subset(Sigres_C_mock72, log2FoldChange > 0)
Sigres_DOWN_C_mock72 <- subset(Sigres_C_mock72, log2FoldChange < 0)
## top 1000 genes based on FC or basemean
topgenesdown72_FC <- head(Sigres_DOWN_C_mock72[order(Sigres_DOWN_C_mock72$log2FoldChange),], 1000)
topgenesup72_FC <- head(Sigres_UP_C_mock72[order(Sigres_UP_C_mock72$log2FoldChange, decreasing = TRUE),], 1000)
# topgenesdown72_baseMean <- head(Sigres_DOWN_C_mock72[order(Sigres_DOWN_C_mock72$baseMean, decreasing = TRUE),], 1000)
# topgenesup72_baseMean <- head(Sigres_UP_C_mock72[order(Sigres_UP_C_mock72$baseMean),], 1000)

### 96hr RNA
Sigres_C_mock96.R <- subset(res_C_mock96.R, padj < padj.cutoff)
Sigres_UP_C_mock96.R <- subset(Sigres_C_mock96.R, log2FoldChange > 0)
Sigres_DOWN_C_mock96.R <- subset(Sigres_C_mock96.R, log2FoldChange < 0)

### 120hr RNA
Sigres_C_mock120.R <- subset(res_C_mock120.R, padj < padj.cutoff)
Sigres_UP_C_mock120.R <- subset(Sigres_C_mock120.R, log2FoldChange >  0)
Sigres_DOWN_C_mock120.R <- subset(Sigres_C_mock120.R, log2FoldChange <  0)

### 96hr TR
Sigres_C_mock96.TR <- subset(res_C_mock96.TR, padj < padj.cutoff)
Sigres_UP_C_mock96.TR <- subset(Sigres_C_mock96.TR, log2FoldChange >  0)
Sigres_DOWN_C_mock96.TR <- subset(Sigres_C_mock96.TR, log2FoldChange <  0)

### 72hr TR
Sigres_C_mock72.TR <- subset(res_C_mock72.TR, padj < padj.cutoff)
Sigres_UP_C_mock72.TR <- subset(Sigres_C_mock72.TR, log2FoldChange >  0)
Sigres_DOWN_C_mock72.TR <- subset(Sigres_C_mock72.TR, log2FoldChange <  0)

#########################################################################################################
############################################ Overlaps ###################################################
#########################################################################################################

### 72 / 96 overlaps
overlap.72.96_DE <- list(all72 = row.names(Sigres_C_mock72), all96 = row.names(Sigres_C_mock96))
list_to_matrix(overlap.72.96_DE)
m1_overlap.72.96_DE = make_comb_mat(overlap.72.96_DE)
# overlap.72.96_DE.genes <- extract_comb(m1_overlap.72.96_DE, "11")

overlap.72.96_DE_UD <- list(up72 = row.names(Sigres_UP_C_mock72), up96 = row.names(Sigres_UP_C_mock96), down72 = row.names(Sigres_DOWN_C_mock72), down96 = row.names(Sigres_DOWN_C_mock96))
list_to_matrix(overlap.72.96_DE_UD)
m1_overlap.72.96_DE_UD = make_comb_mat(overlap.72.96_DE_UD)
overlap.72.96_DE_UP.genes <- extract_comb(m1_overlap.72.96_DE_UD, "1100")
overlap.72.96_DE_DOWN.genes <- extract_comb(m1_overlap.72.96_DE_UD, "0011")

## RNA
overlap.72.96_DE.TR <- list(all72 = row.names(Sigres_C_mock72.TR), all96 = row.names(Sigres_C_mock96.TR))
list_to_matrix(overlap.72.96_DE.TR)
m1_overlap.72.96_DE.TR = make_comb_mat(overlap.72.96_DE.TR)
overlap.72.96_DE.TR.genes <- extract_comb(m1_overlap.72.96_DE.TR, "11")

overlap.72.96_DE_UD.TR <- list(up72 = row.names(Sigres_UP_C_mock72.TR), up96 = row.names(Sigres_UP_C_mock96.TR), down72 = row.names(Sigres_DOWN_C_mock72.TR), down96 = row.names(Sigres_DOWN_C_mock96.TR))
list_to_matrix(overlap.72.96_DE_UD.TR)
m1_overlap.72.96_DE_UD.TR = make_comb_mat(overlap.72.96_DE_UD.TR)
overlap.72.96_DE_UP.TR.genes <- extract_comb(m1_overlap.72.96_DE_UD.TR, "1100")
overlap.72.96_DE_DOWN.TR.genes <- extract_comb(m1_overlap.72.96_DE_UD.TR, "0011")

### TT vs RNA
overlap.72.96_DE_UD.TT.TR <- list(upTT = overlap.72.96_DE_UP.genes, upRNA = overlap.72.96_DE_UP.TR.genes, downTT = overlap.72.96_DE_DOWN.genes, downRNA = overlap.72.96_DE_DOWN.TR.genes)
list_to_matrix(overlap.72.96_DE_UD.TT.TR)
m1_overlap.72.96_DE_UD.TT.TR = make_comb_mat(overlap.72.96_DE_UD.TT.TR)
onlyUPinTT <- extract_comb(m1_overlap.72.96_DE_UD.TT.TR, "1000")
upinBOTH.TTRNA <- extract_comb(m1_overlap.72.96_DE_UD.TT.TR, "1100")
onlyDOWNinTT <- extract_comb(m1_overlap.72.96_DE_UD.TT.TR, "0010")
downinBOTH.TTRNA <- extract_comb(m1_overlap.72.96_DE_UD.TT.TR, "0011")


#### RNA 96/120 different prep -- RNAseq only
overlap.120.96_DE_UD.TR <- list(up120 = row.names(Sigres_UP_C_mock120.R), up96 = row.names(Sigres_UP_C_mock96.R), down120 = row.names(Sigres_DOWN_C_mock120.R), down96 = row.names(Sigres_DOWN_C_mock96.R))
list_to_matrix(overlap.120.96_DE_UD.TR)
m1_overlap.120.96_DE_UD.TR = make_comb_mat(overlap.120.96_DE_UD.TR)
overlap.120.96_DE_UP.TR.genes <- extract_comb(m1_overlap.120.96_DE_UD.TR, "1100")
overlap.120.96_DE_DOWN.TR.genes <- extract_comb(m1_overlap.120.96_DE_UD.TR, "0011")

#########################################################################################################
############################################ KD DE genes ################################################
#########################################################################################################

#### DE with KD
contrast1.TTKDS <- c("condition", "TsiSinf96", "TsiCinf96")
contrast1.TTKDA <- c("condition", "TsiAinf96", "TsiCinf96")
res_C_KDS_96 <- results(dds_TT.all, contrast=contrast1.TTKDS, alpha = alpha_threshold)
res_C_KDA_96 <- results(dds_TT.all, contrast=contrast1.TTKDA, alpha = alpha_threshold)
Sigres_C_KDS_96 <- subset(res_C_KDS_96, pvalue < 0.1)
Sigres_C_KDA_96 <- subset(res_C_KDA_96, pvalue < 0.1)
SigresUP_C_KDS_96 <- subset(Sigres_C_KDS_96, log2FoldChange > 0)
SigresDOWN_C_KDS_96 <- subset(Sigres_C_KDS_96, log2FoldChange < 0)
SigresUP_C_KDA_96 <- subset(Sigres_C_KDA_96, log2FoldChange > 0)
SigresDOWN_C_KDA_96 <- subset(Sigres_C_KDA_96, log2FoldChange < 0)

contrast2.TTKDS <- c("condition", "T72siSinf", "T72siCinf")
contrast2.TTKDA <- c("condition", "T72siAinf", "T72siCinf")
res_C_KDS_72 <- results(dds_TT.all, contrast=contrast2.TTKDS, alpha = alpha_threshold)
res_C_KDA_72 <- results(dds_TT.all, contrast=contrast2.TTKDA, alpha = alpha_threshold)
Sigres_C_KDS_72 <- subset(res_C_KDS_72, pvalue < 0.1)
Sigres_C_KDA_72 <- subset(res_C_KDA_72, pvalue < 0.1)
SigresUP_C_KDS_72 <- subset(Sigres_C_KDS_72, log2FoldChange > 0)
SigresDOWN_C_KDS_72 <- subset(Sigres_C_KDS_72, log2FoldChange < 0)
SigresUP_C_KDA_72 <- subset(Sigres_C_KDA_72, log2FoldChange > 0)
SigresDOWN_C_KDA_72 <- subset(Sigres_C_KDA_72, log2FoldChange < 0)

KDup.TT <- list(siS.72 = rownames(SigresUP_C_KDS_72), siA.72 = rownames(SigresUP_C_KDA_72), siS.96 = rownames(SigresUP_C_KDS_96), siA.96 = rownames(SigresUP_C_KDA_96))
list_to_matrix(KDup.TT)
m1_KDup.TT = make_comb_mat(KDup.TT, mode = "union") 
allKDup.TT <- extract_comb(m1_KDup.TT, "1111")
SUN1.KDup.TT <- extract_comb(m1_KDup.TT, "1010")
ATAT1.KDup.TT <- extract_comb(m1_KDup.TT, "0101")

KDdown.TT <- list(siS.72 = rownames(SigresDOWN_C_KDS_72), siA.72 = rownames(SigresDOWN_C_KDA_72), siS.96 = rownames(SigresDOWN_C_KDS_96), siA.96 = rownames(SigresDOWN_C_KDA_96))
list_to_matrix(KDdown.TT)
m1_KDdown.TT = make_comb_mat(KDdown.TT, mode = "union") 
allKDdown.TT <- extract_comb(m1_KDdown.TT, "1111")
SUN1.KDdown.TT <- extract_comb(m1_KDdown.TT, "1010")
ATAT1.KDdown.TT <- extract_comb(m1_KDdown.TT, "0101")

KD.TT <- list(siSup.72 = rownames(SigresUP_C_KDS_72), siAup.72 = rownames(SigresUP_C_KDA_72), siSup.96 = rownames(SigresUP_C_KDS_96), siAup.96 = rownames(SigresUP_C_KDA_96), siSdown.72 = rownames(SigresDOWN_C_KDS_72), siAdown.72 = rownames(SigresDOWN_C_KDA_72), siSdown.96 = rownames(SigresDOWN_C_KDS_96), siAdown.96 = rownames(SigresDOWN_C_KDA_96))
# list_to_matrix(KD.TT)
# m1_KD.TT = make_comb_mat(KD.TT, mode = "union") 
# allKDall.TT <- extract_comb(m1_KDdown.TT, "11111111")


#########################################################################################################
############################################ Enhanced Volcano ###########################################
#########################################################################################################

vol_res_res_96TTwKD <- EnhancedVolcano(res_C_mock96,
                                    lab = rownames(res_C_mock96),
                                    x = 'log2FoldChange',
                                    y = 'pvalue',
                                    pCutoff = 0.05,
                                    FCcutoff = 2.0,
                                    labSize = 0.0)

#########################################################################################################
############################################ GO terms ###################################################
#########################################################################################################

##----------ClusterProfiler--------#
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html
# https://www.sciencedirect.com/science/article/pii/S2666675821000667#appsec1
# Gene universe file
# All "possible" genes in dataset (after removal of low counts)
# bg_genes.24 = row.names(dds_24control)
bg_genes.TT = row.names(dds_TT.all)
bg_genes.RNA = row.names(dds_RNA)
bg_genes.TR = row.names(dds_TR.all)

######## for single samples
# GO.list <- XXX
# GO.enrich <- enrichGO(gene      = GO.list,
#                   universe      = bg_genes.TT,
#                   OrgDb         = org.Hs.eg.db,
#                   ont           = "ALL",
#                   pAdjustMethod = "BH",
#                   keyType       = "ENSEMBL",
#                   pvalueCutoff  = 0.05,
#                   qvalueCutoff  = 0.01,
#                   readable      = TRUE, #,
#                   minGSSize = 2,
#                   maxGSSize = 500, pool = TRUE)
# head(GO.enrich)
# dotplot(GO.enrich)
# GO.enrich.tree <- pairwise_termsim(GO.enrich)
# treeplot(GO.enrich.tree, hclust_method = "centroid", nWords = 0, label_format = 15, showCategory = 10, group_color = c("#999999", "#483D8B", "#4169E1", "#ADD8E6", "#4169E1"))

## for simplification of terms
# GO.enrich.s <- clusterProfiler::simplify(GO.enrich,
#                                          cutoff = 0.7,
#                                          by = "p.adjust",
#                                          select_fun = min,
#                                          measure = "Wang",
#                                          semData = NULL)
# dotplot(GO.enrich.s)
# GO.enrich.s.tree <- pairwise_termsim(GO.enrich.s)
# treeplot(GO.enrich.s.tree, hclust_method = "centroid", nWords = 0, label_format = 15, showCategory = 10, group_color = c("#999999", "#483D8B", "#4169E1", "#ADD8E6", "#4169E1"))


######## for comparisons
compare.GO.list = KD.TT
compare.GO <- compareCluster(geneCluster = compare.GO.list, fun = enrichGO,
                       readable = TRUE, 
                       OrgDb = org.Hs.eg.db, 
                       keyType="ENSEMBL",
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       universe      = bg_genes.TT,
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05, pool = TRUE)
head(compare.GO)
dotplot(compare.GO, showCategory = 10, size = "Count")

## for simplifying terms
compare.GO.S <- clusterProfiler::simplify(compare.GO,
                    cutoff = 0.7,
                    by = "p.adjust",
                    select_fun = min,
                    measure = "Wang",
                    semData = NULL)


head(compare.GO.S)
dotplot(compare.GO.S)


#########################################################################################################
############################################ Export to BED file #########################################
#########################################################################################################

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
# mart <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl', mirror = "uswest", GRCh = 37)
mart <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl')
# head(listAttributes(mart), 20)

geneList <- listX
geneList.tmp <- biomaRt::getBM(attributes=c('chromosome_name','start_position','end_position','strand','ensembl_gene_id','external_gene_name'),
                               filters='ensembl_gene_id', values = geneList, mart = mart)
#create a table that resembles a BED file
geneList.bed <- data.frame(chr=geneList.tmp$chromosome_name,
                           start=geneList.tmp$start_position,
                           end=geneList.tmp$end_position,
                           id=paste(geneList.tmp$ensembl_gene_id, geneList.tmp$external_gene_name, sep="/"),
                           score=rep(0,nrow(geneList.tmp)),
                           strand=geneList.tmp$strand)
# perform some cosmetics
geneList.bed <- geneList.bed[order(geneList.bed[,'chr'],geneList.bed[,'start']),]
geneList.bed$chr <- paste('chr', geneList.bed$chr, sep = '')
geneList.bed$strand <- gsub(pattern = "-1",replacement = "-",x = geneList.bed$strand)
geneList.bed$strand <- gsub(pattern = "1",replacement = "+",x = geneList.bed$strand)
# write bed file
write.table(x = geneList.bed,
            file = 'listX.bed',
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)


