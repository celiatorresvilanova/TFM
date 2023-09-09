###########################################################################################################################################################
                                                     # Selection of a dataset for further analysis                                                                                                                              
###########################################################################################################################################################

# Load needed packages

library(biomaRt)
library(limma)
library(edgeR)
library(DESeq2)
library(sva)
library(DT)
library(GSA)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GEOquery)
library(org.Hs.eg.db)
library(GSEABase)
library(GOstats)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# 1. Load data and create an Expression Set
## 1.1 Original data 

setwd("C:/Users/Cèlia Torres Vilanov/Documents/SAINT")
load(file="assay_data_orig.RData")
load(file="phenodata.RData")

problematic_samples <- c("S9_2", "S16_1", "S22_1", "S18_1")
assay_data_orig_wp <- assay_data_orig [, !colnames(assay_data_orig) %in% problematic_samples]

rownames(phenodata) <- colnames(assay_data_orig)
phenodata_wp <- phenodata[!(rownames(phenodata) %in% problematic_samples), ]
ExpressionSet_orig_wp <- ExpressionSet(assayData = as.matrix(assay_data_orig_wp), phenoData = AnnotatedDataFrame(phenodata_wp))


## 1.2 Data excluding lowly expressed genes that have high duplication ratios 

load(file="assay_data_filt_genes.RData")
load(file="dm_final.RData")

intersection.wp <- intersect(dm_final$ID, rownames(assay_data_orig_wp))
assay_data_filt_genes_wp <- assay_data_orig_wp[intersection.wp,]
dim(assay_data_filt_genes_wp)

ExpressionSet_filt_genes_wp <- ExpressionSet(assayData = as.matrix(assay_data_filt_genes_wp), phenoData = AnnotatedDataFrame(phenodata_wp))

# 2. Filtration
## 2.1. Filtration by low expressed genes
### 2.1.1 Original data excluding problematic samples

dim(assay_data_orig_wp) # 62,757 genes
keep <- rowSums(assay_data_orig_wp > 10) >= 11 # at least 11 samples have 10 reads per gene
assay_data_orig_wp_f <- assay_data_orig_wp[keep,]
dim(assay_data_orig_wp_f)  # 20,697 genes

### 2.1.2 Data excluding lowly expressed genes that have high duplication ratios and problematic samples

dim(assay_data_filt_genes_wp) #  47,173 genes
keep <- rowSums(assay_data_filt_genes_wp > 10) >= 11 # at least 11 samples have 10 reads per gene
assay_data_filt_genes_wp_f <- assay_data_filt_genes_wp[keep,]
dim(assay_data_filt_genes_wp_f)  # 19,968 genes

## 2.2 Filtration of transcripts that translate to globin domains and chains
### 2.2.1 Original data excluding problematic samples

globin_transcripts <- c("ENSG00000086506", "ENSG00000161544","ENSG00000165553", "ENSG00000188536", "ENSG00000206172", "ENSG00000206177", "ENSG00000213931", "ENSG00000213934", "ENSG00000223609", "ENSG00000244734", "ENSG00000260629","ENSG00000196565")
rownames(assay_data_filt_genes_wp_f) <- gsub ("\\..*", "", rownames(assay_data_filt_genes_wp_f))

assay_data_filt_genes_wp_f <- assay_data_filt_genes_wp_f[!(rownames(assay_data_filt_genes_wp_f) %in% globin_transcripts), ]
dim(assay_data_filt_genes_wp_f) # 19,961 genes

# 3. Differential Expression Analysis using DESeq
## 3.1 Baseline of the ivermectin and placebo groups
### 3.1.1 Original data

load(file="ExpressionSet_orig.RData")
load(file="assay_data_orig_f.RData")

baseline.orig <- ExpressionSet_orig[,ExpressionSet_orig$day==1]
colbase.orig <- pData(baseline.orig)
countsbase.orig <- assay_data_orig_f[,rownames(colbase.orig)]
countsbasem.orig <- as.matrix(countsbase.orig)

colbase.orig$treat <- as.factor(colbase.orig$treat)
DESeq.Orig.1 <- DESeqDataSetFromMatrix(countData = countsbasem.orig,
                              colData = colbase.orig,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.Orig.1 <- DESeq(DESeq.Orig.1)

res.orig.1 <- results(DESeq.Orig.1, contrast = c("treat", '1', '2'))  # treat = treatment (predictor variable)

res.DESeq.orig.1 <- as.data.frame(res.orig.1)
res.DESeq.orig.sign.1 <- res.DESeq.orig.1[!is.na(res.DESeq.orig.1$padj) & res.DESeq.orig.1$padj<0.3 & abs(res.DESeq.orig.1$log2FoldChange) > log2(2), ]
dim(res.DESeq.orig.sign.1)  # 6 Differentially Expressed Genes (DEGs) between ivermectin-treated and placebo groups at baseline from the original data
head(res.DESeq.orig.sign.1) 

### 3.1.2 Original data excluding problematic samples

baseline.orig.wp <- ExpressionSet_orig_wp[,ExpressionSet_orig_wp$day==1]
colbase.orig.wp <- pData(baseline.orig.wp)
countsbase.orig.wp <- assay_data_orig_wp_f[,rownames(colbase.orig.wp)]
countsbasem.orig.wp <- as.matrix(countsbase.orig.wp)

colbase.orig.wp$treat<- as.factor(colbase.orig.wp$treat)
DESeq.Orig.wp.1 <- DESeqDataSetFromMatrix(countData = countsbasem.orig.wp,
                              colData = colbase.orig.wp,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.Orig.wp.1 <- DESeq(DESeq.Orig.wp.1)

res.orig.wp.1 <- results(DESeq.Orig.wp.1, contrast = c("treat", '1', '2')) # treat = treatment (predictor variable)

res.DESeq.orig.wp.1 <- as.data.frame(res.orig.wp.1)
res.DESeq.orig.wp.1.sign <- res.DESeq.orig.wp.1[!is.na(res.DESeq.orig.wp.1$padj) & res.DESeq.orig.wp.1$padj<0.3 & abs(res.DESeq.orig.wp.1$log2FoldChange) > log2(2), ]
dim(res.DESeq.orig.wp.1.sign) # 5 DEGs between ivermectin-treated and placebo groups at baseline from the original data excluding problematic samples
head(res.DESeq.orig.wp.1.sign)

### 3.1.3. Original data excluding lowly expressed genes that have high ratio of duplicates

load(file="ExpressionSet_filt_genes.RData")
load(file="assay_data_filt_genes_f.RData")

baseline.filt.genes <- ExpressionSet_filt_genes[,ExpressionSet_filt_genes$day==1]
colbase.filt.genes <- pData(baseline.filt.genes)
countsbase.filt.genes <- assay_data_filt_genes_f[,rownames(colbase.filt.genes)]
countsbasem.filt.genes <- as.matrix(countsbase.filt.genes)

colbase.filt.genes$treat <- as.factor(colbase.filt.genes$treat) # treat = treatment (predictor variable)
DESeq.filt.genes.1 <- DESeqDataSetFromMatrix(countData = countsbasem.filt.genes,
                              colData = colbase.filt.genes,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.filt.genes.1 <- DESeq(DESeq.filt.genes.1)

res.filt.genes.1 <- results(DESeq.filt.genes.1, contrast = c("treat", '1', '2')) 

res.DESeq.filt.genes.1 <- as.data.frame(res.filt.genes.1)
res.DESeq.filt.genes.sign.1 <- res.DESeq.filt.genes.1[!is.na(res.DESeq.filt.genes.1$padj) & res.DESeq.filt.genes.1$padj<0.3 & abs(res.DESeq.filt.genes.1$log2FoldChange) > log2(2), ]
dim(res.DESeq.filt.genes.sign.1) # 6 DEGs between ivermectin-treated and placebo groups at baseline from the original data excluding lowly expressed genes that have high ratio of duplicates
head(res.DESeq.filt.genes.sign.1)

### 3.1.4 Original data excluding lowly expressed genes that have high ratio of duplicates and problematic samples

baseline.filt.genes.wp <- ExpressionSet_filt_genes_wp[,ExpressionSet_filt_genes_wp$day==1]
colbase.filt.genes.wp <- pData(baseline.filt.genes.wp)
countsbase.filt.genes.wp <- assay_data_filt_genes_wp_f[,rownames(colbase.filt.genes.wp)]
countsbasem.filt.genes.wp <- as.matrix(countsbase.filt.genes.wp)

colbase.filt.genes.wp$treat <- as.factor(colbase.filt.genes.wp$treat) # treat = treatment (predictor variable)
DESeq.filt.genes.wp.1 <- DESeqDataSetFromMatrix(countData = countsbasem.filt.genes.wp,
                              colData = colbase.filt.genes.wp,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.filt.genes.wp.1 <- DESeq(DESeq.filt.genes.wp.1)

res.filt.genes.wp.1 <- results(DESeq.filt.genes.wp.1, contrast = c("treat", '1', '2'))  # treat = treatment (predictor variable)

res.DESeq.filt.genes.wp.1 <- as.data.frame(res.filt.genes.wp.1)
res.DESeq.filt.genes.sign.wp.1 <- res.DESeq.filt.genes.wp.1[!is.na(res.DESeq.filt.genes.wp.1$padj) & res.DESeq.filt.genes.wp.1$padj<0.3 & abs(res.DESeq.filt.genes.wp.1$log2FoldChange) > log2(2), ]
dim(res.DESeq.filt.genes.sign.wp.1) # 5 DEGs between ivermectin-treated and placebo groups  at baseline from the original data excluding lowly expressed genes that have high ratio of duplicates and problematic samples
head(res.DESeq.filt.genes.sign.wp.1)

## 3.2 Differences at timepoints 4 and 7 between ivermectin-treated and placebo groups
### 3.2.1 Original data
# Differences at timepoint 4 between ivermectin-treated and placebo groups

orig.4 <- ExpressionSet_orig[,ExpressionSet_orig$day==4]
col.orig.4 <- pData(orig.4)
counts.orig.4 <- assay_data_orig_f[,rownames(col.orig.4)]
countsm.orig.4 <- as.matrix(counts.orig.4)

col.orig.4$treat <- as.factor(col.orig.4$treat) # treat = treatment (predictor variable)
DESeq.Orig.4 <- DESeqDataSetFromMatrix(countData = countsm.orig.4,
                              colData = col.orig.4,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.Orig.4 <- DESeq(DESeq.Orig.4)

res.orig.4 <- results(DESeq.Orig.4, contrast = c("treat", '1', '2')) # treat = treatment (predictor variable)

res.DESeq.orig.4 <- as.data.frame(res.orig.4)
res.DESeq.orig.sign.4 <- res.DESeq.orig.4[!is.na(res.DESeq.orig.4$padj) & res.DESeq.orig.4$padj<0.3 & abs(res.DESeq.orig.4$log2FoldChange) > log2(2), ]
dim(res.DESeq.orig.sign.4)  # 5 DEGs between ivermectin-treated and placebo groups at timepoint 4 from the original data
head(res.DESeq.orig.sign.4)


# Differences at timepoint 7 between ivermectin-treated and placebo groups

orig.7 <- ExpressionSet_orig[,ExpressionSet_orig$day==7]
col.orig.7 <- pData(orig.7)
counts.orig.7 <- assay_data_orig_f[,rownames(col.orig.7)]
countsm.orig.7 <- as.matrix(counts.orig.7)

col.orig.7$treat <- as.factor(col.orig.7$treat) # treat = treatment (predictor variable)
DESeq.Orig.7 <- DESeqDataSetFromMatrix(countData = countsm.orig.7,
                              colData = col.orig.7,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.Orig.7 <- DESeq(DESeq.Orig.7)

res.orig.7 <- results(DESeq.Orig.7, contrast = c("treat", '1', '2')) # treat = treatment (predictor variable)

res.DESeq.orig.7 <- as.data.frame(res.orig.7)
res.DESeq.orig.sign.7 <- res.DESeq.orig.7[!is.na(res.DESeq.orig.7$padj) & res.DESeq.orig.7$padj<0.3 & abs(res.DESeq.orig.7$log2FoldChange) > log2(2), ]
dim(res.DESeq.orig.sign.7) # 29 DEGs between ivermectin-treated and placebo groups at timepoint 7 from the original data
head(res.DESeq.orig.sign.7)

### 3.2.2 Original data excluding problematic samples
# Differences at timepoint 4 between ivermectin-treated and placebo groups

orig.wp.4 <- ExpressionSet_orig_wp[,ExpressionSet_orig_wp$day==4]
col.orig.wp.4 <- pData(orig.wp.4)
counts.orig.wp.4 <- assay_data_orig_wp_f[,rownames(col.orig.wp.4)]
countsm.orig.wp.4 <- as.matrix(counts.orig.wp.4)

col.orig.wp.4$treat <- as.factor(col.orig.wp.4$treat) # treat = treatment (predictor variable)
DESeq.Orig.wp.4 <- DESeqDataSetFromMatrix(countData = countsm.orig.wp.4,
                              colData = col.orig.wp.4,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.Orig.wp.4 <- DESeq(DESeq.Orig.wp.4)

res.orig.wp.4 <- results(DESeq.Orig.wp.4, contrast = c("treat", '1', '2')) # treat = treatment (predictor variable)

res.DESeq.orig.wp.4 <- as.data.frame(res.orig.wp.4)
res.DESeq.orig.wp.sign.4 <- res.DESeq.orig.wp.4[!is.na(res.DESeq.orig.wp.4$padj) & res.DESeq.orig.wp.4$padj<0.3 & abs(res.DESeq.orig.wp.4$log2FoldChange) > log2(2), ]
dim(res.DESeq.orig.wp.sign.4) # 3 DEGs between ivermectin-treated and placebo groups at timepoint 4 from the original data excluding problematic samples
head(res.DESeq.orig.wp.sign.4)

# Timepoint 7 doesn't have problematic samples, so I won't analyse them, because the results are the same than in the original data

### 3.2.3 Original data excluding lowly expressed genes that have high ratio of duplicates

# Differences at timepoint 4 between ivermectin-treated and placebo groups

filt.genes.4 <- ExpressionSet_filt_genes[,ExpressionSet_filt_genes$day==4]
col.filt.genes.4 <- pData(filt.genes.4)
counts.filt.genes.4 <- assay_data_filt_genes_f[,rownames(col.filt.genes.4)]
countsm.filt.genes.4 <- as.matrix(counts.filt.genes.4)

col.filt.genes.4$treat <- as.factor(col.filt.genes.4$treat) # treat = treatment (predictor variable)
DESeq.filt.genes.4 <- DESeqDataSetFromMatrix(countData = countsm.filt.genes.4,
                              colData = col.filt.genes.4,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.filt.genes.4 <- DESeq(DESeq.filt.genes.4)

res.filt.genes.4 <- results(DESeq.filt.genes.4, contrast = c("treat", '1', '2')) # treat = treatment (predictor variable)

res.DESeq.filt.genes.4 <- as.data.frame(res.filt.genes.4)
res.DESeq.filt.genes.sign.4 <- res.DESeq.filt.genes.4[!is.na(res.DESeq.filt.genes.4$padj) & res.DESeq.filt.genes.4$padj<0.3 & abs(res.DESeq.filt.genes.4$log2FoldChange) > log2(2), ]
dim(res.DESeq.filt.genes.sign.4) # 5 DEGs between ivermectin-treated and placebo groups at timepoint 4 from the original data excluding lowly expressed genes that have high ratio of duplicates
head(res.DESeq.filt.genes.sign.4)

# Differences at timepoint 7 between ivermectin-treated and placebo groups

filt.genes.7 <- ExpressionSet_filt_genes[,ExpressionSet_filt_genes$day==7]
col.filt.genes.7 <- pData(filt.genes.7)
counts.filt.genes.7 <- assay_data_filt_genes_f[,rownames(col.filt.genes.7)]
countsm.filt.genes.7 <- as.matrix(counts.filt.genes.7)

col.filt.genes.7$treat <- as.factor(col.filt.genes.7$treat) # treat = treatment (predictor variable)
DESeq.filt.genes.7 <- DESeqDataSetFromMatrix(countData = countsm.filt.genes.7,
                              colData = col.filt.genes.7,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.filt.genes.7 <- DESeq(DESeq.filt.genes.7)

res.filt.genes.7 <- results(DESeq.filt.genes.7, contrast = c("treat", '1', '2')) # treat = treatment (predictor variable)

res.DESeq.filt.genes.7 <- as.data.frame(res.filt.genes.7)
res.DESeq.filt.genes.sign.7 <- res.DESeq.filt.genes.7[!is.na(res.DESeq.filt.genes.7$padj) & res.DESeq.filt.genes.7$padj<0.3 & abs(res.DESeq.filt.genes.7$log2FoldChange) > log2(2), ]
dim(res.DESeq.filt.genes.sign.7) # 27 DEGs between ivermectin-treated and placebo groups at timepoint 7 from the original data excluding lowly expressed genes that have high ratio of duplicates
head(res.DESeq.filt.genes.sign.7)

### 3.2.4 Original data excluding lowly expressed genes that have high ratio of duplicates and problematic samples

# Differences at timepoint 4 between ivermectin-treated and placebo groups

filt.genes.wp.4 <- ExpressionSet_filt_genes_wp[,ExpressionSet_filt_genes_wp$day==4]
col.filt.genes.wp.4 <- pData(filt.genes.wp.4)
counts.filt.genes.wp.4 <- assay_data_filt_genes_wp_f[,rownames(col.filt.genes.wp.4)]
countsm.filt.genes.wp.4 <- as.matrix(counts.filt.genes.wp.4)

col.filt.genes.wp.4$treat <- as.factor(col.filt.genes.wp.4$treat) # treat = treatment (predictor variable)
DESeq.filt.genes.wp.4 <- DESeqDataSetFromMatrix(countData = countsm.filt.genes.wp.4,
                              colData = col.filt.genes.wp.4,
                              design = ~ treat) # treat = treatment (predictor variable)

DESeq.filt.genes.wp.4 <- DESeq(DESeq.filt.genes.wp.4)

res.filt.genes.wp.4 <- results(DESeq.filt.genes.wp.4, contrast = c("treat", '1', '2')) # treat = treatment (predictor variable)

res.DESeq.filt.genes.wp.4 <- as.data.frame(res.filt.genes.wp.4)
res.DESeq.filt.genes.wp.sign.4 <- res.DESeq.filt.genes.wp.4[!is.na(res.DESeq.filt.genes.wp.4$padj) & res.DESeq.filt.genes.wp.4$padj<0.3 & abs(res.DESeq.filt.genes.wp.4$log2FoldChange) > log2(2), ]
dim(res.DESeq.filt.genes.wp.sign.4) # 3 DEGs between ivermectin-treated and placebo groups at timepoint 4 from the original data excluding lowly expressed genes that have high ratio of duplicates and problematic samples
head(res.DESeq.filt.genes.wp.sign.4)

# Timepoint 7 doesn't have problematic samples, so I won't analyse them, because the results are the same using data excluding lowly expressed genes that have high duplicates ratio

# 4. Differentially Expression Analysis using limma (voom transformation + bayesian adjustment)
## 4.1 Differences at baseline between ivermectin-treated and placebo groups
### 4.1.1 Original data

treat <- as.factor(ifelse(colbase.orig$treat==1, "Ivermectin", "Placebo")) # treat = treatment (predictor variable)
design.treat <- model.matrix(~0 + treat) # treat = treatment (predictor variable)
voom.baseline.orig <- voom(countsbasem.orig, design.treat, plot = F) 

fit.baseline.orig <- lmFit(voom.baseline.orig, design.treat)

contm.treat <- makeContrasts(IvervsPlac=treatIvermectin-treatPlacebo,
                                 levels = design.treat)

fit.baseline.orig <- contrasts.fit(fit.baseline.orig, contm.treat)
fite.baseline.orig <- eBayes(fit.baseline.orig)

summary(decideTests(fite.baseline.orig, methode = "separate"))
summary(decideTests(fite.baseline.orig, methode = "separate", adjust.method = "none"))

table.baseline.orig.adj <- topTable(fite.baseline.orig, number = Inf, adjust = "fdr")

res.limma.orig.sign <- table.baseline.orig.adj[table.baseline.orig.adj$adj.P.Val < 0.3, ]
res.limma.orig.sign # 0 DEGs between ivermectin-treated and placebo groups at baseline from the original data

### 4.1.2 Original data excluding problematic samples

treat.wp <- ifelse(colbase.orig.wp$treat==1, "Ivermectin", "Placebo") # treat = treatment (predictor variable)
design.treat.wp <- model.matrix(~0 + treat.wp) # treat = treatment (predictor variable)
voom.baseline.orig.wp <- voom(countsbasem.orig.wp, design.treat.wp, plot = F) 

fit.baseline.orig.wp <- lmFit(voom.baseline.orig.wp, design.treat.wp)

contm.treat.wp <- makeContrasts(IvervsPlac=treat.wpIvermectin-treat.wpPlacebo,
                                 levels = design.treat.wp)

fit.baseline.orig.wp <- contrasts.fit(fit.baseline.orig.wp, contm.treat.wp)
fite.baseline.orig.wp <- eBayes(fit.baseline.orig.wp)

summary(decideTests(fite.baseline.orig.wp, methode = "separate"))
summary(decideTests(fite.baseline.orig.wp, methode = "separate", adjust.method = "none"))

table.baseline.orig.wp.adj <- topTable(fite.baseline.orig.wp, number = Inf, adjust = "fdr")

res.limma.orig.sign.wp <- table.baseline.orig.wp.adj[table.baseline.orig.wp.adj$adj.P.Val < 0.3, ]
res.limma.orig.sign.wp # 0 DEGs between ivermectin-treated and placebo groups at baseline from the original data excluding problematic samples

### 4.1.3 Original data excluding lowly expressed genes that have high ratio of duplicates

voom.baseline.filt.genes <- voom(countsbasem.filt.genes, design.treat, plot = F) 

fit.baseline.filt.genes <- lmFit(voom.baseline.filt.genes, design.treat)

fit.baseline.filt.genes <- contrasts.fit(fit.baseline.filt.genes, contm.treat)
fite.baseline.filt.genes <- eBayes(fit.baseline.filt.genes)

summary(decideTests(fite.baseline.filt.genes, methode = "separate"))
summary(decideTests(fite.baseline.filt.genes, methode = "separate", adjust.method = "none"))

table.baseline.filt.genes.adj <- topTable(fite.baseline.filt.genes, number = Inf, adjust = "fdr")

res.limma.filt.genes.sign <- table.baseline.filt.genes.adj[table.baseline.filt.genes.adj$adj.P.Val < 0.3,]
res.limma.filt.genes.sign # 0 DEGs between ivermectin-treated and placebo groups at baseline from the data excluding lowly expressed genes that have high duplication ratio

### 4.1.4 Original data excluding lowly expressed genes that have high ratio of duplicates and problematic samples

voom.baseline.filt.genes.wp <- voom(countsbasem.filt.genes.wp, design.treat.wp, plot = F) 
fit.baseline.filt.genes.wp <- lmFit(voom.baseline.filt.genes.wp, design.treat.wp)

fit.baseline.filt.genes.wp <- contrasts.fit(fit.baseline.filt.genes.wp, contm.treat.wp)
fite.baseline.filt.genes.wp <- eBayes(fit.baseline.filt.genes.wp)

summary(decideTests(fite.baseline.filt.genes.wp, methode = "separate"))
summary(decideTests(fite.baseline.filt.genes.wp, methode = "separate", adjust.method = "none"))

table.baseline.filt.genes.wp.adj <- topTable(fite.baseline.filt.genes.wp, number = Inf, adjust = "fdr")

res.limma.filt.genes.sign.wp <- table.baseline.filt.genes.wp.adj[table.baseline.filt.genes.wp.adj$adj.P.Val < 0.3,]
res.limma.filt.genes.sign.wp # 0 DEGs between ivermectin-treated and placebo groups at baseline from the data excluding lowly expressed genes that have high duplication ratio and problematic samples

# No Enrichment Analysis with GO terms was per-formed, as the differentially expressed genes were very few using the filtration criterion described ear-lier. 

# 5 Enrichment Analysis: Limma (voom + camera) using BTMs
## 5.1 Differences at baseline between ivermectin-treated and placebo groups
### 5.1.1 Original data

load(file= "btm.RData")
names(btm$genesets) <- btm$geneset.names

geneIDsbaseorig <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.baseline.orig), columns = c("SYMBOL","GENEID"))
idxbaseorig <- ids2indices(btm$genesets,id=geneIDsbaseorig$SYMBOL)

cam.orig.1 <- camera(countsbasem.orig,idxbaseorig,design.treat,contrast=contm.treat)

res.FDR.orig.1 <- cam.orig.1[cam.orig.1$FDR < 0.05 & cam.orig.1$NGenes > 10,]
dim(res.FDR.orig.1) # 32 Blood Transcriptional Modules (BTMs) were enriched at baseline of ivermectin-treated or placebo groups from the original data 

res.FDR.orig.1<-res.FDR.orig.1%>%dplyr::filter(!str_detect(rownames(res.FDR.orig.1), "TBA"))  # 7 BTMs To Be Annotated (TBA)

dim(res.FDR.orig.1) # 25 BTMs were enriched at baseline of ivermectin-treated or placebo groups from the original data excluding BTMs TBA
head(res.FDR.orig.1)

### 5.1.2 Original data exluding duplicates

row_names_orig.wp <- rownames(voom.baseline.orig.wp$E)
new_row_names_orig.wp <- sub("\\.\\d+", "", row_names_orig.wp)
rownames(voom.baseline.orig.wp) <- new_row_names_orig.wp

geneIDsbaseorig.wp <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.baseline.orig.wp), columns = c("SYMBOL","GENEID"))
idxbaseorig.wp <- ids2indices(btm$genesets,id=geneIDsbaseorig.wp$SYMBOL)

cam.orig.wp.1 <- camera(countsbasem.orig.wp,idxbaseorig.wp,design.treat.wp,contrast=contm.treat.wp)

res.FDR.orig.wp.1 <- cam.orig.wp.1[cam.orig.wp.1$FDR < 0.05 & cam.orig.wp.1$NGenes > 10,]
dim(res.FDR.orig.wp.1) # 5 BTMs were enriched at baseline of ivermectin-treated or placebo groups from the original data excluding problematic samples
res.FDR.orig.wp.1<-res.FDR.orig.wp.1%>%dplyr::filter(!str_detect(rownames(res.FDR.orig.wp.1), "TBA"))  # 0 BTMs TBA

dim(res.FDR.orig.wp.1) # 5 BTMs were enriched at baseline of ivermectin-treated or placebo groups from the original data excluding problematic samples excluding  BTMs TBA
head(res.FDR.orig.wp.1)

### 5.1.3 Original data excluding lowly expressed genes that have high ratio of duplicates

geneIDsfiltgenes.1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.baseline.filt.genes), columns = c("SYMBOL","GENEID"))
idxbaseorig.filtgenes<- ids2indices(btm$genesets,id=geneIDsfiltgenes.1$SYMBOL)

cam.filtgenes.1 <- camera(countsbasem.filt.genes,idxbaseorig.filtgenes,design.treat,contrast=contm.treat)

res.FDR.base.filt.genes.1 <- cam.filtgenes.1[cam.filtgenes.1$FDR < 0.05 & cam.filtgenes.1$NGenes > 10,]
dim(res.FDR.base.filt.genes.1)  # 0 BTMs were enriched at baseline of ivermectin-treated or placebo groups from the data excluding lowly expressed genes that have high ratio of duplicates

### 5.1.4 Original data excluding lowly expressed genes that have high ratio of duplicates and problematic samples

geneIDsbasefiltgenes.wp.1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.baseline.filt.genes.wp), columns = c("SYMBOL","GENEID"))
idxbaseorig.filtgenes.wp.1 <- ids2indices(btm$genesets,id=geneIDsbasefiltgenes.wp.1$SYMBOL)

cam.base.filtgenes.wp.1 <- camera(countsbasem.filt.genes.wp,idxbaseorig.filtgenes.wp.1,design.treat.wp,contrast=contm.treat.wp)

res.FDR.base.filt.genes.wp.1 <- cam.base.filtgenes.wp.1[cam.base.filtgenes.wp.1$FDR < 0.05 & cam.base.filtgenes.wp.1$NGenes > 10,]
dim(res.FDR.base.filt.genes.wp.1) # 0 BTMs were enriched at baseline of ivermectin-treated or placebo groups from the data excluding lowly expressed genes that have high ratio of duplicates and problematic samples

## 5.2 Differences at timepoints 4 and 7 between ivermectin-treated and placebo groups
### 5.2.1 Original data
# Differences at timepoint 4 between ivermectin-treated and placebo groups

voom.orig.4 <- voom(countsm.orig.4, design.treat, plot = F) 

geneIDsorig4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.orig.4), columns = c("SYMBOL","GENEID"))
idxbaseorig.4 <- ids2indices(btm$genesets,id=geneIDsorig4$SYMBOL)

cam.base.orig.4 <- camera(countsm.orig.4,idxbaseorig.4,design.treat, contrast=contm.treat)

res.FDR.base.orig.4 <- cam.base.orig.4[cam.base.orig.4$FDR < 0.05 & cam.base.orig.4$NGenes > 10,]
dim(res.FDR.base.orig.4) # 72 BTMs were enriched at timepoint 4 of ivermectin-treated or placebo groups from the original data 
res.FDR.base.orig.4 <-res.FDR.base.orig.4%>%dplyr::filter(!str_detect(rownames(res.FDR.base.orig.4), "TBA")) # 9 BTMs TBA
dim(res.FDR.base.orig.4) # 63 BTMs were enriched at timepoint 4 of ivermectin-treated or placebo groups from the original data excluding BTMs TBA
head(res.FDR.base.orig.4)

# Differences at timepoint 7 between ivermectin-treated and placebo groups

voom.orig.7 <- voom(countsm.orig.7, design.treat, plot = F) 

geneIDsorig7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.orig.7), columns = c("SYMBOL","GENEID"))
idxbaseorig.7 <- ids2indices(btm$genesets,id=geneIDsorig7$SYMBOL)

cam.base.orig.7 <- camera(countsm.orig.7,idxbaseorig.7,design.treat,contrast=contm.treat)

res.FDR.base.orig.7 <- cam.base.orig.7[cam.base.orig.7$FDR < 0.05 & cam.base.orig.7$NGenes > 10,]
dim(res.FDR.base.orig.7) # 84 BTMs were enriched at timepoint 7 of ivermectin-treated or placebo groups from the original data 
res.FDR.base.orig.7 <-res.FDR.base.orig.7%>%dplyr::filter(!str_detect(rownames(res.FDR.base.orig.7), "TBA")) # 21 BTMs TBA
dim(res.FDR.base.orig.7) # 63 BTMs were enriched at timepoint 7 of ivermectin-treated or placebo groups from the original data excluding BTMs TBA
head(res.FDR.base.orig.7)

### 5.2.2 Original data excluding problematic samples
# Differences at timepoint 4 between ivermectin-treated and placebo groups

treat.wp.4 <- as.factor(ifelse(col.orig.wp.4$treat == 1, "Ivermectin", "Placebo"))
design.treat.wp.4 <- model.matrix(~0 + treat.wp.4)
  
voom.orig.wp.4 <- voom(countsm.orig.wp.4, design.treat.wp.4, plot = F) 

row_names_orig.wp.4 <- rownames(voom.orig.wp.4$E)
new_row_names_orig.wp.4 <- sub("\\.\\d+", "", row_names_orig.wp.4)
rownames(voom.orig.wp.4) <- new_row_names_orig.wp.4

geneIDsorig.wp.4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.orig.wp.4), columns = c("SYMBOL","GENEID"))
idxbaseorig.wp.4 <- ids2indices(btm$genesets,id=geneIDsorig.wp.4$SYMBOL)

contm.treat.wp.4 <- makeContrasts(IvervsPlac=treat.wp.4Ivermectin-treat.wp.4Placebo,
                                 levels = design.treat.wp.4)

cam.base.orig.wp.4 <- camera(counts.orig.wp.4,idxbaseorig.wp.4,design.treat.wp.4,contrast=contm.treat.wp.4)

res.FDR.base.orig.wp.4 <- cam.base.orig.wp.4[cam.base.orig.wp.4$FDR < 0.05 & cam.base.orig.wp.4$NGenes > 10,]
dim(res.FDR.base.orig.wp.4) # 74 BTMs were enriched at timepoint 4 of ivermectin-treated or placebo groups from the original data excluding problematic samples

res.FDR.base.orig.wp.4<-res.FDR.base.orig.wp.4%>%dplyr::filter(!str_detect(rownames(res.FDR.base.orig.wp.4), "TBA")) # 7 BTMs TBA

dim(res.FDR.base.orig.wp.4) # 67 BTMs were enriched at timepoint 4 of ivermectin-treated or placebo groups from the original data excluding problematic samples and BTMs TBA
head(res.FDR.base.orig.wp.4)

### 5.2.3 Data excluding lowly expressed genes that have high ratio of duplicates
# Differences at timepoint 4 between ivermectin-treated and placebo groups

voom.filt.genes.4 <- voom(countsm.filt.genes.4, design.treat, plot = F) 

geneIDsfiltgenes4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.filt.genes.4), columns = c("SYMBOL","GENEID"))
idxbasefiltgenes.4 <- ids2indices(btm$genesets,id=geneIDsfiltgenes4$SYMBOL)

cam.base.filt.genes.4 <- camera(countsm.filt.genes.4,idxbasefiltgenes.4,design.treat,contrast=contm.treat)

res.FDR.base.filtgenes.4 <- cam.base.filt.genes.4[cam.base.filt.genes.4$FDR < 0.05 & cam.base.filt.genes.4$NGenes > 10,]
dim(res.FDR.base.filtgenes.4) # 0 BTMs were enriched at timepoint 4 of ivermectin-treated or placebo groups from the data excluding lowly expressed genes that have high duplication ratio

# Differences at timepoint 7 between ivermectin-treated and placebo groups

voom.filt.genes.7 <- voom(countsm.filt.genes.7, design.treat, plot = F) 

geneIDsfiltgenes7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.filt.genes.7), columns = c("SYMBOL","GENEID"))
idxbasefiltgenes.7 <- ids2indices(btm$genesets,id=geneIDsfiltgenes7$SYMBOL)

cam.base.filt.genes.7 <- camera(countsm.filt.genes.7,idxbasefiltgenes.7,design.treat,contrast=contm.treat)

res.FDR.base.filt.genes.7 <- cam.base.filt.genes.7[cam.base.filt.genes.7$FDR < 0.05 & cam.base.filt.genes.7$NGenes > 10,]
dim(res.FDR.base.filt.genes.7) # 0 BTMs were enriched at timepoint 7 of ivermectin-treated or placebo groups from the data excluding lowly expressed genes that have high duplication ratio

### 5.2.4 Data excluding lowly expressed genes that have high ratio of duplicates and problematic samples
# Differences at timepoint 4 between ivermectin-treated and placebo groups

voom.filt.genes.wp.4 <- voom(countsm.filt.genes.wp.4, design.treat.wp.4, plot = F) 

geneIDsfiltgenes.wp.4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.filt.genes.wp.4), columns = c("SYMBOL","GENEID"))
idxbasefiltgenes.wp.4 <- ids2indices(btm$genesets,id=geneIDsfiltgenes.wp.4$SYMBOL)


cam.base.filt.genes.wp.4 <- camera(countsm.filt.genes.wp.4,idxbasefiltgenes.wp.4,design.treat.wp.4,contrast=contm.treat.wp.4)

res.FDR.base.filtgenes.wp.4 <- cam.base.filt.genes.wp.4[cam.base.filt.genes.wp.4$FDR < 0.05 & cam.base.filt.genes.wp.4$NGenes > 10,]
res.FDR.base.filtgenes.wp.4 # 0 BTMs were enriched at timepoint 4 of ivermectin-treated or placebo groups from the data excluding lowly expressed genes that have high duplication ratio and problematic samples

# 6. Data visualization
## 6.1 Differences at baseline between ivermectin-treated and placebo groups 
### 6.1.1 Volcano plot of DEGs: DESeq from the original data

colors <- c("grey", "red", "blue")
diffDESeq.1 <- res.DESeq.orig.1 %>% mutate(gene = rownames(res.DESeq.orig.1), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Overexpressed in the ivermectin group", ifelse(padj<0.3 & log2FoldChange < (-1), "Underexpressed in the ivermectin group","Non significant")))

diffDESeq.1.Volcano <- ggplot(data=diffDESeq.1, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs in the baseline from the original data") +
     theme(plot.title = element_text(size = 10))

significant_genes.1 <- diffDESeq.1[diffDESeq.1$sig != "Non significant", ]
significant_genes.1 <- distinct(significant_genes.1, rownames(significant_genes.1), .keep_all = TRUE)

diffDESeq.1.Volcano <- diffDESeq.1.Volcano +
  geom_text_repel(data = significant_genes.1, aes(label = c("RAMP3", "PGM5", "CD177", "MAPK8IP1P2", "lncRNA (ENSG00000285668)", "lncRNA (ENSG00000290457)")),
                  box.padding = 0.5, point.padding = 0.2) +
  guides(col = guide_legend(override.aes = list(size = 3)))

print(diffDESeq.1.Volcano)

### 6.1.2 Volcano plot of DEGs: DESeq from the original data excluding problematic samples

diffDESeq.wp.1 <- res.DESeq.orig.wp.1 %>% mutate(gene = rownames(res.DESeq.orig.wp.1), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Upregulated in the ivermectin-treated group", ifelse(padj<0.3 & log2FoldChange < (-1), "Downregulated in the ivermectin-treated group","Non significant")))

diffDESeq.wp.1.Volcano <- ggplot(data=diffDESeq.wp.1, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs in the baseline from the original data excluding problematic samples") +
     theme(plot.title = element_text(size = 10))

significant_genes.wp.1 <- diffDESeq.wp.1[diffDESeq.wp.1$sig != "Non significant", ]
significant_genes.wp.1 <- distinct(significant_genes.wp.1, rownames(significant_genes.wp.1), .keep_all = TRUE)

diffDESeq.wp.1.Volcano <- diffDESeq.wp.1.Volcano +
  geom_text_repel(data = significant_genes.wp.1, aes(label = c("LGR6", "PGM5", "HSPA8P8", "lncRNA (ENSG00000285668)", "lncRNA (ENSG00000290457)"),
                  box.padding = 0.5, point.padding = 0.2))

print(diffDESeq.wp.1.Volcano)

### 6.1.3 Volcano plot of DEGs: DESeq from the data excluding lowly expressed genes that have high ratio of duplicates

diffDESeq.filtgenes.1 <- res.DESeq.filt.genes.1 %>% mutate(gene = rownames(res.DESeq.filt.genes.1), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Upregulated in the ivermectin-treated group", ifelse(padj<0.3 & log2FoldChange < (-1), "Downregulated in the ivermectin-treated group","Non significant")))

diffDESeq.filtgenes.1.Volcano <- ggplot(data=diffDESeq.filtgenes.1, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs in the baseline from the data excluding lowly expressed genes that have high duplication ratios") +
     theme(plot.title = element_text(size = 10))

significant_genes.filtgenes.1 <- diffDESeq.filtgenes.1[diffDESeq.filtgenes.1$sig != "Non significant", ]
significant_genes.filtgenes.1 <- distinct(significant_genes.filtgenes.1, rownames(significant_genes.filtgenes.1), .keep_all = TRUE)

diffDESeq.filtgenes.1.Volcano <- diffDESeq.filtgenes.1.Volcano +
  geom_text_repel(data = significant_genes.filtgenes.1, aes(label = c("RAMP3", "PGM5", "lncRNA (ENSG00000285668)", "MAPK8IP1P2","lncRNA (ENSG00000290457)", "CD177"),
                  box.padding = 0.5, point.padding = 0.2))

print(diffDESeq.filtgenes.1.Volcano)

### 6.1.4 Volcano plot of DEGs: DESeq from the data excluding lowly expressed genes that have high ratio of duplicates and the problematic samples

diffDESeq.filtgenes.wp.1 <- res.DESeq.filt.genes.wp.1 %>% mutate(gene = rownames(res.DESeq.filt.genes.wp.1), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Upregulated in the ivermectin-treated group", ifelse(padj<0.3 & log2FoldChange < (-1), "Downregulated in the ivermectin-treated group","Non significant")))

diffDESeq.filtgenes.wp.1.Volcano <- ggplot(data=diffDESeq.filtgenes.wp.1, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs in the baseline from the data excluding lowly expressed genes that have high duplication ratios and problematic samples") +
     theme(plot.title = element_text(size = 10))

significant_genes.filtgenes.wp.1 <- diffDESeq.filtgenes.wp.1[diffDESeq.filtgenes.wp.1$sig != "Non significant", ]
significant_genes.filtgenes.wp.1 <- distinct(significant_genes.filtgenes.wp.1, rownames(significant_genes.filtgenes.wp.1), .keep_all = TRUE)

diffDESeq.filtgenes.wp.1.Volcano <- diffDESeq.filtgenes.wp.1.Volcano +
  geom_text_repel(data = significant_genes.filtgenes.wp.1, aes(label = c("LGR6", "HSPA8P8", "PGM5", "lncRNA (ENSG00000285668)", "lncRNA (ENSG00000290457)"),
                  box.padding = 0.5, point.padding = 0.2))

print(diffDESeq.filtgenes.wp.1.Volcano)

### 6.1.5 Dot plot of Enrichment Analysis: limma (voom + camera + BTMs) from the original data 

annotationbtms<-read.csv("BTM_annotation_revision.csv", sep=";")
colnames(annotationbtms)[3]<-"Geneset"

res.FDR.orig.1$Geneset <- rownames(res.FDR.orig.1)

baseline.limma.btm.orig <- res.FDR.orig.1 %>%inner_join(annotationbtms, by = "Geneset")
baseline.limma.btm.orig$NGenes<-as.numeric(baseline.limma.btm.orig$NGenes)
baseline.limma.btm.orig$High.level.annotation.group[grep("platelet", baseline.limma.btm.orig$Geneset)]<-"Platelets"
baseline.limma.btm.orig<-baseline.limma.btm.orig[order(baseline.limma.btm.orig$High.level.annotation.group),]
baseline.limma.btm.orig$Geneset<-factor(baseline.limma.btm.orig$Geneset, levels=c(baseline.limma.btm.orig$Geneset))

EA.baseline.orig <- baseline.limma.btm.orig %>% 
  ggplot(aes(x=Direction, y = (Geneset), size = NGenes, color = log10(FDR)), size=4) +
  geom_point()+ scale_x_discrete(labels=c('Paths enriched in placebo', 'Paths enriched in ivermectin'))+ 
  theme(
    axis.text = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(10),  
  ) +
  scale_color_gradient(low = "red2", high = "steelblue2")+
  scale_size_area(max_size = 10) +
  labs(
    title = "Enrichment Analysis using BTMs at baseline from the original data"
  )

print(EA.baseline.orig)

### 6.1.6 Dot plot of Enrichment Analysis: limma (voom + camera + Blood Transcriptional Modules) from the original data excluding the problematic samples

res.FDR.orig.wp.1$Geneset <- rownames(res.FDR.orig.wp.1)

baseline.limma.btm.orig.wp <- res.FDR.orig.wp.1 %>%inner_join(annotationbtms, by = "Geneset")
baseline.limma.btm.orig.wp$NGenes<-as.numeric(baseline.limma.btm.orig.wp$NGenes)
baseline.limma.btm.orig.wp$High.level.annotation.group[grep("platelet", baseline.limma.btm.orig.wp$Geneset)]<-"Platelets"
baseline.limma.btm.orig.wp<-baseline.limma.btm.orig.wp[order(baseline.limma.btm.orig.wp$High.level.annotation.group),]
baseline.limma.btm.orig.wp$Geneset<-factor(baseline.limma.btm.orig.wp$Geneset, levels=c(baseline.limma.btm.orig.wp$Geneset))

EA.baseline.orig.wp <-baseline.limma.btm.orig.wp %>% 
  ggplot(aes(x=Direction, y = (Geneset), size = NGenes, color = log10(FDR)), size=6) +
  geom_point()+ scale_x_discrete(labels=c('Paths enriched in placebo', 'Paths enriched in ivermectin'))+ 
  theme(
    axis.text = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(10),  
  ) +
  scale_color_gradient(low = "red2", high = "steelblue2")+
  scale_size_area(max_size = 10) +
  labs(
    title = "Enrichment Analysis using BTMs at baseline from the original data excluding problematic samples" 
  )

print(EA.baseline.orig.wp)

## 6.2 Differences at timepoints 4 and 7 between the ivermectin-treated and the placebo groups
### 6.2.1 Volcano plot of DEGs: DESeq at timepoint 4 from the original data

diffDESeq.4 <- res.DESeq.orig.4 %>% mutate(gene = rownames(res.DESeq.orig.4), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Overexpressed in the ivermectin group", ifelse(padj<0.3 & log2FoldChange < (-1), "Underexpressed in the ivermectin group","Non significant")))

diffDESeq.4.Volcano <- ggplot(data=diffDESeq.4, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs at timepoint 4 from the original data") +
     theme(plot.title = element_text(size = 8))

significant_genes.4 <- diffDESeq.4[diffDESeq.4$sig != "Non significant", ]
significant_genes.4 <- distinct(significant_genes.4, rownames(significant_genes.4), .keep_all = TRUE)

diffDESeq.4.Volcano <- diffDESeq.4.Volcano +
  geom_text_repel(data = significant_genes.4, aes(label = c("BPI", "TMEM252-DT", "MAPK8IP1P2", "LINC02804","lncRNA (ENSG00000285668)")),
                  box.padding = 0.5, point.padding = 0.2)

print(diffDESeq.4.Volcano)

### 6.2.2 Volcano plot of DEGs: DESeq at timepoint 4 from the original data excluding the problematic samples

diffDESeq.wp.4 <- res.DESeq.orig.wp.4 %>% mutate(gene = rownames(res.DESeq.orig.wp.4), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Overexpressed in the ivermectin group", ifelse(padj<0.3 & log2FoldChange < (-1), "Underexpressed in the ivermectin group","Non significant")))

diffDESeq.wp.4.Volcano <- ggplot(data=diffDESeq.wp.4, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs at timepoint 4 from the original data excluding problematic samples") +
     theme(plot.title = element_text(size = 7))


significant_genes.wp.4 <- diffDESeq.wp.4[diffDESeq.wp.4$sig != "Non significant", ]
significant_genes.wp.4 <- distinct(significant_genes.wp.4, rownames(significant_genes.wp.4), .keep_all = TRUE)

diffDESeq.wp.4.Volcano <- diffDESeq.wp.4.Volcano +
  geom_text_repel(data = significant_genes.wp.4, aes(label = c("BPI", "TMEM252-DT", "lncRNA (ENSG00000285668)"),
                  box.padding = 0.5, point.padding = 0.2))

print(diffDESeq.wp.4.Volcano)

### 6.2.3 Volcano plot of DEGs: DESeq at timepoint 4 from the data excluding lowly expressed genes that have high ratio of duplicates

diffDESeq.filt.genes.4 <- res.DESeq.filt.genes.4 %>% mutate(gene = rownames(res.DESeq.filt.genes.4), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Overexpressed in the ivermectin group", ifelse(padj<0.3 & log2FoldChange < (-1), "Underexpressed in the ivermectin group","Non significant")))

diffDESeq.filt.genes.4.Volcano <- ggplot(data=diffDESeq.filt.genes.4, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs at timepoint 4 from the data excluding lowly expressed genes that have high ratio of duplicates") +
     theme(plot.title = element_text(size = 6))

significant_genes.filt.genes.4 <- diffDESeq.filt.genes.4[diffDESeq.filt.genes.4$sig != "Non significant", ]
significant_genes.filt.genes.4 <- distinct(significant_genes.filt.genes.4, rownames(significant_genes.filt.genes.4), .keep_all = TRUE)

diffDESeq.filt.genes.4.Volcano <- diffDESeq.filt.genes.4.Volcano +
  geom_text_repel(data = significant_genes.filt.genes.4, aes(label = c("LINC02804", "TMEM252-DT", "lncRNA (ENSG00000285668)", "MAPK8IP1P2", "BPI")),
                  box.padding = 0.5, point.padding = 0.2)

print(diffDESeq.filt.genes.4.Volcano)

### 6.2.4 Volcano plot of DEGs: DESeq at timepoint 4 from the data excluding lowly expressed genes that have high ratio of duplicates and the problematic samples

diffDESeq.filt.genes.wp.4 <- res.DESeq.filt.genes.wp.4 %>% mutate(gene = rownames(res.DESeq.filt.genes.wp.4), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Overexpressed in the ivermectin group", ifelse(padj<0.3 & log2FoldChange < (-1), "Underexpressed in the ivermectin group","Non significant")))

diffDESeq.filt.genes.wp.4.Volcano <- ggplot(data=diffDESeq.filt.genes.wp.4, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs in the timepoint 4 from the data excluding lowly expressed genes that have high ratio of duplicates and problematic samples") +
     theme(plot.title = element_text(size = 5))

significant_genes.filt.genes.wp.4 <- diffDESeq.filt.genes.wp.4[diffDESeq.filt.genes.wp.4$sig != "Non significant", ]
significant_genes.filt.genes.wp.4 <- distinct(significant_genes.filt.genes.wp.4, rownames(significant_genes.filt.genes.wp.4), .keep_all = TRUE)

diffDESeq.filt.genes.wp.4.Volcano <- diffDESeq.filt.genes.wp.4.Volcano +
  geom_text_repel(data = significant_genes.filt.genes.wp.4, aes(label = c("TMEM252-DT", "lncRNA (ENSG00000285668)", "BPI"),
                  box.padding = 0.5, point.padding = 0.2))

print(diffDESeq.filt.genes.wp.4.Volcano)

### 6.2.5 Volcano plot of DEGs: DESeq at timepoint 7 from the original data

diffDESeq.7 <- res.DESeq.orig.7 %>% mutate(gene = rownames(res.DESeq.orig.7), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Overexpressed in the ivermectin group", ifelse(padj<0.3 & log2FoldChange < (-1), "Underexpressed in the ivermectin group","Non significant")))

diffDESeq.7 <- na.omit(diffDESeq.7)

diffDESeq.7.Volcano <- ggplot(data=diffDESeq.7, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs in the timepoint 7 from the original data") +
     theme(plot.title = element_text(size = 8))

significant_genes.7 <- diffDESeq.7[diffDESeq.7$sig != "Non significant", ]
significant_genes.7 <- distinct(significant_genes.7, rownames(significant_genes.7), .keep_all = TRUE)

diffDESeq.7.Volcano <- diffDESeq.7.Volcano +
  geom_text_repel(data = significant_genes.7, aes(label = c ("VNN1", "NUP210L", "VWDE", "CLIC6", "RBFOX3", "MYEOV", "BRD7P2", "COL13A1", "CD177", "IGLV1-36", "lncRNA (ENSG00000215765)", "IGHV3-72", "C14orf132", "TRGV5P", "CLEC2L", "TMEM158", "lncRNA (ENSG00000254810)", "LINC01580", "MAPK8IP1P2", "GCATP1", "LOC102723407", "lncRNA (ENSG00000285668)", "Protein coding gene (ENSG00000285837)", "LOC105375423", "lncRNA (ENSG00000289180)", "lncRNA (ENSG00000289320)", "lncRNA (ENSG00000290457)", "lncRNA (ENSG00000291034)", "Uncategorized gene (ENSG00000291194)")),
                  box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf, size = 2)

print(diffDESeq.7.Volcano)

### 6.2.7 Volcano plot of DEGs: DESeq at timepoint 7 from the data excluding lowly expressed genes that have high ratio of duplicates

diffDESeq.filt.genes.7 <- res.DESeq.filt.genes.7 %>% mutate(gene = rownames(res.DESeq.filt.genes.7), logp = -(log10(pvalue)), logadjp = -(log10(padj)),
                          FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(padj<0.3 & log2FoldChange > 1, "Overexpressed in the ivermectin group", ifelse(padj<0.3 & log2FoldChange < (-1), "Underexpressed in the ivermectin group","Non significant")))

diffDESeq.filt.genes.7 <- na.omit(diffDESeq.filt.genes.7)

diffDESeq.filt.genes.7.Volcano <- ggplot(data=diffDESeq.filt.genes.7, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs at timepoint 7 from the data excluding lowly expressed genes that have high ratio of duplicates") +
     theme(plot.title = element_text(size = 6))

significant_genes.filt.genes.7 <- diffDESeq.filt.genes.7[diffDESeq.filt.genes.7$sig != "Non significant", ]
significant_genes.filt.genes.7 <- distinct(significant_genes.filt.genes.7, rownames(significant_genes.filt.genes.7), .keep_all = TRUE)

diffDESeq.filt.genes.7.Volcano <- diffDESeq.filt.genes.7.Volcano +
  geom_text_repel(data = significant_genes.filt.genes.7, aes(label = c ("NUP210L", "TMEM158", "SETP20", "lncRNA (ENSG00000289320)", "VNN1", "VWDE", "TRGV5P", "LOC105375423", "CLEC2L", "	COL13A1", "MYEOV", "lncRNA (ENSG00000254810)", "GCATP1", "Antisense to DCAF5", "C14orf132", "IGHV3-72", "	LINC01580", "lncRNA (ENSG00000285668)", "MAPK8IP1P2", "lncRNA (ENSG00000290457)", "RBFOX3", "c", "CD177", "CLIC6", "IGLV1-36", "lncRNA (ENSG00000291034)","LOC102723407")),
                  box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf, size = 2)

print(diffDESeq.filt.genes.7.Volcano)

# 6.2.8 Dot plot of Enrichment Analysis: limma (voom + camera + BTMs) at timepoint 4 from the original data 

res.FDR.base.orig.4$Geneset <- rownames(res.FDR.base.orig.4)

timepoint.4.limma.btm.orig <- res.FDR.base.orig.4 %>%inner_join(annotationbtms, by = "Geneset")
timepoint.4.limma.btm.orig$NGenes<-as.numeric(timepoint.4.limma.btm.orig$NGenes)
timepoint.4.limma.btm.orig$High.level.annotation.group[grep("platelet", timepoint.4.limma.btm.orig$Geneset)]<-"Platelets"
timepoint.4.limma.btm.orig <- timepoint.4.limma.btm.orig[order(timepoint.4.limma.btm.orig$High.level.annotation.group),]
timepoint.4.limma.btm.orig$Geneset<-factor(timepoint.4.limma.btm.orig$Geneset, levels=c(timepoint.4.limma.btm.orig$Geneset))

EA.timepoint.4.orig <- timepoint.4.limma.btm.orig %>% 
  ggplot(aes(x=Direction, y = (Geneset), size = NGenes, color = log10(FDR)), size=4) +
  geom_point()+ scale_x_discrete(labels=c('Paths enriched in placebo', 'Paths enriched in ivermectin'))+ 
  theme(
    axis.text = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(10),  
  ) +
  scale_color_gradient(low = "red2", high = "steelblue2")+
  scale_size_area(max_size = 10) +
  labs(
    title = "Enrichment Analysis using BTMs at timepoint 4 from the original data"
  )

print(EA.timepoint.4.orig)

# 6.2.9 Dot plot of Enrichment Analysis: limma (voom + camera + BTMs) at timepoint 4 from the original data excluding the problematic samples

res.FDR.base.orig.wp.4$Geneset <- rownames(res.FDR.base.orig.wp.4)

timepoint.4.limma.btm.orig.wp <- res.FDR.base.orig.wp.4 %>%inner_join(annotationbtms, by = "Geneset")
timepoint.4.limma.btm.orig.wp$NGenes<-as.numeric(timepoint.4.limma.btm.orig.wp$NGenes)
timepoint.4.limma.btm.orig.wp$High.level.annotation.group[grep("platelet", timepoint.4.limma.btm.orig.wp$Geneset)]<-"Platelets"
timepoint.4.limma.btm.orig.wp <- timepoint.4.limma.btm.orig.wp[order(timepoint.4.limma.btm.orig.wp$High.level.annotation.group),]
timepoint.4.limma.btm.orig.wp$Geneset<-factor(timepoint.4.limma.btm.orig.wp$Geneset, levels=c(timepoint.4.limma.btm.orig.wp$Geneset))

EA.timepoint.4.orig.wp <- timepoint.4.limma.btm.orig.wp %>% 
  ggplot(aes(x=Direction, y = (Geneset), size = NGenes, color = log10(FDR)), size=4) +
  geom_point()+ scale_x_discrete(labels=c('Enriched paths in placebo', 'Enriched paths in ivermectin'))+ 
  theme(
    axis.text = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(10),  
  ) +
  scale_color_gradient(low = "red2", high = "steelblue2")+
  scale_size_area(max_size = 10) +
  labs(
    title = "Enrichment Analysis using BTMs at timepoint 4 from the original data excluding problematic samples)",  
  )

print(EA.timepoint.4.orig.wp)

# 6.2.10 Dot plot of Enrichment Analysis: limma (voom + camera + BTMs) at timepoint 7 from the original data 

res.FDR.base.orig.7$Geneset <- rownames(res.FDR.base.orig.7)

timepoint.7.limma.btm.orig <- res.FDR.base.orig.7 %>%inner_join(annotationbtms, by = "Geneset")
timepoint.7.limma.btm.orig$NGenes <- as.numeric(timepoint.7.limma.btm.orig$NGenes)
timepoint.7.limma.btm.orig$High.level.annotation.group[grep("platelet", timepoint.7.limma.btm.orig$Geneset)]<-"Platelets"
timepoint.7.limma.btm.orig <- timepoint.7.limma.btm.orig[order(timepoint.7.limma.btm.orig$High.level.annotation.group),]
timepoint.7.limma.btm.orig$Geneset<-factor(timepoint.7.limma.btm.orig$Geneset, levels=c(timepoint.7.limma.btm.orig$Geneset))

EA.timepoint.7.orig <- timepoint.7.limma.btm.orig %>% 
  ggplot(aes(x=Direction, y = (Geneset), size = NGenes, color = log10(FDR)), size=6) +
  geom_point()+ scale_x_discrete(labels=c('Enriched paths in placebo', 'Enriched paths in ivermectin'))+ 
  theme(
    axis.text = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(10),  
  ) +
  scale_color_gradient(low = "red2", high = "steelblue2")+
  scale_size_area(max_size = 10) +
  labs(
    title = "Enrichment Analysis using BTMs at timepoint 7 from the original data"
  )

print(EA.timepoint.7.orig)
