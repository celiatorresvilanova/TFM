###########################################################################################################################################################
                                          # Differential Gene Expression (DGE) analysis and functional analysis                                                                                                                             
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
library(clusterProfiler)

# 1. Differentially Expression Analysis using DESeq
## 1.1 Differences on day 4 between the ivermectin-treated and placebo groups (adjusting for baseline)

setwd("C:/Users/Cèlia Torres Vilanov/Documents/SAINT")
load(file="ExpressionSet_orig.RData")
load(file = "assay_data_orig_f.RData")

d4 <- ExpressionSet_orig[,ExpressionSet_orig$day==1|ExpressionSet_orig$day==4]
colDatad4 <- pData(d4)
counts4 <- assay_data_orig_f[,rownames(colDatad4)]
countData4 <- as.matrix(counts4)

designFormula4<- as.formula( ~ treat+day+treat:day) # treat = treatment

colDatad4$treat <- as.factor(colDatad4$treat) # treat = treatment
colDatad4$day <- as.factor(colDatad4$day)

DESeq.timepoint.4 <- DESeqDataSetFromMatrix(countData = countData4,
                                colData = colDatad4,
                                design = designFormula4)

DESeq.timepoint.4 <- DESeq(DESeq.timepoint.4, test="LRT", reduced = ~ treat+day) # treat = treatment

res.DESeq.timepoint.4 <- results(DESeq.timepoint.4, contrast = c("treat", '1', '2')) # treat = treatment

res.DESeq.timepoint.4 <- as.data.frame(res.DESeq.timepoint.4)

res.sign.DESeq.timepoint.4 <- res.DESeq.timepoint.4[res.DESeq.timepoint.4$pvalue < 0.05  & abs(res.DESeq.timepoint.4$log2FoldChange) > log2(2) & !is.na(res.DESeq.timepoint.4$pvalue),] 

dim(res.sign.DESeq.timepoint.4) # 44 DEGs
head(res.sign.DESeq.timepoint.4)

## 1.2 Differences on day 7 between the ivermectin-treated and placebo groups (adjusting for baseline)

d7 <- ExpressionSet_orig[,ExpressionSet_orig$day==1|ExpressionSet_orig$day==7]
colDatad7 <- pData(d7)
counts7 <- assay_data_orig_f[,rownames(colDatad7)]
countData7 <- as.matrix(counts7)

designFormula7 <- as.formula( ~ treat+day+treat:day) # treat = treatment

colDatad7$treat <- as.factor(colDatad7$treat) # treat = treatment
colDatad7$day <- as.factor(colDatad7$day)

DESeq.timepoint.7 <- DESeqDataSetFromMatrix(countData = countData7,
                                colData = colDatad7,
                                design = designFormula7)

DESeq.timepoint.7 <- DESeq(DESeq.timepoint.7, test="LRT", reduced = ~ treat+day) # treat = treatment

res.DESeq.timepoint.7 <- results(DESeq.timepoint.7, contrast = c("treat", '1', '2')) # treat = treatment
res.DESeq.timepoint.7 <- as.data.frame(res.DESeq.timepoint.7)
res.sign.DESeq.timepoint.7 <- res.DESeq.timepoint.7[res.DESeq.timepoint.7$pvalue < 0.05  & abs(res.DESeq.timepoint.7$log2FoldChange) > log2(2) & !is.na(res.DESeq.timepoint.7$pvalue),] 

dim(res.sign.DESeq.timepoint.7)  # 155 DEGs
head(res.sign.DESeq.timepoint.7)

## 1.3 Disease progression in the control group

control <- ExpressionSet_orig[,ExpressionSet_orig$treat==2]
colDataCtrl <- pData(control)
countsCtrl <- assay_data_orig_f[,rownames(colDataCtrl)]
countsmCtrl <- as.matrix(countsCtrl)

colDataCtrl$day <- as.factor(colDataCtrl$day)

DESeq.progresion.ctrl <- DESeqDataSetFromMatrix(countData = countsmCtrl,
                                colData = colDataCtrl,
                                design = ~ day)

DESeq.progresion.ctrl <- DESeq(DESeq.progresion.ctrl)

res.DESeq.progresion.ctrl.4 <- results(DESeq.progresion.ctrl, contrast = c("day", '4', '1'))
res.DESeq.progresion.ctrl.7 <- results(DESeq.progresion.ctrl, contrast = c("day", '7', '4'))
res.DESeq.progresion.ctrl.1.7 <- results(DESeq.progresion.ctrl, contrast = c("day", '7', '1'))

res.DESeq.progresion.ctrl.4 <- as.data.frame(res.DESeq.progresion.ctrl.4)
res.DESeq.progresion.ctrl.7 <- as.data.frame(res.DESeq.progresion.ctrl.7)
res.DESeq.progresion.ctrl.1.7 <- as.data.frame(res.DESeq.progresion.ctrl.1.7)

res.DESeq.sign.progresion.ctrl.4 <- res.DESeq.progresion.ctrl.4[res.DESeq.progresion.ctrl.4$padj < 0.05 & abs(res.DESeq.progresion.ctrl.4$log2FoldChange) > log2(2) & !is.na(res.DESeq.progresion.ctrl.4$padj), ]
dim(res.DESeq.sign.progresion.ctrl.4) # 350 DEGs
head(res.DESeq.sign.progresion.ctrl.4)

res.DESeq.sign.progresion.ctrl.7 <- res.DESeq.progresion.ctrl.7[res.DESeq.progresion.ctrl.7$padj < 0.05 & abs(res.DESeq.progresion.ctrl.7$log2FoldChange) > log2(2) & !is.na(res.DESeq.progresion.ctrl.7$padj), ]
dim(res.DESeq.sign.progresion.ctrl.7) # 115 DEGs
head(res.DESeq.sign.progresion.ctrl.7)

res.DESeq.sign.progresion.ctrl.1.7 <- res.DESeq.progresion.ctrl.1.7[res.DESeq.progresion.ctrl.1.7$padj < 0.05 & abs(res.DESeq.progresion.ctrl.1.7$log2FoldChange) > log2(2) & !is.na(res.DESeq.progresion.ctrl.1.7$padj), ]
dim(res.DESeq.sign.progresion.ctrl.1.7) # 1140 DEGs
head(res.DESeq.sign.progresion.ctrl.1.7)

immune_gene_list <- read.table("~/SAINT/immune_gene_list.txt", quote="\"", comment.char="")
immune_genes <- as.list(immune_gene_list)

selected_rows_4 <- res.DESeq.sign.progresion.ctrl.4[rownames(res.DESeq.sign.progresion.ctrl.4) %in% immune_gene_list$V1, ]  # 138 immune-related DEGs
selected_rows_4_up <- selected_rows_4[selected_rows_4$log2FoldChange>1.5,]  # 36 immune-related DEGs
selected_rows_4_down <- selected_rows_4[selected_rows_4$log2FoldChange< (-1.5),]  # 16 immune-related DEGs

selected_rows_7 <- res.DESeq.sign.progresion.ctrl.7[rownames(res.DESeq.sign.progresion.ctrl.7) %in% immune_gene_list$V1, ]  # 37 immune-related DEGs
selected_rows_7_up <- selected_rows_7[selected_rows_7$log2FoldChange>1.5,]  # 0 immune-related DEGs
selected_rows_7_down <- selected_rows_7[selected_rows_7$log2FoldChange< (-1.5),]  # 23 immune-related DEGs

selected_rows_1_7 <- res.DESeq.sign.progresion.ctrl.1.7[rownames(res.DESeq.sign.progresion.ctrl.1.7) %in% immune_gene_list$V1, ]  # 338 immune-related DEGs
selected_rows_1.7_up <- selected_rows_1_7[selected_rows_1_7$log2FoldChange>1.5,]  # 48 immune-related DEGs
selected_rows_1.7_down <- selected_rows_1_7[selected_rows_1_7$log2FoldChange< (-1.5),]  # 114 immune-related DEGs

## 1.4 Disease progression in the ivermectin group

ivermectin <- ExpressionSet_orig[,ExpressionSet_orig$treat==1]
colDataTreat <- pData(ivermectin)
countsTreat <- assay_data_orig_f[,rownames(colDataTreat)]
countsmTreat <- as.matrix(countsTreat)

colDataTreat$day <- as.factor(colDataTreat$day)

DESeq.progresion.treat <- DESeqDataSetFromMatrix(countData = countsmTreat,
                                colData = colDataTreat,
                                design = ~ day)

DESeq.progresion.treat <- DESeq(DESeq.progresion.treat)


res.DESeq.progresion.treat.4 <- results(DESeq.progresion.treat, contrast = c("day", '4', '1'))
res.DESeq.progresion.treat.7 <- results(DESeq.progresion.treat, contrast = c("day", '7', '4'))
res.DESeq.progresion.treat.1.7 <- results(DESeq.progresion.treat, contrast = c("day", '7', '1'))

res.DESeq.progresion.treat.4 <- as.data.frame(res.DESeq.progresion.treat.4)
res.DESeq.progresion.treat.7 <- as.data.frame(res.DESeq.progresion.treat.7)
res.DESeq.progresion.treat.1.7 <- as.data.frame(res.DESeq.progresion.treat.1.7)

res.DESeq.sign.progresion.treat.4 <- res.DESeq.progresion.treat.4[res.DESeq.progresion.treat.4$padj < 0.05 & abs(res.DESeq.progresion.treat.4$log2FoldChange) > log2(2) & !is.na(res.DESeq.progresion.treat.4$padj), ]
dim(res.DESeq.sign.progresion.treat.4) # 1092 DEGs
head(res.DESeq.sign.progresion.treat.4)

res.DESeq.sign.progresion.treat.7 <- res.DESeq.progresion.treat.7[res.DESeq.progresion.treat.7$padj < 0.05 & abs(res.DESeq.progresion.treat.7$log2FoldChange) > log2(2) & !is.na(res.DESeq.progresion.treat.7$padj), ]
dim(res.DESeq.sign.progresion.treat.7) # 0 DEGs
head(res.DESeq.sign.progresion.treat.7)

res.DESeq.sign.progresion.treat.1.7 <- res.DESeq.progresion.treat.1.7[res.DESeq.progresion.treat.1.7$padj < 0.05 & abs(res.DESeq.progresion.treat.1.7$log2FoldChange) > log2(2) & !is.na(res.DESeq.progresion.treat.1.7$padj), ]
dim(res.DESeq.sign.progresion.treat.1.7) # 1123 DEGs
head(res.DESeq.sign.progresion.treat.1.7)

selected_rows_t4 <- res.DESeq.sign.progresion.treat.4[rownames(res.DESeq.sign.progresion.treat.4) %in% immune_gene_list$V1, ]  # 325 immune-related DEGs
selected_rows_t4_up <- selected_rows_t4[selected_rows_t4$log2FoldChange>1.5,]  # 113 immune-related DEGs
selected_rows_t4_down <- selected_rows_t4[selected_rows_t4$log2FoldChange< (-1.5),]  # 64 immune-related DEGs

selected_rowst_1_7 <- res.DESeq.sign.progresion.treat.1.7[rownames(res.DESeq.sign.progresion.treat.1.7) %in% immune_gene_list$V1, ]  # 330 immune-related DEGs
selected_rowst_1.7_up <- selected_rowst_1_7[selected_rowst_1_7$log2FoldChange>1.5,]  # 60 immune-related DEGs
selected_rowst_1.7_down <- selected_rowst_1_7[selected_rowst_1_7$log2FoldChange< (-1.5),]  # 99 immune-related DEGs

# 2 Enrichment Analysis: Gene Set Enrichment Analysis (Gene Ontology terms)
## 2.1 Differences on day 4 between the ivermectin-treated and placebo groups (adjusting for baseline)

Overexpressed.timepoint.4 <- res.sign.DESeq.timepoint.4[res.sign.DESeq.timepoint.4$log2FoldChange > 1,]
Underexpressed.timepoint.4 <- res.sign.DESeq.timepoint.4[res.sign.DESeq.timepoint.4$log2FoldChange < 1,]

Overexpressed.timepoint.4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Overexpressed.timepoint.4), columns = c("SYMBOL","GENEID"))
Underexpressed.timepoint.4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Underexpressed.timepoint.4), columns = c("SYMBOL","GENEID"))

res.over.4 <- AnnotationDbi::select(org.Hs.eg.db, 
                          keys =  Overexpressed.timepoint.4$SYMBOL,
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")

res.under.4 <- AnnotationDbi::select(org.Hs.eg.db, 
                             keys =  Underexpressed.timepoint.4$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


res.enrichGO.over.4 <- enrichGO(gene = res.over.4$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE, 
                   pvalueCutoff = 0.05)

res.enrichGO.over.4 <- summary(res.enrichGO.over.4)
dim(res.enrichGO.over.4)  # 0 GO pathways

res.enrichGO.under.4 <- enrichGO(gene = res.under.4$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)

res.enrichGO.under.4 <- summary(res.enrichGO.under.4)
dim(res.enrichGO.under.4)  # 2 GO pathways

## 2.2 Differences on day 7 between the ivermectin-treated and placebo groups (adjusting for baseline)

Overexpressed.timepoint.7 <- res.sign.DESeq.timepoint.7[res.sign.DESeq.timepoint.7$log2FoldChange > 1,]
Underexpressed.timepoint.7 <- res.sign.DESeq.timepoint.7[res.sign.DESeq.timepoint.7$log2FoldChange < 1,]

Overexpressed.timepoint.7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Overexpressed.timepoint.7), columns = c("SYMBOL","GENEID"))
Underexpressed.timepoint.7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Underexpressed.timepoint.7), columns = c("SYMBOL","GENEID"))

res.over.7 <-AnnotationDbi::select(org.Hs.eg.db, 
                          keys =  Overexpressed.timepoint.7$SYMBOL,
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")

res.under.7 <-AnnotationDbi::select(org.Hs.eg.db, 
                             keys =  Underexpressed.timepoint.7$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")

 
res.enrichGO.over.7 <- enrichGO(gene = res.over.7$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE, 
                   pvalueCutoff = 0.05)

res.enrichGO.over.7 <- summary(res.enrichGO.over.7)
dim(res.enrichGO.over.7)  # 14 GO pathways


res.enrichGO.under.7 <- enrichGO(gene = res.under.7$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)

res.enrichGO.under.7 <- summary(res.enrichGO.under.7)
dim(res.enrichGO.under.7)  # 16 GO pathways


## 2.3 Disease progression in the control group

# Progression day 1 to 4
Over.progression.ctrl.4 <- res.DESeq.sign.progresion.ctrl.4[res.DESeq.sign.progresion.ctrl.4$log2FoldChange > 1,]
Under.progression.ctrl.4 <- res.DESeq.sign.progresion.ctrl.4[res.DESeq.sign.progresion.ctrl.4$log2FoldChange < 1,]

Over.progression.ctrl.4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Over.progression.ctrl.4), columns = c("SYMBOL","GENEID"))
Under.progression.ctrl.4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Under.progression.ctrl.4), columns = c("SYMBOL","GENEID"))

res.over.progr.ctrl.4 <-AnnotationDbi::select(org.Hs.eg.db, 
                          keys =  Over.progression.ctrl.4$SYMBOL,
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")

res.under.progr.ctrl.4 <-AnnotationDbi::select(org.Hs.eg.db, 
                             keys =  Under.progression.ctrl.4$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


res.GO.over.progr.ctrl.4 <- enrichGO(gene = res.over.progr.ctrl.4$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE, 
                   pvalueCutoff = 0.05)

res.GO.over.progr.ctrl.4 <- summary(res.GO.over.progr.ctrl.4)
res.GO.over.progr.ctrl.4 <- res.GO.over.progr.ctrl.4[res.GO.over.progr.ctrl.4$Count>20 & res.GO.over.progr.ctrl.4$p.adjust<0.05,]
dim(res.GO.over.progr.ctrl.4)  # 32 GO pathways

res.GO.under.progr.ctrl.4 <- enrichGO(gene = res.under.progr.ctrl.4$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)

res.GO.under.progr.ctrl.4 <- summary(res.GO.under.progr.ctrl.4)
res.GO.under.progr.ctrl.4 <- res.GO.under.progr.ctrl.4[res.GO.under.progr.ctrl.4$Count>20 & res.GO.under.progr.ctrl.4$p.adjust<0.05,]
dim(res.GO.under.progr.ctrl.4)  # 0 GO pathways

# Progression day 4 to 7
Over.progression.ctrl.7 <- res.DESeq.sign.progresion.ctrl.7[res.DESeq.sign.progresion.ctrl.7$log2FoldChange > 1,]
Under.progression.ctrl.7 <- res.DESeq.sign.progresion.ctrl.7[res.DESeq.sign.progresion.ctrl.7$log2FoldChange < 1,]

Under.progression.ctrl.7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Under.progression.ctrl.7), columns = c("SYMBOL","GENEID"))

res.under.progr.ctrl.7 <-AnnotationDbi::select(org.Hs.eg.db, 
                             keys =  Under.progression.ctrl.7$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


res.GO.under.progr.ctrl.7 <- enrichGO(gene = res.under.progr.ctrl.7$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)

res.GO.under.progr.ctrl.7 <- summary(res.GO.under.progr.ctrl.7)
res.GO.under.progr.ctrl.7 <- res.GO.under.progr.ctrl.7[res.GO.under.progr.ctrl.7$Count>20 & res.GO.under.progr.ctrl.7$p.adjust<0.05,]
dim(res.GO.under.progr.ctrl.7) # 5 GO pathways

# Progression day 1 to 7
Over.progression.ctrl.1.7 <- res.DESeq.sign.progresion.ctrl.1.7[res.DESeq.sign.progresion.ctrl.1.7$log2FoldChange > 1,]
Under.progression.ctrl.1.7 <- res.DESeq.sign.progresion.ctrl.1.7[res.DESeq.sign.progresion.ctrl.1.7$log2FoldChange < 1,]

Over.progression.ctrl.1.7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Over.progression.ctrl.1.7), columns = c("SYMBOL","GENEID"))
Under.progression.ctrl.1.7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Under.progression.ctrl.1.7), columns = c("SYMBOL","GENEID"))

res.over.progr.ctrl.1.7 <-AnnotationDbi::select(org.Hs.eg.db, 
                          keys =  Over.progression.ctrl.1.7$SYMBOL,
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")

res.under.progr.ctrl.1.7 <-AnnotationDbi::select(org.Hs.eg.db, 
                             keys =  Under.progression.ctrl.1.7$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


res.GO.over.progr.ctrl.1.7 <- enrichGO(gene = res.over.progr.ctrl.1.7$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE, 
                   pvalueCutoff = 0.05)

res.GO.over.progr.ctrl.1.7 <- summary(res.GO.over.progr.ctrl.1.7)
res.GO.over.progr.ctrl.1.7 <- res.GO.over.progr.ctrl.1.7[res.GO.over.progr.ctrl.1.7$Count>20 & res.GO.over.progr.ctrl.1.7$p.adjust < 0.05,]
dim(res.GO.over.progr.ctrl.1.7)  # 2 GO pathways


res.under.progr.ctrl.1.7 <-AnnotationDbi::select(org.Hs.eg.db, 
                             keys =  Under.progression.ctrl.1.7$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


res.GO.under.progr.ctrl.1.7 <- enrichGO(gene = res.under.progr.ctrl.1.7$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)

res.GO.under.progr.ctrl.1.7 <- summary(res.GO.under.progr.ctrl.1.7)
res.GO.under.progr.ctrl.1.7 <- res.GO.under.progr.ctrl.1.7[res.GO.under.progr.ctrl.1.7$Count>20 & res.GO.under.progr.ctrl.1.7$p.adjust < 0.05,]
dim(res.GO.under.progr.ctrl.1.7)  # 132 GO pathways


## 2.4 Disease progression in the treatment group
# Progression day 1 to 4
Over.progression.treat.4 <- res.DESeq.sign.progresion.treat.4[res.DESeq.sign.progresion.treat.4$log2FoldChange > 1,]
Under.progression.treat.4 <- res.DESeq.sign.progresion.treat.4[res.DESeq.sign.progresion.treat.4$log2FoldChange < 1,]

Over.progression.treat.4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Over.progression.treat.4), columns = c("SYMBOL","GENEID"))
Under.progression.treat.4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Under.progression.treat.4), columns = c("SYMBOL","GENEID"))

res.over.progr.treat.4 <-AnnotationDbi::select(org.Hs.eg.db, 
                          keys =  Over.progression.treat.4$SYMBOL,
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")

res.under.progr.treat.4 <-AnnotationDbi::select(org.Hs.eg.db, 
                             keys =  Under.progression.treat.4$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


res.GO.over.progr.treat.4 <- enrichGO(gene = res.over.progr.treat.4$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE, 
                   pvalueCutoff = 0.05)

res.GO.over.progr.treat.4 <- summary(res.GO.over.progr.treat.4)
res.GO.over.progr.treat.4 <- res.GO.over.progr.treat.4[res.GO.over.progr.treat.4$Count>20 & res.GO.over.progr.treat.4$p.adjust < 0.05,]
dim(res.GO.over.progr.treat.4) # 26 GO pathways

res.GO.under.progr.treat.4 <- enrichGO(gene = res.under.progr.treat.4$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)

res.GO.under.progr.treat.4 <- summary(res.GO.under.progr.treat.4)
res.GO.under.progr.treat.4 <- res.GO.under.progr.treat.4[res.GO.under.progr.treat.4$Count>20 & res.GO.under.progr.treat.4$p.adjust < 0.05,]
dim(res.GO.under.progr.treat.4) # 38 GO pathways

# Progression day 1 to 7
Over.progression.treat.1.7 <- res.DESeq.sign.progresion.treat.1.7[res.DESeq.sign.progresion.treat.1.7$log2FoldChange > 1,]
Under.progression.treat.1.7 <- res.DESeq.sign.progresion.treat.1.7[res.DESeq.sign.progresion.treat.1.7$log2FoldChange < 1,]

Over.progression.treat.1.7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Over.progression.treat.1.7), columns = c("SYMBOL","GENEID"))
Under.progression.treat.1.7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(Under.progression.treat.1.7), columns = c("SYMBOL","GENEID"))

res.over.progr.treat.1.7 <-AnnotationDbi::select(org.Hs.eg.db, 
                          keys =  Over.progression.treat.1.7$SYMBOL,
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")

res.under.progr.treat.1.7 <-AnnotationDbi::select(org.Hs.eg.db, 
                             keys =  Under.progression.treat.1.7$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


res.GO.over.progr.treat.1.7 <- enrichGO(gene = res.over.progr.treat.1.7$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE, 
                   pvalueCutoff = 0.05)

res.GO.over.progr.treat.1.7 <- summary(res.GO.over.progr.treat.1.7)
res.GO.over.progr.treat.1.7 <- res.GO.over.progr.treat.1.7[res.GO.over.progr.treat.1.7$Count>20 & res.GO.over.progr.treat.1.7$p.adjust < 0.05,]
dim(res.GO.over.progr.treat.1.7)  # 2 GO pathways


res.under.progr.treat.1.7 <-AnnotationDbi::select(org.Hs.eg.db, 
                             keys =  Under.progression.treat.1.7$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


res.GO.under.progr.treat.1.7 <- enrichGO(gene = res.under.progr.treat.1.7$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)

res.GO.under.progr.treat.1.7 <- summary(res.GO.under.progr.treat.1.7)
res.GO.under.progr.treat.1.7 <- res.GO.under.progr.treat.1.7[res.GO.under.progr.treat.1.7$Count>20 & res.GO.under.progr.treat.1.7$p.adjust < 0.05,]
dim(res.GO.under.progr.treat.1.7)  # 88 GO pathways

# 3 Enrichment Analysis: Limma (voom + camera) using Blood transcriptional modules (BTM)
## 3.1 Differences on day 4 between the ivermectin-treated and placebo groups (adjusting for baseline)

colDatad4$treat <- ifelse(colDatad4$treat==1, "Iver", "Ctrl") # treat = treatment
colDatad4$treat.timepoint <- paste0(colDatad4$treat, sep=".",colDatad4$day)

design.4 <- model.matrix(~0+treat.timepoint, data=colDatad4)
colnames(design.4) <-c("Ctrl.1", "Ctrl.4", "Iver.1", "Iver.4")

voom.4 <- voom(countData4, design=design.4, plot=TRUE)

geneIDs.4 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.4), columns = c("SYMBOL","GENEID"))

load(file= "btm.RData")
names(btm$genesets) <- btm$geneset.names

idx.4 <- ids2indices(btm$genesets,id=geneIDs.4$SYMBOL)

contr.matrix.4 <- makeContrasts(
  Iver4vsCtrl4 = (Iver.4-Iver.1)-(Ctrl.4-Ctrl.1), 
  levels = colnames(design.4))

cam.iverctrl4 <- camera(countData4,idx.4,design.4,contrast=contr.matrix.4)

sigeneset.4 <- cam.iverctrl4[cam.iverctrl4$FDR<0.05,]
sigeneset.4 <- sigeneset.4[sigeneset.4$NGenes>10,]
dim(sigeneset.4) # 51 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 4

sigs.4 <- sigeneset.4%>%dplyr::filter(!str_detect(rownames(sigeneset.4), "TBA"))  # 12 BTMs to be annotated (TBA)
dim(sigs.4) # 39 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 4 excluding To Be Annotated BTMs
sigs.4$Geneset <- rownames(sigs.4)

## 3.2 Differences on day 7 between the ivermectin-treated and placebo groups (adjusting for baseline)

colDatad7$treat <- ifelse(colDatad7$treat==1, "Iver", "Ctrl") # treat = treatment
colDatad7$treat.timepoint <- paste0(colDatad7$treat, sep=".",colDatad7$day)

design.7 <- model.matrix(~0+treat.timepoint, data=colDatad7)
colnames(design.7) <-c("Ctrl.1", "Ctrl.7", "Iver.1", "Iver.7")

voom.7 <- voom(countData7, design=design.7, plot=TRUE)

geneIDs.7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.7), columns = c("SYMBOL","GENEID"))
idx.7 <- ids2indices(btm$genesets,id=geneIDs.7$SYMBOL)

contr.matrix.7 <- makeContrasts(
  Iver7vsCtrl7 = (Iver.7-Iver.1)-(Ctrl.7-Ctrl.1), 
  levels = colnames(design.7))

cam.iverctrl7 <- camera(countData7,idx.7,design.7,contrast=contr.matrix.7)

sigeneset.7 <- cam.iverctrl7[cam.iverctrl7$FDR<0.05,]
sigeneset.7 <- sigeneset.7[sigeneset.7$NGenes>10,]
dim(sigeneset.7)  # 49 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 4

sigs.7 <- sigeneset.7%>%dplyr::filter(!str_detect(rownames(sigeneset.7), "TBA")) # 21 To Be Annotated (TBA) BTMs
dim(sigs.7) # 28 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 4 excluding to be annotated BTMs
sigs.7$Geneset <- rownames(sigs.7)

## 3.3 Disease progression in the control group

design.ctrl.timepoint <- model.matrix(~0+day, data=colDataCtrl)
voom.ctrl <- voom(countsCtrl, design=design.ctrl.timepoint, plot=TRUE)
geneIDs.ctrl <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.ctrl), columns = c("SYMBOL","GENEID"))
idx.ctrl <- ids2indices(btm$genesets,id=geneIDs.ctrl$SYMBOL)

contr.ctrl.4 <- makeContrasts(
  Ctrl4vsCtrl1 = day4 - day1, 
  levels = colnames(design.ctrl.timepoint))

contr.ctrl.7 <- makeContrasts(
  Ctrl4vsCtrl1 = day7 - day4, 
  levels = colnames(design.ctrl.timepoint))

contr.ctrl.1.7 <- makeContrasts(
  Ctrl7vsCtrl1 = day7 - day1, 
  levels = colnames(design.ctrl.timepoint))

cam.ctrl.4 <- camera(countsmCtrl,idx.ctrl,design.ctrl.timepoint,contrast=contr.ctrl.4)
cam.ctrl.7 <- camera(countsmCtrl,idx.ctrl,design.ctrl.timepoint,contrast=contr.ctrl.7)
cam.ctrl.1.7 <- camera(countsmCtrl,idx.ctrl,design.ctrl.timepoint,contrast=contr.ctrl.1.7)

sign.btm.ctrl.4 <- cam.ctrl.4[cam.ctrl.4$FDR<0.05,]
sign.btm.ctrl.4 <- sign.btm.ctrl.4[sign.btm.ctrl.4$NGenes>20,]
dim(sign.btm.ctrl.4) # 36 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 4 

sign.btm.ctrl.4 <- sign.btm.ctrl.4%>%dplyr::filter(!str_detect(rownames(sign.btm.ctrl.4), "TBA"))
dim(sign.btm.ctrl.4) # 35 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 4 excluding TBA BTMs

sign.btm.ctrl.7 <- cam.ctrl.7[cam.ctrl.7$FDR<0.05,]
sign.btm.ctrl.7 <- sign.btm.ctrl.7[sign.btm.ctrl.7$NGenes>20,]
dim(sign.btm.ctrl.7) # 34 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 7

sign.btm.ctrl.7 <- sign.btm.ctrl.7%>%dplyr::filter(!str_detect(rownames(sign.btm.ctrl.7), "TBA"))
dim(sign.btm.ctrl.7) # 33 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 7 excluding TBA BTMs

sign.btm.ctrl.1.7 <- cam.ctrl.1.7[cam.ctrl.1.7$FDR<0.05,]
sign.btm.ctrl.1.7 <- sign.btm.ctrl.1.7[sign.btm.ctrl.1.7$NGenes>20,]
dim(sign.btm.ctrl.1.7)  # 31 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 7

sign.btm.ctrl.1.7 <- sign.btm.ctrl.1.7%>%dplyr::filter(!str_detect(rownames(sign.btm.ctrl.1.7), "TBA"))
dim(sign.btm.ctrl.1.7)  # 30 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 7 excluding TBA BTMs

## 3.4 Disease progression in the treatment group

design.treat.timepoint <- model.matrix(~0+day, data=colDataTreat)

voom.treat <- voom(countsTreat, design=design.treat.timepoint, plot=TRUE)

geneIDs.treat <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(voom.treat), columns = c("SYMBOL","GENEID"))

idx.treat <- ids2indices(btm$genesets,id=geneIDs.treat$SYMBOL)

contr.treat.4 <- makeContrasts(
  Treat4vsTreat1 = day4 - day1, 
  levels = colnames(design.treat.timepoint))

contr.treat.7 <- makeContrasts(
  Treat4vsTreat1 = day7 - day4, 
  levels = colnames(design.treat.timepoint))

contr.treat.1.7 <- makeContrasts(
  Treat7vsTreat1 = day7 - day1, 
  levels = colnames(design.treat.timepoint))


cam.treat.4 <- camera(countsmTreat,idx.treat,design.treat.timepoint,contrast=contr.treat.4)
cam.treat.7 <- camera(countsmTreat,idx.treat,design.treat.timepoint,contrast=contr.treat.7)
cam.treat.1.7 <- camera(countsmTreat,idx.treat,design.treat.timepoint,contrast=contr.treat.1.7)

sign.btm.treat.4 <- cam.treat.4[cam.treat.4$FDR<0.05,]
sign.btm.treat.4 <- sign.btm.treat.4[sign.btm.treat.4$NGenes>20,]
dim(sign.btm.treat.4) # 37 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 4

sign.btm.treat.4 <- sign.btm.treat.4%>%dplyr::filter(!str_detect(rownames(sign.btm.treat.4), "TBA"))
dim(sign.btm.treat.4) # 35 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 4

sign.btm.treat.7 <- cam.treat.7[cam.treat.7$FDR<0.05,]
sign.btm.treat.7 <- sign.btm.treat.7[sign.btm.treat.7$NGenes>20,]
dim(sign.btm.treat.7) # 21 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 7 

sign.btm.treat.7 <- sign.btm.treat.7%>%dplyr::filter(!str_detect(rownames(sign.btm.treat.7), "TBA"))
dim(sign.btm.treat.7) # 20 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 7 excluding TBA BTMs

sign.btm.treat.1.7 <- cam.treat.1.7[cam.treat.1.7$FDR<0.05,]
sign.btm.treat.1.7 <- sign.btm.treat.1.7[sign.btm.treat.1.7$NGenes>20,]
dim(sign.btm.treat.1.7) # 30 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 7 

sign.btm.treat.1.7 <- sign.btm.treat.1.7%>%dplyr::filter(!str_detect(rownames(sign.btm.treat.1.7), "TBA"))
dim(sign.btm.treat.1.7) # 30 Blood Transcriptional Modules are differentially expressed between ivermectin and placebo groups in the timepoint 7 excluding TBA BTMs

# 4. Data visualization
## 4.1 Differences on day 4 between the ivermectin-treated and placebo groups (adjusting for baseline)
### 4.1.1 Volcano plot of differentially Expressed Genes: DESeq

colors <- c("grey", "red", "blue")

Volcano.4 <- res.DESeq.timepoint.4 %>% mutate(gene = rownames(res.DESeq.timepoint.4), logp = -(log10(pvalue)), logadjp = -(log10(padj)), FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(pvalue<0.05 & log2FoldChange > 1, "Overexpressed in the ivermectin group", ifelse(pvalue<0.05& log2FoldChange < (-1), "Underexpressed in the ivermectin group","Non significant")))

Volcano.4.plot <- ggplot(data=Volcano.4, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs in the timepoint 4 (Adjusting for the baseline)") +
     theme(plot.title = element_text(size = 12))

significant_genes.4 <- Volcano.4[Volcano.4$sig != "Non significant", ]
significant_genes.4 <- distinct(significant_genes.4, rownames(significant_genes.4), .keep_all = TRUE)

diffDESeq.4.Volcano <- Volcano.4.plot +
  geom_text_repel(data = significant_genes.4, aes(label = c("RHBDF1", "CCDC113", "ABHD14A-ACY1", "TTC6", "GALNT13", "MEGF10", "DAAM2", "LGI4", "DNAAF1", "ADAMTS3", "VWA5B1", "B3GALT1", "CCDC184", "C6orf132", "RBM34", "BNIP5", "NRARP", "CD177", "MT-TH", "ADAMTS7P4", "lncRNA(ENSG00000226994)", "lncRNA(ENSG00000228061)", "LOC112268276", "TONSL-AS1", "LINC01765", "lncRNA(ENSG00000242593)", "AADACL2-AS1", "LINC01322", "lncRNA(ENSG00000251031)", "lncRNA(ENSG00000254587)", "lncRNA(ENSG00000256084)", "lncRNA(ENSG00000260126)", "MANEA-DT", "MAFTRR", "TMEM170A-CFDP1", "lncRNA(ENSG00000263990)", "lncRNA(ENSG00000270933)", "lncRNA(ENSG00000272084)", "lncRNA(ENSG00000272720)", "lncRNA(ENSG00000278419)", "IQCJ-SCHIP1", "lncRNA(ENSG00000285216", "lncRNA(ENSG00000290108)", "GBP1P1")),
                  box.padding = 0.5, point.padding = 0.2,
                  max.overlaps = Inf,
                  size = 3) + 
  guides(col = guide_legend(override.aes = list(size = 2)))

print(diffDESeq.4.Volcano)

### 4.1.2 Dot plot of Enrichment Analysis: Gene Set Enrichment Analysis (Gene Ontology) on day 4 (Adjusting for the baseline)

  plot.GO.4 <- ggplot(res.enrichGO.under.4, aes(x = Count, y = Description, size = Count, color = -log10(pvalue))) +
  geom_point() +
  scale_size_continuous(range = c(2), breaks = c(1)) +
  scale_color_gradient(low = "red2", high = "steelblue2") +
  labs(x = "Counts of underexpressed genes in ivermectin group", y = "Gene Ontology Term") +
  theme_minimal() +
  theme( 
    plot.title = element_text(size = 10),  
    axis.title = element_text(size = 10),    
    axis.text = element_text(size = 10),   
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  ) +
  ggtitle("Enrichment Analysis using GO terms on day 4 (Adjusting for the baseline)") +
  scale_x_continuous(breaks = c(2), labels = c("2"))

print(plot.GO.4)

# 4.1.3 Enrichment Analysis: limma (voom + camera + Blood Transcriptional Modules) on day 4 adjusting for the baseline

annotationbtms<-read.csv("BTM_annotation_revision.csv", sep=";")
colnames(annotationbtms)[3]<-"Geneset"

timepoint.4.res <- sigs.4 %>% inner_join(annotationbtms, by = "Geneset")
dim(timepoint.4.res)

timepoint.4.res$NGenes<-as.numeric(timepoint.4.res$NGenes)
timepoint.4.res <- timepoint.4.res[order(timepoint.4.res$High.level.annotation.group),]
timepoint.4.res$High.level.annotation.group[grep("platelet", rownames(timepoint.4.res))]<-"Platelets"
timepoint.4.res$Geneset<-factor(timepoint.4.res$Geneset, levels=c(timepoint.4.res$Geneset))

EA.timepoint.4 <- timepoint.4.res %>% 
  ggplot(aes(x = Direction, y = Geneset, size = NGenes, color = -log10(FDR)), size = 4) +
  geom_point() + 
  scale_x_discrete(labels = c('Enriched paths in placebo', 'Enriched paths in ivermectin')) + 
  theme(
    axis.text = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(7),  
  ) +
  scale_color_gradient(low = "red2", high = "steelblue2") +
  scale_size_area(max_size = 10) +
  scale_size_continuous(breaks = c(50, 100, 150, 200, 250, 300)) +
  labs(
    title = "Enrichment Analysis using BTMs on day 4 adjusting for the baseline"
  )

plot(EA.timepoint.4)

## 4.2 Differences on day 7 between the ivermectin-treated and placebo groups (adjusting for baseline)
### 4.2.1 Volcano plot of differentially Expressed Genes: DESeq

Volcano.7 <- res.DESeq.timepoint.7 %>% mutate(gene = rownames(res.DESeq.timepoint.7), logp = -(log10(pvalue)), logadjp = -(log10(padj)), FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -(2^abs(log2FoldChange)))) %>%
                   mutate(sig = ifelse(pvalue<0.05 & log2FoldChange > 1, "Overexpressed in the ivermectin group", ifelse(pvalue<0.05& log2FoldChange < (-1), "Underexpressed in the ivermectin group","Non significant")))

Volcano.7.plot <- ggplot(data=Volcano.7, aes(x=log2FoldChange, y=logp )) +
     geom_point(alpha = 1, size= 1, aes(col = sig)) + 
     scale_color_manual(values = colors) +
     xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) + labs(col=" ") + 
     geom_vline(xintercept = 1, linetype= "dotted") + geom_vline(xintercept = -1, linetype= "dotted") + 
     geom_hline(yintercept = -log10(0.1), linetype= "dotted")  +  theme_bw() + 
     labs(title = "Volcano Plot of DEGs in the timepoint 7 (Adjusting for the baseline)") +
     theme(plot.title = element_text(size = 12))

significant_genes.7 <- Volcano.7[Volcano.7$sig != "Non significant", ]
significant_genes.7 <- distinct(significant_genes.7, rownames(significant_genes.7), .keep_all = TRUE)

genes.7 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(significant_genes.7), columns = c("SYMBOL","GENEID"))

int <- intersect(genes.7$GENEID, rownames(significant_genes.7))
significant_genes.7 <- significant_genes.7[int,]

genes.7 <- genes.7[genes.7$GENEID %in% int, ]

selected_genes <- significant_genes.7[genes.7$SYMBOL %in% rownames(significant_genes.7), ]

diffDESeq.7.Volcano <- Volcano.7.plot +
  geom_text_repel(data = significant_genes.7, aes(label = genes.7$SYMBOL),
                  box.padding = 0.5, point.padding = 0.2,
                  max.overlaps = Inf,
                  size = 2) + 
  guides(col = guide_legend(override.aes = list(size = 2)))

print(diffDESeq.7.Volcano)

### 4.2.2 Dot plot of Enrichment Analysis: Gene Set Enrichment Analysis (Gene Ontology) on day 7(Adjusting for the baseline)

combined_timepoint.7 <- rbind(res.enrichGO.under.7, res.enrichGO.over.7)
combined_timepoint.7$Direction <- c(rep("Enriched paths in placebo", nrow(res.enrichGO.over.7)), 
                                     rep("Enriched paths in ivermectin", nrow(res.enrichGO.under.7)))

combined_timepoint.7$Direction <- factor(
  combined_timepoint.7$Direction,
  levels = c("Enriched paths in placebo", "Enriched paths in ivermectin")
)

plot.GO.7 <- ggplot(combined_timepoint.7, aes(x = Direction, y = Description, size = Count, color = -log10(pvalue))) +
  geom_point() +
  scale_size_continuous(range = c(1,5), breaks = c(1,3,5,7,9), labels = c("1", "3", "5", "7", "9")) +    
  scale_color_gradient(low = "red2", high = "steelblue2") +
  labs(x = "Enriched pathways", y = "Gene Ontology Term", color = "-log10(P-Value)") +
  theme_minimal() +
  theme( 
    plot.title = element_text(size = 12),  
    axis.title = element_text(size = 10),    
    axis.text = element_text(size = 10),   
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ) +
  ggtitle("Enrichment Analysis using GO terms on day 7 (Adjusting for the baseline)")

  # 4.2.3 Enrichment Analysis: limma (voom + camera + Blood Transcriptional Modules) in the timepoint 7 adjusting for the baseline

  timepoint.7.res <- sigs.7 %>% inner_join(annotationbtms, by = "Geneset")
dim(timepoint.7.res)

timepoint.7.res$NGenes<-as.numeric(timepoint.7.res$NGenes)
timepoint.7.res <- timepoint.7.res[order(timepoint.7.res$High.level.annotation.group),]
timepoint.7.res$High.level.annotation.group[grep("platelet", timepoint.7.res$Geneset)]<-"Platelets"
timepoint.7.res$Geneset<-factor(timepoint.7.res$Geneset, levels=c(timepoint.7.res$Geneset))

EA.timepoint.7 <- timepoint.7.res %>% 
  ggplot(aes(x=Direction, y = (Geneset), size = NGenes, color = -log10(FDR)), size=4) +
  geom_point()+ scale_x_discrete(labels=c('Enriched paths in placebo', 'Enriched paths in ivermectin'))+ 
  theme(
    axis.text = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(8),  
  ) +
  scale_color_gradient(low = "red2", high = "steelblue2")+
  scale_size_area(max_size = 10) +
  scale_size_continuous(breaks = c(50, 100, 150, 200, 250)) +
  labs(
    title = "Enrichment Analysis using BTMs on day 7 (adjusting for the baseline)"

plot(EA.timepoint.7)

## 4.3 Differences in the progression of the ivermectin and placebo groups
### 4.3.1. Dot plot of Enrichment Analysis: Gene Set Enrichment Analysis (Gene Ontology) in the disease progression

intersection <- intersect(rownames(combined_timepoint.progr.ctrl), rownames(combined_timepoint.progr.treat))

combined_timepoint.progr.ctrl_d <- combined_timepoint.progr.ctrl[!(row.names(combined_timepoint.progr.ctrl) %in% intersection), ]

combined_timepoint.progr.treat_d <- combined_timepoint.progr.treat[!(row.names(combined_timepoint.progr.treat) %in% intersection), ]


combined_ctrl_treat <- rbind(combined_timepoint.progr.ctrl_d, combined_timepoint.progr.treat_d)

dim(combined_ctrl_treat)

combined_ctrl_treat_47 <- combined_ctrl_treat[combined_ctrl_treat$Direction=="Enriched paths on day 4 (vs day 1)" | combined_ctrl_treat$Direction=="Enriched paths on day 7 (vs day 1)",]

plot_comb_47 <- ggplot(combined_ctrl_treat_47, aes(x = Direction, y = Description, size = Count, color = -log10(pvalue))) +
  geom_point() +
  scale_size_continuous(range = c(1, 4), breaks = c(21, 22, 35), labels = c("21", "22", "30")) +
  scale_color_gradient(low = "red2", high = "steelblue2") +
  labs(x = "Enriched pathways", y = "Gene Ontology Term", color = "-log10(P-Value)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

print(plot_comb_47)


### 4.3.2. Enrichment Analysis: limma (voom + camera + Blood Transcriptional Modules) in the disease progression

intersection_BTMs <- intersect(merged_control_BTM_order$Geneset, merged_treat_BTM_order$Geneset)

merged_control_BTM_order_f <- merged_control_BTM_order %>% filter(!merged_control_BTM_order$Geneset %in% intersection_BTMs) %>% distinct(Geneset, .keep_all = TRUE)
merged_control_BTM_order_f$Group <- "Placebo"
merged_treat_BTM_order_f <- merged_treat_BTM_order %>% filter(!merged_treat_BTM_order$Geneset %in% intersection_BTMs) %>% distinct(Geneset, .keep_all = TRUE)
merged_treat_BTM_order_f$Group <- "Ivermectin"

ctrl_treat_BTMs <- rbind(merged_control_BTM_order_f, merged_treat_BTM_order_f)

dim(ctrl_treat_BTMs)

ctrl_treat_BTMs_47 <- ctrl_treat_BTMs[-c(1,2),]

plot_BTMs_47 <- ggplot(ctrl_treat_BTMs_47, aes(x = Direction, y = Geneset, size = NGenes, color = -log10(FDR), shape = Group)) +
  geom_point() +
  scale_size_continuous(range = c(1, 5), breaks = c(22, 24, 25, 34, 43), labels = c("22", "24", "25", "34", "43")) +
  scale_color_gradient(low = "red2", high = "steelblue2") +
  labs(x = "Enriched pathways", y = "Gene Ontology Term", color = "-log10(FDR)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

print(plot_BTMs_47)

## 4.4 Pathways in common in the progression of the ivermectin and placebo groups
### 4.4.1. Dot plot of Enrichment Analysis: Gene Set Enrichment Analysis (Gene Ontology) in the disease progression

combined_ctrl_treat_common <- combined_timepoint.progr.ctrl[intersection,]
dim(combined_ctrl_treat_common)

combined_ctrl_treat_common_f <- combined_ctrl_treat_common %>%
  filter(str_detect(Description, "(virus|immune|viral|cytokine|interferon|inflammatory|leukocyte|myeloid|mononuclear|lymphocyte|T cell|immunoglobulin)"))

plot_comb_47_common <- ggplot(combined_ctrl_treat_common_f, aes(x = Direction, y = Description, size = Count, color = -log10(pvalue))) +
  geom_point() +
  scale_size_continuous(range = c(1, 5), breaks = c(30, 40, 50, 60, 70), labels = c("30", "40", "50", "60", "70")) +
  scale_color_gradient(low = "red2", high = "steelblue2") +
  labs(x = "Enriched pathways", y = "Gene Ontology Term", color = "-log10(P-Value)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

print(plot_comb_47_common)

### 4.4.2. Enrichment Analysis: limma (voom + camera + Blood Transcriptional Modules) in the disease progression

intersection_BTMs <- intersect(merged_control_BTM_order$Geneset, merged_treat_BTM_order$Geneset)
merged_treat_BTM_order_common_f <- merged_treat_BTM_order %>%
  filter(Geneset %in% intersection_BTMs)

merged_treat_BTM_order_common_f$NGenes<-as.numeric(merged_treat_BTM_order_common_f$NGenes)
merged_treat_BTM_order_common_f <- merged_treat_BTM_order_common_f[order(merged_treat_BTM_order_common_f$High.level.annotation.group),]
merged_treat_BTM_order_common_f$High.level.annotation.group[grep("platelet", rownames(merged_treat_BTM_order_common_f))]<-"Platelets"

plot_BTMs_47_common <- ggplot(merged_treat_BTM_order_common_f, aes(x = Direction, y = Geneset, size = NGenes, color = -log10(FDR))) +
  geom_point() +
  scale_size_continuous(range = c(1, 6), breaks = c(50, 100, 150, 200, 250, 300), labels = c("50", "100", "150", "200", "250", "300")) +
  scale_color_gradient(low = "red2", high = "steelblue2") +
  labs(x = "Enriched pathways", y = "Gene Ontology Term", color = "-log10(FDR)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

print(plot_BTMs_47_common)
