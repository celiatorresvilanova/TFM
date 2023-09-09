###########################################################################################################################################################
                                                             # Preprocessing and quality control                                                                                                                              
###########################################################################################################################################################

# Load needed packages

library(stringr)
library(gtools)
library(edgeR)
library(Rsubread)
library(dplyr)
library(Biobase)
library(ggplot2)
library(plotly)

# 1. Merge counts files into a single counts matrix
## 1.1 Original data

setwd("/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Quantifications")

counts_files_orig <- list.files(path = "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Quantifications", pattern = ".out.bam.txt")
file_names <- str_extract(counts_files_orig, "orig(.*)Aligned")
extract_names <- str_replace(file_names, "orig", "")
extract_names <- str_replace(extract_names, "Aligned", "")

merged_counts_orig <- read.table(counts_files_orig[1], col.names = c("gene_id", extract_names[1]))

for (i in 2:length(counts_files_orig)) {
  file <- read.table(counts_files_orig[i], col.names = c("gene_id", extract_names[i]))
  merged_counts_orig <- merge(merged_counts_orig, file, by = "gene_id", all = TRUE)
}

merged_counts_orig <- merged_counts_orig[-c(1:5),]
column_order <- c("gene_id", mixedsort(extract_names))
merged_counts_orig <- merged_counts_orig[,column_order]
dim(merged_counts_orig)

write.table(merged_counts_orig, file = "Complete_matrix/complete_counts_matrix_orig.txt", sep = "\t", row.names = FALSE)

## 1.2 Data excluding all the duplicates

setwd("/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Quantifications")

counts_files <- list.files(path = "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Quantifications", pattern = ".dup.bam.txt")
file_names <- str_extract(counts_files, "dup(.*)Aligned")
extract_names <- str_replace(file_names, "dup", "")
extract_names <- str_replace(extract_names, "Aligned", "")

merged_counts_dup <- read.table(counts_files[1], col.names = c("gene_id", extract_names[1]))

for (i in 2:length(counts_files)) {
  file <- read.table(counts_files[i], col.names = c("gene_id", extract_names[i]))
  merged_counts_dup <- merge(merged_counts_dup, file, by = "gene_id", all = TRUE)
}

merged_counts_dup <- merged_counts_dup[-c(1:5),]
column_order <- c("gene_id", mixedsort(extract_names))
merged_counts_dup <- merged_counts_dup[,column_order]
dim(merged_counts_dup)

write.table(merged_counts_dup, file = "Complete_matrix/complete_counts_matrix_dup.txt", sep = "\t", row.names = FALSE)

# 2. Load data and create an Expression Set
## 2.1 Original data 

assay_data_orig <- read.delim("/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Quantifications/Complete_matrix/complete_counts_matrix_orig.txt")
rownames(assay_data_orig) <- assay_data_orig$gene_id
assay_data_orig <- assay_data_orig[, -1]

phenodata <- read.csv("/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/saint.csv")
phenodata<- phenodata[phenodata$day %in% c(1,4,7),]
phenodata$extraction.batch <- rep(c(1, 2, 3, 4, 5, 6), times = c(12, 12, 12, 12, 12, 12))
phenodata$available.ng <- c(684.79, 483.20, 690.37, 1000, 690.52, 777.42,
                            1000, 480.00, 1000, 1000, 1000, 910.00, 
                            1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 1000, 343.85, 1000,
                            328.64, 931.68, 877.92, 1000, 1000, 615.22,
                            669.12, 1000, 1000, 1000, 1000, 598.58,
                            1000, 589.22, 833.35, 1000, 988.77, 1000,
                            522.58, 369.58, 447.78, 1000, 942.00, 891.33, 
                            294.15, 339.48, 1000, 1000, 1000, 1000, 
                            1000, 784.72, 647.02, 1000, 919.02, 999.26,
                            251.26, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 703.44, 1000, 1000, 1000)

phenodata$quality.rna <- c(9.2, 9.3, 9.1, 8.5, 9.1, 9.4, 9.1, 9.7, 9.2, 9.1,
                           8.8, 9.1, 8.3, 8.9, 9.3, 8.8, 9.4, 9.2, 8.2, 8.8,
                           8.8, 8.9, 9, 9.2, 9, 9.2, 9.1, 8.7, 9.5, 9.2, 8.6,
                           8.7, 8.5, 7.5, 8.4, 9.1, 8.7, 9.2, 9.2, 9, 9, 9.2,
                           9.1, 9.5, 9.4, 8, 7.9, 8.8, 8.4, 8.4, 8.9, 9, 9.4,
                           9.2, 8.6, 9.2, 9.3, 8.7, 9, 8.8, 8.9, 9.2, 8.9, 8.8,
                           8.7, 7.9, 8.7, 8.7, 9.1, 8.7, 8.9, 9)

phenodata$nano260230 <- c(0.18, 0.12, 0.17, 0.27, 0.17, 0.09, NA, NA, NA, NA, 
                          NA, NA, 0.22, 0.14, 0.33, 0.21, 0.33, 0.27, 0.32, 
                          0.33, 0.30, 0.07, 0.03, 0.17, 0.09, 0.07, 0.19, 0.43,
                          0.25, 0.08, 0.06, 0.18, 0.15, 0.47, 0.14, 0.06, 0.13,
                          0.06, 0.08, 0.14, 0.06, 0.07, 0.10, 0.06, 0.07, 0.19,
                          0.08, 0.13, 0.06, 0.05, 0.30, 0.29, 0.10, 0.16, 0.13,
                          0.07, 0.06, 0.24, 0.16, 0.12, 0.02, 0.26, 0.28, 0.39,
                          0.26, 0.34, 0.16, 0.27, 0.04, 0.14, 0.15, 0.40)

phenodata$nano260280 <- c(2.46, 2.41, 2.31, 2.28, 2.65, 2.26, NA, NA, NA, NA, 
                          NA, NA, 2.18, 2.36, 2.14, 2.14, 2.22, 2.17, 2.23,
                          2.21, 2.15, 2.36, 3.22, 2.26, 2.80, 2.54, 2.36, 2.18,
                          2.36, 2.53, 2.56, 2.19, 2.18, 2.17, 2.29, 2.40, 2.33,
                          2.64, 2.47, 2.26, 2.54, 2.43, 2.85, 3.17, 2.88, 2.31,
                          2.76, 2.36, 3.67, 3.35, 2.21, 2.15, 2.40, 2.24, 2.28,
                          2.35, 2.45, 2.21, 2.42, 2.33, 4.90, 2.30, 2.17, 2.16,
                          2.17, 2.23, 2.35, 2.21, 2.68, 2.33, 2.16, 2.13)

phenodata$library.batch <- rep(c(1, 2), times = c(24,48))

rownames(phenodata) <- colnames(assay_data_orig)
ExpressionSet_orig <- ExpressionSet(assayData = as.matrix(assay_data_orig), phenoData = AnnotatedDataFrame(phenodata))

## 2.2 Data excluding lowly expressed genes that have high duplication ratios 

load(file="dm_final.RData")
intersection <- intersect(dm_final$ID, rownames(assay_data_orig))
length(intersection)
assay_data_filt_genes <- assay_data_orig[intersection,] 
dim(assay_data_filt_genes) # 47,173 genes

ExpressionSet_filt_genes <- ExpressionSet(assayData = as.matrix(assay_data_filt_genes), phenoData = AnnotatedDataFrame(phenodata))

## 2.3 Data excluding all the duplicates 

assay_data_dup <- read.delim("/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Quantifications/Complete_matrix/complete_counts_matrix_dup.txt")
rownames(assay_data_dup) <- assay_data_dup$gene_id
assay_data_dup <- assay_data_dup[, -1]

ExpressionSet_dup <- ExpressionSet(assayData = as.matrix(assay_data_dup), phenoData = AnnotatedDataFrame(phenodata))

# 3. Filtration
## 3.1. Filtration by lowly expressed genes
### 3.1.1 Original data

dim(assay_data_orig) # 62,757 genes
keep <- rowSums(assay_data_orig > 10) >= 12 # at least 12 samples have 10 reads per gene
assay_data_orig_f <- assay_data_orig[keep,]
dim(assay_data_orig_f)  # 20,729 genes

### 3.1.2 Data excluding lowly expressed genes that have high duplication ratios

dim(assay_data_filt_genes) #  47,173 genes
keep <- rowSums(assay_data_filt_genes > 10) >= 12 # at least 12 samples have 10 reads per gene
assay_data_filt_genes_f <- assay_data_filt_genes[keep,]
dim(assay_data_filt_genes_f)  # 19,980 genes

### 3.1.3 Data excluding all the duplicates

dim(assay_data_dup) # 62,757 genes
keep <- rowSums(assay_data_dup > 10) >= 12 # at least 12 samples have 10 reads per gene
assay_data_dup_f <- assay_data_dup[keep,]
dim(assay_data_dup_f)  # 17,295 genes

## 3.2 Filtration of transcripts that translate to globin domains and chains
### 3.2.1 Original data

rownames(assay_data_orig_f) <- gsub ("\\..*", "", rownames(assay_data_orig_f))

globin_transcripts <- c("ENSG00000086506", "ENSG00000161544","ENSG00000165553", "ENSG00000188536", "ENSG00000206172", "ENSG00000206177", "ENSG00000213931", "ENSG00000213934", "ENSG00000223609", "ENSG00000244734", "ENSG00000260629","ENSG00000196565")

assay_data_orig_f <- assay_data_orig_f[!(rownames(assay_data_orig_f) %in% globin_transcripts), ]
dim(assay_data_orig_f) # 20,722 genes

### 3.2.2 Data excluding lowly expressed genes that have high duplication ratios 

rownames(assay_data_filt_genes_f) <- gsub ("\\..*", "", rownames(assay_data_filt_genes_f))

assay_data_filt_genes_f <- assay_data_filt_genes_f[!(rownames(assay_data_filt_genes_f) %in% globin_transcripts), ]
dim(assay_data_filt_genes_f) # 19,973 genes

### 3.2.3 Data excluding all the duplicates

rownames(assay_data_dup_f) <- gsub ("\\..*", "", rownames(assay_data_dup_f))

assay_data_dup_f <- assay_data_dup_f[!(rownames(assay_data_dup_f) %in% globin_transcripts), ]
dim(assay_data_dup_f) # 17,288 genes

# 4. Data Visualization
## 4.1 Total reads per sample 
### 4.1.1 Histogram of the distribution of total reads per sample
#### 4.1.1.1 Original data

orig_palette <- c("#2C7DA0", "#013A63", "#2A6F97", rep("#2C7DA0", 2), "#2A6F97",
                  "#61A5C2", "#2C7DA0", rep("#2A6F97", 4), "#2C7DA0", 
                  rep("#2A6F97", 2), "#89C2D9", rep("#2A6F97", 2), "#89C2D9", 
                  rep("#2C7DA0", 3), rep("#2A6F97", 4), "#2C7DA0",
                  "#2A6F97", "#013A63", "#2C7DA0", "#61A5C2", "#2C7DA0",
                  "#2A6F97", "#2A6F97", rep("#2C7DA0", 2), "#014F86", "#61A5C2",
                  rep("#2A6F97", 4), "#014F86", "#2A6F97", "#014F86", "#2C7DA0",
                  rep("#2A6F97", 2), rep("#014F86", 2), "#013A63",
                  rep("#2C7DA0", 3), "#2A6F97", rep("#014F86", 2), "#61A5C2",
                  "#2C7DA0", "#014F86", "#2C7DA0", rep("#2A6F97", 2), "#2C7DA0",
                  "#2A6F97", "#014F86", "#2C7DA0", "#014F86", "#013A63", 
                  "#2A6F97", "#012A4A","#014F86")
 
sampleT_orig <- apply(assay_data_orig,2,sum)/10^6
sampleT_orig

orig_df <- data.frame(Sample = names(sampleT_orig), Value = sampleT_orig)

theme_set(theme_minimal(base_size = 15))

hist_orig <- ggplot(orig_df, aes(x = Value)) +
  geom_histogram(aes(fill = Sample), binwidth = 1, color = "black") +
  scale_fill_manual(values = orig_palette) +
  labs(title = "Histogram of reads per sample in the original data",
       x = "Reads per sample", y = "Frequency")

print(hist_orig)

#### 4.1.1.2 Data excluding lowly expressed genes that have high duplication ratios

filt_genes_palette <- c("#2C7DA0", "#012A4A", "#2A6F97", rep("#2C7DA0", 2), "#2A6F97",
                  "#61A5C2", "#2C7DA0", rep("#2A6F97", 4), "#2C7DA0", 
                  rep("#2A6F97", 2), "#89C2D9", rep("#2A6F97", 2), "#89C2D9", 
                  rep("#2C7DA0", 3), rep("#2A6F97", 4), "#2C7DA0",
                  "#014F86", "#013A63", "#2C7DA0", "#61A5C2", "#2C7DA0",
                  "#2A6F97", "#2A6F97", rep("#2C7DA0", 2), "#014F86", "#2C7DA0",
                  rep("#2A6F97", 4), "#014F86", "#2A6F97", "#014F86", "#2C7DA0",
                  rep("#2A6F97", 2), rep("#014F86", 2), "#013A63",
                  rep("#2C7DA0", 3), "#2A6F97", rep("#014F86", 2), "#61A5C2",
                  "#2C7DA0", "#014F86", "#2C7DA0", "#014F86", "#2A6F97", "#2C7DA0",
                  "#2A6F97", "#014F86", "#2C7DA0", "#014F86", "#013A63", 
                  "#2A6F97", "#012A4A","#014F86")
 

sampleT_filt_genes <- apply(assay_data_filt_genes,2,sum)/10^6
sampleT_filt_genes

filt_genes_df <- data.frame(Sample = names(sampleT_filt_genes), Value = sampleT_filt_genes)

theme_set(theme_minimal(base_size = 15))

hist_filt_genes <- ggplot(filt_genes_df, aes(x = Value)) +
  geom_histogram(aes(fill = Sample), binwidth = 1, color = "black") +
  scale_fill_manual(values = filt_genes_palette) +
  labs(title = "Histogram of reads per sample in the data excluding lowly expressed genes that have high ratio of duplicates",
       x = "Reads per sample", y = "Frequency")

print(hist_filt_genes)

#### 4.1.1.3 Data excluding all the duplicates 

sampleT_dup <- apply(assay_data_dup, 2, sum) / 10^6
sampleT_dup

dup_palette <- c("#2C7DA0", "#014F86", "#2C7DA0", "#2A6F97", "#61A5C2",
                 "#89C2D9", "#61A5C2", "#89C2D9", "#61A5C2", "#2C7DA0",
                 "#468FAF", rep("#2C7DA0", 2), "#61A5C2", "#2C7DA0",
                 "#89C2D9", "#468FAF", "#61A5C2", "#89C2D9", 
                 rep("#A9D6E5",2), "#2C7DA0", "#01497C", "#2C7DA0",  
                 "#468FAF", "#2C7DA0", "#61A5C2", "#013A63", "#01497C", 
                 rep("#61A5C2", 2),  rep("#89C2D9", 2), rep("#468FAF", 2),
                 rep("#61A5C2", 2), "#468FAF", "#2C7DA0", "#468FAF", "#89C2D9",
                 "#61A5C2", "#013A63", "#89C2D9", "#A9D6E5", "#89C2D9",
                 "#2C7DA0", "#468FAF", "#2C7DA0", rep("#2A6F97", 2),
                 "#89C2D9", "#468FAF", rep("#61A5C2", 2), "#468FAF", "#61A5C2",
                 rep("#89C2D9", 3), "#468FAF", rep("#61A5C2", 4), "#2C7DA0",
                 rep("#61A5C2", 2), "#01497C", "#468FAF", "#012A4A", "#89C2D9")

dup_df <- data.frame(Sample = names(sampleT_dup), Value = sampleT_dup)

theme_set(theme_minimal(base_size = 15))

hist_dup <- ggplot(dup_df, aes(x = Value)) +
  geom_histogram(aes(fill = Sample), binwidth = 1, color = "black") +
  scale_fill_manual(values = dup_palette) +
  labs(title = "Histogram of reads per sample in the data excluding the duplicates",
       x = "Reads per sample", y = "Frequency")

print(hist_dup)

### 4.1.2 Boxplots of the total reads per sample
#### 4.1.2.1 Original data

dataframe.orig <- data.frame(Value = sampleT_orig)

boxplot.orig <- ggplot(dataframe.orig, aes(y = Value)) +
  geom_boxplot(fill = "#4590caff", color = "black") +
  labs(title = "Boxplot of the reads per sample in the original data",
       x = "Reads per sample", y = "Frequency")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# Add labels for outliers
outliers.orig <- boxplot.stats(sampleT_orig)$out
outliers_orig_df <- data.frame(Value = outliers.orig, Label = names(outliers.orig))
boxplot.orig <- boxplot.orig + geom_text(data = outliers_orig_df, aes(x = 1, 
                                         y = Value, label = Label),
                                         color = "black", size = 7, 
                                         vjust = -0.5, nudge_y = -0.1,
                                         nudge_x = -0.96)

print(boxplot.orig)

#### 4.1.2.2 Data excluding lowly expressed genes that have high duplication ratios
dataframe.filt.genes <- data.frame(Value = sampleT_filt_genes)

boxplot.filt.genes <- ggplot(dataframe.filt.genes, aes(y = Value)) +
  geom_boxplot(fill = "#4590caff", color = "black") +
  labs(title = "Boxplot of the reads per sample in the data excluding lowly expressed genes that have high duplication ratios",
       x = "Reads per sample", y = "Frequency")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))


# Add labels for outliers

outliers.filt.genes <- boxplot.stats(sampleT_filt_genes)$out
outliers_filt_genes_df <- data.frame(Value = outliers.filt.genes, Label = names(outliers.filt.genes))
boxplot.filt.genes <- boxplot.filt.genes + geom_text(data = outliers_filt_genes_df, aes(x = 1, 
                                         y = Value, label = Label),
                                         color = "black", size = 7, 
                                         vjust = -0.5, nudge_y = -0.05,
                                         nudge_x = -0.96)

print(boxplot.filt.genes)

#### 4.1.2.3 Data excluding all the duplicates

dataframe.dup <- data.frame(Value = sampleT_dup)

boxplot.dup <- ggplot(dataframe.dup, aes(y = Value)) +
  geom_boxplot(fill = "#4590caff", color = "black") +
  labs(title = "Boxplot of the reads per sample in the data without duplicates",
       x = "Reads per sample", y = "Frequency")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# Add labels for outliers
outliers.dup <- boxplot.stats(sampleT_dup)$out
outliers_dup_df <- data.frame(Value = outliers.dup, Label = names(outliers.dup))
boxplot.dup <- boxplot.dup + geom_text(data = outliers_dup_df, aes(x = 1, 
                                         y = Value, label = Label),
                                         color = "black", size = 7, 
                                         vjust = -0.5, nudge_y = -0.05,
                                         nudge_x = -0.96)

print(boxplot.dup)

# 5. TMM Normalization
## 5.1 Original data

DGEListOrig <- DGEList(counts = assay_data_orig_f)
Norm.FactorOrig <- calcNormFactors(DGEListOrig, method = "TMM")
countsTMMOrig <- cpm(Norm.FactorOrig, log = T)

palette_TMM <- c("#A9D6E5", "#9ECDDF", "#94C8DC", "#89C2D9", "#66A2BD", "#5799b6", "#2C7DA0", "#2B769C", "#2A6F97", "#206793", "#165F8F", "#0C578B", "#075389",  "#01497C",  "#014270", "#013A63", "#013257", "#012A4A")
              
hist(countsTMMOrig[,1], xlab="log2-ratio", col = palette_TMM, main="TMM normalization of original data") 

## 5.2 Data excluding lowly expressed genes that have high duplication ratios 

DGEListfiltgenes <- DGEList(counts = assay_data_filt_genes_f)
Norm.FactorFiltGenes <- calcNormFactors(DGEListfiltgenes, method = "TMM")
countsTMMFiltGenes <- cpm(Norm.FactorFiltGenes, log = T)
              
hist(countsTMMFiltGenes[,1], xlab="log2-ratio", col = palette_TMM, main="TMM normalization of the data excluding lowly expressed genes that have high duplication ratios") 

### 5.3 Data excluding all the duplicates

DGEListDup <- DGEList(counts = assay_data_dup_f)
Norm.FactorDup <- calcNormFactors(DGEListDup, method = "TMM")
countsTMMDup <- cpm(Norm.FactorDup, log = T)

hist(countsTMMDup[,1], xlab="log2-ratio", col = palette_TMM, main="TMM normalization of data excluding duplicates")

# 6. Unsupervised Analysis: Principal Component Analysis
## 6.1. Duplications
### 6.1.1. Original data

sampleT_orig_f <- apply(assay_data_orig_f,2,sum)/10^6
sampleTDF_orig_f <- data.frame(sample=names(sampleT_orig_f), total=sampleT_orig_f)

sampleT_dup_f <- apply(assay_data_dup_f,2,sum)/10^6
sampleTDF_dup_f <- data.frame(sample=names(sampleT_dup_f), total=sampleT_dup_f)

duplicates.n <- 1 - (sampleTDF_dup_f$total/sampleTDF_orig_f$total)
duplicates_cat <- as.factor(ntile(duplicates.n, 3))
duplicates_cat <- ifelse(duplicates_cat == 1, "Low duplicate rate", ifelse(duplicates_cat == 2, "Medium duplicate rate", "High duplicate rate"))
duplicates_cat <- factor(duplicates_cat, levels = c("High duplicate rate", "Medium duplicate rate", "Low duplicate rate"))

norm.expres.orig <- as.data.frame(apply(countsTMMOrig, 1, function(x) (x/max(x, na.rm = TRUE))))
norm.expres.orig[is.na(norm.expres.orig)] <- 0
 
PCA.orig <- prcomp(norm.expres.orig,retx=TRUE, center=TRUE, scale=FALSE)
summary(PCA.orig)$importance[,c(1:3)] # PC1=0.26, PC2=0.14 and PC3=0.09

scores.data.orig <- data.frame(PCA.orig$x[, c("PC1", "PC2")])
scores.data.orig$duplicates <- as.factor(duplicates_cat)

scores.plot.orig <- ggplot(data = scores.data.orig, aes(x = PC1, y = PC2, colour = duplicates)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the duplicates in the the normalized original data")+
  theme_bw()
ggplotly(scores.plot.orig)

### 6.1.2 Data excluding lowly expressed genes that have high duplication ratios

norm.expres.filt.genes <- as.data.frame(apply(countsTMMFiltGenes, 1, function(x) (x/max(x, na.rm = TRUE))))
norm.expres.filt.genes[is.na(norm.expres.filt.genes)] <- 0
 
PCA.filt.genes <- prcomp(norm.expres.filt.genes,retx=TRUE, center=TRUE, scale=FALSE)
summary(PCA.filt.genes)$importance[,c(1:3)] # PC1=0.49, PC2=0.24 and PC3=0.09

scores.data.filt.genes <- data.frame(PCA.filt.genes$x[, c("PC1", "PC2")])
scores.data.filt.genes$duplicates <- as.factor(duplicates_cat)

scores.plot.filt.genes <- ggplot(data = scores.data.filt.genes, aes(x = PC1, y = PC2, colour = duplicates)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the duplicates in the the normalized data exlcluding lowly expressed genes that have high duplication ratio")+
  theme(plot.title = element_text(size = 8))
ggplotly(scores.plot.filt.genes)

### 6.1.3 Data excluding all the duplicates

norm.expres.dup <- as.data.frame(apply(countsTMMDup, 1, function(x) (x/max(x, na.rm = TRUE))))
norm.expres.dup[is.na(norm.expres.dup)] <- 0
 
PCA.dup <- prcomp(norm.expres.dup,retx=TRUE, center=TRUE, scale=FALSE)
summary(PCA.dup)$importance[,c(1:3)] # PC1=0.15, PC2=0.24 and PC3=0.06

scores.data.dup <- data.frame(PCA.dup$x[, c("PC1", "PC2")])
scores.data.dup$duplicates <- duplicates_cat

scores.plot.dup <- ggplot(data = scores.data.dup, aes(x = PC1, y = PC2, colour = duplicates)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the duplicates in the the normalized data exlcluding all the duplicates")+
  theme(plot.title = element_text(size = 12))
ggplotly(scores.plot.dup)

## 6.2 Treatment
### 6.2.1 Original data
treatment <- as.factor(pData(ExpressionSet_orig)$treat)

scores.plot.orig <- ggplot(data = scores.data.orig, aes(x = PC1, y = PC2, colour = treatment)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the treatment in the normalized original data")+
  theme_bw()
ggplotly(scores.plot.orig)

### 6.2.2 Data excluding lowly expressed genes that have high duplication ratios 

scores.data.filt.genes$treatment <-ifelse(treatment==1, "Treatment", "Control")
  
scores.plot.filt.genes <- ggplot(data = scores.data.filt.genes, aes(x = PC1, y = PC2, colour = treatment)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the treatment from the normalized data excluding lowly expressed genes that have high duplication ratios")+
  theme(plot.title = element_text(size = 7))
ggplotly(scores.plot.filt.genes)

### 6.2.3 Data excluding all the duplicates

scores.data.dup$treatment <-ifelse(treatment==1, "Treatment", "Control")

scores.plot.dup <- ggplot(data = scores.data.dup, aes(x = PC1, y = PC2, colour = treatment)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the treatment from the normalized data excluding all the duplicates")+
  theme(plot.title = element_text(size = 12))
ggplotly(scores.plot.dup)

## 6.3 Treatment and timepoints
### 6.3.1 Original data

timepoints <- as.factor(pData(ExpressionSet_orig)$day)

scores.plot.orig <- ggplot(data = scores.data.orig, aes(x = PC1, y = PC2, colour = treatment, shape = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the treatment and the timepoints from the normalized original data")+
  theme(plot.title = element_text(size = 12))
ggplotly(scores.plot.orig)

### 6.3.2 Data excluding lowly expressed genes that have high duplication ratios

scores.plot.filt.genes <- ggplot(data = scores.data.filt.genes, aes(x = PC1, y = PC2, colour = treatment, shape = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the treatment and timepoints from the normalized data excluding lowly expressed genes that have high duplication ratios")+
  theme(plot.title = element_text(size = 7))
ggplotly(scores.plot.filt.genes)

### 6.3.3 Data excluding all the duplicates

scores.plot.dup <- ggplot(data = scores.data.dup, aes(x = PC1, y = PC2, colour = treatment, shape = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the treatment and the timepoints from the normalized data excluding all the duplicates")+
  theme(plot.title = element_text(size = 10))
ggplotly(scores.plot.dup)

## 6.4 Timepoints
### 6.4.1 Original data

scores.plot.orig <- ggplot(data = scores.data.orig, aes(x = PC1, y = PC2, colour = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the timepoints from the normalized original data")+
  theme_bw()
ggplotly(scores.plot.orig)

### 6.4.2 Data excluding lowly expressed genes that have high duplication ratios

scores.plot.filt.genes <- ggplot(data = scores.data.filt.genes, aes(x = PC1, y = PC2, colour = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the timepoints from the normalized data excluding lowly expressed genes that have high duplication ratios")+
  theme(plot.title = element_text(size = 7))
ggplotly(scores.plot.filt.genes)

### 6.4.3 Data excluding lowly expressed genes that have high duplication ratios 

scores.plot.dup <- ggplot(data = scores.data.dup, aes(x = PC1, y = PC2, colour = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the timepoints from the normalized data excluding the duplicates")+
  theme(plot.title = element_text(size = 12))
ggplotly(scores.plot.dup)

## 6.5 Sex and timepoints
### 6.5.1 Original data

sex <- as.factor(pData(ExpressionSet_orig)$sex)
scores.data.orig$sex <-ifelse(sex==1, "Female", "Male")
  
scores.plot.orig <- ggplot(data = scores.data.orig, aes(x = PC1, y = PC2, colour = sex, shape = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the sex and the timepoints from the normalized original data")+
  theme_bw()
ggplotly(scores.plot.orig)

### 6.5.2 Data excluding lowly expressed genes that have high duplication ratios 

scores.data.filt.genes$sex <- ifelse(sex==1, "Female", "Male")

scores.plot.filt.genes <- ggplot(data = scores.data.filt.genes, aes(x = PC1, y = PC2, colour = sex, shape=timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the sex and the timepoints from the normalized data excluding lowly expressed genes that have high duplication ratios")+
  theme(plot.title = element_text(size = 6))
ggplotly(scores.plot.filt.genes)

### 6.5.3 Data excluding all the duplicates

scores.data.dup$sex <-ifelse(sex==1, "Female", "Male")

scores.plot.dup <- ggplot(data = scores.data.dup, aes(x = PC1, y = PC2, colour = sex, shape = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the sex and the timepoints of the normalized data excluding all the duplicates")+
  theme(plot.title = element_text(size = 11))
ggplotly(scores.plot.dup)

## 6.6 Patients and timepoints
### 6.6.1 Original data

patients <- as.factor(pData(ExpressionSet_orig)$studyno)

scores.plot.orig <- ggplot(data = scores.data.orig, aes(x = PC1, y = PC2, colour = patients, shape = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the patients and the timepoints from the normalized original data")+
  theme(plot.title = element_text(size = 11))
ggplotly(scores.plot.orig)

### 6.6.2 Data excluding lowly expressed genes that have high duplication ratios

scores.plot.filt.genes <- ggplot(data = scores.data.filt.genes, aes(x = PC1, y = PC2, colour = patients, shape=timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the patients and the timepoints from the normalized data excluding lowly expressed genes that have high duplication ratios")+
  theme(plot.title = element_text(size = 7))
ggplotly(scores.plot.filt.genes)

### 6.6.3 Data excluding all the duplicates

scores.plot.dup <- ggplot(data = scores.data.dup, aes(x = PC1, y = PC2, colour = patients, shape = timepoints)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the patients and the timepoints from the normalized data excluding all the duplicates")+
  theme(plot.title = element_text(size = 9))
ggplotly(scores.plot.dup)

## 6.7 Extraction batch
### 6.7.1 Original data

extraction.batch <- as.factor(pData(ExpressionSet_orig)$extraction.batch)
scores.plot.orig <- ggplot(data = scores.data.orig, aes(x = PC1, y = PC2, colour = extraction.batch)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the extraction batch from the normalized original data")+
  theme_bw()
ggplotly(scores.plot.orig)

### 6.7.2 Data excluding lowly expressed genes that have high duplication ratios

scores.plot.filt.genes <- ggplot(data = scores.data.filt.genes, aes(x = PC1, y = PC2, colour = extraction.batch)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the extraction batch from the normalized data excluding lowly expressed genes that have high duplication ratios")+
  theme(plot.title = element_text(size = 8))
ggplotly(scores.plot.filt.genes)

### 6.7.3 Data excluding all the duplicates

scores.plot.dup <- ggplot(data = scores.data.dup, aes(x = PC1, y = PC2, colour = extraction.batch)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the extraction batch of the normalized data excluding all the duplicates")+
  theme(plot.title = element_text(size = 12))

ggplotly(scores.plot.dup)

## 6.8 Ages
### 6.8.1 Original data

ages <- as.factor(pData(ExpressionSet_orig)$age)
scores.data.orig$ages <- as.factor(ntile(ages, 3))

scores.plot.orig <- ggplot(data = scores.data.orig, aes(x = PC1, y = PC2, colour = ages)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the ages from the normalized original data")+
  theme_bw()
ggplotly(scores.plot.orig)

### 6.8.2 Data excluding lowly expressed genes that have high duplication ratios

scores.data.filt.genes$ages <- as.factor(ntile(ages, 3))

scores.plot.filt.genes <- ggplot(data = scores.data.filt.genes, aes(x = PC1, y = PC2, colour = ages)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the ages from the normalized data excluding lowly expressed genes that have high duplication ratios")+
  theme(plot.title = element_text(size = 7))
ggplotly(scores.plot.filt.genes)

### 6.8.3 Data excluding all the duplicates

scores.data.dup$ages <-  as.factor(ntile(ages, 3))
scores.plot.dup <- ggplot(data = scores.data.dup, aes(x = PC1, y = PC2, colour = ages)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "PCA to study sample aggregation based on the ages of the normalized data excluding all the duplicates")+
  theme_bw()
ggplotly(scores.plot.dup)

