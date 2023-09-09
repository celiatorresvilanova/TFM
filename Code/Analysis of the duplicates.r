###########################################################################################################################################################
                                                             # Analysis of the duplicates                                                                                                                             
###########################################################################################################################################################

# Load needed packages

library(dupRadar)
library(ggplot2)
library(plotly)
library(dplyr)

# 1. DupRadar 
## 1.1 Comprehensive Duplication Analysis of all combined Samples

bamDuprm <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Duplicates/Aligned.sortedByCoord.out.mark.bam"	# the duplicate marked and aligned bam file
gtf <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/GRCh38.gtf"	# Gene model
stranded <- 1		# Stranded
paired   <- FALSE	#  Single-end reads
threads  <- 32		

dm.all.combined.samples <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
save(dm.all.combined.samples, file = "dm.all.RData")

par(mfrow=c(1,2))
duprateExpDensPlot(DupMat=dm.all.combined.samples)

## 1.1.1 Filtering criteria to remove lowly expressed genes that have high duplicate rates

dim(dm.all.combined.samples)
dm_f1 <- dm.all.combined.samples [!(dm.all.combined.samples$dupRateMulti>0.5 & dm.all.combined.samples$RPKMulti<10),]
dim(dm_f1)
dm_f2 <- dm_f1[!(dm_f1$dupRateMulti>0.7 & dm_f1$RPKMulti>=10 & dm_f1$RPKMulti<=100),]
dim(dm_f2)
dm_f3 <- dm_f2[!(dm_f2$dupRateMulti>0.75 & dm_f2$RPKMulti>=100 & dm_f2$RPKMulti<=500),]
dim(dm_f3)
dm_final <- dm_f3[!(dm_f3$dupRateMulti>0.8 & dm_f3$RPKMulti>=500 & dm_f3$RPKMulti<=1000),]
dim(dm_final)

## 1.2. Separate Duplication Analysis of all individual samples

bamDir <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Duplicates/Samples/"
gtf <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/GRCh38.gtf"
stranded <- 1
paired <- FALSE
threads <- 32

bamFiles <- list.files(bamDir, pattern = "\\.mark\\.bam$", full.names = TRUE)
pdf("dupRadar_plots.pdf", width = 12, height = 8)

for (bamFile in bamFiles) {
  dm.separated.samples <- analyzeDuprates(bamFile, gtf, stranded, paired, threads)
  sampleName <- tools::file_path_sans_ext(basename(bamFile))
  par(mfrow = c(1, 2))
  duprateExpDensPlot(DupMat = dm.separated.samples)
  title(main = sampleName)
}

dev.off()

## 1.3 Separate Duplication Analysis of problematic samples
### 1.3.1 First sequencing run

bamDir <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Duplicates/Problematic_samples/1st_read"
gtf <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/GRCh38.gtf"
stranded <- 1
paired <- FALSE
threads <- 32

bamFiles <- list.files(bamDir, full.names = TRUE)

pdf("problematic_samples_1st_read_dupRadar_plots.pdf", width = 12, height = 8)

for (bamFile in bamFiles) {
  dm.separated.samples <- analyzeDuprates(bamFile, gtf, stranded, paired, threads)
  sampleName <- tools::file_path_sans_ext(basename(bamFile))
  par(mfrow = c(1, 2))
  duprateExpDensPlot(DupMat = dm.separated.samples)
  title(main = sampleName)
}

dev.off()

### 1.3.2 Second sequencing run

bamDir <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Duplicates/Problematic_samples/2nd_read"
gtf <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/GRCh38.gtf"
stranded <- 1
paired <- FALSE
threads <- 32

bamFiles <- list.files(bamDir, full.names = TRUE)

pdf("problematic_samples_2nd_read_dupRadar_plots.pdf", width = 12, height = 8)

for (bamFile in bamFiles) {
  dm.separated.samples <- analyzeDuprates(bamFile, gtf, stranded, paired, threads)
  sampleName <- tools::file_path_sans_ext(basename(bamFile))
  par(mfrow = c(1, 2))
  duprateExpDensPlot(DupMat = dm.separated.samples)
  title(main = sampleName)
}

dev.off()

### 1.3.3 Third sequencing run

bamDir <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Duplicates/Problematic_samples/3rd_read"
gtf <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/GRCh38.gtf"
stranded <- 1
paired <- FALSE
threads <- 32

bamFiles <- list.files(bamDir, full.names = TRUE)

pdf("problematic_samples_3rd_read_dupRadar_plots.pdf", width = 12, height = 8)

for (bamFile in bamFiles) {
  dm.separated.samples <- analyzeDuprates(bamFile, gtf, stranded, paired, threads)
  sampleName <- tools::file_path_sans_ext(basename(bamFile))
  par(mfrow = c(1, 2))
  duprateExpDensPlot(DupMat = dm.separated.samples)
  title(main = sampleName)
}

dev.off()

# 2. Unsupervised Analysis: Principal Component Analysis
## 2.1 Principal Component Analysis to study sample aggregation based on the duplicates

load(file="assay_data_orig_f.RData")
load(file="assay_data_filt_genes_f.RData")
load(file="assay_data_dup_f.RData")

sampleT_orig_f <- apply(assay_data_orig_f,2,sum)/10^6
sampleTDF_orig_f <- data.frame(sample=names(sampleT_orig_f), total=sampleT_orig_f)
sampleT_filt_genes_f <- apply(assay_data_filt_genes_f,2,sum)/10^6
sampleT_dup_f <- apply(assay_data_dup_f,2,sum)/10^6
sampleTDF_dup_f <- data.frame(sample=names(sampleT_dup_f), total=sampleT_dup_f)

duplicates.n <- 1 - (sampleTDF_dup_f$total/sampleTDF_orig_f$total)
duplicates_cat <- as.factor(ntile(duplicates.n, 3))
duplicates_cat <- ifelse(duplicates_cat == 1, "Low duplicate rate", ifelse(duplicates_cat == 2, "Medium duplicate rate", "High duplicate rate"))
duplicates_cat <- factor(duplicates_cat, levels = c("High duplicate rate", "Medium duplicate rate", "Low duplicate rate"))

### 2.1.1 Original data

norm.expres.orig.d <- as.data.frame(apply(assay_data_orig_f, 1, function(x) (x/max(x, na.rm = TRUE))))
norm.expres.orig.d[is.na(norm.expres.orig.d)] <- 0
 
PCA.orig.d <- prcomp(norm.expres.orig.d,retx=TRUE, center=TRUE, scale=FALSE)
summary(PCA.orig.d)$importance[,c(1:3)] # PC1=0.24, PC2=0.09 and PC3=0.06

scores.data.orig.d <- data.frame(PCA.orig.d$x[, c("PC1", "PC2")])
scores.data.orig.d$duplicates <- as.factor(duplicates_cat)

scores.plot.orig.d <- ggplot(data = scores.data.orig.d, aes(x = PC1, y = PC2, colour = duplicates)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig.d)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig.d)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "Principal Component Analysis to study sample aggregation based on the duplication in the original data")+
  theme_bw()
ggplotly(scores.plot.orig.d)

### 2.1.2 Data excluding lowly expressed genes that have high duplication ratios

norm.expres.filt.genes.d <- as.data.frame(apply(assay_data_filt_genes_f, 1, function(x) (x/max(x, na.rm = TRUE))))
norm.expres.filt.genes.d[is.na(norm.expres.filt.genes.d)] <- 0
 
PCA.filt.genes.d <- prcomp(norm.expres.filt.genes.d,retx=TRUE, center=TRUE, scale=FALSE)
summary(PCA.filt.genes.d)$importance[,c(1:3)] # PC1=0.24, PC2=0.10 and PC3=0.06

scores.data.filt.genes.d <- data.frame(PCA.filt.genes.d$x[, c("PC1", "PC2")])
scores.data.filt.genes.d$duplicates <- as.factor(duplicates_cat)

scores.plot.filt.genes.d <- ggplot(data = scores.data.filt.genes.d, aes(x = PC1, y = PC2, colour = duplicates)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes.d)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes.d)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "Principal Component Analysis to study sample aggregation based on the duplication in the data excluding lowly expressed genes that have high duplication ratio")+
  theme(plot.title = element_text(size = 9))
ggplotly(scores.plot.filt.genes.d)

### 2.1.3 Data excluding all the duplicates

norm.expres.dup.d <- as.data.frame(apply(assay_data_dup_f, 1, function(x) (x/max(x, na.rm = TRUE))))
norm.expres.dup.d[is.na(norm.expres.dup.d)] <- 0
 
PCA.dup.d <- prcomp(norm.expres.dup.d,retx=TRUE, center=TRUE, scale=FALSE)
summary(PCA.dup.d)$importance[,c(1:3)] # PC1=0.54, PC2=0.10 and PC3=0.05

scores.data.dup.d <- data.frame(PCA.dup.d$x[, c("PC1", "PC2")])
scores.data.dup.d$duplicates <- duplicates_cat

scores.plot.dup.d <- ggplot(data = scores.data.dup.d, aes(x = PC1, y = PC2, colour = duplicates)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup.d)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup.d)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "Principal Component Analysis to study sample aggregation based on the duplication in the data excluding all the duplicates")+
  theme_bw()
ggplotly(scores.plot.dup.d)

## 2.2 Principal Component Analysis to study sample aggregation based on the duplicates and the total reads per sample
### 2.2.1 Original data 

scores.data.orig.d$reads_per_sample <- sampleT_orig_f

scores.plot.orig.d <- ggplot(data = scores.data.orig.d, aes(x = PC1, y = PC2, colour = reads_per_sample, shape = duplicates)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.orig.d)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.orig.d)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "Principal Component Analysis to study sample aggregation based on the duplicates and the total reads per sample in the original data")+
  theme(plot.title = element_text(size = 12))
ggplotly(scores.plot.orig.d)

shapiro.test(scores.data.orig.d$reads_per_sample) # The reads per sample in original data was normally distributed

anova_reads_original <- aov(reads_per_sample ~ duplicates, data = scores.data.orig.d)
summary(anova_reads_original)

posthoc_reads_original <- TukeyHSD(anova_reads_original)
print(posthoc_reads_original)

#### 2.2.2 Data excluding lowly expressed genes that have high duplication ratios  

scores.data.filt.genes.d$reads_per_sample <- sampleT_filt_genes_f

scores.plot.filt.genes.d <- ggplot(data = scores.data.filt.genes.d, aes(x = PC1, y = PC2, colour = reads_per_sample, shape = duplicates)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.filt.genes.d)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.filt.genes.d)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "Principal Component Analysis to study sample aggregation based on the duplicates and the total reads per sample in the data excluding lowly expressed genes that have high duplication ratio")+
  theme(plot.title = element_text(size = 7))
ggplotly(scores.plot.filt.genes.d)

shapiro.test(scores.data.filt.genes.d$reads_per_sample) # The reads per sample in original data was normally distributed

anova_reads_filt_genes <- aov(reads_per_sample ~ duplicates, data = scores.data.filt.genes.d)
summary(anova_reads_filt_genes)

posthoc_reads_filt_genes_orig <- TukeyHSD(anova_reads_filt_genes)
print(posthoc_reads_filt_genes_orig)

#### 2.2.3 Data excluding all the duplicates

scores.data.dup.d$reads_per_sample <- sampleT_dup_f

scores.plot.dup.d <- ggplot(data = scores.data.dup.d, aes(x = PC1, y = PC2, colour = reads_per_sample, shape = duplicates)) +
  geom_point(alpha = I(0.7), size = 4) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab(paste("PC1 (", round(summary(PCA.dup.d)$importance[2,1], 2) * 100, "%)"))+
  ylab(paste("PC2 (", round(summary(PCA.dup.d)$importance[2,2], 2) * 100, "%)"))+
  stat_ellipse() +
  labs(title = "Principal Component Analysis to study sample aggregation based on the duplicates and the total reads per sample in the data excluding all the duplicates")+
  theme(plot.title = element_text(size = 11))
ggplotly(scores.plot.dup.d)

shapiro.test(scores.data.dup.d$reads_per_sample) # The reads per sample in original data was not normally distributed

kruskal.test(scores.data.dup.d$reads_per_sample, scores.data.dup.d$duplicates)

group_high <- scores.data.dup.d$reads_per_sample[scores.data.dup.d$duplicates == "High duplicate rate"]
group_medium <- scores.data.dup.d$reads_per_sample[scores.data.dup.d$duplicates == "Medium duplicate rate"]
group_low <- scores.data.dup.d$reads_per_sample[scores.data.dup.d$duplicates == "Low duplicate rate"]

result_reads_high_low_dup <- wilcox.test(group_high, group_low)
print(result_reads_high_low_dup)

result_reads_high_medium_dup <- wilcox.test(group_high, group_medium)
print(result_reads_high_medium_dup)

result_reads_medium_low_dup <- wilcox.test(group_medium, group_low)
print(result_reads_medium_low_dup)

## 3. Data visualization of the difference in the total reads per sample after filtering
### 3.1 Original data  

load(file="assay_data_orig.RData")
sampleT_orig <- apply(assay_data_orig,2,sum)/10^6

sampleT_orig_f
range(sampleT_orig_f) # 33.4-39.1
differences_orig <- sampleT_orig - sampleT_orig_f
range(differences_orig) # 0.0196-0.173

palette <- scale_fill_gradient(low = "#A9D6E5", high="#012A4A")
sampleTDF_orig_f<- data.frame(sample=names(sampleT_orig_f), total=differences_orig)

ggplot(sampleTDF_orig_f, aes(x = sample, y = total, fill = total)) +
  geom_bar(stat = "identity") +
  palette +
  labs(title = "Differences in the total reads per sample in the filtered original data", x = "Sample names", y = "Differences in the total reads per sample after filtering") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

### 3.2 Data excluding lowly expressed genes that have high duplication ratios  

load(file="assay_data_filt_genes.RData")
sampleT_filt_genes <- apply(assay_data_filt_genes,2,sum)/10^6

sampleT_filt_genes_f
range(sampleT_filt_genes_f) # 33.3-39.1
differences_filt_genes <- sampleT_filt_genes - sampleT_filt_genes_f
range(differences_filt_genes) # 0.0178-0.169

sampleTDF_filt_genes_f<- data.frame(sample=names(sampleT_filt_genes_f), total=differences_filt_genes)

ggplot(sampleTDF_filt_genes_f, aes(x = sample, y = total, fill = total)) +
  geom_bar(stat = "identity") +
  palette +
  labs(title = "Differences in the reads per sample in the data excluding lowly expressed genes that have high duplication ratios", x = "Sample names", y = "Differences in the total reads per sample after filtering") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        plot.title = element_text(size = 15))

## 3.3 Data excluding all the duplicates 

load(file="assay_data_dup.RData")
sampleT_dup <- apply(assay_data_dup, 2, sum) / 10^6

range(sampleT_dup_f) # 1.17-9.92
differences_dup <- sampleT_dup - sampleT_dup_f
range(differences_dup) # 0.00960-0.0689

sampleTDF_dup_f<- data.frame(sample=names(sampleT_dup), total=differences_dup)

ggplot(sampleTDF_dup_f, aes(x = sample, y = total, fill = total)) +
  geom_bar(stat = "identity") +
  palette +
  labs(title = "Differences in the total reads per sample in the filtered data excluding all the duplicates", x = "Sample names", y = "Differences in the total reads per sample after filtering") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# 4. Descriptive analysis of the duplicates
## 4.1 Available RNA

load(file="ExpressionSet_orig.RData")
available.RNA <- pData(ExpressionSet_orig)$available.ng
available.RNA <- as.factor(ifelse(available.RNA < 1000, "Limited RNA availability", "Non-limited RNA availability"))

fisher.availableRNA <- fisher.test(available.RNA, duplicates_cat)
print(fisher.availableRNA) # P-value = 0.0022

colors <- c("#013A63", "#2A6F97", "#61A5C2")

availableRNA.dup.df <- data.frame(available.RNA, duplicates_cat)

ass.av.RNA <- ggplot(availableRNA.dup.df, aes(x = available.RNA, fill = duplicates_cat)) +
  geom_bar(position = "dodge") +
  xlab("Available RNA") +
  ylab("Frequency") +
  ggtitle("Association between RNA availability and duplicate rate") +
  scale_fill_manual(values = colors)

print(ass.av.RNA)

## 4.1.2 Available RNA and timepoints

timepoints <- pData(ExpressionSet_orig)$day

fisher.availableRNA.timepoints <- fisher.test(available.RNA, timepoints)
print(fisher.availableRNA.timepoints)  # P-value = 0.15

colors <- c("#013A63", "#61A5C2")

 available.RNA.timepoints.df <- data.frame(timepoints, available.RNA)

available.RNA.tmp <-ggplot(available.RNA.timepoints.df, aes(x = factor(timepoints), fill = available.RNA)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = colors) +
  labs(x = "Timepoints", y = "Count", title = "Association between available RNA and timepoints",
       fill = "Available RNA") +
  theme_minimal()

print(available.RNA.tmp)

## 4.2 RNA quality

quality.rna <- pData(ExpressionSet_orig)$quality.rna

shapiro.test(quality.rna) # The reads per sample in original data was not normally distributed

kruskal.test(quality.rna, duplicates_cat)  # P-value = 0.14

RNAquality.dup.df <- data.frame(quality.rna, duplicates_cat)

ggplot(RNAquality.dup.df, aes(x = duplicates_cat, y = quality.rna, fill = duplicates_cat)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#013A63", "#2A6F97", "#61A5C2")) +
  xlab("Duplicate Rate") +
  ylab("RNA quality") +
  ggtitle("Association between duplicate rate and RNA quality")


## 4.3 Contaminations (Nanodrop 260/280)

nano260280 <- pData(ExpressionSet_orig)$nano260280

shapiro.test(nano260280) # The reads per sample in original data were not normally distributed

kruskal.test(nano260280, duplicates_cat)  # P-value = 0.041

group_high_260280 <- nano260280[duplicates_cat == "High duplicate rate"]
group_medium_260280 <-  nano260280[duplicates_cat == "Medium duplicate rate"]
group_low_260280 <-  nano260280[duplicates_cat == "Low duplicate rate"]

result_cont260180_high_low_dup <- wilcox.test(group_high_260280, group_low_260280)
print(result_cont260180_high_low_dup)   # P-value=0.02074

result_cont260280_high_med_dup <- wilcox.test(group_high_260280, group_medium_260280)
print(result_cont260280_high_med_dup)   # P-value = 0.048

result_cont260280_med_low_dup <- wilcox.test(group_medium_260280, group_low_260280)
print(result_cont260280_med_low_dup)    # P-value = 0.85

nano260280.dup.df <- data.frame(nano260280, duplicates_cat)

nano260280 <- ggplot(nano260280.dup.df, aes(x = duplicates_cat, y = nano260280, fill = duplicates_cat)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#013A63", "#2A6F97", "#61A5C2")) +
  xlab("Duplicate Rate") +
  ylab("nano260280") +
  ggtitle("Association between duplicate rate and nanodrop 260/280")

print(nano260280)


## 4.4 Contaminations (Nanodrop 260/230)

nano260230 <- pData(ExpressionSet_orig)$nano260230

shapiro.test(nano260230) # The reads per sample in original data were not normally distributed

kruskal.test(nano260230, duplicates_cat)  # P-value = 0.029

group_high_260230 <- nano260230[duplicates_cat == "High duplicate rate"]
group_medium_260230 <-  nano260230[duplicates_cat == "Medium duplicate rate"]
group_low_260230 <-  nano260230[duplicates_cat == "Low duplicate rate"]

result_cont260230_high_low_dup <- wilcox.test(group_high_260230, group_low_260230)
print(result_cont260230_high_low_dup)  # P-value = 0.017

result_cont260230_high_medium_dup <- wilcox.test(group_high_260230, group_medium_260230)
print(result_cont260230_high_medium_dup)  # P-value = 0.039

result_cont260230_medium_low_dup <- wilcox.test(group_medium_260230, group_low_260230)
print(result_cont260230_medium_low_dup)  # P-value = 0.51

nano260230.dup.df <- data.frame(nano260230, duplicates_cat)

nano260230 <- ggplot(nano260230.dup.df, aes(x = duplicates_cat, y = nano260230, fill = duplicates_cat)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#013A63", "#2A6F97", "#61A5C2")) +
  xlab("Duplicate Rate") +
  ylab("nano260230") +
  ggtitle("Association between duplicate rate and nanodrop 260/230")

print(nano260230)

## 4.5 Library preparation batch

library.batch <- pData(ExpressionSet_orig)$library.batch

fisher.library.batch <- fisher.test(library.batch, duplicates_cat)
print(fisher.library.batch) # P-value = 0.38

## 4.6 Treatment

treatment <- as.factor(pData(ExpressionSet_orig)$treat)

fisher.treatment<- fisher.test(treatment, duplicates_cat)
print(fisher.treatment) # P-value = 0.15

## 4.7 Timepoints

timepoints <- as.factor(pData(ExpressionSet_orig)$day)

fisher.timepoints<- chisq.test(timepoints, duplicates_cat)
print(fisher.timepoints) # P-value = 0.74 

## 4.8 RNA extraction batch

extraction.batch <- as.factor(pData(ExpressionSet_orig)$extraction.batch)

fisher.extr.batch<- fisher.test(extraction.batch, duplicates_cat)
print(fisher.extr.batch) # P-value = 0.82

## 4.9 Sex

sex <- as.factor(pData(ExpressionSet_orig)$sex)

fisher.sex <- fisher.test(sex, duplicates_cat)
print(fisher.sex)  # P-value = 0.42

## 4.10 Ages

ages <- as.numeric(pData(ExpressionSet_orig)$age)

age_groups <- ifelse(ages<=26, "Young", "Adult")

fisher.age <- fisher.test(age_groups, duplicates_cat)
print(fisher.age)  # P-value = 0.42

## 4.11 Patients

patients <- as.factor(pData(ExpressionSet_orig)$studyno)
fisher.patients <- fisher.test(patients, duplicates_cat)
print(fisher.patients)

contingency_table <- table(patients, duplicates_cat)

patients.dup <- mosaicplot(contingency_table, main = "Association between patients and duplicates",
           xlab = "Patients", ylab = "Duplicates",
           color = c("#013A63", "#2A6F97", "#61A5C2"))

print(patients.dup)  # P-value = 3.31e-06