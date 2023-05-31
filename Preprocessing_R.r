
###################################################################################################################################################################
                                                                     # R PREPROCESSING                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
###################################################################################################################################################################


# Three readings were performed for each sample with the aim of achieving 40 million reads per sample. The data was merged using R version 4.3.0. The files 
# were read using the readFastq function from the ShortRead R package, and they were merged using the WriteFastq function from the same package.

########################################################################### MERGE FILES ###########################################################################

library(ShortRead) 

# The directory path containing the three read folders was set as the working directory
dir_path <- setwd("/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data")

# A list of filenames in the "1st read" folder with the suffix "001.fastq" was obtained
filenames <- list.files(file.path(dir_path, "1st read"), pattern = "001.fastq", full.names = FALSE)

for (filename in filenames) {
  # For each filename, the file paths for the corresponding files in the three read folders ("1st read", "2nd read", "3rd read") were constructed
  read1_file <- file.path(dir_path, "1st read", filename)
  read2_file <- file.path(dir_path, "2nd read", filename)
  read3_file <- file.path(dir_path, "3rd read", filename)
  read_total <- c(read1_file, read2_file, read3_file)
   
  # The existence of the corresponding files in "read1", "read2", and "read3" was checked
  if (file.exists(read1_file) && file.exists(read2_file) && file.exists(read3_file)) {
    merged_file_path = file.path(dir_path, "merged", filename)

# f the files existed, the merged file path was determined in the "merged" folder
for (fastq_files in read_total) {
    fq = readFastq(fastq_files)
    writeFastq(fq, merged_file_path, mode="a", compress=TRUE)
}
 
  }
}

########################################################################### DUPLICATIONS ###########################################################################

# After the files were merged, fastqc reports, trimming, and alignment were carried out using Linux 4.18.0-425.13.1.el8_7.x86_64 (github). The fastqc reports were
# examined, revealing a significant proportion of duplicate reads. In RNA-Seq, the distinction between artificial reads and normal read duplication caused by
# over-sequencing highly expressed genes poses a challenge. This is because duplicate reads are frequently observed in RNA-Seq data and cannot be exclusively 
# attributed to technical errors. Highly expressed genes often exceed the threshold of one read per base pair of the exon model, resulting in inevitable read
# duplication. To address this issue, the decision was made to mark these duplicates using Sambamba in Linux 4.18.0-425.13.1.el8_7.x86_64 (github). Subsequently, 
# the analysis was performed using dupRadar, a visualization tool designed for investigating duplications and gene expression patterns. DupRadar facilitates the
# examination of the relationship between the percentage of duplications and gene expression. A significant proportion of duplications in genes with low expression 
# levels was observed when Sambamba was used. Although the level of duplication increased with gene expression, it remained unreasonably high in genes with lower
# expression. Therefore, the duplicates were removed, and the analysis proceeded with two datasets: the original data and the data with duplicates eliminated.

# The dupRadar was conducted with the following parameters:
    # The BAM and GTF files were specified in the variables bamDuprm and gtf, respectively. The BAM file needs to have the duplicates marked.
    # Stranded was set to 1 to account for the alignment of reads in a strand-specific manner
    # Paired was set to FALSE since we were working with single-end reads

library(dupRadar)
library(Rsubread)

bamDuprm <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/merged_quality/Duplicates/Aligned.sortedByCoord.out.mark.bam" 
gtf <- "/PROJECTES/MALARIA_IMMUNO/22.SAINT/Transcriptomics/Data/GRCh38.gtf" 
stranded <- 1 
paired   <- FALSE 
threads  <- 32 

dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
par(mfrow=c(1,2))

duprateExpDensPlot(DupMat=dm)