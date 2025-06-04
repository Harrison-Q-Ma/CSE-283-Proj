library(dplyr)
library(tidyverse)
library(ggplot2)
metadata <- read.csv("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/optional/pnas_patient_info.csv")

# 1) Load TPM table
#    - adjust the file path as needed
tpm <- read.delim("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/pnas_tpm_96_nodup.txt",
                  header = FALSE, sep = "\t", stringsAsFactors = FALSE)
readcount <- read.delim("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/pnas_readcounts_96_nodup.txt",
                  header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# First column is gene IDs
rownames(tpm) <- tpm[,1]
tpm <- tpm[,-1]

# 2) Call detected where TPM > 5
detected <- tpm > 5

# 3a) Compute detection frequency per gene
detection_freq <- rowSums(detected) / ncol(tpm)

# Build QC table
qc <- data.frame(
  GeneID = rownames(tpm),
  DetectionFrequency = detection_freq,
  row.names = NULL,
  stringsAsFactors = FALSE
)

# 3b) Compute detection frequency per sample
detection_freq_sample <- colSums(detected) / nrow(tpm)

# Build QC table
qc_sample <- data.frame(
  SampleID = colnames(tpm),
  DetectionFrequency = detection_freq_sample,
  row.names = NULL,
  stringsAsFactors = FALSE
)


# 4a) Filter genes with freq ≥ 0.1
min_freq <- 0.1
qc_filtered <- subset(qc, DetectionFrequency >= min_freq)

# 4b) Write out results
write.csv(qc, "./2_recurrence.classifier/output/2_QC/exRNA_detection_frequency_all_genes.csv", row.names = T)
write.csv(qc_filtered, "./2_recurrence.classifier/output/2_QC/exRNA_detection_frequency_filtered.csv", row.names = T)
write.csv(qc_sample, "./2_recurrence.classifier/output/2_QC/exRNA_detection_frequency_all_samples.csv", row.names = T)

log_tpm <- log2(tpm + 1)
log_tpm_filtered <- log_tpm%>%rownames_to_column("gene_id")%>%filter(gene_id %in% qc_filtered$GeneID)%>%column_to_rownames("gene_id")
readcount_filtered <- readcount%>%filter(V1 %in% qc_filtered$GeneID)%>%column_to_rownames("V1")
write.csv(log_tpm, "./2_recurrence.classifier/output/2_QC/pnas_tpm_96_nodup_normalized.csv", row.names = T)
write.csv(log_tpm_filtered, "./2_recurrence.classifier/output/2_QC/pnas_tpm_96_nodup_normalized_filtered.csv", row.names = T)
write.csv(readcount_filtered, "./2_recurrence.classifier/output/2_QC/pnas_readcounts_96_nodup_filtered.csv", row.names = T)
