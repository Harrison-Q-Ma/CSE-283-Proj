library(dplyr)
library(tidyverse)
library(ggplot2)

# 1) Load TPM table
#    - adjust the file path as needed
tpm <- read.delim("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/pnas_normal_tpm.txt",
                  header = T, sep = "\t", stringsAsFactors = FALSE)
readcount <- read.delim("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/pnas_normal_readcounts.txt",
                        header = T, sep = "\t", stringsAsFactors = FALSE)

# 2) Call detected where TPM > 5
detected <- tpm > 5

# 3) Compute detection frequency per gene
detection_freq <- rowSums(detected) / ncol(tpm)

# Build QC table
qc <- data.frame(
  GeneID = rownames(tpm),
  DetectionFrequency = detection_freq,
  row.names = NULL,
  stringsAsFactors = FALSE
)

# 4a) Filter genes with freq ≥ 0.1
min_freq <- 0.1
qc_filtered <- subset(qc, DetectionFrequency >= min_freq)

# 4b) Write out results
write.csv(qc, "./1_bc.classifier/output/1_QC/exRNA_detection_frequency_all_genes.csv", row.names = T)
write.csv(qc_filtered, "./1_bc.classifier/output/1_QC/exRNA_detection_frequency_filtered.csv", row.names = T)

log_tpm <- log2(tpm + 1)
log_tpm_filtered <- log_tpm%>%rownames_to_column("gene_id")%>%filter(gene_id %in% qc_filtered$GeneID)%>%column_to_rownames("gene_id")
readcount_filtered <- readcount%>%rownames_to_column("gene_id")%>%filter(gene_id %in% qc_filtered$GeneID)%>%column_to_rownames("gene_id")
write.csv(log_tpm, "./1_bc.classifier/output/1_QC/pnas_normal_tpm_normalized.csv", row.names = T)
write.csv(log_tpm_filtered, "./1_bc.classifier/output/1_QC/pnas_normal_tpm_normalized_filtered.csv", row.names = T)
write.csv(readcount_filtered, "./1_bc.classifier/output/1_QC/pnas_normal_readcounts_filtered.csv", row.names = T)
