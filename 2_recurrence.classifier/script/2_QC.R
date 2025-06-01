library(dplyr)
library(tidyverse)
library(ggplot2)
metadata <- read.csv("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/optional/pnas_patient_info.csv")

# 1) Load TPM table
#    - adjust the file path as needed
tpm <- read.delim("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/pnas_tpm_96_nodup.txt",
                  header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# First column is gene IDs
rownames(tpm) <- tpm[,1]
tpm <- tpm[,-1]

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
write.csv(qc, "./2_recurrence.classifier/output/exRNA_detection_frequency_all_genes.csv", row.names = FALSE)
write.csv(qc_filtered, "./2_recurrence.classifier/output/exRNA_detection_frequency_filtered.csv", row.names = FALSE)

log_tpm <- log2(tpm + 1)
write.csv(log_tpm, "./2_recurrence.classifier/output/pnas_tpm_96_nodup_normalized.csv", row.names = FALSE)
