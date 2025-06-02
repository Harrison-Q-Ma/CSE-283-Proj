library(tidyverse)    
library(magrittr)     
library(tibble)     

# ── 1) Load TPM matrix ────────────────────────────
tpm <- read.delim("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/pnas_normal_tpm.txt",
                  header = T, sep = "\t", stringsAsFactors = FALSE)

# ── 2) Compute mean & variance per gene ─────────────────────────────────────
hvg_stats <- tpm %>%
  rownames_to_column("gene_id")%>%
  # select HVG only in those that have present expression greater than 0.1
  filter(gene_id %in% qc_filtered$GeneID) %>% 
  rowwise() %>%
  mutate(
    avg_expr = mean(c_across(-gene_id)),  # mean TPM across all sample‐columns
    var_expr = var(c_across(-gene_id))    # variance of TPM across all sample‐columns
  ) %>%
  ungroup() %>%
  # compute CV² = var / (mean²)
  mutate(cv2 = var_expr / (avg_expr^2)) %>%
  select(gene_id, avg_expr, var_expr, cv2)

# ── 3) Rank genes by dispersion (cv2) ────────────────────────────────────────
hvg_ranked <- hvg_stats %>%
  arrange(desc(cv2))  # descending → highest CV² first

# ── 4) Choose top N HVGs ─────────────────────────────────────────
top_n <- 2000
top_hvgs <- hvg_ranked %>%
  slice(1:top_n) %>%
  pull(gene_id)

# ── 5) Save tables ────────────────────────────────────────
write_csv(hvg_ranked,"./1_bc.classifier/output/1_QC/HVG_stats_normal_all_genes_filtered.csv")
write_csv(hvg_ranked %>% slice(1:top_n),
          glue::glue("./1_bc.classifier/output/1_QC/normal_top_{top_n}_HVGs.csv"))

