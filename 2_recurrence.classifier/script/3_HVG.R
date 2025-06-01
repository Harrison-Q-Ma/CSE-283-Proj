library(tidyverse)    
library(magrittr)     
library(tibble)     

# ── 1) Load TPM matrix ────────────────────────────
tpm <- read_tsv("data/tpm_96_nodup.tsv", col_names = TRUE)

# ── 2) Compute mean & variance per gene ─────────────────────────────────────
hvg_stats <- tpm %>%
  rownames_to_column(var = "gene_id") %>%   
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

# ── 4) (Optional) Choose top N HVGs ─────────────────────────────────────────
top_n <- 2000
top_hvgs <- hvg_ranked %>%
  slice(1:top_n) %>%
  pull(gene_id)

# ── 5) (Optional) Save tables to disk ────────────────────────────────────────
write_csv(hvg_ranked,"./2_recurrence.classifier/output/HVG_stats_all_genes.csv")
write_csv(hvg_ranked %>% slice(1:top_n),
          glue::glue("top_{top_n}_HVGs.csv"))

