library(tidyverse)

# ─────────────────────────────────────────────────────────────────────────────
# 1) LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
# 1a) Read in log2TPM + 1 matrix
#     - Assume file has header: first column = gene_id, then one column per sample
#     - Each row = one gene
log2tpm_df <- read.csv("./2_recurrence.classifier/output/2_QC/pnas_tpm_96_nodup_normalized_filtered.csv",sep=',',header = T)

# 1b) Convert gene_id column into rownames, then drop it
#     (so that `expr_mat` is a numeric matrix with rownames = gene IDs)
expr_mat <- log2tpm_df %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix()         # now each row is a gene, each col is a sample

# 1c) Read in metadata: one row per sample, with at least these two columns:
#       • sampleID    (must match exactly the column names of expr_mat)
#       • recurrence  (e.g. "R" for recurrence, "N" for non‐recurrence)
metadata <- read_csv("data/metadata.csv", col_names = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 2) MAKE SURE SAMPLE ORDER MATCHES BETWEEN expr_mat AND metadata
# ─────────────────────────────────────────────────────────────────────────────
# 2a) Subset or reorder metadata so that row order = column order of expr_mat
#     (drop any samples in metadata that aren’t in expr_mat, and vice versa)
common_samples <- intersect(colnames(expr_mat), metadata$sampleID)

# Subset expression matrix to common samples
expr_mat <- expr_mat[, common_samples]

# Subset metadata to common samples, then reorder rows to match expr_mat columns
metadata <- metadata %>%
  filter(sampleID %in% common_samples) %>%
  slice(match(common_samples, sampleID))

# At this point, 
#   colnames(expr_mat)[i] == metadata$sampleID[i]  for all i  
#   metadata$recurrence is in the same order as the columns of expr_mat

# Extract a simple grouping vector (“R” vs “N”) in the same order as columns
group_vec <- metadata$recurrence
# (If your column is named differently—e.g. “recurStatus”—change above accordingly.)

# ─────────────────────────────────────────────────────────────────────────────
# 3) LOOP OVER GENES TO COMPUTE log2FC and p‐VALUE
# ─────────────────────────────────────────────────────────────────────────────
gene_ids <- rownames(expr_mat)
n_genes  <- length(gene_ids)

# Pre-allocate vectors
log2FC   <- numeric(n_genes)   # will store mean(R) – mean(N)
p_value  <- numeric(n_genes)   # will store raw p‐values from t.test

for (i in seq_len(n_genes)) {
  # Extract the i‐th gene’s expression across all samples
  gene_expr <- expr_mat[i, ]    # a numeric vector of length = ncol(expr_mat)
  
  # Split into two groups by recurrence status
  expr_R <- gene_expr[group_vec == "R"]
  expr_N <- gene_expr[group_vec == "N"]
  
  # 3a) Compute log2FC = mean_R – mean_N
  log2FC[i] <- mean(expr_R, na.rm = TRUE) - mean(expr_N, na.rm = TRUE)
  
  # 3b) Perform two‐sample t‐test (unpaired, equal variance assumed by default)
  #     We catch warnings or errors (e.g. if a group has constant values).
  t_test_res <- try(t.test(expr_R, expr_N), silent = TRUE)
  if (inherits(t_test_res, "try-error")) {
    p_value[i] <- NA
  } else {
    p_value[i] <- t_test_res$p.value
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 4) ADJUST P-VALUES (BONFERRONI)
# ─────────────────────────────────────────────────────────────────────────────
p_adj <- p.adjust(p_value, method = "bonferroni")

# ─────────────────────────────────────────────────────────────────────────────
# 5) ASSEMBLE RESULTS INTO A DATA FRAME
# ─────────────────────────────────────────────────────────────────────────────
de_results <- tibble(
  gene     = gene_ids,
  log2FC   = log2FC,
  p_value  = p_value,
  p_adj    = p_adj
) %>%
  arrange(p_adj)      # sorted by smallest adjusted p‐value first

# ─────────────────────────────────────────────────────────────────────────────
# 6) OPTIONAL: SAVE TO CSV
# ─────────────────────────────────────────────────────────────────────────────
write_csv(de_results, "DE_results_log2TPM_ttest.csv")

# ─────────────────────────────────────────────────────────────────────────────
# 7) QUICK SUMMARY OF OUTPUT
# ─────────────────────────────────────────────────────────────────────────────
#   • de_results$log2FC:  
#       – positive means “higher in Recurrence (R)”  
#       – negative means “higher in Non‐recurrence (N)”  
#   • de_results$p_value:  raw t-test p-value  
#   • de_results$p_adj:    Bonferroni‐adjusted p-value  
# ─────────────────────────────────────────────────────────────────────────────
