library(tidyverse)

# ─────────────────────────────────────────────────────────────────────────────
# 1) LOAD DATA: breast cancer vs. normal healthy samples
# ─────────────────────────────────────────────────────────────────────────────

# 1a) Read in log2(TPM + 1) for breast cancer samples
#     - assume file has header: first column = gene_id, then one column per cancer sample
log2tpm_cancer <- read.table(
  "./2_recurrence.classifier/output/2_QC/pnas_tpm_96_nodup_normalized_filtered.csv",
  sep = ",", header = TRUE, row.names = 1
) %>%
  as.matrix()   

# 1b) Read in log2(TPM + 1) for normal healthy samples
#     - same format: first column = gene_id, then one column per normal sample
log2tpm_normal <- read.table(
  "./1_bc.classifier/output/1_QC/pnas_normal_tpm_normalized_filtered.csv",
  sep = ",", header = TRUE, row.names = 1
) %>%
  as.matrix()   

# ─────────────────────────────────────────────────────────────────────────────
# 2) ALIGN GENE SETS BETWEEN CANCER AND NORMAL
# ─────────────────────────────────────────────────────────────────────────────

# Find genes present in both matrices
common_genes <- intersect(rownames(log2tpm_cancer), rownames(log2tpm_normal))

# Subset each matrix to the common genes only, in the same order
log2tpm_cancer <- log2tpm_cancer[common_genes, , drop = FALSE]
log2tpm_normal <- log2tpm_normal[common_genes, , drop = FALSE]

# ─────────────────────────────────────────────────────────────────────────────
# 3) COMBINE MATRICES AND DEFINE GROUP VECTOR
# ─────────────────────────────────────────────────────────────────────────────

# Bind columns: first all cancer samples, then all normal samples
expr_mat <- cbind(
  log2tpm_cancer,
  log2tpm_normal
)

group_vec <- c(
  rep("Cancer",  ncol(log2tpm_cancer)),
  rep("Normal",  ncol(log2tpm_normal))
)

# ─────────────────────────────────────────────────────────────────────────────
# 4) DE ANALYSIS LOOP: log₂ fold-change + t-test p-values
# ─────────────────────────────────────────────────────────────────────────────

gene_ids <- rownames(expr_mat)
n_genes  <- length(gene_ids)

# Pre-allocate vectors
log2FC  <- numeric(n_genes)   # will store mean(Cancer) – mean(Normal)
p_value <- numeric(n_genes)   # will store raw p-values from t.test

for (i in seq_len(n_genes)) {
  # Extract the i-th gene’s log2(TPM+1) across all samples
  gene_expr <- expr_mat[i, ]     # numeric vector of length = ncol(expr_mat)
  
  # Split into two groups by “Cancer” vs “Normal”
  expr_C <- gene_expr[group_vec == "Cancer"]
  expr_N <- gene_expr[group_vec == "Normal"]
  
  # 4a) Compute log2FC = mean(Cancer) – mean(Normal)
  log2FC[i] <- mean(expr_C, na.rm = TRUE) - mean(expr_N, na.rm = TRUE)
  
  # 4b) Perform two-sample t-test (unpaired, equal variance assumed)
  t_test_res <- try(t.test(expr_C, expr_N), silent = TRUE)
  if (inherits(t_test_res, "try-error")) {
    p_value[i] <- NA
  } else {
    p_value[i] <- t_test_res$p.value
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 5) ADJUST P-VALUES (BONFERRONI)
# ─────────────────────────────────────────────────────────────────────────────

p_adj <- p.adjust(p_value, method = "bonferroni")

# ─────────────────────────────────────────────────────────────────────────────
# 6) ASSEMBLE RESULTS INTO A DATA FRAME
# ─────────────────────────────────────────────────────────────────────────────

de_results_cn <- tibble(
  gene    = gene_ids,
  log2FC  = log2FC,
  p_value = p_value,
  p_adj   = p_adj
) %>%
  arrange(p_adj)   # sorted by smallest adjusted p-value first

# ─────────────────────────────────────────────────────────────────────────────
# 7) SAVE TO CSV
# ─────────────────────────────────────────────────────────────────────────────

write_csv(
  de_results_cn,
  "./1_bc.classifier/output/3_bcvsnormal/DE_CancervsNormal_log2TPM_ttest.csv"
)

de_results_cn<- read.csv("./1_bc.classifier/output/3_bcvsnormal/DE_CancervsNormal_log2TPM_ttest.csv")

ensembl_ids<- de_results_cn%>%.[,1]


library(biomaRt)

# Connect to Ensembl (human GRCh38)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# Query Ensembl for HGNC symbol (and biotype, if desired)
annotations <- getBM(
  filters    = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  values     = ensembl_ids,
  mart       = ensembl
)

de_results_cn <-merge(x=de_results_cn,y=annotations,by.x="gene",by.y="ensembl_gene_id",all.x = T)
library(EnhancedVolcano)
EnhancedVolcano(de_results_cn, 
                #NA,
                de_results_cn$hgnc_symbol,
                x ="log2FC", 
                y ="p_adj",
                FCcutoff=1,
                pCutoff = 0.05)+ coord_flip()


ggsave(paste0("./1_bc.classifier/output/3_bcvsnormal/volcano.png"),width=10, height=6, limitsize = FALSE)
de_results_cn%>%filter(log2FC>1 & p_adj<0.05)%>%dim()

