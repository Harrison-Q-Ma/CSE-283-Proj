library(tidyverse)    
library(magrittr)     
library(tibble)     

# ── 1) Load TPM matrix ────────────────────────────
tpm <- read.delim("~/Desktop/BISB program/25SP/BENG 203:CSE 283/CSE-283-Proj/2_recurrence.classifier/output/2_QC/pnas_tpm_96_nodup_normalized_filtered.csv",
                  header = T, sep = ",", stringsAsFactors = FALSE)

# ── 2) Compute mean & variance per gene ─────────────────────────────────────
hvg_stats <- tpm %>%
  # select HVG only in those that have present expression greater than 0.1
  #filter(X %in% qc_filtered$GeneID) %>% 
  rowwise() %>%
  mutate(
    avg_expr = mean(c_across(-X)),  # mean TPM across all sample‐columns
    var_expr = var(c_across(-X))    # variance of TPM across all sample‐columns
  ) %>%
  ungroup() %>%
  # compute CV² = var / (mean²)
  mutate(cv2 = var_expr / (avg_expr^2)) %>%
  dplyr::select(X, avg_expr, var_expr, cv2)

# ── 3) Rank genes by dispersion (cv2) ────────────────────────────────────────
hvg_ranked <- hvg_stats %>%
  arrange(desc(cv2))  # descending → highest CV² first

# ── 4) Choose top N HVGs ─────────────────────────────────────────
top_n <- 2000
top_hvgs <- hvg_ranked %>%
  slice(1:top_n) %>%
  pull(X)

# ── 5) Save tables ────────────────────────────────────────
write_csv(hvg_ranked,"./2_recurrence.classifier/output/HVG_stats_all_genes_filtered.csv")
write_csv(hvg_ranked %>% slice(1:top_n),
          glue::glue("./2_recurrence.classifier/output/top_{top_n}_HVGs.csv"))
hvg_ranked$category<-case_when(hvg_ranked$cv2>=7.176915 ~ "Highly Variable Genes",
                               TRUE ~ "other genes")

library(ggplot2)

ggplot(hvg_ranked, aes(x = avg_expr, y = cv2)) +
  geom_point(alpha = 0.5, size = 1.5,aes(color = category)) +
  # Draw a horizontal line at log2(14.32987) to show the cutoff
  #geom_hline(yintercept = log2(14.32987), linetype = "dashed", color = "red") +
  # Force specific colors for the two categories
  #scale_color_manual(
  #  values = c("Other Genes" = "grey70",
  #             "Highly Variable Genes" = "black")
  #) +
  labs(
    x     = "Average expression", 
    y     = "dispersions of genes", 
    color = ""
  ) +
  #xlim(0,500)+
  theme_bw(base_size = 12) +
  theme(
    #legend.position = c(0.85, 0.85),
    #legend.background = element_rect(fill = "white", color = "grey80"),
    #panel.grid.major = element_line(color = "grey90", size = 0.3),
    #panel.grid.minor = element_blank()
  )
