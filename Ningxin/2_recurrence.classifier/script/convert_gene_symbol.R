ensembl_ids <- c(
  # HVG‐only (20)
  "ENSG00000210195", "ENSG00000231875", "ENSG00000225394", "ENSG00000260811",
  "ENSG00000222017", "ENSG00000232054", "ENSG00000254993", "ENSG00000255167",
  "ENSG00000158639", "ENSG00000249729", "ENSG00000256837", "ENSG00000186678",
  "ENSG00000243601", "ENSG00000257141", "ENSG00000215313", "ENSG00000244541",
  "ENSG00000256101", "ENSG00000235646", "ENSG00000267399", "ENSG00000250516",
  # HVG+DEG (20)
  "ENSG00000252316", "ENSG00000210077", "ENSG00000114739", "ENSG00000201861",
  "ENSG00000201098", "ENSG00000202354", "ENSG00000061455", "ENSG00000199839",
  "ENSG00000201413", "ENSG00000163359", "ENSG00000102287", "ENSG00000172992",
  "ENSG00000070808", "ENSG00000130032", "ENSG00000232035", "ENSG00000240869",
  "ENSG00000067066", "ENSG00000188536", "ENSG00000163736", "ENSG00000237649",
  # LR (20)
  "ENSG00000272240", "ENSG00000241223", "ENSG00000264824", "ENSG00000234898",
  "ENSG00000277466", "ENSG00000258646", "ENSG00000258595", "ENSG00000237632",
  "ENSG00000258427", "ENSG00000269877", "ENSG00000265520", "ENSG00000248416",
  "ENSG00000274770", "ENSG00000267317", "ENSG00000272085", "ENSG00000268896",
  "ENSG00000255336", "ENSG00000206746", "ENSG00000252749", "ENSG00000260687",
  # SVM (20)
  "ENSG00000272240", "ENSG00000241223", "ENSG00000229213", "ENSG00000259910",
  "ENSG00000258646", "ENSG00000234898", "ENSG00000264824", "ENSG00000277466",
  "ENSG00000267317", "ENSG00000248785", "ENSG00000228955", "ENSG00000281485",
  "ENSG00000212187", "ENSG00000268896", "ENSG00000226565", "ENSG00000264982",
  "ENSG00000274770", "ENSG00000234521", "ENSG00000258595", "ENSG00000255336"
)

ensembl_ids<- de_results_cn%>%filter(log2FC>1 & p_adj<0.05)%>%.[,1]
# Remove duplicates to get a unique set
ensembl_ids <- unique(ensembl_ids)
ensembl_ids

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

library(tidyverse)

table(annotations$gene_biotype) %>% 
  as.data.frame() %>% 
  rename(biotype = Var1, count = Freq) %>% 
  mutate(
    percent = count / sum(count) * 100,
    full_label = paste0(biotype, "\n", round(percent, 1), "%"),
    # only keep label if slice ≥ 3%
    plot_label = ifelse(percent >= 3, full_label, NA_character_)
  ) %>% 
  ggplot(aes(x = "", y = count, fill = biotype)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(label = plot_label),
            position      = position_stack(vjust = 0.5),
            size          = 3,
            check_overlap = F) +
  theme_void(base_size = 14) +
  labs(
    title = "Gene biotype distribution for 2000 HVGs found in cancer samples",
    fill  = "Gene Biotype"
  ) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "right"
  )
