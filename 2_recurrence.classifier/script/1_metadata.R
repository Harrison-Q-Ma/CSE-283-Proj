library(dplyr)
library(tidyverse)
library(ggplot2)
metadata <- read.csv("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/optional/pnas_patient_info.csv")
metadata$cancer.subtype<-
  case_when(metadata$er_cat==1 & metadata$pr_cat==1 & metadata$her2_cat==0 ~ "Luminal_A",
            metadata$er_cat==1  & metadata$her2_cat==1 ~ "Luminal_B",
            metadata$er_cat==0 & metadata$pr_cat==0 & metadata$her2_cat==1 ~ "HER2_enriched",
            metadata$er_cat==0 & metadata$pr_cat==0 & metadata$her2_cat==0 ~ "Triple_negative")
donor_summary<-metadata[,c(1,21:31)]
donor_summary<-donor_summary%>%dplyr::distinct(donor_summary$poiseid,.keep_all = TRUE)%>%.[,c(-13)]



# Load required libraries
library(tidyverse)
library(lubridate)

# Load the metadata CSV
metadata <- read.csv("~/Desktop/BISB program/25SP/BENG 203:CSE 283/data/optional/pnas_patient_info.csv")

# Convert date columns
metadata$datespecimen <- dmy(metadata$datespecimen)
metadata$daterecurrence <- dmy(metadata$daterecurrence)
metadata$datechemostart <- dmy(metadata$datechemostart)
metadata$datechemoend <- dmy(metadata$datechemoend)

# Remove recurrence date if recurStatus == "N"
metadata$daterecurrence[metadata$recurStatus == "N"] <- NA

# Construct donor label: poiseid + recurStatus
metadata$donor_label <- paste0(metadata$poiseid, " (", metadata$recurStatus, ")")

# Order y-axis: R first, then N
metadata$recur_order <- ifelse(metadata$recurStatus == "R", 1, 2)
metadata <- metadata %>%
  arrange(recur_order, poiseid) %>%
  mutate(donor_label = factor(donor_label, levels = unique(donor_label)))

# Timeline data: one row per sample
timeline_data <- metadata %>%
  select(donor_label, poiseid, datespecimen, recurStatus) %>%
  distinct()

# Recurrence data: one per donor
recurrence_data <- metadata %>%
  select(donor_label, poiseid, daterecurrence) %>%
  distinct() %>%
  group_by(donor_label) %>%
  summarise(daterecurrence = first(na.omit(daterecurrence))) %>%
  filter(!is.na(daterecurrence))

# Chemotherapy period data
chemo_data <- metadata %>%
  select(donor_label, datechemostart, datechemoend) %>%
  distinct() %>%
  filter(!is.na(datechemostart) & !is.na(datechemoend))

# Final Plot
ggplot(timeline_data, aes(x = datespecimen, y = donor_label)) +
  # Chemotherapy duration bars
  geom_segment(data = chemo_data,
               aes(x = datechemostart, xend = datechemoend, y = donor_label, yend = donor_label),
               color = "black", size = 1.2) +
  
  # Timeline lines
  geom_line(aes(group = donor_label), color = "gray60") +
  
  # Sample timepoints
  geom_point(aes(color = recurStatus), size = 2) +
  
  # Recurrence events
  geom_point(data = recurrence_data,
             aes(x = daterecurrence, y = donor_label),
             shape = 4, color = "red", size = 3, stroke = 1.5) +
  
  # Color coding
  scale_color_manual(
    values = c("R" = "steelblue", "N" = "darkgreen"),
    labels = c("R" = "Recurrence", "N" = "No Recurrence"),
    name = "Sample Type"
  ) +
  
  # Labels and themes
  labs(
    x = "Date",
    y = "Donor",
    title = "Breast cancer sample timeline",
    subtitle = "Red cross = recurrence event, Black bar = chemo period"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

