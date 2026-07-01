#collapse 18S ASV table to OTU table
library(readr)
library(dplyr)

# Input files
asv_table_file <- "DADA2_counts_as_matrix.tsv"
mapping_file <- "ASV_to_OTU99_map.tsv"

# Output file
otu_table_file <- "otu99_table_18S.tsv"

# Read ASV table and mapping file
asv_table <- read_tsv(asv_table_file, show_col_types = FALSE)
mapping <- read_tsv(mapping_file, show_col_types = FALSE)

# Check required columns
stopifnot("ASV" %in% colnames(asv_table))
stopifnot(all(c("ASV", "OTU99") %in% colnames(mapping)))

# Join ASV table to OTU mapping, then collapse counts by OTU99
otu_table <- asv_table %>%
  left_join(mapping, by = "ASV") %>%
  mutate(
    OTU99 = if_else(is.na(OTU99), ASV, OTU99)
  ) %>%
  select(-ASV) %>%
  group_by(OTU99) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(OTU = OTU99)

# Write OTU table
write_tsv(otu_table, otu_table_file)

# Quick checks
cat("Number of ASVs in original table:", nrow(asv_table), "\n")
cat("Number of OTUs in collapsed table:", nrow(otu_table), "\n")
cat("Output written to:", otu_table_file, "\n")
