#remove 18S contaminated samples from the dataset (just mammals)

library(readr)
library(dplyr)

# -------------------------------
# Input files
# -------------------------------

taxonomy_file <- "combined_PR2_with_SILVA_fungi.tsv"
otu_table_file <- "otu99_table_18S.tsv"

# -------------------------------
# Output files
# -------------------------------

filtered_taxonomy_file <- "combined_PR2_with_SILVA_fungi_no_mammals.tsv"
filtered_otu_table_file <- "otu99_table_no_mammals.tsv"
mammal_otus_file <- "mammal_OTUs_removed.tsv"

# -------------------------------
# Read files
# -------------------------------

taxonomy <- read_tsv(taxonomy_file, show_col_types = FALSE, na = c("", "NA"))
otu_table <- read_tsv(otu_table_file, show_col_types = FALSE)

# -------------------------------
# Standardize OTU ID column names
# -------------------------------

# Taxonomy table should have OTU_ID
stopifnot("OTU_ID" %in% colnames(taxonomy))

# OTU table may have first column named OTU, OTU99, or OTU_ID
otu_id_col <- colnames(otu_table)[1]

otu_table <- otu_table %>%
  rename(OTU_ID = all_of(otu_id_col))

# -------------------------------
# Identify mammal OTUs from combined taxonomy
# Looks across taxonomy columns for Mammalia
# -------------------------------

taxonomy_cols <- c(
  "Domain",
  "Supergroup",
  "Division",
  "Subdivision/Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus"
)

mammal_otus <- taxonomy %>%
  filter(
    if_any(
      all_of(taxonomy_cols),
      ~ grepl("Mammalia", .x, ignore.case = TRUE)
    )
  ) %>%
  distinct(OTU_ID)

# -------------------------------
# Calculate percent of total reads that are mammal reads
# -------------------------------

count_cols <- setdiff(colnames(otu_table), "OTU_ID")

total_reads <- otu_table %>%
  summarise(total = sum(across(all_of(count_cols)), na.rm = TRUE)) %>%
  pull(total)

mammal_reads <- otu_table %>%
  semi_join(mammal_otus, by = "OTU_ID") %>%
  summarise(total = sum(across(all_of(count_cols)), na.rm = TRUE)) %>%
  pull(total)

percent_mammal_reads <- 100 * mammal_reads / total_reads

cat("Number of mammal OTUs found:", nrow(mammal_otus), "\n")
cat("Total reads:", total_reads, "\n")
cat("Mammal reads:", mammal_reads, "\n")
cat("Percent mammal reads:", round(percent_mammal_reads, 4), "%\n")

# -------------------------------
# Remove mammal OTUs from both taxonomy and OTU table
# -------------------------------

taxonomy_no_mammals <- taxonomy %>%
  anti_join(mammal_otus, by = "OTU_ID")

otu_table_no_mammals <- otu_table %>%
  anti_join(mammal_otus, by = "OTU_ID")

# -------------------------------
# Write outputs
# -------------------------------

write_tsv(taxonomy_no_mammals, filtered_taxonomy_file, na = "NA")
write_tsv(otu_table_no_mammals, filtered_otu_table_file)
write_tsv(mammal_otus, mammal_otus_file)

# -------------------------------
# Final checks
# -------------------------------

cat("Original taxonomy rows:", nrow(taxonomy), "\n")
cat("Filtered taxonomy rows:", nrow(taxonomy_no_mammals), "\n")
cat("Original OTU table rows:", nrow(otu_table), "\n")
cat("Filtered OTU table rows:", nrow(otu_table_no_mammals), "\n")
cat("Removed mammal OTUs written to:", mammal_otus_file, "\n")
cat("Filtered taxonomy written to:", filtered_taxonomy_file, "\n")

# Identify count columns
count_cols <- setdiff(colnames(otu_table), "OTU_ID")

# Total reads across all OTUs
total_reads <- otu_table %>%
  summarise(total_reads = sum(across(all_of(count_cols)), na.rm = TRUE)) %>%
  pull(total_reads)

# Total mammal reads
mammal_reads <- otu_table %>%
  semi_join(mammal_otus, by = "OTU_ID") %>%
  summarise(total_mammal_reads = sum(across(all_of(count_cols)), na.rm = TRUE)) %>%
  pull(total_mammal_reads)

# Percent mammal reads
percent_mammal_reads <- 100 * mammal_reads / total_reads

cat("Total reads:", total_reads, "\n")
cat("Total mammal reads:", mammal_reads, "\n")
cat("Percent mammal reads:", round(percent_mammal_reads, 4), "%\n")