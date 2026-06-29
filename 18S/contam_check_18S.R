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






# Remove 18S contaminant OTUs identified from DNA or RNA blanks
#
# Logic:
# 1. Use metadata column "JMF sample ID" to identify blank samples.
# 2. Match metadata JMF sample IDs to OTU table columns.
#    - Ignores one trailing letter in OTU table names, e.g.
#      JMF-2301-09-0061 matches JMF-2301-09-0061C
# 3. Flag any OTU with >0 reads in ANY DNA or RNA blank.
# 4. Remove those flagged OTUs from the entire OTU table.
# 5. Filter taxonomy to match remaining OTUs.
# 6. Print:
#    - total number of OTUs
#    - number of OTUs flagged
#    - percent of OTUs flagged
#    - total reads
#    - reads from flagged OTUs
#    - percent of reads from flagged OTUs

# -------------------------------
# Input files
# -------------------------------

taxonomy_file <- "combined_PR2_with_SILVA_fungi_no_mammals.tsv"
otu_table_file <- "otu99_table_no_mammals.tsv"
metadata_file <- "../jmf_esw_d15_with_all_metadata_and_isotopes.csv"

# -------------------------------
# Output files
# -------------------------------

filtered_otu_table_file <- "otu99_table_no_mammals_no_blank_contaminant_OTUs.tsv"
filtered_taxonomy_file <- "combined_PR2_with_SILVA_fungi_no_mammals_no_blank_contaminant_OTUs.tsv"

flagged_contaminant_otus_file <- "blank_flagged_contaminant_OTUs.tsv"
blank_samples_used_file <- "blank_samples_used_for_contaminant_detection.tsv"
contaminant_summary_file <- "blank_contaminant_OTU_summary.tsv"

# -------------------------------
# Helper: remove one trailing letter from sample IDs
#
# Example:
# JMF-2301-09-0061C -> JMF-2301-09-0061
# JMF-2301-09-0061  -> JMF-2301-09-0061
# -------------------------------

sample_core_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("[A-Za-z]$", "", x)
  return(x)
}

# -------------------------------
# Read files
# -------------------------------

taxonomy <- read_tsv(taxonomy_file, show_col_types = FALSE, na = c("", "NA"))
otu_table <- read_tsv(otu_table_file, show_col_types = FALSE)
metadata <- read_csv(metadata_file, show_col_types = FALSE, na = c("", "NA"))

# Clean column names
colnames(taxonomy) <- trimws(colnames(taxonomy))
colnames(otu_table) <- trimws(colnames(otu_table))
colnames(metadata) <- trimws(colnames(metadata))

# -------------------------------
# Standardize OTU ID column
# -------------------------------

stopifnot("OTU_ID" %in% colnames(taxonomy))

otu_id_col <- colnames(otu_table)[1]

otu_table <- otu_table %>%
  rename(OTU_ID = all_of(otu_id_col))

# -------------------------------
# Identify OTU table sample columns
# -------------------------------

count_cols <- setdiff(colnames(otu_table), "OTU_ID")

otu_table <- otu_table %>%
  mutate(
    across(
      all_of(count_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

# -------------------------------
# Check metadata columns
# -------------------------------

stopifnot("JMF sample ID" %in% colnames(metadata))

if (!"act" %in% colnames(metadata)) {
  metadata$act <- NA_character_
}

if (!"site" %in% colnames(metadata)) {
  metadata$site <- NA_character_
}

if (!"Sample description" %in% colnames(metadata)) {
  metadata$`Sample description` <- NA_character_
}

# -------------------------------
# Build OTU table sample lookup
# -------------------------------

otu_sample_lookup <- tibble(
  otu_sample_col = count_cols,
  sample_core = sample_core_id(count_cols)
)

# -------------------------------
# Prepare metadata and identify blanks
#
# A metadata row is treated as a blank if:
# - site == "Blank"
# OR
# - Sample description contains "blank"
# OR
# - JMF sample ID contains "blank"
# -------------------------------

metadata_for_matching <- metadata %>%
  mutate(
    jmf_sample_id = trimws(as.character(`JMF sample ID`)),
    sample_core = sample_core_id(jmf_sample_id),
    act = tolower(trimws(as.character(act))),
    site = trimws(as.character(site)),
    sample_description_clean = trimws(as.character(`Sample description`)),
    
    is_blank = site == "Blank" |
      grepl("blank", sample_description_clean, ignore.case = TRUE) |
      grepl("blank", jmf_sample_id, ignore.case = TRUE)
  )

# -------------------------------
# Match metadata samples to OTU table columns
# using JMF sample ID core
# -------------------------------

metadata_otu_matching <- metadata_for_matching %>%
  filter(
    !is.na(jmf_sample_id),
    jmf_sample_id != "",
    !is.na(sample_core),
    sample_core != ""
  ) %>%
  left_join(
    otu_sample_lookup,
    by = "sample_core"
  ) %>%
  distinct(
    jmf_sample_id,
    sample_core,
    otu_sample_col,
    act,
    site,
    sample_description_clean,
    is_blank
  )

# -------------------------------
# Get blank sample columns in OTU table
#
# Includes both DNA and RNA blanks.
# -------------------------------

blank_samples_used <- metadata_otu_matching %>%
  filter(
    is_blank,
    !is.na(otu_sample_col)
  ) %>%
  distinct(
    jmf_sample_id,
    sample_core,
    otu_sample_col,
    act,
    site,
    sample_description_clean
  ) %>%
  arrange(act, jmf_sample_id, otu_sample_col)

blank_sample_cols <- blank_samples_used %>%
  pull(otu_sample_col) %>%
  unique()

cat("Blank sample columns found in OTU table:", length(blank_sample_cols), "\n")
print(blank_samples_used, n = Inf)

# -------------------------------
# Stop if no blanks were found
# -------------------------------

if (length(blank_sample_cols) == 0) {
  stop("No blank samples were matched to OTU table columns. Check JMF sample IDs and OTU table column names.")
}

# -------------------------------
# Identify non-blank sample columns
# -------------------------------

nonblank_sample_cols <- setdiff(count_cols, blank_sample_cols)

# -------------------------------
# Flag contaminant OTUs
#
# An OTU is flagged if it has >0 reads in ANY DNA or RNA blank.
# -------------------------------

flagged_contaminant_otus <- otu_table %>%
  mutate(
    reads_in_blanks = rowSums(across(all_of(blank_sample_cols)), na.rm = TRUE),
    reads_in_nonblank_samples = if (length(nonblank_sample_cols) > 0) {
      rowSums(across(all_of(nonblank_sample_cols)), na.rm = TRUE)
    } else {
      0
    },
    total_reads_for_OTU = rowSums(across(all_of(count_cols)), na.rm = TRUE)
  ) %>%
  filter(reads_in_blanks > 0) %>%
  select(
    OTU_ID,
    reads_in_blanks,
    reads_in_nonblank_samples,
    total_reads_for_OTU
  )

# -------------------------------
# Calculate summary statistics
# -------------------------------

total_otus <- nrow(otu_table)
flagged_otus <- nrow(flagged_contaminant_otus)
percent_otus_flagged <- 100 * flagged_otus / total_otus

total_reads_all_samples <- otu_table %>%
  summarise(
    total = sum(across(all_of(count_cols)), na.rm = TRUE)
  ) %>%
  pull(total)

total_reads_nonblank_samples <- otu_table %>%
  summarise(
    total = if (length(nonblank_sample_cols) > 0) {
      sum(across(all_of(nonblank_sample_cols)), na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(total)

flagged_reads_all_samples <- otu_table %>%
  semi_join(flagged_contaminant_otus, by = "OTU_ID") %>%
  summarise(
    total = sum(across(all_of(count_cols)), na.rm = TRUE)
  ) %>%
  pull(total)

flagged_reads_nonblank_samples <- otu_table %>%
  semi_join(flagged_contaminant_otus, by = "OTU_ID") %>%
  summarise(
    total = if (length(nonblank_sample_cols) > 0) {
      sum(across(all_of(nonblank_sample_cols)), na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(total)

percent_reads_flagged_all_samples <- 100 * flagged_reads_all_samples / total_reads_all_samples

percent_reads_flagged_nonblank_samples <- 100 * flagged_reads_nonblank_samples / total_reads_nonblank_samples

# -------------------------------
# Print summary
# -------------------------------

cat("\n===== Blank contaminant OTU summary =====\n")
cat("Total OTUs in OTU table:", total_otus, "\n")
cat("OTUs flagged from DNA or RNA blanks:", flagged_otus, "\n")
cat("Percent OTUs flagged:", round(percent_otus_flagged, 4), "%\n")

cat("\nTotal reads across all samples:", total_reads_all_samples, "\n")
cat("Reads belonging to flagged OTUs across all samples:", flagged_reads_all_samples, "\n")
cat("Percent of all reads belonging to flagged OTUs:", round(percent_reads_flagged_all_samples, 4), "%\n")

cat("\nTotal reads across non-blank samples:", total_reads_nonblank_samples, "\n")
cat("Reads belonging to flagged OTUs in non-blank samples:", flagged_reads_nonblank_samples, "\n")
cat("Percent of non-blank reads belonging to flagged OTUs:", round(percent_reads_flagged_nonblank_samples, 4), "%\n")

# -------------------------------
# Remove flagged OTUs from entire OTU table
# -------------------------------

otu_table_no_contaminants <- otu_table %>%
  anti_join(
    flagged_contaminant_otus %>% select(OTU_ID),
    by = "OTU_ID"
  )

# Optional: remove blank sample columns from final OTU table
# Comment this out if you want to keep blank columns.
otu_table_no_contaminants <- otu_table_no_contaminants %>%
  select(-any_of(blank_sample_cols))

# -------------------------------
# Filter taxonomy to remaining OTUs
# -------------------------------

taxonomy_no_contaminants <- taxonomy %>%
  semi_join(
    otu_table_no_contaminants %>% select(OTU_ID),
    by = "OTU_ID"
  )

# -------------------------------
# Add taxonomy to flagged contaminant OTU table
# -------------------------------

flagged_contaminant_otus_with_taxonomy <- flagged_contaminant_otus %>%
  left_join(taxonomy, by = "OTU_ID") %>%
  arrange(desc(reads_in_blanks), desc(reads_in_nonblank_samples))

# -------------------------------
# Save summary table
# -------------------------------

contaminant_summary <- tibble(
  total_otus = total_otus,
  flagged_otus = flagged_otus,
  percent_otus_flagged = percent_otus_flagged,
  total_reads_all_samples = total_reads_all_samples,
  flagged_reads_all_samples = flagged_reads_all_samples,
  percent_reads_flagged_all_samples = percent_reads_flagged_all_samples,
  total_reads_nonblank_samples = total_reads_nonblank_samples,
  flagged_reads_nonblank_samples = flagged_reads_nonblank_samples,
  percent_reads_flagged_nonblank_samples = percent_reads_flagged_nonblank_samples,
  blank_sample_columns_used = length(blank_sample_cols)
)

# -------------------------------
# Write outputs
# -------------------------------

write_tsv(otu_table_no_contaminants, filtered_otu_table_file)
write_tsv(taxonomy_no_contaminants, filtered_taxonomy_file, na = "NA")
write_tsv(flagged_contaminant_otus_with_taxonomy, flagged_contaminant_otus_file, na = "NA")
write_tsv(blank_samples_used, blank_samples_used_file, na = "NA")
write_tsv(contaminant_summary, contaminant_summary_file, na = "NA")

# -------------------------------
# Final checks
# -------------------------------

cat("\n===== Final output checks =====\n")
cat("Original OTU table rows:", nrow(otu_table), "\n")
cat("Filtered OTU table rows:", nrow(otu_table_no_contaminants), "\n")
cat("Original taxonomy rows:", nrow(taxonomy), "\n")
cat("Filtered taxonomy rows:", nrow(taxonomy_no_contaminants), "\n")

cat("\nFiltered OTU table written to:", filtered_otu_table_file, "\n")
cat("Filtered taxonomy written to:", filtered_taxonomy_file, "\n")
cat("Flagged contaminant OTUs written to:", flagged_contaminant_otus_file, "\n")
cat("Blank samples used written to:", blank_samples_used_file, "\n")
cat("Contaminant summary written to:", contaminant_summary_file, "\n")
