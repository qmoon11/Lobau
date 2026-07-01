
# ============================================================
# 16S/SSU ASV blank-contaminant filtering
#
# Inputs:
# - ASV count table with blanks still included:
#   JMF-2301-09__all__rRNA_SSU_515_806__JMFR_MSG4_KP48D/DADA2_counts_as_matrix.tsv
#
# - ASV taxonomy table:
#   DADA2_ASVs.rRNA_SSU.SILVA_reference.DADA2_classified.tsv
#
# Logic:
# 1. Read ASV table and taxonomy table.
# 2. Use metadata to identify DNA and RNA blank sample columns.
# 3. Flag ASVs present in DNA blanks and/or RNA blanks.
# 4. Remove those ASVs from the full ASV table.
# 5. Remove blank sample columns from final ASV table.
# 6. Filter taxonomy to match remaining ASVs.
# 7. Report ASV/read statistics overall and for non-blank samples.
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# -------------------------------
# Input files
# -------------------------------

asv_table_file <- "JMF-2301-09__all__rRNA_SSU_515_806__JMFR_MSG4_KP48D/DADA2_counts_as_matrix.tsv"
taxonomy_file <- "JMF-2301-09__all__rRNA_SSU_515_806__JMFR_MSG4_KP48D/DADA2_ASVs.rRNA_SSU.SILVA_reference.DADA2_classified.tsv"
metadata_file <- "jmf_esw_d15_with_all_metadata_and_isotopes.csv"

# -------------------------------
# Output files
# -------------------------------

filtered_asv_table_file <- "SSU_515_806_ASV_table_no_blank_contaminants.tsv"
filtered_taxonomy_file <- "SSU_515_806_ASV_taxonomy_no_blank_contaminants.tsv"

blank_samples_used_file <- "SSU_515_806_blank_samples_used.tsv"
flagged_asvs_file <- "SSU_515_806_blank_flagged_ASVs.tsv"
blank_summary_file <- "SSU_515_806_blank_contaminant_summary.tsv"
reads_removed_per_sample_file <- "SSU_515_806_reads_removed_per_sample.tsv"

# -------------------------------
# Helper: clean sample IDs
#
# Handles:
# JMF-2301-09-0152B-16S   -> JMF-2301-09-0152
# JMF-2301-09-0152B-16S.F -> JMF-2301-09-0152
# JMF-2301-09-0152B       -> JMF-2301-09-0152
# -------------------------------

sample_core_id <- function(x) {
  x <- trimws(as.character(x))
  
  # Remove FASTQ-like suffixes if present
  x <- sub("\\.[12](\\.filtered)?\\.fastq\\.gz$", "", x, ignore.case = TRUE)
  
  # Remove marker/orientation suffixes
  x <- sub("[-.](ITS|18S|16S|SSU|rRNA_SSU)(\\.[FR])?$", "", x, ignore.case = TRUE)
  
  # Remove one trailing letter from JMF IDs only
  x <- ifelse(
    grepl("^JMF-", x),
    sub("[A-Za-z]$", "", x),
    x
  )
  
  return(x)
}

# -------------------------------
# Read files
# -------------------------------

asv_table <- read_tsv(asv_table_file, show_col_types = FALSE)
taxonomy <- read_tsv(taxonomy_file, show_col_types = FALSE, na = c("", "NA"))
metadata <- read_csv(metadata_file, show_col_types = FALSE, na = c("", "NA"))

# Clean column names
colnames(asv_table) <- trimws(colnames(asv_table))
colnames(taxonomy) <- trimws(colnames(taxonomy))
colnames(metadata) <- trimws(colnames(metadata))

colnames(asv_table)[1] <- sub("^\ufeff", "", colnames(asv_table)[1])
colnames(taxonomy)[1] <- sub("^\ufeff", "", colnames(taxonomy)[1])
colnames(metadata)[1] <- sub("^\ufeff", "", colnames(metadata)[1])

# -------------------------------
# Standardize ASV table
# -------------------------------

# First column should be ASV ID
colnames(asv_table)[1] <- "ASV"

asv_table <- asv_table %>%
  mutate(ASV = as.character(ASV))

asv_count_cols <- setdiff(colnames(asv_table), "ASV")

asv_table <- asv_table %>%
  mutate(
    across(
      all_of(asv_count_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

# -------------------------------
# Standardize taxonomy table
# -------------------------------

# Try to find ASV ID column in taxonomy
tax_asv_col <- case_when(
  "ASV" %in% colnames(taxonomy) ~ "ASV",
  "qseqid" %in% colnames(taxonomy) ~ "qseqid",
  "Sequence_ID" %in% colnames(taxonomy) ~ "Sequence_ID",
  "Feature.ID" %in% colnames(taxonomy) ~ "Feature.ID",
  "FeatureID" %in% colnames(taxonomy) ~ "FeatureID",
  TRUE ~ NA_character_
)

if (is.na(tax_asv_col)) {
  stop(
    "Could not find ASV ID column in taxonomy table. Columns found: ",
    paste(colnames(taxonomy), collapse = ", ")
  )
}

colnames(taxonomy)[colnames(taxonomy) == tax_asv_col] <- "ASV"

taxonomy <- taxonomy %>%
  mutate(ASV = as.character(ASV))

# -------------------------------
# Prepare metadata for blank detection
# -------------------------------

if (!"JMF sample ID" %in% colnames(metadata)) {
  stop("Metadata file must contain column: JMF sample ID")
}

if (!"act" %in% colnames(metadata)) {
  metadata$act <- NA_character_
}

if (!"site" %in% colnames(metadata)) {
  metadata$site <- NA_character_
}

if (!"Sample description" %in% colnames(metadata)) {
  metadata$`Sample description` <- NA_character_
}

metadata_for_matching <- metadata %>%
  mutate(
    jmf_sample_id = trimws(as.character(`JMF sample ID`)),
    sample_core = sample_core_id(jmf_sample_id),
    act = tolower(trimws(as.character(act))),
    site = trimws(as.character(site)),
    sample_description_clean = trimws(as.character(`Sample description`)),
    
    row_contains_blank = if_any(
      everything(),
      ~ !is.na(.x) & grepl("blank", as.character(.x), ignore.case = TRUE)
    ),
    
    is_blank = site == "Blank" |
      grepl("blank", sample_description_clean, ignore.case = TRUE) |
      grepl("blank", jmf_sample_id, ignore.case = TRUE) |
      row_contains_blank
  )

# -------------------------------
# Match metadata samples to ASV table columns
# -------------------------------

asv_sample_lookup <- tibble(
  sample_col = as.character(asv_count_cols),
  sample_core = as.character(sample_core_id(asv_count_cols))
)

metadata_asv_matching <- metadata_for_matching %>%
  filter(
    !is.na(jmf_sample_id),
    jmf_sample_id != "",
    !is.na(sample_core),
    sample_core != ""
  ) %>%
  left_join(
    asv_sample_lookup,
    by = "sample_core"
  ) %>%
  distinct(
    jmf_sample_id,
    sample_core,
    sample_col,
    act,
    site,
    sample_description_clean,
    is_blank
  )

# -------------------------------
# Identify blank sample columns
# -------------------------------

blank_samples_from_metadata <- metadata_asv_matching %>%
  filter(
    is_blank,
    !is.na(sample_col)
  ) %>%
  transmute(
    jmf_sample_id = as.character(jmf_sample_id),
    sample_core = as.character(sample_core),
    sample_col = as.character(sample_col),
    act = as.character(act),
    site = as.character(site),
    sample_description_clean = as.character(sample_description_clean)
  ) %>%
  distinct()

# Also catch columns that directly contain "blank"
blank_cols_direct <- asv_count_cols[grepl("blank", asv_count_cols, ignore.case = TRUE)]

blank_cols_by_name <- tibble(
  jmf_sample_id = rep(NA_character_, length(blank_cols_direct)),
  sample_core = as.character(sample_core_id(blank_cols_direct)),
  sample_col = as.character(blank_cols_direct),
  act = case_when(
    grepl("rna", blank_cols_direct, ignore.case = TRUE) ~ "rna",
    grepl("dna", blank_cols_direct, ignore.case = TRUE) ~ "dna",
    TRUE ~ NA_character_
  ),
  site = rep("Blank", length(blank_cols_direct)),
  sample_description_clean = rep(
    "Detected as blank from ASV table column name",
    length(blank_cols_direct)
  )
)

blank_samples_used <- bind_rows(
  blank_samples_from_metadata,
  blank_cols_by_name
) %>%
  filter(!is.na(sample_col)) %>%
  distinct(sample_col, .keep_all = TRUE) %>%
  arrange(act, sample_col)

blank_sample_cols <- blank_samples_used %>%
  pull(sample_col) %>%
  unique()

cat("\nBlank sample columns found in ASV table:", length(blank_sample_cols), "\n")
print(blank_samples_used, n = Inf)

if (length(blank_sample_cols) == 0) {
  cat("\nMetadata rows identified as blanks:\n")
  print(
    metadata_for_matching %>%
      filter(is_blank) %>%
      select(jmf_sample_id, sample_core, act, site, sample_description_clean),
    n = Inf
  )
  
  cat("\nFirst 30 ASV table sample columns and cleaned cores:\n")
  print(
    asv_sample_lookup %>%
      head(30),
    n = Inf
  )
  
  stop("No blank sample columns were found in the ASV table.")
}

write_tsv(blank_samples_used, blank_samples_used_file, na = "NA")

# -------------------------------
# Separate DNA and RNA blank columns
# -------------------------------

dna_blank_cols <- blank_samples_used %>%
  filter(act == "dna") %>%
  pull(sample_col) %>%
  unique()

rna_blank_cols <- blank_samples_used %>%
  filter(act == "rna") %>%
  pull(sample_col) %>%
  unique()

cat("\nDNA blank columns:", length(dna_blank_cols), "\n")
print(dna_blank_cols)

cat("\nRNA blank columns:", length(rna_blank_cols), "\n")
print(rna_blank_cols)

# -------------------------------
# Flag ASVs found in DNA or RNA blanks
# -------------------------------

nonblank_sample_cols <- setdiff(asv_count_cols, blank_sample_cols)

flagged_asvs <- asv_table %>%
  mutate(
    reads_in_DNA_blanks = if (length(dna_blank_cols) > 0) {
      rowSums(across(all_of(dna_blank_cols)), na.rm = TRUE)
    } else {
      0
    },
    reads_in_RNA_blanks = if (length(rna_blank_cols) > 0) {
      rowSums(across(all_of(rna_blank_cols)), na.rm = TRUE)
    } else {
      0
    },
    reads_in_any_blank = reads_in_DNA_blanks + reads_in_RNA_blanks,
    reads_in_nonblank_samples = if (length(nonblank_sample_cols) > 0) {
      rowSums(across(all_of(nonblank_sample_cols)), na.rm = TRUE)
    } else {
      0
    },
    total_reads_for_ASV_all_samples = rowSums(across(all_of(asv_count_cols)), na.rm = TRUE),
    contaminant_type = case_when(
      reads_in_DNA_blanks > 0 & reads_in_RNA_blanks > 0 ~ "DNA_and_RNA_blank",
      reads_in_DNA_blanks > 0 & reads_in_RNA_blanks == 0 ~ "DNA_blank_only",
      reads_in_DNA_blanks == 0 & reads_in_RNA_blanks > 0 ~ "RNA_blank_only",
      TRUE ~ "not_blank_detected"
    )
  ) %>%
  filter(reads_in_any_blank > 0) %>%
  select(
    ASV,
    contaminant_type,
    reads_in_DNA_blanks,
    reads_in_RNA_blanks,
    reads_in_any_blank,
    reads_in_nonblank_samples,
    total_reads_for_ASV_all_samples
  ) %>%
  arrange(desc(reads_in_any_blank), desc(reads_in_nonblank_samples))

flagged_asv_ids <- flagged_asvs$ASV

# -------------------------------
# Summary before removal
# -------------------------------

total_asvs <- nrow(asv_table)
flagged_asv_count <- length(flagged_asv_ids)
percent_asvs_flagged <- ifelse(
  total_asvs > 0,
  100 * flagged_asv_count / total_asvs,
  NA_real_
)

total_reads_all_samples <- asv_table %>%
  summarise(total = sum(across(all_of(asv_count_cols)), na.rm = TRUE)) %>%
  pull(total)

total_reads_nonblank_samples <- asv_table %>%
  summarise(
    total = if (length(nonblank_sample_cols) > 0) {
      sum(across(all_of(nonblank_sample_cols)), na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(total)

flagged_reads_all_samples <- asv_table %>%
  semi_join(flagged_asvs %>% select(ASV), by = "ASV") %>%
  summarise(total = sum(across(all_of(asv_count_cols)), na.rm = TRUE)) %>%
  pull(total)

flagged_reads_nonblank_samples <- asv_table %>%
  semi_join(flagged_asvs %>% select(ASV), by = "ASV") %>%
  summarise(
    total = if (length(nonblank_sample_cols) > 0) {
      sum(across(all_of(nonblank_sample_cols)), na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(total)

percent_flagged_reads_all_samples <- ifelse(
  total_reads_all_samples > 0,
  100 * flagged_reads_all_samples / total_reads_all_samples,
  NA_real_
)

percent_flagged_reads_nonblank_samples <- ifelse(
  total_reads_nonblank_samples > 0,
  100 * flagged_reads_nonblank_samples / total_reads_nonblank_samples,
  NA_real_
)

# -------------------------------
# Print summary
# -------------------------------

cat("\n===== SSU 515/806 blank-contaminant ASV summary =====\n")
cat("Total ASVs:", total_asvs, "\n")
cat("ASVs flagged from DNA or RNA blanks:", flagged_asv_count, "\n")
cat("Percent ASVs flagged:", round(percent_asvs_flagged, 4), "%\n")

cat("\nTotal reads across all samples:", total_reads_all_samples, "\n")
cat("Reads belonging to flagged ASVs across all samples:", flagged_reads_all_samples, "\n")
cat("Percent of all reads belonging to flagged ASVs:", round(percent_flagged_reads_all_samples, 4), "%\n")

cat("\nTotal reads across non-blank samples:", total_reads_nonblank_samples, "\n")
cat("Reads belonging to flagged ASVs in non-blank samples:", flagged_reads_nonblank_samples, "\n")
cat("Percent of non-blank reads belonging to flagged ASVs:", round(percent_flagged_reads_nonblank_samples, 4), "%\n")

cat("\nFlagged ASVs by contaminant type:\n")
print(flagged_asvs %>% count(contaminant_type), n = Inf)

# -------------------------------
# Remove flagged ASVs and blank sample columns
# -------------------------------

asv_table_clean <- asv_table %>%
  anti_join(flagged_asvs %>% select(ASV), by = "ASV") %>%
  select(-any_of(blank_sample_cols))

clean_count_cols <- setdiff(colnames(asv_table_clean), "ASV")

# Drop ASVs with zero reads after blank columns and contaminant ASVs are removed
asv_table_clean <- asv_table_clean %>%
  mutate(
    total_reads_after_filtering = rowSums(across(all_of(clean_count_cols)), na.rm = TRUE)
  ) %>%
  filter(total_reads_after_filtering > 0) %>%
  select(-total_reads_after_filtering)

# -------------------------------
# Filter taxonomy to remaining ASVs
# -------------------------------

taxonomy_clean <- taxonomy %>%
  semi_join(asv_table_clean %>% select(ASV), by = "ASV")

# Add taxonomy to flagged ASV file
flagged_asvs_with_taxonomy <- flagged_asvs %>%
  left_join(taxonomy, by = "ASV")

# -------------------------------
# Per-sample read removal summary
# -------------------------------

total_reads_per_sample <- asv_table %>%
  summarise(
    across(
      all_of(asv_count_cols),
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "total_reads_original"
  )

removed_reads_per_sample <- flagged_asvs %>%
  select(ASV) %>%
  left_join(asv_table, by = "ASV") %>%
  summarise(
    across(
      all_of(asv_count_cols),
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "removed_reads_from_flagged_ASVs"
  )

final_reads_per_sample <- asv_table_clean %>%
  summarise(
    across(
      all_of(clean_count_cols),
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "final_reads"
  )

reads_removed_per_sample <- total_reads_per_sample %>%
  left_join(removed_reads_per_sample, by = "sample") %>%
  left_join(final_reads_per_sample, by = "sample") %>%
  mutate(
    removed_reads_from_flagged_ASVs = if_else(
      is.na(removed_reads_from_flagged_ASVs),
      0,
      removed_reads_from_flagged_ASVs
    ),
    final_reads = if_else(is.na(final_reads), 0, final_reads),
    was_blank_sample = sample %in% blank_sample_cols,
    percent_reads_removed_from_flagged_ASVs = if_else(
      total_reads_original > 0,
      100 * removed_reads_from_flagged_ASVs / total_reads_original,
      NA_real_
    )
  ) %>%
  arrange(desc(percent_reads_removed_from_flagged_ASVs))

cat("\nReads removed per sample:\n")
print(reads_removed_per_sample, n = Inf)

# -------------------------------
# Final summary after removal
# -------------------------------

total_reads_final <- asv_table_clean %>%
  summarise(total = sum(across(all_of(clean_count_cols)), na.rm = TRUE)) %>%
  pull(total)

summary_df <- tibble(
  total_asvs_original = total_asvs,
  asvs_flagged_from_blanks = flagged_asv_count,
  percent_asvs_flagged = percent_asvs_flagged,
  asvs_remaining_after_filtering = nrow(asv_table_clean),
  taxonomy_rows_original = nrow(taxonomy),
  taxonomy_rows_after_filtering = nrow(taxonomy_clean),
  blank_sample_columns_used = length(blank_sample_cols),
  dna_blank_columns_used = length(dna_blank_cols),
  rna_blank_columns_used = length(rna_blank_cols),
  total_reads_all_samples_original = total_reads_all_samples,
  flagged_reads_all_samples = flagged_reads_all_samples,
  percent_flagged_reads_all_samples = percent_flagged_reads_all_samples,
  total_reads_nonblank_samples_original = total_reads_nonblank_samples,
  flagged_reads_nonblank_samples = flagged_reads_nonblank_samples,
  percent_flagged_reads_nonblank_samples = percent_flagged_reads_nonblank_samples,
  total_reads_final_after_removing_flagged_ASVs_and_blank_columns = total_reads_final
)

# -------------------------------
# Write outputs
# -------------------------------

write_tsv(asv_table_clean, filtered_asv_table_file)
write_tsv(taxonomy_clean, filtered_taxonomy_file, na = "NA")
write_tsv(blank_samples_used, blank_samples_used_file, na = "NA")
write_tsv(flagged_asvs_with_taxonomy, flagged_asvs_file, na = "NA")
write_tsv(summary_df, blank_summary_file, na = "NA")
write_tsv(reads_removed_per_sample, reads_removed_per_sample_file, na = "NA")

# -------------------------------
# Final messages
# -------------------------------

cat("\n===== Final output checks =====\n")
cat("Original ASV table rows:", nrow(asv_table), "\n")
cat("Filtered ASV table rows:", nrow(asv_table_clean), "\n")
cat("Original taxonomy rows:", nrow(taxonomy), "\n")
cat("Filtered taxonomy rows:", nrow(taxonomy_clean), "\n")
cat("Final reads after filtering:", total_reads_final, "\n")

cat("\nFiltered ASV table written to:", filtered_asv_table_file, "\n")
cat("Filtered taxonomy written to:", filtered_taxonomy_file, "\n")
cat("Blank samples used written to:", blank_samples_used_file, "\n")
cat("Flagged ASVs written to:", flagged_asvs_file, "\n")
cat("Summary written to:", blank_summary_file, "\n")
cat("Reads removed per sample written to:", reads_removed_per_sample_file, "\n")
