
# ============================================================
# 18S trimming workflow
#
# This script:
# 1. Removes mammal OTUs based on taxonomy.
# 2. Uses DNA/RNA blanks to flag contaminant OTUs.
# 3. Removes any OTU detected in DNA or RNA blanks from the full table.
# 4. Removes blank sample columns from the final OTU table.
# 5. Removes anything assigned to Bacteria at the end.
# 6. Reports OTU and read totals at each major filtering step.
# ============================================================

library(readr)
library(dplyr)
library(tidyr)

# -------------------------------
# Input files
# -------------------------------

taxonomy_file <- "18S/combined_PR2_with_SILVA_fungi.tsv"
otu_table_file <- "18S/otu99_table_18S.tsv"
metadata_file <- "jmf_esw_d15_with_all_metadata_and_isotopes.csv"

# -------------------------------
# Output files
# -------------------------------

mammal_otus_file <- "18S/mammal_OTUs_removed.tsv"
taxonomy_no_mammals_file <- "18S/combined_PR2_with_SILVA_fungi_no_mammals.tsv"
otu_table_no_mammals_file <- "18S/otu99_table_no_mammals.tsv"

filtered_otu_table_file <- "18S/otu99_table_no_mammals_no_blank_contaminant_OTUs_no_bacteria.tsv"
filtered_taxonomy_file <- "18S/combined_PR2_with_SILVA_fungi_no_mammals_no_blank_contaminant_OTUs_no_bacteria.tsv"

flagged_contaminant_otus_file <- "18S/blank_flagged_contaminant_OTUs.tsv"
blank_samples_used_file <- "18S/blank_samples_used_for_contaminant_detection.tsv"
contaminant_summary_file <- "18S/blank_contaminant_OTU_summary.tsv"

bacteria_otus_file <- "18S/bacteria_OTUs_removed.tsv"
bacteria_summary_file <- "18S/bacteria_removal_summary.tsv"

# -------------------------------
# Helper: clean sample IDs
#
# Handles:
# JMF-2301-09-0061C       -> JMF-2301-09-0061
# JMF-2301-09-0061C-18S   -> JMF-2301-09-0061
# JMF-2301-09-0061C-18S.F -> JMF-2301-09-0061
# -------------------------------

sample_core_id <- function(x) {
  x <- trimws(as.character(x))
  
  # Remove FASTQ-like suffixes if present
  x <- sub("\\.[12](\\.filtered)?\\.fastq\\.gz$", "", x, ignore.case = TRUE)
  
  # Remove marker/orientation suffixes if present
  x <- sub("[-.](ITS|18S|16S)(\\.[FR])?$", "", x, ignore.case = TRUE)
  
  # Remove one trailing letter from JMF IDs only
  x <- ifelse(
    grepl("^JMF-", x),
    sub("[A-Za-z]$", "", x),
    x
  )
  
  return(x)
}

# ============================================================
# 1. Read files
# ============================================================

taxonomy <- read_tsv(taxonomy_file, show_col_types = FALSE, na = c("", "NA"))
otu_table <- read_tsv(otu_table_file, show_col_types = FALSE)

# Clean column names
colnames(taxonomy) <- trimws(colnames(taxonomy))
colnames(otu_table) <- trimws(colnames(otu_table))

colnames(taxonomy)[1] <- sub("^\ufeff", "", colnames(taxonomy)[1])
colnames(otu_table)[1] <- sub("^\ufeff", "", colnames(otu_table)[1])

# -------------------------------
# Standardize OTU ID column names
# -------------------------------

stopifnot("OTU_ID" %in% colnames(taxonomy))

# Avoid rename/all_of issues by directly setting first column name
colnames(otu_table)[1] <- "OTU_ID"

taxonomy <- taxonomy %>%
  mutate(OTU_ID = as.character(OTU_ID))

otu_table <- otu_table %>%
  mutate(OTU_ID = as.character(OTU_ID))

# Identify count columns
count_cols <- setdiff(colnames(otu_table), "OTU_ID")

otu_table <- otu_table %>%
  mutate(
    across(
      all_of(count_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

# ============================================================
# 2. Remove mammal OTUs
# ============================================================

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

taxonomy_cols <- intersect(taxonomy_cols, colnames(taxonomy))

mammal_otus <- taxonomy %>%
  filter(
    if_any(
      all_of(taxonomy_cols),
      ~ grepl("Mammalia", as.character(.x), ignore.case = TRUE)
    )
  ) %>%
  distinct(OTU_ID)

total_reads_before_mammal <- otu_table %>%
  summarise(total = sum(across(all_of(count_cols)), na.rm = TRUE)) %>%
  pull(total)

mammal_reads <- otu_table %>%
  semi_join(mammal_otus, by = "OTU_ID") %>%
  summarise(total = sum(across(all_of(count_cols)), na.rm = TRUE)) %>%
  pull(total)

percent_mammal_reads <- ifelse(
  total_reads_before_mammal > 0,
  100 * mammal_reads / total_reads_before_mammal,
  NA_real_
)

cat("\n===== Mammal removal summary =====\n")
cat("Original taxonomy rows:", nrow(taxonomy), "\n")
cat("Original OTU table rows:", nrow(otu_table), "\n")
cat("Number of mammal OTUs found:", nrow(mammal_otus), "\n")
cat("Total reads before mammal removal:", total_reads_before_mammal, "\n")
cat("Mammal reads:", mammal_reads, "\n")
cat("Percent mammal reads:", round(percent_mammal_reads, 4), "%\n")

taxonomy_no_mammals <- taxonomy %>%
  anti_join(mammal_otus, by = "OTU_ID")

otu_table_no_mammals <- otu_table %>%
  anti_join(mammal_otus, by = "OTU_ID")

write_tsv(mammal_otus, mammal_otus_file)
write_tsv(taxonomy_no_mammals, taxonomy_no_mammals_file, na = "NA")
write_tsv(otu_table_no_mammals, otu_table_no_mammals_file)

cat("Filtered taxonomy rows after mammal removal:", nrow(taxonomy_no_mammals), "\n")
cat("Filtered OTU rows after mammal removal:", nrow(otu_table_no_mammals), "\n")

# Replace working objects for next step
taxonomy <- taxonomy_no_mammals
otu_table <- otu_table_no_mammals

count_cols <- setdiff(colnames(otu_table), "OTU_ID")

# ============================================================
# 3. Remove blank-associated contaminant OTUs
# ============================================================

metadata <- read_csv(metadata_file, show_col_types = FALSE, na = c("", "NA"))

colnames(metadata) <- trimws(colnames(metadata))
colnames(metadata)[1] <- sub("^\ufeff", "", colnames(metadata)[1])

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

otu_sample_lookup <- tibble(
  otu_sample_col = as.character(count_cols),
  sample_core = as.character(sample_core_id(count_cols))
)

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

cat("\n===== Blank sample matching =====\n")
cat("Blank sample columns found in OTU table:", length(blank_sample_cols), "\n")
print(blank_samples_used, n = Inf)

if (length(blank_sample_cols) == 0) {
  cat("\nMetadata rows identified as blanks:\n")
  print(
    metadata_for_matching %>%
      filter(is_blank) %>%
      select(jmf_sample_id, sample_core, act, site, sample_description_clean),
    n = Inf
  )
  
  cat("\nFirst 30 OTU table sample columns and cleaned cores:\n")
  print(
    otu_sample_lookup %>%
      head(30),
    n = Inf
  )
  
  stop("No blank samples were matched to OTU table columns. Check JMF sample IDs and OTU table column names.")
}

write_tsv(blank_samples_used, blank_samples_used_file, na = "NA")

nonblank_sample_cols <- setdiff(count_cols, blank_sample_cols)

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

total_otus_before_blank <- nrow(otu_table)
flagged_otus <- nrow(flagged_contaminant_otus)
percent_otus_flagged <- ifelse(
  total_otus_before_blank > 0,
  100 * flagged_otus / total_otus_before_blank,
  NA_real_
)

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

percent_reads_flagged_all_samples <- ifelse(
  total_reads_all_samples > 0,
  100 * flagged_reads_all_samples / total_reads_all_samples,
  NA_real_
)

percent_reads_flagged_nonblank_samples <- ifelse(
  total_reads_nonblank_samples > 0,
  100 * flagged_reads_nonblank_samples / total_reads_nonblank_samples,
  NA_real_
)

cat("\n===== Blank contaminant OTU summary =====\n")
cat("Total OTUs in OTU table after mammal removal:", total_otus_before_blank, "\n")
cat("OTUs flagged from DNA or RNA blanks:", flagged_otus, "\n")
cat("Percent OTUs flagged:", round(percent_otus_flagged, 4), "%\n")

cat("\nTotal reads across all samples:", total_reads_all_samples, "\n")
cat("Reads belonging to flagged OTUs across all samples:", flagged_reads_all_samples, "\n")
cat("Percent of all reads belonging to flagged OTUs:", round(percent_reads_flagged_all_samples, 4), "%\n")

cat("\nTotal reads across non-blank samples:", total_reads_nonblank_samples, "\n")
cat("Reads belonging to flagged OTUs in non-blank samples:", flagged_reads_nonblank_samples, "\n")
cat("Percent of non-blank reads belonging to flagged OTUs:", round(percent_reads_flagged_nonblank_samples, 4), "%\n")

# Remove flagged OTUs from entire OTU table
otu_table_no_contaminants <- otu_table %>%
  anti_join(
    flagged_contaminant_otus %>% select(OTU_ID),
    by = "OTU_ID"
  )

# Remove blank sample columns from final table
otu_table_no_contaminants <- otu_table_no_contaminants %>%
  select(-any_of(blank_sample_cols))

# Filter taxonomy to remaining OTUs
taxonomy_no_contaminants <- taxonomy %>%
  semi_join(
    otu_table_no_contaminants %>% select(OTU_ID),
    by = "OTU_ID"
  )

flagged_contaminant_otus_with_taxonomy <- flagged_contaminant_otus %>%
  left_join(taxonomy, by = "OTU_ID") %>%
  arrange(desc(reads_in_blanks), desc(reads_in_nonblank_samples))

contaminant_summary <- tibble(
  total_otus_after_mammal_removal = total_otus_before_blank,
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

write_tsv(flagged_contaminant_otus_with_taxonomy, flagged_contaminant_otus_file, na = "NA")
write_tsv(contaminant_summary, contaminant_summary_file, na = "NA")

cat("\nOTU rows after blank contaminant removal:", nrow(otu_table_no_contaminants), "\n")
cat("Taxonomy rows after blank contaminant removal:", nrow(taxonomy_no_contaminants), "\n")

# ============================================================
# 4. Remove anything assigned to Bacteria at the end
# ============================================================

bacteria_taxonomy_cols <- c(
  "Domain",
  "Kingdom",
  "kingdom",
  "Supergroup",
  "Division",
  "Subdivision/Kingdom",
  "Phylum",
  "phylum",
  "Class",
  "class",
  "Order",
  "order",
  "Family",
  "family",
  "Genus",
  "genus"
)

bacteria_taxonomy_cols <- intersect(
  bacteria_taxonomy_cols,
  colnames(taxonomy_no_contaminants)
)

if (length(bacteria_taxonomy_cols) == 0) {
  stop("No taxonomy columns found for bacterial filtering.")
}

bacteria_otus <- taxonomy_no_contaminants %>%
  filter(
    if_any(
      all_of(bacteria_taxonomy_cols),
      ~ grepl("Bacteria", as.character(.x), ignore.case = TRUE)
    )
  ) %>%
  distinct(OTU_ID)

count_cols_after_blank <- setdiff(colnames(otu_table_no_contaminants), "OTU_ID")

otu_table_no_contaminants <- otu_table_no_contaminants %>%
  mutate(
    across(
      all_of(count_cols_after_blank),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

total_reads_before_bacteria <- otu_table_no_contaminants %>%
  summarise(
    total = sum(across(all_of(count_cols_after_blank)), na.rm = TRUE)
  ) %>%
  pull(total)

bacteria_reads <- otu_table_no_contaminants %>%
  semi_join(bacteria_otus, by = "OTU_ID") %>%
  summarise(
    total = sum(across(all_of(count_cols_after_blank)), na.rm = TRUE)
  ) %>%
  pull(total)

percent_bacteria_reads <- ifelse(
  total_reads_before_bacteria > 0,
  100 * bacteria_reads / total_reads_before_bacteria,
  NA_real_
)

taxonomy_no_bacteria <- taxonomy_no_contaminants %>%
  anti_join(bacteria_otus, by = "OTU_ID")

otu_table_no_bacteria <- otu_table_no_contaminants %>%
  anti_join(bacteria_otus, by = "OTU_ID")

total_reads_after_bacteria <- otu_table_no_bacteria %>%
  summarise(
    total = sum(across(all_of(count_cols_after_blank)), na.rm = TRUE)
  ) %>%
  pull(total)

bacteria_otus_with_taxonomy_and_reads <- taxonomy_no_contaminants %>%
  semi_join(bacteria_otus, by = "OTU_ID") %>%
  left_join(
    otu_table_no_contaminants %>%
      semi_join(bacteria_otus, by = "OTU_ID") %>%
      mutate(
        total_reads_for_OTU = rowSums(across(all_of(count_cols_after_blank)), na.rm = TRUE)
      ) %>%
      select(OTU_ID, total_reads_for_OTU),
    by = "OTU_ID"
  ) %>%
  arrange(desc(total_reads_for_OTU))

bacteria_summary <- tibble(
  total_otus_before_bacteria_removal = nrow(otu_table_no_contaminants),
  bacteria_otus_removed = nrow(bacteria_otus),
  percent_otus_removed_as_bacteria = ifelse(
    nrow(otu_table_no_contaminants) > 0,
    100 * nrow(bacteria_otus) / nrow(otu_table_no_contaminants),
    NA_real_
  ),
  total_reads_before_bacteria_removal = total_reads_before_bacteria,
  bacteria_reads_removed = bacteria_reads,
  percent_reads_removed_as_bacteria = percent_bacteria_reads,
  total_reads_after_bacteria_removal = total_reads_after_bacteria,
  final_otus_after_bacteria_removal = nrow(otu_table_no_bacteria)
)

cat("\n===== Bacteria removal summary =====\n")
cat("OTUs before bacterial removal:", nrow(otu_table_no_contaminants), "\n")
cat("Bacterial OTUs removed:", nrow(bacteria_otus), "\n")
cat(
  "Percent OTUs removed as bacteria:",
  round(bacteria_summary$percent_otus_removed_as_bacteria, 4),
  "%\n"
)

cat("\nTotal reads before bacterial removal:", total_reads_before_bacteria, "\n")
cat("Bacterial reads removed:", bacteria_reads, "\n")
cat("Percent reads removed as bacteria:", round(percent_bacteria_reads, 4), "%\n")
cat("Total reads after bacterial removal:", total_reads_after_bacteria, "\n")
cat("Final OTUs after bacterial removal:", nrow(otu_table_no_bacteria), "\n")

write_tsv(bacteria_otus_with_taxonomy_and_reads, bacteria_otus_file, na = "NA")
write_tsv(bacteria_summary, bacteria_summary_file, na = "NA")

# ============================================================
# 5. Write final outputs
# ============================================================

write_tsv(otu_table_no_bacteria, filtered_otu_table_file)
write_tsv(taxonomy_no_bacteria, filtered_taxonomy_file, na = "NA")

cat("\n===== Final output checks =====\n")
cat("Original OTU table rows:", nrow(otu_table_no_mammals) + nrow(mammal_otus), "\n")
cat("Final OTU table rows:", nrow(otu_table_no_bacteria), "\n")
cat("Original taxonomy rows:", nrow(taxonomy_no_mammals) + nrow(mammal_otus), "\n")
cat("Final taxonomy rows:", nrow(taxonomy_no_bacteria), "\n")

cat("\nMammal OTUs written to:", mammal_otus_file, "\n")
cat("Blank contaminant OTUs written to:", flagged_contaminant_otus_file, "\n")
cat("Bacterial OTUs written to:", bacteria_otus_file, "\n")
cat("Final filtered OTU table written to:", filtered_otu_table_file, "\n")
cat("Final filtered taxonomy written to:", filtered_taxonomy_file, "\n")
cat("Blank samples used written to:", blank_samples_used_file, "\n")
cat("Blank contaminant summary written to:", contaminant_summary_file, "\n")
cat("Bacteria removal summary written to:", bacteria_summary_file, "\n")
