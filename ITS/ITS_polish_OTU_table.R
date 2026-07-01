
# ============================================================
# Update ITS OTU table using blank contaminants detected at ASV level
# AND remove specified human contaminant OTUs
#
# Logic:
# 1. Read ASV count table that still contains blank samples.
# 2. Use metadata to identify DNA/RNA blank sample columns in the ASV table.
# 3. Flag any ASV with >0 reads in any DNA or RNA blank.
# 4. Read already-collapsed OTU table.
# 5. Remove an entire OTU if ANY representative ASV inside that OTU
#    was flagged from a blank.
# 6. Also remove specified human contaminant OTUs:
#      OTU98_78, OTU98_861, OTU98_281
# 7. Final OTU table has both:
#      - blank-associated OTUs removed
#      - human contaminant OTUs removed
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# -------------------------------
# Input files
# -------------------------------

ASV_table_file <- "JMF-2301-09__all__rRNA_ITS2_E__JMFR_MSG4_KP48D/DADA2_counts_as_matrix.tsv"
OTU_table_file <- "ITS/ITSA_OTU98_table_renamed.tsv"
metadata_file <- "jmf_esw_d15_with_all_metadata_and_isotopes.csv"

# -------------------------------
# Human contaminant OTUs to remove
# -------------------------------

human_contaminant_otus <- c(
  "OTU98_78",
  "OTU98_861",
  "OTU98_281"
)

# -------------------------------
# Output files
# -------------------------------

flagged_asvs_file <- "ITS_blank_flagged_ASVs.tsv"
flagged_otus_file <- "ITS_OTUs_removed_due_to_blank_ASVs.tsv"
human_contaminant_otus_file <- "ITS_human_contaminant_OTUs_removed.tsv"

updated_otu_table_file <- "ITS/ITSA_OTU98_table_no_blank_ASV_contaminant_OTUs_no_human_contaminants.tsv"

blank_samples_used_file <- "ITS_ASV_blank_samples_used.tsv"
summary_file <- "ITS_blank_ASV_and_human_contaminant_filtering_summary.tsv"
reads_removed_per_sample_file <- "ITS_reads_removed_per_sample_blank_and_human.tsv"

# -------------------------------
# Helper: clean sample IDs
#
# Handles:
# JMF-2301-09-0152B-ITS   -> JMF-2301-09-0152
# JMF-2301-09-0152B-ITS.F -> JMF-2301-09-0152
# JMF-2301-09-0152B-ITS.R -> JMF-2301-09-0152
# JMF-2301-09-0152B       -> JMF-2301-09-0152
# ESW_20210406_rna        -> ESW_20210406_rna
# -------------------------------

sample_core_id <- function(x) {
  x <- trimws(as.character(x))
  
  # Remove FASTQ-like suffixes if present
  x <- sub("\\.[12](\\.filtered)?\\.fastq\\.gz$", "", x, ignore.case = TRUE)
  
  # Remove marker/orientation suffixes
  x <- sub("[-.](ITS|18S|16S)(\\.[FR])?$", "", x, ignore.case = TRUE)
  
  # Remove one trailing letter from JMF IDs only
  x <- ifelse(
    grepl("^JMF-", x),
    sub("[A-Za-z]$", "", x),
    x
  )
  
  return(x)
}

# -------------------------------
# Helper: split representative ASV strings
# -------------------------------

split_representatives <- function(x) {
  reps <- unlist(strsplit(as.character(x), ";"))
  reps <- trimws(reps)
  reps <- reps[!is.na(reps) & reps != ""]
  reps
}

# -------------------------------
# Read ASV table, OTU table, and metadata
# -------------------------------

asv_table <- read_tsv(ASV_table_file, show_col_types = FALSE)
otu_table <- read_tsv(OTU_table_file, show_col_types = FALSE)
metadata <- read_csv(metadata_file, show_col_types = FALSE, na = c("", "NA"))

# Clean column names
colnames(asv_table) <- trimws(colnames(asv_table))
colnames(otu_table) <- trimws(colnames(otu_table))
colnames(metadata) <- trimws(colnames(metadata))

colnames(asv_table)[1] <- sub("^\ufeff", "", colnames(asv_table)[1])
colnames(otu_table)[1] <- sub("^\ufeff", "", colnames(otu_table)[1])
colnames(metadata)[1] <- sub("^\ufeff", "", colnames(metadata)[1])

# -------------------------------
# Standardize ASV table
# -------------------------------

# Assume first column contains ASV IDs
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
# Standardize OTU table
# -------------------------------

# Expected:
# OTU    representatives    sample1    sample2 ...
if (!"OTU" %in% colnames(otu_table)) {
  colnames(otu_table)[1] <- "OTU"
}

if (!"representatives" %in% colnames(otu_table)) {
  stop("OTU table must contain a 'representatives' column listing ASVs in each OTU.")
}

otu_table <- otu_table %>%
  mutate(
    OTU = as.character(OTU),
    representatives = as.character(representatives)
  )

otu_count_cols <- setdiff(colnames(otu_table), c("OTU", "representatives"))

otu_table <- otu_table %>%
  mutate(
    across(
      all_of(otu_count_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

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
# Match metadata blanks to ASV table sample columns
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
# Identify blank columns in ASV table
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

# Also catch any ASV table columns that directly contain "blank"
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
# Flag ASVs found in any DNA or RNA blank
# -------------------------------

nonblank_asv_sample_cols <- setdiff(asv_count_cols, blank_sample_cols)

flagged_asvs <- asv_table %>%
  mutate(
    reads_in_blanks = rowSums(across(all_of(blank_sample_cols)), na.rm = TRUE),
    reads_in_nonblank_samples = rowSums(across(all_of(nonblank_asv_sample_cols)), na.rm = TRUE),
    total_reads_for_ASV = rowSums(across(all_of(asv_count_cols)), na.rm = TRUE)
  ) %>%
  filter(reads_in_blanks > 0) %>%
  select(
    ASV,
    reads_in_blanks,
    reads_in_nonblank_samples,
    total_reads_for_ASV
  ) %>%
  arrange(desc(reads_in_blanks), desc(reads_in_nonblank_samples))

flagged_asv_ids <- flagged_asvs$ASV

cat("\nTotal ASVs in ASV table:", nrow(asv_table), "\n")
cat("ASVs flagged from DNA or RNA blanks:", length(flagged_asv_ids), "\n")
cat(
  "Percent ASVs flagged:",
  round(100 * length(flagged_asv_ids) / nrow(asv_table), 4),
  "%\n"
)

write_tsv(flagged_asvs, flagged_asvs_file, na = "NA")

# -------------------------------
# Flag OTUs whose representatives include any flagged ASV
# -------------------------------

otu_representative_check <- otu_table %>%
  mutate(
    representative_list = map(representatives, split_representatives),
    flagged_representative_ASVs = map_chr(
      representative_list,
      ~ paste(intersect(.x, flagged_asv_ids), collapse = ";")
    ),
    n_flagged_representative_ASVs = if_else(
      flagged_representative_ASVs == "",
      0L,
      str_count(flagged_representative_ASVs, ";") + 1L
    ),
    remove_due_to_blank_ASV = n_flagged_representative_ASVs > 0
  )

flagged_otus <- otu_representative_check %>%
  filter(remove_due_to_blank_ASV) %>%
  select(
    OTU,
    representatives,
    flagged_representative_ASVs,
    n_flagged_representative_ASVs,
    all_of(otu_count_cols)
  )

cat("\nTotal OTUs in OTU table:", nrow(otu_table), "\n")
cat("OTUs removed because at least one representative ASV was in a blank:", nrow(flagged_otus), "\n")
cat(
  "Percent OTUs removed due to blank ASVs:",
  round(100 * nrow(flagged_otus) / nrow(otu_table), 4),
  "%\n"
)

write_tsv(flagged_otus, flagged_otus_file, na = "NA")

# -------------------------------
# Calculate reads removed from OTU table due to blank ASV representative OTUs
# -------------------------------

total_reads_otu_table <- otu_table %>%
  summarise(
    total_reads = sum(across(all_of(otu_count_cols)), na.rm = TRUE)
  ) %>%
  pull(total_reads)

reads_in_blank_removed_otus <- flagged_otus %>%
  summarise(
    removed_reads = sum(across(all_of(otu_count_cols)), na.rm = TRUE)
  ) %>%
  pull(removed_reads)

percent_reads_removed_blank <- ifelse(
  total_reads_otu_table > 0,
  100 * reads_in_blank_removed_otus / total_reads_otu_table,
  NA_real_
)

cat("\nTotal reads in original OTU table:", total_reads_otu_table, "\n")
cat("Reads in OTUs removed due to blank ASVs:", reads_in_blank_removed_otus, "\n")
cat("Percent reads removed due to blank ASVs:", round(percent_reads_removed_blank, 4), "%\n")

# -------------------------------
# Remove blank-associated OTUs from OTU table
# -------------------------------

otu_table_no_blank_contam <- otu_table %>%
  anti_join(
    flagged_otus %>% select(OTU),
    by = "OTU"
  )

# ============================================================
# Remove human contaminant OTUs
# ============================================================

human_contaminant_otus_present <- otu_table_no_blank_contam %>%
  filter(OTU %in% human_contaminant_otus)

human_contaminant_otus_missing <- setdiff(
  human_contaminant_otus,
  human_contaminant_otus_present$OTU
)

cat("\nHuman contaminant OTUs requested for removal:\n")
print(human_contaminant_otus)

cat("\nHuman contaminant OTUs present after blank filtering:", nrow(human_contaminant_otus_present), "\n")
print(human_contaminant_otus_present %>% select(OTU, representatives), n = Inf)

if (length(human_contaminant_otus_missing) > 0) {
  cat("\nHuman contaminant OTUs not present after blank filtering, maybe already removed or absent:\n")
  print(human_contaminant_otus_missing)
}

# Reads in human contaminant OTUs after blank filtering
total_reads_after_blank_filter <- otu_table_no_blank_contam %>%
  summarise(
    total_reads = sum(across(all_of(otu_count_cols)), na.rm = TRUE)
  ) %>%
  pull(total_reads)

reads_in_human_contaminant_otus <- human_contaminant_otus_present %>%
  summarise(
    human_contaminant_reads = sum(across(all_of(otu_count_cols)), na.rm = TRUE)
  ) %>%
  pull(human_contaminant_reads)

percent_reads_removed_human_after_blank <- ifelse(
  total_reads_after_blank_filter > 0,
  100 * reads_in_human_contaminant_otus / total_reads_after_blank_filter,
  NA_real_
)

percent_reads_removed_human_original_total <- ifelse(
  total_reads_otu_table > 0,
  100 * reads_in_human_contaminant_otus / total_reads_otu_table,
  NA_real_
)

cat("\nTotal reads after blank-ASV contaminant OTU removal:", total_reads_after_blank_filter, "\n")
cat("Reads in human contaminant OTUs to remove:", reads_in_human_contaminant_otus, "\n")
cat(
  "Percent of post-blank-filter reads removed as human contaminants:",
  round(percent_reads_removed_human_after_blank, 4),
  "%\n"
)
cat(
  "Percent of original OTU-table reads removed as human contaminants:",
  round(percent_reads_removed_human_original_total, 4),
  "%\n"
)

# Per-human-OTU read totals
human_contaminant_otu_summary <- human_contaminant_otus_present %>%
  rowwise() %>%
  mutate(
    total_reads_for_OTU = sum(c_across(all_of(otu_count_cols)), na.rm = TRUE),
    percent_of_post_blank_reads = ifelse(
      total_reads_after_blank_filter > 0,
      100 * total_reads_for_OTU / total_reads_after_blank_filter,
      NA_real_
    ),
    percent_of_original_reads = ifelse(
      total_reads_otu_table > 0,
      100 * total_reads_for_OTU / total_reads_otu_table,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  select(
    OTU,
    representatives,
    total_reads_for_OTU,
    percent_of_post_blank_reads,
    percent_of_original_reads
  )

cat("\nHuman contaminant OTU read breakdown:\n")
print(human_contaminant_otu_summary, n = Inf)

write_tsv(human_contaminant_otu_summary, human_contaminant_otus_file, na = "NA")

# -------------------------------
# Remove human contaminant OTUs from the blank-filtered table
# -------------------------------

otu_table_clean <- otu_table_no_blank_contam %>%
  filter(!OTU %in% human_contaminant_otus)

write_tsv(otu_table_clean, updated_otu_table_file, na = "NA")

# -------------------------------
# Final read totals after both filters
# -------------------------------

total_reads_final <- otu_table_clean %>%
  summarise(
    total_reads = sum(across(all_of(otu_count_cols)), na.rm = TRUE)
  ) %>%
  pull(total_reads)

total_reads_removed_all_filters <- total_reads_otu_table - total_reads_final

percent_reads_removed_all_filters <- ifelse(
  total_reads_otu_table > 0,
  100 * total_reads_removed_all_filters / total_reads_otu_table,
  NA_real_
)

cat("\nTotal reads in final OTU table after blank + human removal:", total_reads_final, "\n")
cat("Total reads removed by all filters:", total_reads_removed_all_filters, "\n")
cat("Percent reads removed by all filters:", round(percent_reads_removed_all_filters, 4), "%\n")

# -------------------------------
# Per-sample read summary
# -------------------------------

total_reads_per_sample <- otu_table %>%
  summarise(
    across(
      all_of(otu_count_cols),
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "total_reads_original"
  )

blank_removed_reads_per_sample <- flagged_otus %>%
  summarise(
    across(
      all_of(otu_count_cols),
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "blank_removed_reads"
  )

human_removed_reads_per_sample <- human_contaminant_otus_present %>%
  summarise(
    across(
      all_of(otu_count_cols),
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "human_removed_reads"
  )

final_reads_per_sample <- otu_table_clean %>%
  summarise(
    across(
      all_of(otu_count_cols),
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "final_reads"
  )

reads_removed_per_sample <- total_reads_per_sample %>%
  left_join(blank_removed_reads_per_sample, by = "sample") %>%
  left_join(human_removed_reads_per_sample, by = "sample") %>%
  left_join(final_reads_per_sample, by = "sample") %>%
  mutate(
    blank_removed_reads = if_else(is.na(blank_removed_reads), 0, blank_removed_reads),
    human_removed_reads = if_else(is.na(human_removed_reads), 0, human_removed_reads),
    final_reads = if_else(is.na(final_reads), 0, final_reads),
    total_removed_reads = blank_removed_reads + human_removed_reads,
    percent_reads_removed = if_else(
      total_reads_original > 0,
      100 * total_removed_reads / total_reads_original,
      NA_real_
    ),
    percent_blank_removed = if_else(
      total_reads_original > 0,
      100 * blank_removed_reads / total_reads_original,
      NA_real_
    ),
    percent_human_removed = if_else(
      total_reads_original > 0,
      100 * human_removed_reads / total_reads_original,
      NA_real_
    )
  ) %>%
  arrange(desc(percent_reads_removed))

cat("\nReads removed per OTU-table sample:\n")
print(reads_removed_per_sample, n = Inf)

write_tsv(reads_removed_per_sample, reads_removed_per_sample_file, na = "NA")

# -------------------------------
# Summary output
# -------------------------------

summary_df <- tibble(
  total_asvs = nrow(asv_table),
  flagged_asvs_from_blanks = length(flagged_asv_ids),
  percent_asvs_flagged = 100 * length(flagged_asv_ids) / nrow(asv_table),
  
  total_otus_original = nrow(otu_table),
  otus_removed_due_to_blank_asv_representatives = nrow(flagged_otus),
  percent_otus_removed_due_to_blank_asvs = 100 * nrow(flagged_otus) / nrow(otu_table),
  
  human_contaminant_otus_requested = length(human_contaminant_otus),
  human_contaminant_otus_present_after_blank_filter = nrow(human_contaminant_otus_present),
  human_contaminant_otus_missing_after_blank_filter = length(human_contaminant_otus_missing),
  
  total_otus_after_blank_filter = nrow(otu_table_no_blank_contam),
  total_otus_final_after_blank_and_human_filter = nrow(otu_table_clean),
  
  total_reads_original_otu_table = total_reads_otu_table,
  reads_removed_due_to_blank_asv_otus = reads_in_blank_removed_otus,
  percent_reads_removed_due_to_blank_asv_otus = percent_reads_removed_blank,
  
  total_reads_after_blank_filter = total_reads_after_blank_filter,
  reads_removed_due_to_human_contaminant_otus = reads_in_human_contaminant_otus,
  percent_post_blank_reads_removed_due_to_human_contaminants = percent_reads_removed_human_after_blank,
  percent_original_reads_removed_due_to_human_contaminants = percent_reads_removed_human_original_total,
  
  total_reads_final = total_reads_final,
  total_reads_removed_all_filters = total_reads_removed_all_filters,
  percent_reads_removed_all_filters = percent_reads_removed_all_filters,
  
  blank_sample_columns_used_from_asv_table = length(blank_sample_cols)
)

write_tsv(summary_df, summary_file, na = "NA")

cat("\n===== Final summary =====\n")
print(summary_df)

cat("\nFlagged ASVs written to:", flagged_asvs_file, "\n")
cat("Blank-associated OTUs written to:", flagged_otus_file, "\n")
cat("Human contaminant OTUs written to:", human_contaminant_otus_file, "\n")
cat("Final OTU table written to:", updated_otu_table_file, "\n")
cat("Blank samples used written to:", blank_samples_used_file, "\n")
cat("Reads removed per sample written to:", reads_removed_per_sample_file, "\n")
cat("Summary written to:", summary_file, "\n")


# -------------------------------
# Trim taxonomy table to only remaining OTUs
# Add this after the final cleaned OTU table has been created:
# otu_table_clean
# -------------------------------

taxonomy_file <- "ClassifyITSoutputs/initial_assignments.csv"
updated_taxonomy_file <- "ClassifyITSoutputs/initial_assingments_no_blank_no_human_contaminants.csv"

taxonomy <- read_csv(
  taxonomy_file,
  show_col_types = FALSE,
  na = c("", "NA"),
  trim_ws = TRUE
)

# Clean taxonomy column names
colnames(taxonomy) <- trimws(colnames(taxonomy))
colnames(taxonomy)[1] <- sub("^\ufeff", "", colnames(taxonomy)[1])

# Check taxonomy ID column
if (!"qseqid" %in% colnames(taxonomy)) {
  stop(
    "Taxonomy file must contain column 'qseqid'. Columns found: ",
    paste(colnames(taxonomy), collapse = ", ")
  )
}

taxonomy <- taxonomy %>%
  mutate(qseqid = as.character(qseqid))

remaining_otus <- otu_table_clean %>%
  select(OTU) %>%
  mutate(OTU = as.character(OTU)) %>%
  distinct()

taxonomy_clean <- taxonomy %>%
  semi_join(
    remaining_otus,
    by = c("qseqid" = "OTU")
  )

cat("\nOriginal taxonomy rows:", nrow(taxonomy), "\n")
cat("Cleaned taxonomy rows matching final OTU table:", nrow(taxonomy_clean), "\n")
cat("Taxonomy rows removed:", nrow(taxonomy) - nrow(taxonomy_clean), "\n")
cat("Remaining OTUs in final OTU table:", nrow(remaining_otus), "\n")
cat("Remaining OTUs with taxonomy:", nrow(taxonomy_clean), "\n")
cat("Remaining OTUs missing taxonomy:",
    nrow(remaining_otus) - nrow(taxonomy_clean),
    "\n")

write_csv(taxonomy_clean, updated_taxonomy_file, na = "NA")

cat("Final cleaned taxonomy written to:", updated_taxonomy_file, "\n")
