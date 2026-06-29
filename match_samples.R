
library(readr)
library(readxl)
library(dplyr)
library(lubridate)
library(tidyr)

## -------------------------
## Set directory
## -------------------------

setwd("/Users/quinnmoon/Downloads/Lobau_18S_16S_ITS/")


## -------------------------
## Helper: parse mixed date formats
## -------------------------

parse_mixed_date <- function(x) {
  
  if (inherits(x, "Date")) {
    return(x)
  }
  
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) {
    return(as.Date(x))
  }
  
  if (is.numeric(x)) {
    return(as.Date(x, origin = "1899-12-30"))
  }
  
  parsed <- suppressWarnings(ymd(trimws(as.character(x))))
  
  return(as.Date(parsed))
}


## -------------------------
## Helper: read JMF file
## Handles .csv and .xlsx
## -------------------------

read_jmf_file <- function(file) {
  
  if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
    read_excel(file)
  } else if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    read_csv(file, show_col_types = FALSE)
  } else {
    stop("File must be .csv or .xlsx: ", file)
  }
}


## -------------------------
## Helper: read isotope file
## Handles .tsv, .txt, .csv, and .xlsx
## -------------------------

read_isotope_file <- function(file) {
  
  if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
    read_excel(file)
  } else if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    read_csv(file, show_col_types = FALSE)
  } else if (grepl("\\.tsv$|\\.txt$", file, ignore.case = TRUE)) {
    read_tsv(file, show_col_types = FALSE)
  } else {
    stop("Isotope file must be .tsv, .txt, .csv, or .xlsx: ", file)
  }
}


## -------------------------
## Helper: process one JMF file
##
## Blank detection:
## A row is marked as Blank if ANY cell contains "blank",
## case-insensitive.
##
## Keeps only:
## - ESW
## - D15
## - Blank
##
## Primer columns are intentionally dropped.
## -------------------------

process_jmf_file <- function(file, act_vector, dataset_name) {
  
  raw_df <- read_jmf_file(file)
  
  ## Remove rows that are completely blank
  raw_df <- raw_df %>%
    filter(!if_all(everything(), is.na))
  
  ## Detect blank rows if ANY cell contains "blank"
  raw_df <- raw_df %>%
    mutate(
      row_contains_blank = if_any(
        everything(),
        ~ !is.na(.x) & grepl("blank", as.character(.x), ignore.case = TRUE)
      )
    )
  
  ## Make sure dna/rna assignment vector matches number of rows BEFORE site filtering
  if (nrow(raw_df) != length(act_vector)) {
    stop(
      paste0(
        "Row count mismatch for ", file, "\n",
        "Rows in file after removing completely blank rows: ", nrow(raw_df), "\n",
        "Length of act_vector: ", length(act_vector), "\n"
      )
    )
  }
  
  processed_df <- raw_df %>%
    select(
      `JMF sample ID`,
      `Sample description`,
      Date,
      row_contains_blank
    ) %>%
    mutate(
      Date = parse_mixed_date(Date),
      
      Sample_desc_trimmed = trimws(as.character(`Sample description`)),
      
      site = case_when(
        row_contains_blank ~ "Blank",
        startsWith(Sample_desc_trimmed, "ESW") ~ "ESW",
        startsWith(Sample_desc_trimmed, "D15") ~ "D15",
        TRUE ~ NA_character_
      ),
      
      act = act_vector,
      dataset = dataset_name
    ) %>%
    filter(site %in% c("ESW", "D15", "Blank")) %>%
    select(
      dataset,
      `JMF sample ID`,
      `Sample description`,
      Date,
      site,
      act
    )
  
  return(processed_df)
}


## -------------------------
## Process ESW/D15 JMF file
##
## This can be .xlsx or .csv
## -------------------------

jmf_esw_d15 <- process_jmf_file(
  file = "jmf_2301_09_esw_d15.xlsx",
  act_vector = c(
    rep("rna", 78),
    rep("dna", 75)
  ),
  dataset_name = "jmf_2301_09_esw_d15"
)


## -------------------------
## Main JMF dataframe
## -------------------------

jmf_big_df <- jmf_esw_d15


## -------------------------
## JMF sanity checks
## -------------------------

cat("Rows in jmf_esw_d15:", nrow(jmf_esw_d15), "\n")
cat("Rows in combined jmf_big_df:", nrow(jmf_big_df), "\n\n")

site_act_counts <- jmf_big_df %>%
  count(dataset, site, act) %>%
  arrange(dataset, site, act)

print(site_act_counts)

missing_site <- jmf_big_df %>%
  filter(is.na(site))

cat("\nRows with missing site:", nrow(missing_site), "\n")
print(missing_site, n = Inf)

failed_dates <- jmf_big_df %>%
  filter(is.na(Date))

cat("\nRows with failed Date parsing:", nrow(failed_dates), "\n")
print(failed_dates, n = Inf)


## -------------------------
## Check ESW/D15 site x act combinations
## -------------------------

presence_df <- jmf_big_df %>%
  filter(site != "Blank") %>%
  mutate(
    act = tolower(act),
    site_act = paste(site, toupper(act))
  ) %>%
  filter(
    site %in% c("ESW", "D15"),
    act %in% c("dna", "rna")
  ) %>%
  distinct(Date, site_act) %>%
  mutate(present = TRUE) %>%
  complete(
    Date,
    site_act = c(
      "ESW DNA", "ESW RNA",
      "D15 DNA", "D15 RNA"
    ),
    fill = list(present = FALSE)
  ) %>%
  pivot_wider(
    names_from = site_act,
    values_from = present
  ) %>%
  arrange(Date)

missing_presence <- presence_df %>%
  pivot_longer(
    cols = -Date,
    names_to = "site_act",
    values_to = "present"
  ) %>%
  filter(present == FALSE) %>%
  arrange(Date, site_act)

cat("\nMissing Date/site/act combinations:", nrow(missing_presence), "\n")
print(missing_presence, n = Inf)


## -------------------------
## Save JMF-only dataframe
## -------------------------

write.csv(
  jmf_big_df,
  "jmf_esw_d15_combined.csv",
  row.names = FALSE
)


## ============================================================
## Environmental metadata matching
##
## IMPORTANT:
## - Metadata may contain sites other than ESW/D15.
## - We keep ALL metadata rows after optional well_water filtering.
## - If a metadata row does not match a JMF sample,
##   JMF columns will be NA.
## - This is done with full_join().
## ============================================================

metadata_raw <- read_excel("lobau_data_combined.xlsx")


## -------------------------
## Remove well_water metadata rows
##
## If you truly want every metadata row including well_water,
## comment out this block and set:
## metadata_filtered <- metadata_raw
## -------------------------

metadata_filtered <- metadata_raw %>%
  filter(is.na(sample_type) | sample_type != "well_water")


## -------------------------
## Standardize metadata Date and site
##
## Keeps all metadata sites.
## -------------------------

metadata_for_join <- metadata_filtered %>%
  mutate(
    Date = parse_mixed_date(date_yyyy_mm_dd),
    
    well_id_trimmed = trimws(as.character(well_id)),
    
    site = case_when(
      grepl("blank", well_id_trimmed, ignore.case = TRUE) ~ "Blank",
      startsWith(well_id_trimmed, "ESW") ~ "ESW",
      startsWith(well_id_trimmed, "D05") ~ "D05",
      startsWith(well_id_trimmed, "D5")  ~ "D05",
      startsWith(well_id_trimmed, "D10") ~ "D10",
      startsWith(well_id_trimmed, "D15") ~ "D15",
      TRUE ~ well_id_trimmed
    )
  )


## -------------------------
## Check metadata duplicate Date/site keys
## -------------------------

metadata_duplicate_keys <- metadata_for_join %>%
  count(Date, site) %>%
  filter(n > 1)

cat("\nDuplicate metadata Date/site keys:", nrow(metadata_duplicate_keys), "\n")
print(metadata_duplicate_keys, n = Inf)

if (nrow(metadata_duplicate_keys) > 0) {
  cat("\nWarning: Some metadata Date/site keys are duplicated.\n")
  cat("The join may create duplicate rows.\n")
}


## -------------------------
## FULL JOIN JMF samples to environmental metadata
##
## full_join keeps:
## - JMF rows with matching metadata
## - JMF rows without metadata
## - metadata rows without JMF samples
## -------------------------

jmf_with_metadata <- full_join(
  jmf_big_df,
  metadata_for_join,
  by = c("Date", "site"),
  suffix = c("", "_metadata")
)


## ============================================================
## Isotope metadata matching
##
## IMPORTANT:
## - Isotope file may contain sites other than ESW/D15.
## - Keep all isotope rows.
## - If isotope rows do not match JMF/environmental metadata,
##   JMF/environmental columns will be NA.
## ============================================================


## -------------------------
## Read isotope file
##
## Update filename if needed.
## -------------------------

isotopes_raw <- read_isotope_file("lobau_isotopes_july2020_to_july2021.xlsx")


## -------------------------
## Standardize isotope Date and site
##
## Keeps all isotope sites.
## -------------------------

isotopes_for_join <- isotopes_raw %>%
  mutate(
    Date = parse_mixed_date(sampling_date_mm_dd_yyyy),
    
    sample_id_trimmed = trimws(as.character(sample_id)),
    
    site = case_when(
      startsWith(sample_id_trimmed, "ESW") ~ "ESW",
      startsWith(sample_id_trimmed, "D05") ~ "D05",
      startsWith(sample_id_trimmed, "D5")  ~ "D05",
      startsWith(sample_id_trimmed, "D10") ~ "D10",
      startsWith(sample_id_trimmed, "D15") ~ "D15",
      TRUE ~ sub("_.*$", "", sample_id_trimmed)
    )
  ) %>%
  rename(
    isotope_sample_id = sample_id,
    isotope_sampling_date = sampling_date_mm_dd_yyyy,
    isotope_sample_type = sample_type
  ) %>%
  select(
    Date,
    site,
    isotope_sample_id,
    isotope_sampling_date,
    isotope_sample_type,
    dO18_permil,
    stdev_dO18_permil,
    dH2_permil,
    stdev_dH2_permil,
    d_excess_permil
  )


## -------------------------
## Check isotope duplicate Date/site keys
## -------------------------

isotope_duplicate_keys <- isotopes_for_join %>%
  count(Date, site) %>%
  filter(n > 1)

cat("\nDuplicate isotope Date/site keys:", nrow(isotope_duplicate_keys), "\n")
print(isotope_duplicate_keys, n = Inf)

if (nrow(isotope_duplicate_keys) > 0) {
  cat("\nWarning: Some isotope Date/site keys are duplicated.\n")
  cat("The isotope join may create duplicate rows.\n")
}


## -------------------------
## FULL JOIN isotope metadata by Date + site
##
## This keeps isotope-only rows too.
## -------------------------

jmf_with_metadata <- full_join(
  jmf_with_metadata,
  isotopes_for_join,
  by = c("Date", "site")
)


## -------------------------
## Post-join sanity checks
## -------------------------

cat("\nRows in JMF-only table:", nrow(jmf_big_df), "\n")
cat("Rows in environmental metadata table:", nrow(metadata_for_join), "\n")
cat("Rows in isotope metadata table:", nrow(isotopes_for_join), "\n")
cat("Rows after full metadata + isotope joins:", nrow(jmf_with_metadata), "\n")


## Rows with JMF samples but no environmental metadata
jmf_without_env_metadata <- jmf_with_metadata %>%
  filter(!is.na(`JMF sample ID`), site != "Blank", is.na(well_id))

cat("\nJMF rows without environmental metadata:", nrow(jmf_without_env_metadata), "\n")

print(
  jmf_without_env_metadata %>%
    select(
      dataset,
      `JMF sample ID`,
      `Sample description`,
      Date,
      site,
      act
    ),
  n = Inf
)


## Environmental metadata rows without JMF samples
env_metadata_without_jmf <- jmf_with_metadata %>%
  filter(is.na(`JMF sample ID`), !is.na(well_id))

cat("\nEnvironmental metadata rows without JMF samples:", nrow(env_metadata_without_jmf), "\n")

print(
  env_metadata_without_jmf %>%
    select(
      Date,
      site,
      well_id,
      sample_type
    ),
  n = Inf
)


## JMF rows without isotope metadata
jmf_without_isotopes <- jmf_with_metadata %>%
  filter(!is.na(`JMF sample ID`), site != "Blank", is.na(isotope_sample_id))

cat("\nJMF rows without isotope metadata:", nrow(jmf_without_isotopes), "\n")

print(
  jmf_without_isotopes %>%
    select(
      dataset,
      `JMF sample ID`,
      `Sample description`,
      Date,
      site,
      act
    ),
  n = Inf
)


## Isotope rows without JMF samples
isotopes_without_jmf <- jmf_with_metadata %>%
  filter(is.na(`JMF sample ID`), !is.na(isotope_sample_id))

cat("\nIsotope rows without JMF samples:", nrow(isotopes_without_jmf), "\n")

print(
  isotopes_without_jmf %>%
    select(
      Date,
      site,
      isotope_sample_id,
      isotope_sample_type,
      dO18_permil,
      dH2_permil,
      d_excess_permil
    ),
  n = Inf
)


## Joined sample type counts
joined_sample_types <- jmf_with_metadata %>%
  count(sample_type)

cat("\nJoined environmental sample_type counts:\n")
print(joined_sample_types)

joined_isotope_sample_types <- jmf_with_metadata %>%
  count(isotope_sample_type)

cat("\nJoined isotope sample_type counts:\n")
print(joined_isotope_sample_types)


## -------------------------
## View and save final matched table
## -------------------------

head(jmf_with_metadata)

View(jmf_with_metadata)

write.csv(
  jmf_with_metadata,
  "jmf_esw_d15_with_all_metadata_and_isotopes.csv",
  row.names = FALSE
)
