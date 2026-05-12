## ============================================================
## LoBAU JMF sample metadata + environmental metadata workflow
##
## This script:
## 1. Reads two JMF Excel files
## 2. Keeps only selected JMF columns
## 3. Adds:
##      - dataset
##      - site
##      - act: dna/rna
## 4. Handles mixed date formats:
##      - 2021/04/06
##      - 2022-03-07
## 5. Combines both JMF files into one dataframe
## 6. Reads environmental metadata from lobau_data_combined.xlsx
## 7. Removes metadata rows where sample_type == "well_water"
## 8. Joins metadata to JMF samples by Date + site
## 9. Saves final outputs
## ============================================================


## -------------------------
## 0) Set directory and load packages
## -------------------------

setwd("/Users/quinnmoon/Downloads/Lobau_18S_16S_ITS/")

library(readxl)
library(dplyr)
library(lubridate)
library(tidyr)


## -------------------------
## 1) General helper: parse mixed date formats
##
## Handles:
## - Date objects
## - POSIX date-times
## - Excel numeric dates
## - character dates like "2021/04/06" or "2022-03-07"
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
## 2) Helper: process one JMF file
##
## Keeps:
## - JMF sample ID
## - Sample description
## - Date
## - Forward primer
## - Reverse primer
##
## Adds:
## - dataset
## - site
## - act
## -------------------------

process_jmf <- function(file, act_vector, dataset_name) {
  
  raw_df <- read_excel(file)
  
  ## Remove rows that are completely blank
  raw_df <- raw_df %>%
    filter(!if_all(everything(), is.na))
  
  ## Make sure dna/rna assignment vector matches number of rows
  if (nrow(raw_df) != length(act_vector)) {
    stop(
      paste0(
        "Row count mismatch for ", file, "\n",
        "Rows in file: ", nrow(raw_df), "\n",
        "Length of act_vector: ", length(act_vector), "\n"
      )
    )
  }
  
  processed_df <- raw_df %>%
    select(
      `JMF sample ID`,
      `Sample description`,
      Date,
      `Forward primer`,
      `Reverse primer`
    ) %>%
    mutate(
      ## Standardize date format
      Date = parse_mixed_date(Date),
      
      ## Trim sample description before checking prefixes
      Sample_desc_trimmed = trimws(`Sample description`),
      
      ## Assign site from beginning of Sample description
      site = case_when(
        startsWith(Sample_desc_trimmed, "Blank") ~ "Blank",
        startsWith(Sample_desc_trimmed, "D05")   ~ "D05",
        startsWith(Sample_desc_trimmed, "D5")    ~ "D05",
        startsWith(Sample_desc_trimmed, "D10")   ~ "D10",
        startsWith(Sample_desc_trimmed, "D15")   ~ "D15",
        startsWith(Sample_desc_trimmed, "ESW")   ~ "ESW",
        TRUE                                     ~ NA_character_
      ),
      
      ## Add dna/rna assignment and dataset name
      act = act_vector,
      dataset = dataset_name
    ) %>%
    select(
      dataset,
      `JMF sample ID`,
      `Sample description`,
      Date,
      site,
      act,
      `Forward primer`,
      `Reverse primer`
    )
  
  return(processed_df)
}


## -------------------------
## 3) Process first JMF file
##
## File:
## jmf_2302_05_d05_d10.xlsx
##
## Assignment:
## rows 1-50   = dna
## rows 51-102 = rna
## -------------------------

jmf_2302_05_d05_d10 <- process_jmf(
  file = "jmf_2302_05_d05_d10.xlsx",
  act_vector = c(
    rep("dna", 50),
    rep("rna", 52)
  ),
  dataset_name = "jmf_2302_05_d05_d10"
)


## -------------------------
## 4) Process second JMF file
##
## File:
## jmf_2301_09_esw_d15.xlsx
##
## Assignment:
## rows 1-78   = rna
## rows 79-153 = dna
## -------------------------

jmf_2301_09_esw_d15 <- process_jmf(
  file = "jmf_2301_09_esw_d15.xlsx",
  act_vector = c(
    rep("rna", 78),
    rep("dna", 75)
  ),
  dataset_name = "jmf_2301_09_esw_d15"
)


## -------------------------
## 5) Combine both JMF dataframes
## -------------------------

jmf_big_df <- bind_rows(
  jmf_2302_05_d05_d10,
  jmf_2301_09_esw_d15
)


## -------------------------
## 6) JMF sanity checks
## -------------------------

cat("Rows in jmf_2302_05_d05_d10:", nrow(jmf_2302_05_d05_d10), "\n")
cat("Rows in jmf_2301_09_esw_d15:", nrow(jmf_2301_09_esw_d15), "\n")
cat("Rows in combined jmf_big_df:", nrow(jmf_big_df), "\n\n")

## Count samples by dataset, site, and act
site_act_counts <- jmf_big_df %>%
  count(dataset, site, act) %>%
  arrange(dataset, site, act)

print(site_act_counts)

## Check rows where site could not be assigned
missing_site <- jmf_big_df %>%
  filter(is.na(site))

cat("\nRows with missing site:", nrow(missing_site), "\n")
print(missing_site, n = Inf)

## Check rows where Date could not be parsed
failed_dates <- jmf_big_df %>%
  filter(is.na(Date))

cat("\nRows with failed Date parsing:", nrow(failed_dates), "\n")
print(failed_dates, n = Inf)


## -------------------------
## 7) Optional check:
##    For each Date, do we have all site x act combinations?
##
## Expected non-blank combinations:
## - ESW dna
## - ESW rna
## - D05 dna
## - D05 rna
## - D10 dna
## - D10 rna
## - D15 dna
## - D15 rna
## -------------------------

presence_df <- jmf_big_df %>%
  filter(site != "Blank") %>%
  mutate(
    act = tolower(act),
    site_act = paste(site, toupper(act))
  ) %>%
  filter(
    site %in% c("ESW", "D05", "D10", "D15"),
    act %in% c("dna", "rna")
  ) %>%
  distinct(Date, site_act) %>%
  mutate(present = TRUE) %>%
  complete(
    Date,
    site_act = c(
      "ESW DNA", "ESW RNA",
      "D05 DNA", "D05 RNA",
      "D10 DNA", "D10 RNA",
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
## 8) Save combined JMF-only dataframe
## -------------------------

write.csv(
  jmf_big_df,
  "jmf_big_df_combined.csv",
  row.names = FALSE
)


## ============================================================
## Environmental metadata matching
##
## Metadata file:
## lobau_data_combined.xlsx
##
## Matching columns:
## - JMF Date  -> metadata date_yyyy_mm_dd
## - JMF site  -> metadata well_id, standardized to ESW/D05/D10/D15
##
## Metadata exclusion:
## - remove rows where sample_type == "well_water"
##
## Blanks:
## - kept in final table
## - usually do not match metadata, so metadata columns remain NA
## ============================================================


## -------------------------
## 9) Read metadata
## -------------------------

metadata_raw <- read_excel("lobau_data_combined.xlsx")


## -------------------------
## 10) Remove well_water metadata rows
##
## This keeps pumped groundwater and any other non-well_water rows.
## If sample_type is NA, this keeps those rows too.
## -------------------------

metadata_filtered <- metadata_raw %>%
  filter(is.na(sample_type) | sample_type != "well_water")


## -------------------------
## 11) Standardize metadata Date and site
## -------------------------

metadata_for_join <- metadata_filtered %>%
  mutate(
    ## Standardize metadata date to same Date class as JMF Date
    Date = parse_mixed_date(date_yyyy_mm_dd),
    
    ## Standardize metadata well/site IDs
    well_id_trimmed = trimws(as.character(well_id)),
    
    site = case_when(
      startsWith(well_id_trimmed, "ESW") ~ "ESW",
      startsWith(well_id_trimmed, "D05") ~ "D05",
      startsWith(well_id_trimmed, "D5")  ~ "D05",
      startsWith(well_id_trimmed, "D10") ~ "D10",
      startsWith(well_id_trimmed, "D15") ~ "D15",
      TRUE                               ~ well_id_trimmed
    )
  )


## -------------------------
## 12) Check metadata join keys
##
## Ideally there should be only one metadata row per Date + site
## after removing well_water.
##
## If duplicates exist, left_join may duplicate JMF rows.
## -------------------------

metadata_duplicate_keys <- metadata_for_join %>%
  count(Date, site) %>%
  filter(n > 1)

cat("\nDuplicate metadata Date/site keys:", nrow(metadata_duplicate_keys), "\n")
print(metadata_duplicate_keys, n = Inf)

if (nrow(metadata_duplicate_keys) > 0) {
  cat("\nWarning: Some metadata Date/site keys are duplicated.\n")
  cat("The join may create duplicate JMF sample rows.\n")
}


## -------------------------
## 13) Join JMF samples to metadata by Date + site
##
## left_join keeps all JMF samples.
## Metadata columns are appended where Date + site match.
## -------------------------

jmf_with_metadata <- jmf_big_df %>%
  left_join(
    metadata_for_join,
    by = c("Date", "site"),
    suffix = c("", "_metadata")
  )


## -------------------------
## 14) Post-join sanity checks
## -------------------------

cat("\nRows in jmf_big_df:", nrow(jmf_big_df), "\n")
cat("Rows after metadata join:", nrow(jmf_with_metadata), "\n")

if (nrow(jmf_with_metadata) > nrow(jmf_big_df)) {
  cat("\nWarning: joined dataframe has more rows than jmf_big_df.\n")
  cat("This usually means duplicate Date/site rows in metadata.\n")
}

## Check unmatched non-blank samples
## Blank samples are expected not to match metadata.
unmatched_samples <- jmf_with_metadata %>%
  filter(site != "Blank", is.na(well_id))

cat("\nUnmatched non-blank JMF samples:", nrow(unmatched_samples), "\n")

print(
  unmatched_samples %>%
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

## Confirm no well_water metadata was joined
joined_sample_types <- jmf_with_metadata %>%
  count(sample_type)

cat("\nJoined sample_type counts:\n")
print(joined_sample_types)


## -------------------------
## 15) View and save final matched table
## -------------------------

head(jmf_with_metadata)

View(jmf_with_metadata)

write.csv(
  jmf_with_metadata,
  "jmf_with_metadata.csv",
  row.names = FALSE
)

