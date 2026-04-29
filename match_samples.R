## -------------------------
## 0) set directory, packages
## -------------------------
setwd("/Users/quinnmoon/Downloads/Lobau_18S_16S_ITS/")
library(readxl)
library(dplyr)


## -------------------------
## 1) Read JMF sheet (includes blanks)
## -------------------------
jmf_blank_info <- read_excel("jmf_2301_09_esw_d15.xlsx")

## If act is assigned by row order, guard against sheet size changes
stopifnot(nrow(jmf_blank_info) == 78 + 75)

jmf_blank_info <- jmf_blank_info %>%
  mutate(act = c(rep("rna", 78), rep("dna", 75))) %>%
  select(`JMF sample ID`, `Sample description`, Date, act) %>%
  rename(
    Sample_desc = `Sample description`,
    JMF.ID      = `JMF sample ID`
  ) %>%
  mutate(blank_sample = if_else(startsWith(Sample_desc, "Blank"),
                                "blank", "sample")) %>%
  mutate(
    well_id = case_when(
      startsWith(Sample_desc, "ESW") ~ "ESW_",
      startsWith(Sample_desc, "D05") ~ "D05p_",
      startsWith(Sample_desc, "D10") ~ "D10p_",
      startsWith(Sample_desc, "D15") ~ "D15p_",
      TRUE ~ NA_character_
    )
  ) %>%
  ## keep all 4 wells (drop this line if you truly want *everything*)
  filter(!is.na(well_id)) %>%
  mutate(
    date_tag      = format(as.Date(Date), "%Y%m%d"),
    sample_id     = paste0(well_id, date_tag),
    sample_id_act = paste0(sample_id, "_", act)
  ) %>%
  select(JMF.ID, Sample_desc, Date, blank_sample, well_id, sample_id_act, act)

## -------------------------
## 2) Read environmental metadata and duplicate for RNA/DNA
## -------------------------
env_data <- read_excel("lobau_data_combined.xlsx") %>%
  bind_rows(
    mutate(., sample_id_act = paste0(sample_id, "_rna")),
    mutate(., sample_id_act = paste0(sample_id, "_dna"))
  )

## -------------------------
## 3) LEFT JOIN: keep ALL JMF rows (env columns become NA if no match)
## -------------------------
env_lobau_final <- left_join(jmf_blank_info, env_data, by = "sample_id_act") %>%
  as.data.frame()

rownames(env_lobau_final) <- env_lobau_final$JMF.ID

## -------------------------
## 4) Sanity checks
## -------------------------
cat("JMF rows:", nrow(jmf_blank_info), "\n")
cat("Rows after left_join:", nrow(env_lobau_final), "\n")

unmatched <- sum(is.na(env_lobau_final$sample_id))  # assumes env_data has sample_id column
cat("Unmatched JMF rows (env fields NA):", unmatched, "\n")

head(env_lobau_final)
View(env_lobau_final)

# Save final table
write.csv(env_lobau_final, "env_lobau_final.csv", row.names = FALSE)
