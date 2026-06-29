#Description: combine PR2 and SILVA taxonomy assignments for 18S rRNA sequences, replacing PR2 Fungi assignments with SILVA assignments.
library(readr)
library(dplyr)

# Input files
pr2_file <- "18s_taxonomic_pr2_bootstrap.csv"
silva_file <- "18s_taxonomic_silva_bootstrap.csv"

# Output file
out_file <- "combined_PR2_with_SILVA_fungi.tsv"

# Read taxonomy files
pr2 <- read_csv(pr2_file, show_col_types = FALSE, na = c("", "NA"))
silva <- read_csv(silva_file, show_col_types = FALSE, na = c("", "NA"))

# Prepare SILVA table for joining
silva_for_join <- silva %>%
  select(
    ASV,
    Sequence_ID_silva = Sequence_ID,
    Kingdom_silva = Kingdom,
    Phylum_silva = Phylum,
    Class_silva = Class,
    Order_silva = Order,
    Family_silva = Family,
    Genus_silva = Genus
  )

# Start with PR2 assignments; replace taxa where PR2 Subdivision == "Fungi" with SILVA assignments
combined <- pr2 %>%
  left_join(silva_for_join, by = "ASV") %>%
  mutate(
    is_fungi_pr2 = if_else(is.na(Subdivision), FALSE, Subdivision == "Fungi"),
    
    final_sequence = Sequence_ID,
    
    final_Domain = if_else(
      is_fungi_pr2,
      "Eukaryota",
      Domain
    ),
    
    final_Supergroup = if_else(
      is_fungi_pr2,
      "Obazoa",
      Supergroup
    ),
    
    final_Division = if_else(
      is_fungi_pr2,
      "Opisthokonta",
      Division
    ),
    
    `final_Subdivision/Kingdom` = if_else(
      is_fungi_pr2,
      Kingdom_silva,
      Subdivision
    ),
    
    final_Phylum = if_else(
      is_fungi_pr2,
      Phylum_silva,
      NA_character_
    ),
    
    final_Class = if_else(
      is_fungi_pr2,
      Class_silva,
      Class
    ),
    
    final_Order = if_else(
      is_fungi_pr2,
      Order_silva,
      Order
    ),
    
    final_Family = if_else(
      is_fungi_pr2,
      Family_silva,
      Family
    ),
    
    final_Genus = if_else(
      is_fungi_pr2,
      Genus_silva,
      Genus
    ),
    
    final_OTU_ID = ASV
  ) %>%
  transmute(
    sequence = final_sequence,
    Domain = final_Domain,
    Supergroup = final_Supergroup,
    Division = final_Division,
    `Subdivision/Kingdom` = `final_Subdivision/Kingdom`,
    Phylum = final_Phylum,
    Class = final_Class,
    Order = final_Order,
    Family = final_Family,
    Genus = final_Genus,
    OTU_ID = final_OTU_ID
  ) %>%
  mutate(
    across(
      c(
        Domain,
        Supergroup,
        Division,
        `Subdivision/Kingdom`,
        Phylum,
        Class,
        Order,
        Family,
        Genus
      ),
      ~ if_else(
        grepl("_X|Incertae_Sedis", .x),
        NA_character_,
        .x
      )
    )
  )

View(combined)

# sanity checks
cat("Number of PR2 rows:", nrow(pr2), "\n")
cat("Number of combined rows:", nrow(combined), "\n")
cat("Number of PR2 Fungi replaced with SILVA:", sum(pr2$Subdivision == "Fungi", na.rm = TRUE), "\n")

combined %>%
  filter(`Subdivision/Kingdom` == "Fungi") %>%
  count(Supergroup, Division)

# Write combined taxonomy table
write_tsv(combined, out_file, na = "Unclassified")
