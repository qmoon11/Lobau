
# ============================================================
# Prune ITS T-BAS Newick tree to OTUs in your final ITS OTU table
# and optionally map collapsed duplicate OTUs to representative tips
#
# Goal:
#   Trim the ITS Newick tree so it contains only tips corresponding
#   to OTUs in your final ITS OTU table.
#
# Handles:
#   1. Direct matches between OTU table feature IDs and tree tip labels
#   2. Flexible first-token and numeric-suffix matches
#   3. Optional T-BAS duplicate summary mapping:
#        duplicate/member OTU -> representative tree tip
#
# Outputs:
#   ITS_OTU_to_tree_tip_mapping.tsv
#   ITS_OTUs_not_found_in_tree_or_duplicates.tsv
#   ITS_tree_pruned_to_final_OTU_table.treefile
#   ITS_tree_pruned_to_final_OTU_table.pdf
#
# Notes:
#   - This script avoids dplyr::slice(), which can conflict with
#     Bioconductor/S4Vectors and cause Rle/list errors.
#   - If you have a T-BAS duplicates summary file, set duplicates_file
#     to its actual path.
# ============================================================


# -------------------------------
# Packages
# -------------------------------

library(ape)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)


# -------------------------------
# Input files
# -------------------------------

# T-BAS Newick tree
tree_file <- "treeTU3OSOG3.nwk"

# Your final ITS OTU table
otu_table_file <- "final_microbiome_tables/final_ITS_OTU_table_ESW_D15_DNA800_RNA100.tsv"

# Optional T-BAS duplicate summary file.
# Change this to the actual path/name if downloaded.
# If it does not exist, the script still runs using direct tree-tip matching only.
duplicates_file <- "duplicates_summary.tsv"


# -------------------------------
# Output files
# -------------------------------

pruned_tree_file <- "ITS_tree_pruned_to_final_OTU_table.treefile"

otu_tip_mapping_file <- "ITS_OTU_to_tree_tip_mapping.tsv"

missing_otus_file <- "ITS_OTUs_not_found_in_tree_or_duplicates.tsv"

pruned_tree_pdf_file <- "ITS_tree_pruned_to_final_OTU_table.pdf"


# ============================================================
# Safety checks
# ============================================================

if (!file.exists(tree_file)) {
  stop("Could not find tree file: ", tree_file)
}

if (!file.exists(otu_table_file)) {
  stop("Could not find ITS OTU table file: ", otu_table_file)
}


# ============================================================
# Read tree
# ============================================================

tree <- read.tree(tree_file)

tree_tips <- tibble(
  tree_tip_label = as.character(tree$tip.label)
) %>%
  mutate(
    tree_tip_first_token = str_split_fixed(tree_tip_label, "\\s+", 2)[, 1],
    tree_tip_numeric_suffix = str_extract(tree_tip_first_token, "[0-9]+$")
  )

cat("\n============================================================\n")
cat("Original tree summary\n")
cat("============================================================\n")
cat("Tree file:", tree_file, "\n")
cat("Tree tips:", nrow(tree_tips), "\n")
cat("Internal nodes:", tree$Nnode, "\n")


# ============================================================
# Read OTU table and get OTU IDs
# ============================================================

otu_table <- read_tsv(
  otu_table_file,
  show_col_types = FALSE
)

colnames(otu_table) <- trimws(colnames(otu_table))
colnames(otu_table)[1] <- sub("^\ufeff", "", colnames(otu_table)[1])
colnames(otu_table)[1] <- "otu_id"

otu_ids <- otu_table %>%
  transmute(
    otu_id = as.character(otu_id),
    otu_first_token = str_split_fixed(otu_id, "\\s+", 2)[, 1],
    otu_numeric_suffix = str_extract(otu_first_token, "[0-9]+$")
  ) %>%
  distinct()

cat("\n============================================================\n")
cat("OTU table summary\n")
cat("============================================================\n")
cat("OTU table:", otu_table_file, "\n")
cat("OTUs in OTU table:", nrow(otu_ids), "\n")


# ============================================================
# Direct matches between OTU IDs and tree tips
# ============================================================

# -------------------------------
# Named lookup vectors
# -------------------------------

tree_tip_by_label <- tree_tips$tree_tip_label
names(tree_tip_by_label) <- tree_tips$tree_tip_label

tree_tip_by_first_token <- tree_tips$tree_tip_label
names(tree_tip_by_first_token) <- tree_tips$tree_tip_first_token

tree_tip_by_numeric_suffix_tbl <- tree_tips %>%
  filter(
    !is.na(tree_tip_numeric_suffix)
  ) %>%
  distinct(
    tree_tip_numeric_suffix,
    .keep_all = TRUE
  )

tree_tip_numeric_lookup <- tree_tip_by_numeric_suffix_tbl$tree_tip_label
names(tree_tip_numeric_lookup) <- tree_tip_by_numeric_suffix_tbl$tree_tip_numeric_suffix


# -------------------------------
# Create direct mapping
# -------------------------------

direct_mapping <- otu_ids %>%
  mutate(
    tree_tip_exact = unname(tree_tip_by_label[otu_id]),
    tree_tip_first_token_match = unname(tree_tip_by_first_token[otu_first_token]),
    tree_tip_numeric_suffix_match = unname(tree_tip_numeric_lookup[otu_numeric_suffix]),
    
    tree_tip_label = case_when(
      !is.na(tree_tip_exact) ~ tree_tip_exact,
      !is.na(tree_tip_first_token_match) ~ tree_tip_first_token_match,
      !is.na(tree_tip_numeric_suffix_match) ~ tree_tip_numeric_suffix_match,
      TRUE ~ NA_character_
    ),
    
    match_type = case_when(
      !is.na(tree_tip_exact) ~ "exact",
      !is.na(tree_tip_first_token_match) ~ "first_token",
      !is.na(tree_tip_numeric_suffix_match) ~ "numeric_suffix",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(tree_tip_label)
  ) %>%
  transmute(
    otu_id = as.character(otu_id),
    otu_first_token = as.character(otu_first_token),
    otu_numeric_suffix = as.character(otu_numeric_suffix),
    tree_tip_label = as.character(tree_tip_label),
    match_type = as.character(match_type)
  ) %>%
  distinct()

cat("\n============================================================\n")
cat("Direct OTU-to-tree-tip matches\n")
cat("============================================================\n")

if (nrow(direct_mapping) > 0) {
  direct_mapping %>%
    count(match_type, name = "n_matches") %>%
    print(n = Inf)
} else {
  cat("No direct matches found.\n")
}


# ============================================================
# Optional: read T-BAS duplicate summary to recover collapsed tips
# ============================================================

duplicate_mapping <- tibble(
  duplicate_id = character(),
  representative_id = character()
)

if (file.exists(duplicates_file)) {
  
  cat("\n============================================================\n")
  cat("Reading T-BAS duplicate summary\n")
  cat("============================================================\n")
  cat("Duplicates file:", duplicates_file, "\n")
  
  dup_raw <- read_tsv(
    duplicates_file,
    show_col_types = FALSE,
    guess_max = 100000
  )
  
  cat("\nDuplicates file columns:\n")
  print(colnames(dup_raw))
  
  cn <- tolower(colnames(dup_raw))
  
  # Try to guess representative and duplicate/member columns.
  rep_idx <- which(
    str_detect(cn, "representative|rep|centroid")
  )
  
  dup_idx <- which(
    str_detect(cn, "duplicate|member|sequence|query|strain|taxon|otu")
  )
  
  rep_col <- if (length(rep_idx) > 0) {
    colnames(dup_raw)[rep_idx[1]]
  } else {
    NA_character_
  }
  
  dup_col <- if (length(dup_idx) > 0) {
    colnames(dup_raw)[dup_idx[1]]
  } else {
    NA_character_
  }
  
  cat("\nGuessed representative column:", rep_col, "\n")
  cat("Guessed duplicate/member column:", dup_col, "\n")
  
  if (!is.na(rep_col) && !is.na(dup_col)) {
    
    duplicate_mapping <- dup_raw %>%
      transmute(
        duplicate_id = as.character(.data[[dup_col]]),
        representative_id = as.character(.data[[rep_col]])
      ) %>%
      filter(
        !is.na(duplicate_id),
        !is.na(representative_id),
        duplicate_id != "",
        representative_id != ""
      ) %>%
      distinct()
    
    cat("\nDuplicate mapping rows:", nrow(duplicate_mapping), "\n")
    
  } else {
    
    warning(
      "Could not automatically identify representative and duplicate columns in duplicates file. ",
      "Inspect colnames(dup_raw) and set rep_col/dup_col manually if needed."
    )
  }
  
} else {
  
  warning(
    "Duplicates file not found: ",
    duplicates_file,
    ". Proceeding with direct tree-tip matching only."
  )
}


# ============================================================
# Use duplicate mapping to map collapsed OTUs to representative tips
# ============================================================

collapsed_mapping <- tibble(
  otu_id = character(),
  otu_first_token = character(),
  otu_numeric_suffix = character(),
  tree_tip_label = character(),
  match_type = character(),
  representative_id = character()
)

if (nrow(duplicate_mapping) > 0) {
  
  duplicate_mapping2 <- duplicate_mapping %>%
    mutate(
      duplicate_id = as.character(duplicate_id),
      representative_id = as.character(representative_id),
      duplicate_first_token = str_split_fixed(duplicate_id, "\\s+", 2)[, 1],
      representative_first_token = str_split_fixed(representative_id, "\\s+", 2)[, 1],
      duplicate_numeric_suffix = str_extract(duplicate_first_token, "[0-9]+$"),
      representative_numeric_suffix = str_extract(representative_first_token, "[0-9]+$")
    )
  
  # -------------------------------
  # Named duplicate -> representative lookups
  # -------------------------------
  
  rep_by_duplicate_id <- duplicate_mapping2$representative_id
  names(rep_by_duplicate_id) <- duplicate_mapping2$duplicate_id
  
  rep_by_duplicate_first_token <- duplicate_mapping2$representative_id
  names(rep_by_duplicate_first_token) <- duplicate_mapping2$duplicate_first_token
  
  dup_numeric_unique <- duplicate_mapping2 %>%
    filter(
      !is.na(duplicate_numeric_suffix)
    ) %>%
    distinct(
      duplicate_numeric_suffix,
      .keep_all = TRUE
    )
  
  rep_by_duplicate_numeric <- dup_numeric_unique$representative_id
  names(rep_by_duplicate_numeric) <- dup_numeric_unique$duplicate_numeric_suffix
  
  # -------------------------------
  # Map OTU IDs to representative IDs
  # -------------------------------
  
  otu_to_rep <- otu_ids %>%
    mutate(
      representative_exact = unname(rep_by_duplicate_id[otu_id]),
      representative_first_token_match = unname(rep_by_duplicate_first_token[otu_first_token]),
      representative_numeric_match = unname(rep_by_duplicate_numeric[otu_numeric_suffix]),
      
      representative_id = case_when(
        !is.na(representative_exact) ~ representative_exact,
        !is.na(representative_first_token_match) ~ representative_first_token_match,
        !is.na(representative_numeric_match) ~ representative_numeric_match,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(representative_id)
    ) %>%
    mutate(
      representative_id = as.character(representative_id),
      representative_first_token = str_split_fixed(representative_id, "\\s+", 2)[, 1],
      representative_numeric_suffix = str_extract(representative_first_token, "[0-9]+$")
    )
  
  # -------------------------------
  # Map representative IDs to tree tips
  # -------------------------------
  
  collapsed_mapping <- otu_to_rep %>%
    mutate(
      tree_tip_exact = unname(tree_tip_by_label[representative_id]),
      tree_tip_first_token_match = unname(tree_tip_by_first_token[representative_first_token]),
      tree_tip_numeric_suffix_match = unname(tree_tip_numeric_lookup[representative_numeric_suffix]),
      
      tree_tip_label = case_when(
        !is.na(tree_tip_exact) ~ tree_tip_exact,
        !is.na(tree_tip_first_token_match) ~ tree_tip_first_token_match,
        !is.na(tree_tip_numeric_suffix_match) ~ tree_tip_numeric_suffix_match,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(tree_tip_label)
    ) %>%
    transmute(
      otu_id = as.character(otu_id),
      otu_first_token = as.character(otu_first_token),
      otu_numeric_suffix = as.character(otu_numeric_suffix),
      tree_tip_label = as.character(tree_tip_label),
      match_type = "collapsed_duplicate_to_representative",
      representative_id = as.character(representative_id)
    ) %>%
    distinct()
}

cat("\n============================================================\n")
cat("Collapsed duplicate-to-representative matches\n")
cat("============================================================\n")

if (nrow(collapsed_mapping) > 0) {
  collapsed_mapping %>%
    count(match_type, name = "n_matches") %>%
    print(n = Inf)
} else {
  cat("No collapsed duplicate mappings found or duplicates file unavailable.\n")
}


# ============================================================
# Combine direct and collapsed mappings
#
# Avoids dplyr::slice(), which can conflict with Bioconductor
# S4Vectors/Rle classes in some R sessions.
# ============================================================

direct_mapping_clean <- direct_mapping %>%
  mutate(
    representative_id = NA_character_
  ) %>%
  mutate(
    otu_id = as.character(otu_id),
    otu_first_token = as.character(otu_first_token),
    otu_numeric_suffix = as.character(otu_numeric_suffix),
    tree_tip_label = as.character(tree_tip_label),
    match_type = as.character(match_type),
    representative_id = as.character(representative_id)
  )

collapsed_mapping_clean <- collapsed_mapping %>%
  mutate(
    otu_id = as.character(otu_id),
    otu_first_token = as.character(otu_first_token),
    otu_numeric_suffix = as.character(otu_numeric_suffix),
    tree_tip_label = as.character(tree_tip_label),
    match_type = as.character(match_type),
    representative_id = as.character(representative_id)
  )

otu_tip_mapping_all <- bind_rows(
  direct_mapping_clean,
  collapsed_mapping_clean
) %>%
  mutate(
    match_priority = case_when(
      match_type == "exact" ~ 1,
      match_type == "first_token" ~ 2,
      match_type == "numeric_suffix" ~ 3,
      match_type == "collapsed_duplicate_to_representative" ~ 4,
      TRUE ~ 99
    )
  ) %>%
  arrange(
    otu_id,
    match_priority,
    tree_tip_label
  )

# Keep the best match per OTU without using slice()
otu_tip_mapping <- otu_tip_mapping_all %>%
  group_by(
    otu_id
  ) %>%
  mutate(
    row_number_within_otu = dplyr::row_number()
  ) %>%
  ungroup() %>%
  filter(
    row_number_within_otu == 1
  ) %>%
  select(
    -row_number_within_otu,
    -match_priority
  )

write_tsv(
  otu_tip_mapping,
  otu_tip_mapping_file,
  na = "NA"
)

cat("\n============================================================\n")
cat("Final OTU-to-tree-tip mapping\n")
cat("============================================================\n")
cat("Mapping file:", otu_tip_mapping_file, "\n")
cat("Mapped OTUs:", nrow(otu_tip_mapping), "\n")

if (nrow(otu_tip_mapping) > 0) {
  otu_tip_mapping %>%
    count(match_type, name = "n_otus") %>%
    print(n = Inf)
}


# ============================================================
# Report missing OTUs
# ============================================================

missing_otus <- otu_ids %>%
  anti_join(
    otu_tip_mapping,
    by = "otu_id"
  )

write_tsv(
  missing_otus,
  missing_otus_file,
  na = "NA"
)

cat("\n============================================================\n")
cat("Missing OTUs\n")
cat("============================================================\n")
cat("Missing OTUs not found in tree or duplicate mapping:", nrow(missing_otus), "\n")
cat("Missing OTUs file:", missing_otus_file, "\n")

if (nrow(missing_otus) > 0) {
  missing_otus %>%
    print(n = Inf)
}


# ============================================================
# Prune tree to mapped tips
# ============================================================

tips_to_keep <- otu_tip_mapping %>%
  pull(tree_tip_label) %>%
  unique() %>%
  as.character()

cat("\n============================================================\n")
cat("Pruning tree\n")
cat("============================================================\n")
cat("Unique tree tips to keep:", length(tips_to_keep), "\n")

if (length(tips_to_keep) == 0) {
  stop("No tree tips to keep. Check mapping.")
}

tips_to_keep_missing_from_tree <- setdiff(
  tips_to_keep,
  tree$tip.label
)

if (length(tips_to_keep_missing_from_tree) > 0) {
  stop(
    "Some mapped tips are not present in the tree: ",
    paste(tips_to_keep_missing_from_tree, collapse = ", ")
  )
}

pruned_tree <- keep.tip(
  phy = tree,
  tip = tips_to_keep,
  trim.internal = TRUE,
  subtree = FALSE,
  root.edge = 0,
  collapse.singles = FALSE
)

write.tree(
  pruned_tree,
  file = pruned_tree_file
)

cat("\nPruned tree written to:\n")
cat(pruned_tree_file, "\n")

cat("\nPruned tree tips:", length(pruned_tree$tip.label), "\n")
cat("Pruned tree internal nodes:", pruned_tree$Nnode, "\n")


# ============================================================
# Optional quick PDF plot
# ============================================================

pdf(
  pruned_tree_pdf_file,
  width = 10,
  height = 12
)

plot(
  pruned_tree,
  cex = 0.35,
  no.margin = FALSE,
  main = "ITS tree pruned to final OTU table tips"
)

axisPhylo()

dev.off()

cat("\nPruned tree PDF written to:\n")
cat(pruned_tree_pdf_file, "\n")


# ============================================================
# Final summary
# ============================================================

cat("\n============================================================\n")
cat("Done\n")
cat("============================================================\n")
cat("Original tree tips:", length(tree$tip.label), "\n")
cat("Final OTU table OTUs:", nrow(otu_ids), "\n")
cat("Mapped OTUs:", nrow(otu_tip_mapping), "\n")
cat("Missing OTUs:", nrow(missing_otus), "\n")
cat("Unique tree tips retained:", length(pruned_tree$tip.label), "\n")
cat("Pruned tree:", pruned_tree_file, "\n")
cat("Mapping table:", otu_tip_mapping_file, "\n")
cat("Missing OTUs table:", missing_otus_file, "\n")
cat("PDF:", pruned_tree_pdf_file, "\n")
