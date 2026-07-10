
# ============================================================
# Prune ITS T-BAS Newick tree to study OTUs, then fungi-only OTUs,
# and make iTOL annotation files
#
# Main goals:
#   1. Read full ITS T-BAS Newick tree.
#   2. Read final ITS OTU table from this study.
#   3. Map study OTUs to tree tips.
#   4. Optionally use T-BAS duplicate summary to map collapsed OTUs
#      to representative tree tips.
#   5. Prune the tree to all study OTU tips.
#   6. Use taxonomy to identify confirmed fungal OTUs.
#   7. Prune the study tree to confirmed fungal tips only.
#   8. Make iTOL color strip/ring for fungal phylum.
#   9. Print total numbers at each filtering level.
#
# Inputs:
#   treeTU3OSOG3.nwk
#   final_microbiome_tables/final_ITS_OTU_table_ESW_D15_DNA800_RNA100.tsv
#   final_microbiome_tables/final_ITS_taxonomy_ESW_D15_DNA800_RNA100.csv
#   duplicates_summary.tsv, optional
#
# Outputs:
#   ITS_OTU_to_tree_tip_mapping.tsv
#   ITS_OTUs_not_found_in_tree_or_duplicates.tsv
#   ITS_tree_pruned_to_study_OTUs.treefile
#   ITS_study_tree_tip_taxonomy.tsv
#   ITS_study_tree_tips_EXCLUDED_not_confirmed_fungi.tsv
#   ITS_tree_pruned_to_study_OTUs_FUNGI_ONLY.treefile
#   ITS_fungi_only_tip_fungal_phylum_table.tsv
#   iTOL_ITS_fungi_only_phylum_colorstrip.txt
#   ITS_fungi_only_tree_colored_by_phylum.pdf
#   ITS_tree_pruning_summary.tsv
# ============================================================


# -------------------------------
# Packages
# -------------------------------

library(ape)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)


# ============================================================
# Input files
# ============================================================

tree_file <- "treeTU3OSOG3.nwk"

otu_table_file <- "final_microbiome_tables/final_ITS_OTU_table_ESW_D15_DNA800_RNA100.tsv"

taxonomy_file <- "final_microbiome_tables/final_ITS_taxonomy_ESW_D15_DNA800_RNA100.csv"

# Optional. If this file is not present, script still runs with direct matching.
duplicates_file <- "duplicates_summary.tsv"


# ============================================================
# Taxonomy parsing settings
# ============================================================

# Your ITS taxonomy file should have:
# qseqid, kingdom, phylum, class, order, family, genus, species, notes

tax_id_column <- "qseqid"
kingdom_column <- "kingdom"
phylum_column <- "phylum"

# If your taxonomy file has one full taxonomy string column instead,
# put its name here. Otherwise leave NA.
taxonomy_column <- NA_character_


# ============================================================
# Optional duplicate-summary column settings
#
# Leave as NA to auto-detect. If auto-detection fails, inspect:
#   colnames(read_tsv(duplicates_file))
# and set these manually.
# ============================================================

duplicates_representative_column <- NA_character_
duplicates_duplicate_column <- NA_character_


# ============================================================
# Output files
# ============================================================

otu_tip_mapping_file <- "ITS_OTU_to_tree_tip_mapping.tsv"

missing_otus_file <- "ITS_OTUs_not_found_in_tree_or_duplicates.tsv"

study_pruned_tree_file <- "ITS_tree_pruned_to_study_OTUs.treefile"

study_tip_taxonomy_file <- "ITS_study_tree_tip_taxonomy.tsv"

excluded_not_fungi_file <- "ITS_study_tree_tips_EXCLUDED_not_confirmed_fungi.tsv"

fungi_only_tree_file <- "ITS_tree_pruned_to_study_OTUs_FUNGI_ONLY.treefile"

fungi_only_tip_phylum_table_file <- "ITS_fungi_only_tip_fungal_phylum_table.tsv"

itol_fungi_only_phylum_colorstrip_file <- "iTOL_ITS_fungi_only_phylum_colorstrip.txt"

fungi_only_tree_pdf_file <- "ITS_fungi_only_tree_colored_by_phylum.pdf"

tree_pruning_summary_file <- "ITS_tree_pruning_summary.tsv"


# ============================================================
# Helper functions
# ============================================================

unknown_taxon_regex <- regex(
  paste(
    c(
      "^\\s*$",
      "^NA$",
      "^N/A$",
      "^none$",
      "^null$",
      "unclassified",
      "unassigned",
      "unknown",
      "no blast hit",
      "not assigned",
      "incertae sedis",
      "environmental sample",
      "uncultured"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

known_fungal_phyla <- c(
  "Ascomycota",
  "Basidiomycota",
  "Mortierellomycota",
  "Mucoromycota",
  "Chytridiomycota",
  "Glomeromycota",
  "Rozellomycota",
  "Zoopagomycota",
  "Kickxellomycota",
  "Blastocladiomycota",
  "Neocallimastigomycota",
  "Entomophthoromycota",
  "Olpidiomycota",
  "Aphelidiomycota",
  "Calcarisporiellomycota",
  "Monoblepharomycota",
  "Basidiobolomycota"
)

collapse_unique <- function(x) {
  y <- unique(na.omit(as.character(x)))
  y <- y[y != ""]
  
  if (length(y) == 0) {
    return(NA_character_)
  }
  
  paste(sort(y), collapse = ";")
}

clean_kingdom <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  x <- str_replace(x, regex("^k__", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^k:", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^kingdom[:=]", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^D_0__", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^d__0__", ignore_case = TRUE), "")
  
  x <- trimws(x)
  x <- str_replace_all(x, "^['\"]|['\"]$", "")
  
  unknown_idx <- !is.na(x) & str_detect(x, unknown_taxon_regex)
  x[unknown_idx] <- "Unclassified"
  x[is.na(x) | x == ""] <- "Unclassified"
  
  x[str_detect(x, regex("^fungi$", ignore_case = TRUE))] <- "Fungi"
  x[str_detect(x, regex("^bacteria$", ignore_case = TRUE))] <- "Bacteria"
  x[str_detect(x, regex("^archaea$", ignore_case = TRUE))] <- "Archaea"
  x[str_detect(x, regex("^animalia$", ignore_case = TRUE))] <- "Animalia"
  x[str_detect(x, regex("^plantae$", ignore_case = TRUE))] <- "Plantae"
  x[str_detect(x, regex("^protista$", ignore_case = TRUE))] <- "Protista"
  
  x
}

clean_phylum <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  x <- str_replace(x, regex("^p__", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^p:", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^phylum[:=]", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^D_1__", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^d__1__", ignore_case = TRUE), "")
  
  x <- trimws(x)
  x <- str_replace_all(x, "^['\"]|['\"]$", "")
  
  unknown_idx <- !is.na(x) & str_detect(x, unknown_taxon_regex)
  x[unknown_idx] <- "Unclassified"
  x[is.na(x) | x == ""] <- "Unclassified"
  
  x
}

split_taxonomy <- function(tax_string) {
  tax_string <- as.character(tax_string)
  
  if (is.na(tax_string) || trimws(tax_string) == "") {
    return(character())
  }
  
  parts <- unlist(
    str_split(
      tax_string,
      "\\s*;\\s*|\\s*\\|\\s*|\\s*,\\s*"
    )
  )
  
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  
  parts
}

extract_kingdom_one <- function(tax_string) {
  parts <- split_taxonomy(tax_string)
  
  if (length(parts) == 0) {
    return("Unclassified")
  }
  
  kingdom_token <- parts[
    str_detect(
      parts,
      regex("^\\s*(k__|k:|kingdom[:=]|D_0__|d__0__)", ignore_case = TRUE)
    )
  ]
  
  if (length(kingdom_token) > 0) {
    return(clean_kingdom(kingdom_token[1]))
  }
  
  clean_kingdom(parts[1])
}

extract_phylum_one <- function(tax_string) {
  parts <- split_taxonomy(tax_string)
  
  if (length(parts) == 0) {
    return("Unclassified")
  }
  
  phylum_token <- parts[
    str_detect(
      parts,
      regex("^\\s*(p__|p:|phylum[:=]|D_1__|d__1__)", ignore_case = TRUE)
    )
  ]
  
  if (length(phylum_token) > 0) {
    return(clean_phylum(phylum_token[1]))
  }
  
  if (length(parts) >= 2) {
    return(clean_phylum(parts[2]))
  }
  
  clean_phylum(parts[1])
}

extract_kingdom <- function(x) {
  vapply(x, extract_kingdom_one, character(1))
}

extract_phylum <- function(x) {
  vapply(x, extract_phylum_one, character(1))
}

classified_phyla <- function(x) {
  y <- clean_phylum(x)
  
  y <- y[
    !is.na(y) &
      y != "" &
      !str_detect(y, unknown_taxon_regex) &
      !y %in% c("Unclassified", "Mixed", "No_consensus")
  ]
  
  y
}

consensus_phylum <- function(x, min_fraction = 0.50) {
  y <- classified_phyla(x)
  
  if (length(y) == 0) {
    return("Unclassified")
  }
  
  tab <- sort(table(y), decreasing = TRUE)
  top_n <- as.integer(tab[1])
  total_n <- sum(tab)
  top_fraction <- top_n / total_n
  
  if (length(tab) > 1 && as.integer(tab[1]) == as.integer(tab[2])) {
    return("Mixed")
  }
  
  if (top_fraction < min_fraction) {
    return("Mixed")
  }
  
  names(tab)[1]
}

consensus_fraction <- function(x) {
  y <- classified_phyla(x)
  
  if (length(y) == 0) {
    return(NA_real_)
  }
  
  tab <- sort(table(y), decreasing = TRUE)
  as.integer(tab[1]) / sum(tab)
}

read_table_auto <- function(file) {
  header_line <- readLines(file, n = 1, warn = FALSE)
  
  if (str_count(header_line, ",") > str_count(header_line, "\t")) {
    read_csv(file, show_col_types = FALSE, guess_max = 100000)
  } else {
    read_tsv(file, show_col_types = FALSE, guess_max = 100000)
  }
}


# ============================================================
# Safety checks
# ============================================================

if (!file.exists(tree_file)) {
  stop("Could not find tree file: ", tree_file)
}

if (!file.exists(otu_table_file)) {
  stop("Could not find OTU table file: ", otu_table_file)
}

if (!file.exists(taxonomy_file)) {
  stop("Could not find taxonomy file: ", taxonomy_file)
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
cat("Original T-BAS tree summary\n")
cat("============================================================\n")
cat("Tree file:", tree_file, "\n")
cat("Original tree tips:", nrow(tree_tips), "\n")
cat("Original internal nodes:", tree$Nnode, "\n")


# ============================================================
# Read OTU table
# ============================================================

otu_table <- read_tsv(
  otu_table_file,
  show_col_types = FALSE,
  guess_max = 100000
)

colnames(otu_table) <- trimws(colnames(otu_table))
colnames(otu_table)[1] <- sub("^\ufeff", "", colnames(otu_table)[1])
colnames(otu_table)[1] <- "otu_id"

otu_table_colnames_lower <- tolower(colnames(otu_table))

representatives_column <- NA_character_

rep_hits <- which(
  otu_table_colnames_lower %in% c(
    "representative",
    "representatives",
    "rep",
    "centroid"
  )
)

if (length(rep_hits) > 0) {
  representatives_column <- colnames(otu_table)[rep_hits[1]]
}

otu_ids <- otu_table %>%
  transmute(
    otu_id = as.character(otu_id),
    
    otu_representative = if (!is.na(representatives_column) &&
                             representatives_column %in% colnames(otu_table)) {
      as.character(.data[[representatives_column]])
    } else {
      NA_character_
    },
    
    otu_first_token = str_split_fixed(otu_id, "\\s+", 2)[, 1],
    otu_numeric_suffix = str_extract(otu_first_token, "[0-9]+$"),
    
    representative_first_token = str_split_fixed(otu_representative, "\\s+", 2)[, 1],
    representative_numeric_suffix = str_extract(representative_first_token, "[0-9]+$")
  ) %>%
  distinct()

cat("\n============================================================\n")
cat("Final ITS OTU table summary\n")
cat("============================================================\n")
cat("OTU table:", otu_table_file, "\n")
cat("OTUs in study table:", nrow(otu_ids), "\n")
cat("Representative column detected:", representatives_column, "\n")


# ============================================================
# Read taxonomy
# ============================================================

tax_raw <- read_table_auto(taxonomy_file)

colnames(tax_raw) <- trimws(colnames(tax_raw))
colnames(tax_raw)[1] <- sub("^\ufeff", "", colnames(tax_raw)[1])

cat("\n============================================================\n")
cat("Taxonomy file summary\n")
cat("============================================================\n")
cat("Taxonomy file:", taxonomy_file, "\n")
cat("Taxonomy columns:\n")
print(colnames(tax_raw))

if (!tax_id_column %in% colnames(tax_raw)) {
  stop("tax_id_column not found in taxonomy file: ", tax_id_column)
}

if (!kingdom_column %in% colnames(tax_raw)) {
  warning("kingdom_column not found; kingdom will be parsed/fallback as Unclassified: ", kingdom_column)
  kingdom_column <- NA_character_
}

if (!phylum_column %in% colnames(tax_raw)) {
  warning("phylum_column not found; phylum will be parsed/fallback as Unclassified: ", phylum_column)
  phylum_column <- NA_character_
}

if (!is.na(taxonomy_column) && !taxonomy_column %in% colnames(tax_raw)) {
  warning("taxonomy_column was set but not found: ", taxonomy_column)
  taxonomy_column <- NA_character_
}

taxonomy_info <- tax_raw %>%
  transmute(
    tax_id = as.character(.data[[tax_id_column]]),
    taxonomy_record_present = TRUE,
    
    taxonomy_raw = if (!is.na(taxonomy_column) &&
                       taxonomy_column %in% colnames(tax_raw)) {
      as.character(.data[[taxonomy_column]])
    } else {
      NA_character_
    },
    
    kingdom = if (!is.na(kingdom_column) &&
                  kingdom_column %in% colnames(tax_raw)) {
      clean_kingdom(.data[[kingdom_column]])
    } else if (!is.na(taxonomy_column) &&
               taxonomy_column %in% colnames(tax_raw)) {
      extract_kingdom(.data[[taxonomy_column]])
    } else {
      "Unclassified"
    },
    
    phylum = if (!is.na(phylum_column) &&
                 phylum_column %in% colnames(tax_raw)) {
      clean_phylum(.data[[phylum_column]])
    } else if (!is.na(taxonomy_column) &&
               taxonomy_column %in% colnames(tax_raw)) {
      extract_phylum(.data[[taxonomy_column]])
    } else {
      "Unclassified"
    }
  ) %>%
  filter(
    !is.na(tax_id),
    tax_id != ""
  ) %>%
  mutate(
    tax_first_token = str_split_fixed(tax_id, "\\s+", 2)[, 1],
    tax_numeric_suffix = str_extract(tax_first_token, "[0-9]+$"),
    
    raw_mentions_fungi = if_else(
      !is.na(taxonomy_raw) &
        str_detect(taxonomy_raw, regex("\\bFungi\\b", ignore_case = TRUE)),
      TRUE,
      FALSE
    ),
    
    kingdom_classified = !is.na(kingdom) &
      !kingdom %in% c("Unclassified", ""),
    
    is_nonfungal_kingdom = kingdom_classified & kingdom != "Fungi",
    
    is_fungi = case_when(
      kingdom == "Fungi" ~ TRUE,
      is_nonfungal_kingdom ~ FALSE,
      kingdom %in% c("Unclassified", "") &
        phylum %in% known_fungal_phyla ~ TRUE,
      kingdom %in% c("Unclassified", "") &
        raw_mentions_fungi ~ TRUE,
      TRUE ~ NA
    )
  ) %>%
  distinct()

cat("\nKingdom summary in taxonomy file:\n")

taxonomy_info %>%
  count(kingdom, sort = TRUE, name = "n_records") %>%
  print(n = Inf)

cat("\nFungal-status summary in taxonomy file:\n")

taxonomy_info %>%
  mutate(
    fungal_status = case_when(
      is_fungi == TRUE ~ "confirmed_fungi",
      is_fungi == FALSE ~ "confirmed_non_fungi",
      TRUE ~ "unknown_or_unclassified"
    )
  ) %>%
  count(fungal_status, sort = TRUE, name = "n_records") %>%
  print(n = Inf)


# ============================================================
# Map study OTUs to taxonomy
# ============================================================

tax_by_id <- taxonomy_info$tax_id
names(tax_by_id) <- taxonomy_info$tax_id

tax_by_first_token <- taxonomy_info$tax_id
names(tax_by_first_token) <- taxonomy_info$tax_first_token

tax_numeric_unique <- taxonomy_info %>%
  filter(!is.na(tax_numeric_suffix)) %>%
  distinct(tax_numeric_suffix, .keep_all = TRUE)

tax_by_numeric_suffix <- tax_numeric_unique$tax_id
names(tax_by_numeric_suffix) <- tax_numeric_unique$tax_numeric_suffix

otu_tax_mapping <- otu_ids %>%
  mutate(
    tax_id_exact = unname(tax_by_id[otu_id]),
    tax_id_first_token = unname(tax_by_first_token[otu_first_token]),
    tax_id_numeric_suffix = unname(tax_by_numeric_suffix[otu_numeric_suffix]),
    
    tax_id_rep_exact = unname(tax_by_id[otu_representative]),
    tax_id_rep_first_token = unname(tax_by_first_token[representative_first_token]),
    tax_id_rep_numeric_suffix = unname(tax_by_numeric_suffix[representative_numeric_suffix]),
    
    tax_id = case_when(
      !is.na(tax_id_exact) ~ tax_id_exact,
      !is.na(tax_id_first_token) ~ tax_id_first_token,
      !is.na(tax_id_numeric_suffix) ~ tax_id_numeric_suffix,
      !is.na(tax_id_rep_exact) ~ tax_id_rep_exact,
      !is.na(tax_id_rep_first_token) ~ tax_id_rep_first_token,
      !is.na(tax_id_rep_numeric_suffix) ~ tax_id_rep_numeric_suffix,
      TRUE ~ NA_character_
    ),
    
    tax_match_type = case_when(
      !is.na(tax_id_exact) ~ "otu_exact",
      !is.na(tax_id_first_token) ~ "otu_first_token",
      !is.na(tax_id_numeric_suffix) ~ "otu_numeric_suffix",
      !is.na(tax_id_rep_exact) ~ "representative_exact",
      !is.na(tax_id_rep_first_token) ~ "representative_first_token",
      !is.na(tax_id_rep_numeric_suffix) ~ "representative_numeric_suffix",
      TRUE ~ "not_found"
    )
  ) %>%
  select(
    otu_id,
    otu_representative,
    otu_first_token,
    otu_numeric_suffix,
    representative_first_token,
    representative_numeric_suffix,
    tax_id,
    tax_match_type
  ) %>%
  left_join(
    taxonomy_info,
    by = "tax_id"
  ) %>%
  mutate(
    taxonomy_record_present = if_else(
      is.na(taxonomy_record_present),
      FALSE,
      taxonomy_record_present
    ),
    kingdom = if_else(
      is.na(kingdom) | kingdom == "",
      "Unclassified",
      kingdom
    ),
    phylum = if_else(
      is.na(phylum) | phylum == "",
      "Unclassified",
      phylum
    )
  )

cat("\n============================================================\n")
cat("OTU-to-taxonomy mapping summary\n")
cat("============================================================\n")

otu_tax_mapping %>%
  count(tax_match_type, name = "n_otus") %>%
  print(n = Inf)

cat("\nOTU fungal-status summary:\n")

otu_tax_mapping %>%
  mutate(
    fungal_status = case_when(
      is_fungi == TRUE ~ "confirmed_fungi",
      is_fungi == FALSE ~ "confirmed_non_fungi",
      TRUE ~ "unknown_or_unclassified"
    )
  ) %>%
  count(fungal_status, sort = TRUE, name = "n_otus") %>%
  print(n = Inf)


# ============================================================
# Direct match OTUs to tree tips
# ============================================================

tree_tip_by_label <- tree_tips$tree_tip_label
names(tree_tip_by_label) <- tree_tips$tree_tip_label

tree_tip_by_first_token <- tree_tips$tree_tip_label
names(tree_tip_by_first_token) <- tree_tips$tree_tip_first_token

tree_tip_by_numeric_suffix_tbl <- tree_tips %>%
  filter(!is.na(tree_tip_numeric_suffix)) %>%
  distinct(tree_tip_numeric_suffix, .keep_all = TRUE)

tree_tip_numeric_lookup <- tree_tip_by_numeric_suffix_tbl$tree_tip_label
names(tree_tip_numeric_lookup) <- tree_tip_by_numeric_suffix_tbl$tree_tip_numeric_suffix

direct_mapping <- otu_ids %>%
  mutate(
    tree_tip_exact = unname(tree_tip_by_label[otu_id]),
    tree_tip_first_token_match = unname(tree_tip_by_first_token[otu_first_token]),
    tree_tip_numeric_suffix_match = unname(tree_tip_numeric_lookup[otu_numeric_suffix]),
    
    tree_tip_rep_exact = unname(tree_tip_by_label[otu_representative]),
    tree_tip_rep_first_token_match = unname(tree_tip_by_first_token[representative_first_token]),
    tree_tip_rep_numeric_suffix_match = unname(tree_tip_numeric_lookup[representative_numeric_suffix]),
    
    tree_tip_label = case_when(
      !is.na(tree_tip_exact) ~ tree_tip_exact,
      !is.na(tree_tip_first_token_match) ~ tree_tip_first_token_match,
      !is.na(tree_tip_numeric_suffix_match) ~ tree_tip_numeric_suffix_match,
      !is.na(tree_tip_rep_exact) ~ tree_tip_rep_exact,
      !is.na(tree_tip_rep_first_token_match) ~ tree_tip_rep_first_token_match,
      !is.na(tree_tip_rep_numeric_suffix_match) ~ tree_tip_rep_numeric_suffix_match,
      TRUE ~ NA_character_
    ),
    
    match_type = case_when(
      !is.na(tree_tip_exact) ~ "otu_exact",
      !is.na(tree_tip_first_token_match) ~ "otu_first_token",
      !is.na(tree_tip_numeric_suffix_match) ~ "otu_numeric_suffix",
      !is.na(tree_tip_rep_exact) ~ "representative_exact",
      !is.na(tree_tip_rep_first_token_match) ~ "representative_first_token",
      !is.na(tree_tip_rep_numeric_suffix_match) ~ "representative_numeric_suffix",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(tree_tip_label)) %>%
  transmute(
    otu_id = as.character(otu_id),
    otu_representative = as.character(otu_representative),
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
# Optional duplicate mapping
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
  
  dup_raw <- read_table_auto(duplicates_file)
  
  cat("\nDuplicates file columns:\n")
  print(colnames(dup_raw))
  
  cn <- tolower(colnames(dup_raw))
  
  if (!is.na(duplicates_representative_column)) {
    rep_col <- duplicates_representative_column
  } else {
    rep_idx <- which(str_detect(cn, "representative|rep|centroid"))
    rep_col <- if (length(rep_idx) > 0) colnames(dup_raw)[rep_idx[1]] else NA_character_
  }
  
  if (!is.na(duplicates_duplicate_column)) {
    dup_col <- duplicates_duplicate_column
  } else {
    dup_idx <- which(str_detect(cn, "duplicate|member|sequence|query|strain|taxon|otu"))
    dup_col <- if (length(dup_idx) > 0) colnames(dup_raw)[dup_idx[1]] else NA_character_
  }
  
  cat("\nRepresentative column:", rep_col, "\n")
  cat("Duplicate/member column:", dup_col, "\n")
  
  if (!is.na(rep_col) && !is.na(dup_col) &&
      rep_col %in% colnames(dup_raw) &&
      dup_col %in% colnames(dup_raw)) {
    
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
    warning("Could not identify duplicate summary columns. Continuing without duplicate mapping.")
  }
  
} else {
  warning("Duplicates file not found: ", duplicates_file, ". Continuing without duplicate mapping.")
}


# ============================================================
# Collapsed duplicate-to-representative mapping
# ============================================================

collapsed_mapping <- tibble(
  otu_id = character(),
  otu_representative = character(),
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
  
  rep_by_duplicate_id <- duplicate_mapping2$representative_id
  names(rep_by_duplicate_id) <- duplicate_mapping2$duplicate_id
  
  rep_by_duplicate_first_token <- duplicate_mapping2$representative_id
  names(rep_by_duplicate_first_token) <- duplicate_mapping2$duplicate_first_token
  
  dup_numeric_unique <- duplicate_mapping2 %>%
    filter(!is.na(duplicate_numeric_suffix)) %>%
    distinct(duplicate_numeric_suffix, .keep_all = TRUE)
  
  rep_by_duplicate_numeric <- dup_numeric_unique$representative_id
  names(rep_by_duplicate_numeric) <- dup_numeric_unique$duplicate_numeric_suffix
  
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
    filter(!is.na(representative_id)) %>%
    mutate(
      representative_id = as.character(representative_id),
      representative_first_token = str_split_fixed(representative_id, "\\s+", 2)[, 1],
      representative_numeric_suffix = str_extract(representative_first_token, "[0-9]+$")
    )
  
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
    filter(!is.na(tree_tip_label)) %>%
    transmute(
      otu_id = as.character(otu_id),
      otu_representative = as.character(otu_representative),
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
# Combine mapping
# ============================================================

direct_mapping_clean <- direct_mapping %>%
  mutate(
    representative_id = NA_character_,
    otu_id = as.character(otu_id),
    otu_representative = as.character(otu_representative),
    otu_first_token = as.character(otu_first_token),
    otu_numeric_suffix = as.character(otu_numeric_suffix),
    tree_tip_label = as.character(tree_tip_label),
    match_type = as.character(match_type),
    representative_id = as.character(representative_id)
  )

collapsed_mapping_clean <- collapsed_mapping %>%
  mutate(
    otu_id = as.character(otu_id),
    otu_representative = as.character(otu_representative),
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
      match_type == "otu_exact" ~ 1,
      match_type == "otu_first_token" ~ 2,
      match_type == "otu_numeric_suffix" ~ 3,
      match_type == "representative_exact" ~ 4,
      match_type == "representative_first_token" ~ 5,
      match_type == "representative_numeric_suffix" ~ 6,
      match_type == "collapsed_duplicate_to_representative" ~ 7,
      TRUE ~ 99
    )
  ) %>%
  arrange(
    otu_id,
    match_priority,
    tree_tip_label
  )

otu_tip_mapping <- otu_tip_mapping_all %>%
  group_by(otu_id) %>%
  mutate(row_number_within_otu = dplyr::row_number()) %>%
  ungroup() %>%
  filter(row_number_within_otu == 1) %>%
  select(
    -row_number_within_otu,
    -match_priority
  ) %>%
  left_join(
    otu_tax_mapping %>%
      select(
        otu_id,
        kingdom,
        phylum,
        is_fungi,
        taxonomy_record_present,
        tax_id,
        tax_match_type
      ),
    by = "otu_id"
  )

write_tsv(
  otu_tip_mapping,
  otu_tip_mapping_file,
  na = "NA"
)


# ============================================================
# Missing OTUs
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


# ============================================================
# Prune tree to study OTUs
# ============================================================

study_tips_to_keep <- otu_tip_mapping %>%
  pull(tree_tip_label) %>%
  unique() %>%
  as.character()

if (length(study_tips_to_keep) == 0) {
  stop("No mapped tree tips found for study OTUs.")
}

study_pruned_tree <- keep.tip(
  phy = tree,
  tip = study_tips_to_keep,
  trim.internal = TRUE,
  subtree = FALSE,
  root.edge = 0,
  collapse.singles = FALSE
)

write.tree(
  study_pruned_tree,
  file = study_pruned_tree_file
)

cat("\n============================================================\n")
cat("Study-pruned tree summary\n")
cat("============================================================\n")
cat("Study tree file:", study_pruned_tree_file, "\n")
cat("Study OTUs mapped:", nrow(otu_tip_mapping), "\n")
cat("Unique tree tips retained in study tree:", length(study_pruned_tree$tip.label), "\n")
cat("Missing OTUs:", nrow(missing_otus), "\n")


# ============================================================
# Study tree tip taxonomy summary
# ============================================================

study_tip_taxonomy <- otu_tip_mapping %>%
  group_by(tree_tip_label) %>%
  summarise(
    n_mapped_otus = n_distinct(otu_id),
    otu_ids_mapped_to_tip = collapse_unique(otu_id),
    kingdoms_observed = collapse_unique(kingdom),
    phyla_observed = collapse_unique(phylum),
    n_taxonomy_records = sum(taxonomy_record_present == TRUE, na.rm = TRUE),
    n_missing_taxonomy_records = sum(
      is.na(taxonomy_record_present) | taxonomy_record_present == FALSE,
      na.rm = TRUE
    ),
    n_confirmed_fungal_otus = sum(is_fungi == TRUE, na.rm = TRUE),
    n_confirmed_nonfungal_otus = sum(is_fungi == FALSE, na.rm = TRUE),
    n_unknown_fungal_status_otus = sum(is.na(is_fungi), na.rm = TRUE),
    all_mapped_otus_confirmed_fungi = (
      n_mapped_otus > 0 &
        n_missing_taxonomy_records == 0 &
        n_confirmed_nonfungal_otus == 0 &
        n_unknown_fungal_status_otus == 0 &
        n_confirmed_fungal_otus == n_mapped_otus
    ),
    .groups = "drop"
  )

write_tsv(
  study_tip_taxonomy,
  study_tip_taxonomy_file,
  na = "NA"
)

excluded_not_fungi <- study_tip_taxonomy %>%
  filter(!all_mapped_otus_confirmed_fungi)

write_tsv(
  excluded_not_fungi,
  excluded_not_fungi_file,
  na = "NA"
)


# ============================================================
# Prune study tree to confirmed fungi only
# ============================================================

fungal_tips_to_keep <- study_tip_taxonomy %>%
  filter(all_mapped_otus_confirmed_fungi) %>%
  pull(tree_tip_label) %>%
  unique() %>%
  as.character()

if (length(fungal_tips_to_keep) == 0) {
  stop("No confirmed fungal tree tips remain after filtering.")
}

fungi_only_tree <- keep.tip(
  phy = study_pruned_tree,
  tip = fungal_tips_to_keep,
  trim.internal = TRUE,
  subtree = FALSE,
  root.edge = 0,
  collapse.singles = FALSE
)

write.tree(
  fungi_only_tree,
  file = fungi_only_tree_file
)


# ============================================================
# Fungi-only tip phylum consensus
# ============================================================

fungi_only_tip_phylum <- otu_tip_mapping %>%
  filter(tree_tip_label %in% fungi_only_tree$tip.label) %>%
  group_by(tree_tip_label) %>%
  summarise(
    n_mapped_otus = n_distinct(otu_id),
    otu_ids_mapped_to_tip = collapse_unique(otu_id),
    kingdom_observed = collapse_unique(kingdom),
    phyla_observed = collapse_unique(phylum),
    fungal_phylum = consensus_phylum(phylum),
    consensus_fraction = consensus_fraction(phylum),
    .groups = "drop"
  )

write_tsv(
  fungi_only_tip_phylum,
  fungi_only_tip_phylum_table_file,
  na = "NA"
)

cat("\n============================================================\n")
cat("Fungi-only tree summary\n")
cat("============================================================\n")
cat("Fungi-only tree file:", fungi_only_tree_file, "\n")
cat("Fungi-only tree tips:", length(fungi_only_tree$tip.label), "\n")
cat("Excluded non-fungal/unknown study tree tips:", nrow(excluded_not_fungi), "\n")

cat("\nFungal phylum counts in fungi-only tree:\n")

fungi_only_tip_phylum %>%
  count(fungal_phylum, name = "n_tips") %>%
  arrange(desc(n_tips), fungal_phylum) %>%
  print(n = Inf)


# ============================================================
# iTOL color strip for fungi-only phyla
# ============================================================

fungal_phyla <- sort(unique(fungi_only_tip_phylum$fungal_phylum))

phylum_colors <- grDevices::hcl.colors(
  n = length(fungal_phyla),
  palette = "Dark 3"
)

names(phylum_colors) <- fungal_phyla

if ("Unclassified" %in% names(phylum_colors)) {
  phylum_colors["Unclassified"] <- "#969696"
}

if ("Mixed" %in% names(phylum_colors)) {
  phylum_colors["Mixed"] <- "#000000"
}

fungi_only_tip_phylum <- fungi_only_tip_phylum %>%
  mutate(
    color = phylum_colors[fungal_phylum]
  )

legend_shapes <- rep(1, length(fungal_phyla))
legend_colors <- unname(phylum_colors[fungal_phyla])
legend_labels <- fungal_phyla

itol_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tFungal phylum",
  "COLOR\t#000000",
  "STRIP_WIDTH\t35",
  "MARGIN\t10",
  "BORDER_WIDTH\t1",
  "BORDER_COLOR\t#000000",
  "LEGEND_TITLE\tFungal phylum",
  paste0("LEGEND_SHAPES\t", paste(legend_shapes, collapse = "\t")),
  paste0("LEGEND_COLORS\t", paste(legend_colors, collapse = "\t")),
  paste0("LEGEND_LABELS\t", paste(legend_labels, collapse = "\t")),
  "DATA"
)

itol_data_lines <- fungi_only_tip_phylum %>%
  transmute(
    line = paste(
      tree_tip_label,
      color,
      fungal_phylum,
      sep = "\t"
    )
  ) %>%
  pull(line)

writeLines(
  c(itol_lines, itol_data_lines),
  con = itol_fungi_only_phylum_colorstrip_file
)

cat("\niTOL fungi-only phylum color strip written to:\n")
cat(itol_fungi_only_phylum_colorstrip_file, "\n")


# ============================================================
# Quick local PDF
# ============================================================

tip_colors <- fungi_only_tip_phylum$color[
  match(
    fungi_only_tree$tip.label,
    fungi_only_tip_phylum$tree_tip_label
  )
]

pdf(
  fungi_only_tree_pdf_file,
  width = 11,
  height = 14
)

plot(
  fungi_only_tree,
  tip.color = tip_colors,
  cex = 0.35,
  no.margin = FALSE,
  main = "ITS fungi-only tree colored by fungal phylum"
)

legend(
  "topleft",
  legend = names(phylum_colors),
  col = phylum_colors,
  pch = 19,
  cex = 0.7,
  bty = "n"
)

dev.off()


# ============================================================
# Final summary table
# ============================================================

tree_pruning_summary <- tibble(
  level = c(
    "original_TBAS_tree",
    "final_study_ITS_OTU_table",
    "study_OTUs_mapped_to_tree_or_representative",
    "study_OTUs_missing_from_tree_mapping",
    "study_tree_unique_tips_retained",
    "study_tree_tips_confirmed_fungi",
    "study_tree_tips_excluded_not_confirmed_fungi",
    "fungi_only_tree_final_tips"
  ),
  n = c(
    length(tree$tip.label),
    nrow(otu_ids),
    nrow(otu_tip_mapping),
    nrow(missing_otus),
    length(study_pruned_tree$tip.label),
    length(fungal_tips_to_keep),
    nrow(excluded_not_fungi),
    length(fungi_only_tree$tip.label)
  )
)

write_tsv(
  tree_pruning_summary,
  tree_pruning_summary_file,
  na = "NA"
)

cat("\n============================================================\n")
cat("Final pruning totals\n")
cat("============================================================\n")

tree_pruning_summary %>%
  print(n = Inf)

cat("\n============================================================\n")
cat("Files written\n")
cat("============================================================\n")
cat("OTU-to-tree-tip mapping:", otu_tip_mapping_file, "\n")
cat("Missing OTUs:", missing_otus_file, "\n")
cat("Study-pruned tree:", study_pruned_tree_file, "\n")
cat("Study tree tip taxonomy:", study_tip_taxonomy_file, "\n")
cat("Excluded non-fungal/unknown tips:", excluded_not_fungi_file, "\n")
cat("Fungi-only tree:", fungi_only_tree_file, "\n")
cat("Fungi-only tip phylum table:", fungi_only_tip_phylum_table_file, "\n")
cat("iTOL fungi-only phylum color strip:", itol_fungi_only_phylum_colorstrip_file, "\n")
cat("Fungi-only local PDF:", fungi_only_tree_pdf_file, "\n")
cat("Summary table:", tree_pruning_summary_file, "\n")

cat("\nFor iTOL:\n")
cat("1. Upload tree:", fungi_only_tree_file, "\n")
cat("2. Upload annotation:", itol_fungi_only_phylum_colorstrip_file, "\n")

