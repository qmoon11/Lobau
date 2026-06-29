library(dada2)
library(Biostrings)
library(R.utils)

# Databases:
#   1) PR2 version_5.1.1_SSU
#   2) SILVA Fungi v138

setwd("/scratch/tyjames_root/tyjames0/qmoon/Lobau")

fasta_file <- "Lobau_18S_16S_ITS/JMF-2301-09__all__rRNA_SSU_Euk_D__JMFR_MSHI_KP237/OTU99_centroids.fna"

pr2_ref <- "databases/pr2_version_5.1.1_SSU_dada2.fasta"
silva_ref <- "databases/SILVA_SSUfungi_nr99_v138_2_toGenus_trainset.fasta"

sequences_18s <- readDNAStringSet(fasta_file)
headers <- names(sequences_18s)

tax_levels <- c(
  "Domain",
  "Supergroup",
  "Division",
  "Subdivision",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

# -----------------------------
# Run PR2
# -----------------------------

taxa_pr2 <- assignTaxonomy(
  sequences_18s,
  pr2_ref,
  multithread = TRUE,
  minBoot = 80,
  outputBootstraps = TRUE,
  taxLevels = tax_levels,
  verbose = TRUE
)

pr2_tax_df <- as.data.frame(taxa_pr2$tax)
pr2_boot_df <- as.data.frame(taxa_pr2$boot)

names(pr2_boot_df) <- paste0(names(pr2_boot_df), "_bootstrap")

pr2_tax_df$Sequence_ID <- rownames(pr2_tax_df)
pr2_boot_df$Sequence_ID <- rownames(pr2_boot_df)

pr2_df <- merge(
  pr2_tax_df,
  pr2_boot_df,
  by = "Sequence_ID",
  all = TRUE
)

pr2_df$ASV <- headers

saveRDS(pr2_df, "18s_taxonomic_pr2_bootstrap.rds")
write.csv(pr2_df, "18s_taxonomic_pr2_bootstrap.csv", row.names = FALSE)


# -----------------------------
# Run SILVA
# -----------------------------

taxa_silva <- assignTaxonomy(
  sequences_18s,
  silva_ref,
  multithread = TRUE,
  minBoot = 80,
  outputBootstraps = TRUE,
  verbose = TRUE
)

silva_tax_df <- as.data.frame(taxa_silva$tax)
silva_boot_df <- as.data.frame(taxa_silva$boot)

names(silva_boot_df) <- paste0(names(silva_boot_df), "_bootstrap")

silva_tax_df$Sequence_ID <- rownames(silva_tax_df)
silva_boot_df$Sequence_ID <- rownames(silva_boot_df)

silva_df <- merge(
  silva_tax_df,
  silva_boot_df,
  by = "Sequence_ID",
  all = TRUE
)

silva_df$ASV <- headers

saveRDS(silva_df, "18s_taxonomic_silva_bootstrap.rds")
write.csv(silva_df, "18s_taxonomic_silva_bootstrap.csv", row.names = FALSE)
