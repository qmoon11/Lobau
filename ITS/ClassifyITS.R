# Load the package
library(ClassifyITS)


ITS_taxonomy <- ITS_assignment(
  blast_file       = "ITS/megablast_ITS2_ClassifyITS.tsv",
  rep_fasta        = "ITS/OTU98_centroids_ITS.FASTA",
  cutoff_fraction  = 0.3,        # Optional: fraction for alignment length QC
  n_cutoff         = 1,          # Optional: percent N cutoff
  outdir           = "ClassifyITSoutputs",  # Optional: output directory (writes CSV/PDF when provided)
  verbose          = FALSE       # Optional: print progress messages
)



#Manual Inspection:
# Read in BLAST results and assignments
blast <- read.table("ITS/megablast_ITS2_ClassifyITS.tsv", sep = "\t", header = TRUE) #user provided BLAST results
assignments <- read.csv("ClassifyITSoutputs/initial_assignments.csv", stringsAsFactors = FALSE) #preliminary taxonomy assignments from ClassifyITS


# Count total fungal OTUs
total_fungal_otus <- assignments[
  assignments$kingdom == "Fungi",
  "qseqid"
]


