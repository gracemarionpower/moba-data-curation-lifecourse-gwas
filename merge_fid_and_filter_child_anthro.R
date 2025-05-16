# --------------------------------------------------------------------------------------
# Script Name : merge_fid_and_filter_child_anthro.R
# Purpose     : Merge FID into cleaned child anthropometric data and remove 18/19y data
# Date created: 15-05-2025
# Last created: 15-05-2025
# Author      : Grace M. Power
# --------------------------------------------------------------------------------------

# ----------------------------- SETUP ----------------------------------

# Input files
cleaned_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_anthro_cleaned_new.txt"
psam_file    <- "/home/grace.power/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.psam"
# Output
output_file  <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_8_anthro_filtered_fid.txt"

# ----------------------------- READ FILES -----------------------------

# Read cleaned child anthropometric data
child <- read.delim(cleaned_file, stringsAsFactors = FALSE)

# Read .psam to get FID-IID mapping
psam <- read.table(psam_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "")
colnames(psam)[1] <- "FID" 

# ----------------------------- MERGE FID -----------------------------

# Merge to add FID based on IID
child <- merge(psam[, c("FID", "IID")], child, by = "IID")

# ----------------------------- DROP 18Y & 19Y -----------------------------

# Identify and remove columns related to 18y and 19y
drop_cols <- grep("_(18y|19y)$", colnames(child), value = TRUE)
child <- child[, !colnames(child) %in% drop_cols]

# ----------------------------- REORDER COLUMNS -----------------------------

# Move FID and IID to the front
other_cols <- setdiff(colnames(child), c("FID", "IID"))
child <- child[, c("FID", "IID", other_cols)]

# ----------------------------- EXPORT -----------------------------

write.table(child, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

