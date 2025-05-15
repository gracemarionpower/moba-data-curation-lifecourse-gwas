
# --------------------------------------------------------------------------------------
# Script Name: child_anthro_filtered_fid_no13_14.R
# Purpose      : Drop 13/14-year variables from filtered child dataset with FID/IID
# Date created : 15-05-2025
# Last updated : 15-05-2025
# Author       : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# ----------------------------- SETUP ----------------------------------

# Input file (previous output with FID + all timepoints)
input_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_anthro_filtered_fid.txt"

# Output file (new version without 13y/14y)
output_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_anthro_filtered_fid_no13_14.txt"

# ----------------------------- READ DATA -----------------------------

child <- read.delim(input_file, stringsAsFactors = FALSE)

# ----------------------------- DROP 13/14Y VARIABLES -----------------------------

# Identify columns containing "13" or "14" (e.g., weight_13y, age_14c)
drop_cols <- grep("13|14", colnames(child), value = TRUE)
cat("Dropping", length(drop_cols), "columns matching '13' or '14'.\n")

# Drop those columns
child <- child[, !colnames(child) %in% drop_cols]

# ----------------------------- EXPORT -----------------------------

write.table(child, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

cat("\n Dataset without 13/14y timepoints saved to:\n", output_file, "\n")
