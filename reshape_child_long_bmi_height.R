# --------------------------------------------------------------------------------------
# Script Name : reshape_child_long_bmi_height.R
# Purpose     : Create long-format BMI and height files
# Date created : 16-05-2025
# Last updated : 16-05-2025
# Author       : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# ----------------------------- SETUP ----------------------------------

input_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_8_anthro_filtered_fid.txt"

# Final output
output_bmi_final <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/bmi.txt"               # complete case only
output_height_final <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/height.txt"         # allows missing

# Full long-format outputs (all columns)
output_bmi_full <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/bmi_8_long.txt"
output_height_full <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/height_8_long.txt"

# ----------------------------- READ DATA -----------------------------

child <- read.delim(input_file, stringsAsFactors = FALSE)

# ----------------------------- LONG FORMAT: BMI -----------------------------

bmi_cols <- grep("^bmi_", names(child), value = TRUE)
bmi_long <- do.call(rbind, lapply(bmi_cols, function(col) {
  tp <- sub("bmi_", "", col)
  age_col <- paste0("age_", tp)
  data.frame(
    FID = child$FID,
    IID = child$IID,
    sex = child$sex,
    timepoint = tp,
    age = if (age_col %in% names(child)) child[[age_col]] else NA,
    bmi = child[[col]],
    stringsAsFactors = FALSE
  )
}))

# ----------------------------- LONG FORMAT: HEIGHT -----------------------------

height_cols <- grep("^height_", names(child), value = TRUE)
height_long <- do.call(rbind, lapply(height_cols, function(col) {
  tp <- sub("height_", "", col)
  age_col <- paste0("age_", tp)
  data.frame(
    FID = child$FID,
    IID = child$IID,
    sex = child$sex,
    timepoint = tp,
    age = if (age_col %in% names(child)) child[[age_col]] else NA,
    height = child[[col]],
    stringsAsFactors = FALSE
  )
}))

# ----------------------------- FINAL VERSIONS (narrowed columns) -----------------------------

# BMI: narrow to key columns
bmi_narrow <- bmi_long[, c("FID", "IID", "age", "bmi")]
colnames(bmi_narrow)[colnames(bmi_narrow) == "bmi"] <- "value"

# HEIGHT: narrow to key columns
height_narrow <- height_long[, c("FID", "IID", "age", "height")]
colnames(height_narrow)[colnames(height_narrow) == "height"] <- "value"

# ----------------------------- HANDLE MISSINGNESS STRATEGY -----------------------------

# BMI: keep only complete cases
bmi_final <- bmi_narrow[!is.na(bmi_narrow$value) & bmi_narrow$value != ".", ]

# Height: keep only complete cases
height_final <- height_narrow[!is.na(height_narrow$value) & height_narrow$value != ".", ]


# ----------------------------- EXPORT -----------------------------

# Final narrow outputs
write.table(bmi_final, output_bmi_final, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
write.table(height_final, output_height_final, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

# Full long-format outputs with all metadata columns
write.table(bmi_long, output_bmi_full, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
write.table(height_long, output_height_full, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

