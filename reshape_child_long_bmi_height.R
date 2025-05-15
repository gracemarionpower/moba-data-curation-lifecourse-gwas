# --------------------------------------------------------------------------------------
# Script Name : reshape_child_long_bmi_height.R
# Purpose     : Create BMI and height long-format files (base R only)
# Date created : 13-05-2025
# Last updated : 13-05-2025
# Author       : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# ----------------------------- SETUP ----------------------------------

input_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_anthro_filtered_fid_no13_14.txt"
output_bmi <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_bmi_long.txt"
output_height <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_height_long.txt"

# ----------------------------- READ DATA -----------------------------

child <- read.delim(input_file, stringsAsFactors = FALSE)

# ----------------------------- LONG FORMAT: BMI -----------------------------

bmi_cols <- grep("^bmi_", names(child), value = TRUE)

bmi_long <- NULL
for (col in bmi_cols) {
  tp <- sub("bmi_", "", col)
  age_col <- paste0("age_", tp)
  
  df <- data.frame(
    FID = child$FID,
    IID = child$IID,
    sex = child$sex,
    timepoint = tp,
    age = if (age_col %in% colnames(child)) child[[age_col]] else NA,
    bmi = child[[col]],
    stringsAsFactors = FALSE
  )
  
  bmi_long <- rbind(bmi_long, df)
}

# ----------------------------- LONG FORMAT: HEIGHT -----------------------------

height_cols <- grep("^height_", names(child), value = TRUE)

height_long <- NULL
for (col in height_cols) {
  tp <- sub("height_", "", col)
  age_col <- paste0("age_", tp)

  df <- data.frame(
    FID = child$FID,
    IID = child$IID,
    sex = child$sex,
    timepoint = tp,
    age = if (age_col %in% colnames(child)) child[[age_col]] else NA,
    height = child[[col]],
    stringsAsFactors = FALSE
  )
  
  height_long <- rbind(height_long, df)
}

# ----------------------------- EXPORT -----------------------------

write.table(bmi_long, output_bmi, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
write.table(height_long, output_height, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

