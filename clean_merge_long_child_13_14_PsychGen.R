# --------------------------------------------------------------------------------------
# Script Name : clean_merge_long_child_13_14_PsychGen.R
# Purpose     : (1) Clean raw phenotypes to PsychGen-compatible format
#               (2) Merge with FID from PsychGen .fam
#               (3) Output both snapshot and long-format BMI/height (13y & 14c)
# Date created: 02-06-2026
# Author      : Grace Power
# --------------------------------------------------------------------------------------

# Packages
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
library(data.table)

# ----------------------------- PATHS ---------------------------------

raw_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-06-26_PsychGen_compatible/child_anthropometrics_raw.gz"

fam_file <- "/home/grace.power/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam"

out_dir <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/adol/psychgen_link"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cleaned_file       <- file.path(out_dir, "child_anthro_all_timepoints_cleaned-adol-PsychGen_compatible.txt")
merged_out         <- file.path(out_dir, "child_adol_13_14_fid.txt")
output_bmi_long    <- file.path(out_dir, "bmi.txt")
output_height_long <- file.path(out_dir, "height.txt")

# ----------------------------- STEP 1: CLEAN RAW ---------------------

cols_to_keep <- c(
  "child_sentrix_id",
  "height_13", "weight_13", "age_answering_q_13",
  "height_14c", "weight_14c", "age_answering_q_14c"
)

dt <- fread(raw_file, select = cols_to_keep)

num_only <- function(x) {
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(gsub("[^0-9.]+", "", x)))
}

dt[, height_13_cm := num_only(height_13)]
dt[, weight_13_kg := num_only(weight_13)]
dt[, height_14c_cm := num_only(height_14c)]
dt[, weight_14c_kg := num_only(weight_14c)]

dt[, age_13_years := as.numeric(age_answering_q_13) / 12]
dt[, age_14c_years := as.numeric(age_answering_q_14c) / 365.25]

dt[, bmi_13 := weight_13_kg / ((height_13_cm / 100)^2)]
dt[, bmi_14c := weight_14c_kg / ((height_14c_cm / 100)^2)]

# Basic sanity filters
dt[height_13_cm < 100 | height_13_cm > 220, height_13_cm := NA_real_]
dt[height_14c_cm < 100 | height_14c_cm > 220, height_14c_cm := NA_real_]

dt[weight_13_kg < 20 | weight_13_kg > 200, weight_13_kg := NA_real_]
dt[weight_14c_kg < 20 | weight_14c_kg > 200, weight_14c_kg := NA_real_]

dt[bmi_13 < 8 | bmi_13 > 60, bmi_13 := NA_real_]
dt[bmi_14c < 8 | bmi_14c > 60, bmi_14c := NA_real_]

dt[age_13_years < 10 | age_13_years > 18, age_13_years := NA_real_]
dt[age_14c_years < 10 | age_14c_years > 18, age_14c_years := NA_real_]

cleaned <- dt[, .(
  child_sentrix_id,
  height_13_cm, weight_13_kg, bmi_13, age_13_years,
  height_14c_cm, weight_14c_kg, bmi_14c, age_14c_years
)]

fwrite(cleaned, cleaned_file, sep = "\t", na = ".", quote = FALSE)

# ----------------------------- STEP 2: MERGE FID ---------------------

child <- fread(cleaned_file, na.strings = ".", sep = "\t", check.names = FALSE)

fam <- fread(fam_file, header = FALSE)
setnames(fam, c("FID", "IID", "father", "mother", "sex", "phenotype"))
fam <- fam[, .(FID, IID)]

child[, IID := as.character(child_sentrix_id)]

child <- merge(fam, child, by = "IID", all = FALSE)

cat("Rows after PsychGen .fam merge:", nrow(child), "\n")
cat("Unique merged IIDs:", uniqueN(child$IID), "\n")
cat("Duplicated merged IIDs:", sum(duplicated(child$IID)), "\n")

# ----------------------------- STEP 3: SNAPSHOT ----------------------

snapshot <- child[, .(
  FID, IID,
  age_13_years, height_13_cm, weight_13_kg, bmi_13,
  age_14c_years, height_14c_cm, weight_14c_kg, bmi_14c
)]

fwrite(snapshot, merged_out, sep = "\t", na = ".", quote = FALSE)

# ----------------------------- STEP 4: LONG FILES --------------------

bmi_13 <- child[, .(FID, IID, value = bmi_13, age = age_13_years)]
bmi_14 <- child[, .(FID, IID, value = bmi_14c, age = age_14c_years)]

bmi_long <- rbind(bmi_13, bmi_14, use.names = TRUE)
bmi_long <- bmi_long[!is.na(value) & !is.na(age)]

fwrite(bmi_long, output_bmi_long, sep = "\t", na = ".", quote = FALSE)

height_13 <- child[, .(FID, IID, value = height_13_cm, age = age_13_years)]
height_14 <- child[, .(FID, IID, value = height_14c_cm, age = age_14c_years)]

height_long <- rbind(height_13, height_14, use.names = TRUE)
height_long <- height_long[!is.na(value) & !is.na(age)]

fwrite(height_long, output_height_long, sep = "\t", na = ".", quote = FALSE)

# ----------------------------- DONE ---------------------------------

cat("Pipeline complete.\n",
    "- Cleaned file:", cleaned_file, "\n",
    "- Snapshot with FID/IID:", merged_out, "\n",
    "- BMI long:", output_bmi_long, "\n",
    "- Height long:", output_height_long, "\n")
