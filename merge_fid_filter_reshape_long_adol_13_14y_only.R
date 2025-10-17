# --------------------------------------------------------------------------------------
# Script Name : clean_merge_long_child_13_14_HDGB.R
# Purpose     : (1) Clean raw phenotypes to HDGB-compatible format
#               (2) Merge with FID from .psam
#               (3) Output both snapshot and long-format BMI/height (13y & 14c)
# Date created: 17-10-2025
# Author      : Grace Power
# --------------------------------------------------------------------------------------

# Packages
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
library(data.table)

# ----------------------------- PATHS ---------------------------------

# Input raw phenotype file
raw_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-09-04_HDGB_compatible/child_anthropometrics_raw.gz"

# PSAM with FID/IID
psam_file <- "/home/grace.power/archive/moba/geno/HDGB-MoBaGenetics/2025.01.30_beta/moba_genotypes_2025.01.30_common.psam"

# Outputs
cleaned_file        <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_anthro_all_timepoints_cleaned-adol-HDGB_compatible.txt"
merged_out          <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_adol_13_14_fid.txt"
output_bmi_long     <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/bmi_13_14_long.txt"
output_height_long  <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/height_13_14_long.txt"

dir.create(dirname(cleaned_file), recursive = TRUE, showWarnings = FALSE)

# ----------------------------- STEP 1: CLEAN RAW ---------------------

# Columns to keep from raw
cols_to_keep <- c(
  "child_sentrix_id",
  "height_13", "weight_13", "age_answering_q_13",
  "height_14c", "weight_14c", "age_answering_q_14c"
)

# Read only needed cols
dt <- fread(raw_file, select = cols_to_keep)

# Helper to extract numbers only
num_only <- function(x) as.numeric(gsub("[^0-9.]+", "", x))

# Clean height/weight units
dt[, height_13_cm := num_only(height_13)]
dt[, weight_13_kg := num_only(weight_13)]
dt[, height_14c_cm := num_only(height_14c)]
dt[, weight_14c_kg := num_only(weight_14c)]

# Compute BMI
dt[, bmi_13  := ifelse(!is.na(height_13_cm)  & height_13_cm  > 0,
                       weight_13_kg  / ((height_13_cm/100)^2), NA_real_)]
dt[, bmi_14c := ifelse(!is.na(height_14c_cm) & height_14c_cm > 0,
                       weight_14c_kg / ((height_14c_cm/100)^2), NA_real_)]

# Age conversions
dt[, age_answering_q_13_years  := as.numeric(age_answering_q_13) / 12]            # months -> years
dt[, age_answering_q_14c_months := as.numeric(age_answering_q_14c) / 365.25 * 12] # days -> months

# Final cleaned dataset
cleaned <- dt[, .(
  child_sentrix_id,
  height_13_cm,  weight_13_kg,  bmi_13,  age_answering_q_13_years,
  height_14c_cm, weight_14c_kg, bmi_14c, age_answering_q_14c_months
)]

# Save cleaned file (NA as ".")
fwrite(cleaned, cleaned_file, sep = "\t", na = ".", quote = FALSE)

# ----------------------------- STEP 2: MERGE FID ---------------------

# Reload cleaned file to be sure "." are NA internally
child <- fread(cleaned_file, na.strings = ".", sep = "\t", check.names = FALSE)

# Read PSAM FID/IID
psam <- fread(psam_file, header = TRUE)
setnames(psam, old = names(psam)[1:2], new = c("FID","IID"))
psam <- psam[, .(FID, IID)]

# Add IID from child_sentrix_id
child[, IID := as.character(child_sentrix_id)]

# Merge
child <- merge(psam, child, by = "IID", all = FALSE)

# ----------------------------- STEP 3: SNAPSHOT ----------------------

snapshot_cols <- intersect(
  c("FID","IID",
    "age_answering_q_13_years","height_13_cm","weight_13_kg","bmi_13",
    "age_answering_q_14c_months","height_14c_cm","weight_14c_kg","bmi_14c"),
  names(child)
)
merged_snapshot <- child[, ..snapshot_cols]
fwrite(merged_snapshot, merged_out, sep = "\t", na = ".", quote = FALSE)

# ----------------------------- STEP 4: LONG FILES --------------------

# BMI long
bmi_13 <- child[, .(FID, IID, value = bmi_13, age = age_answering_q_13_years)]
bmi_14 <- child[, .(FID, IID, value = bmi_14c, age = age_answering_q_14c_months)]
bmi_long <- rbind(bmi_13, bmi_14, use.names = TRUE)
bmi_long <- bmi_long[!is.na(value) & !is.na(age)]
fwrite(bmi_long, output_bmi_long, sep = "\t", na = ".", quote = FALSE)

# Height long
height_13 <- child[, .(FID, IID, value = height_13_cm, age = age_answering_q_13_years)]
height_14 <- child[, .(FID, IID, value = height_14c_cm, age = age_answering_q_14c_months)]
height_long <- rbind(height_13, height_14, use.names = TRUE)
height_long <- height_long[!is.na(value) & !is.na(age)]
fwrite(height_long, output_height_long, sep = "\t", na = ".", quote = FALSE)

# ----------------------------- DONE ---------------------------------

cat("Pipeline complete.\n",
    "- Cleaned file:", cleaned_file, "\n",
    "- Snapshot with FID/IID:", merged_out, "\n",
    "- BMI long:", output_bmi_long, "\n",
    "- Height long:", output_height_long, "\n")


