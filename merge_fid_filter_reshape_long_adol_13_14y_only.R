# --------------------------------------------------------------------------------------
# Script Name : merge_and_reshape_child_13_14.R
# Purpose     : Merge FID, derive BMI from height/weight, and output long BMI/height (13y & 14y only)
# Date created: 25-09-2025
# Author      : Grace Power
# --------------------------------------------------------------------------------------

# ----------------------------- SETUP ----------------------------------

# Inputs
cleaned_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_anthro_all_timepoints_cleaned-adol.txt"
psam_file    <- "/home/grace.power/archive/moba/geno/HDGB-MoBaGenetics/2025.01.30_beta/moba_genotypes_2025.01.30_common.psam"

merged_out   <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_adol_13_14_fid.txt"

# Final long outputs
output_bmi_long    <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/bmi_13_14_long.txt"
output_height_long <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/height_13_14_long.txt"

# ----------------------------- READ FILES -----------------------------

child <- read.delim(cleaned_file, stringsAsFactors = FALSE, check.names = FALSE)

psam  <- read.table(psam_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "")
colnames(psam)[1] <- "FID"  # first column is family ID in .psam
# Keep only FID/IID for merge
psam   <- psam[, c("FID", "IID")]

# ----------------------------- MERGE FID -----------------------------

# Ensure IID exists in both
if (!"IID" %in% names(child)) stop("IID column not found in child data.")
child <- merge(psam, child, by = "IID")  # inner join keeps only IIDs present in .psam

# ----------------------------- KEEP 13y & 14y -----------------------------

keep_cols <- c(
  "FID", "IID", "sex",
  "weight_13y", "height_13y", "age_13y",
  "weight_14y", "height_14y", "age_14y"
)
keep_cols <- keep_cols[keep_cols %in% names(child)]
child <- child[, keep_cols]

# Coerce numeric columns (gracefully)
num_cols <- intersect(c("weight_13y","height_13y","age_13y","weight_14y","height_14y","age_14y"), names(child))
for (cc in num_cols) child[[cc]] <- suppressWarnings(as.numeric(child[[cc]]))

# ----------------------------- DERIVE BMI (robust units) --------------

# Helper to standardize height to meters when computing BMI.
# If a height value > 3, treat as cm and divide by 100.
to_meters <- function(h) ifelse(is.na(h), NA_real_, ifelse(h > 3, h/100, h))

h13_m <- to_meters(child$height_13y)
h14_m <- to_meters(child$height_14y)

child$bmi_13y <- ifelse(!is.na(child$weight_13y) & !is.na(h13_m),
                        child$weight_13y / (h13_m^2), NA_real_)
child$bmi_14y <- ifelse(!is.na(child$weight_14y) & !is.na(h14_m),
                        child$weight_14y / (h14_m^2), NA_real_)

# ----------------------------- LONG FORMAT: BMI -----------------------

bmi_long <- rbind(
  data.frame(
    FID = child$FID, IID = child$IID, sex = child$sex,
    timepoint = "13y",
    age = child$age_13y,
    value = child$bmi_13y,
    stringsAsFactors = FALSE
  ),
  data.frame(
    FID = child$FID, IID = child$IID, sex = child$sex,
    timepoint = "14y",
    age = child$age_14y,
    value = child$bmi_14y,
    stringsAsFactors = FALSE
  )
)

# Remove rows with missing BMI or age
bmi_long <- bmi_long[!is.na(bmi_long$value) & !is.na(bmi_long$age), ]

# ----------------------------- LONG FORMAT: HEIGHT (in cm) -----------

# Convert original heights to centimeters in output
to_cm <- function(h) ifelse(is.na(h), NA_real_, ifelse(h > 3, h, h*100))

height_long <- rbind(
  data.frame(
    FID = child$FID, IID = child$IID, sex = child$sex,
    timepoint = "13y",
    age = child$age_13y,
    value = to_cm(child$height_13y),
    stringsAsFactors = FALSE
  ),
  data.frame(
    FID = child$FID, IID = child$IID, sex = child$sex,
    timepoint = "14y",
    age = child$age_14y,
    value = to_cm(child$height_14y),
    stringsAsFactors = FALSE
  )
)

# Remove rows with missing height or age
height_long <- height_long[!is.na(height_long$value) & !is.na(height_long$age), ]

# ----------------------------- ORDERING & SAVE ------------------------

# Save a compact, merged snapshot
merged_snapshot <- child[, c("FID","IID","sex",
                             "age_13y","height_13y","weight_13y","bmi_13y",
                             "age_14y","height_14y","weight_14y","bmi_14y")]
write.table(merged_snapshot, file = merged_out, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

# Final long files
write.table(bmi_long,    file = output_bmi_long,    sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
write.table(height_long, file = output_height_long, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")


