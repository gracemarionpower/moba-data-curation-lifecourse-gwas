# --------------------------------------------------------------------------------------
# Script Name : child_data_wide.R
# Purpose     : Generate timepoint-specific child datasets
# Date created: 15-05-2025
# Last updated: 15-05-2025
# Author      : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# ----------------------------- SETUP ----------------------------------

child_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-05-07_Mikko/child_anthropometrics_raw.gz"
child <- read.delim(gzfile(child_file, "rt"), stringsAsFactors = FALSE)

# ----------------------------- FUNCTIONS -------------------------------

clean_numeric <- function(x) as.numeric(gsub(",", ".", x))

convert_height <- function(h) {
  h_clean <- clean_numeric(h)
  h_clean[h_clean < 0.3 | h_clean > 2.2] <- NA
  return(h_clean)
}

convert_weight <- function(w) {
  w_clean <- clean_numeric(w)
  w_clean[w_clean < 2 | w_clean > 200] <- NA
  return(w_clean)
}

# ----------------------------- TIMEPOINT DEFINITIONS -------------------------------

timepoints <- list(
  "6m"  = list(weight = "weight_6m",  height = "length_6m",    age = "age_6m"),
  "1y"  = list(weight = "weight_1y",  height = "length_1y",    age = "age_1y"),
  "2y"  = list(weight = "weight_2y",  height = "length_2y",    age = "age_2y"),
  "3y"  = list(weight = "weight_3y",  height = "length_3y",    age = "age_3y"),
  "5y"  = list(weight = "weight_5y",  height = "length_5y",    age = "age_5y"),
  "7y"  = list(weight = "weight_7y",  height = "height_7y_m",  age = "age_7y"),
  "8y"  = list(weight = "weight_8y",  height = "length_8y",    age = "age_8y"),
  "13y" = list(weight = "weight_13", height = "height_13",    age = "age_answering_q_13"),
  "14y" = list(weight = "weight_14c",height = "height_14c",   age = "age_answering_q_14c"),
  "16m" = list(weight = "weight_16m",height = "length_16m",   age = "age_16m"),
  "18y" = list(weight = "weight_18", height = "height_18",    age = "age_answering_q_18"),
  "19y" = list(weight = "weight_19", height = "height_19",    age = "age_answering_q_19")
)

# ----------------------------- CLEANING LOOP -------------------------------

child$IID <- rownames(child)
final_df <- data.frame(IID = child$IID, sex = child$sex)

for (tp in names(timepoints)) {
  vars <- timepoints[[tp]]

  weight <- convert_weight(child[[vars$weight]])
  height <- convert_height(child[[vars$height]])
  age <- if (!is.null(vars$age)) clean_numeric(child[[vars$age]]) else NA
  bmi <- weight / (height^2)

  df_tp <- data.frame(
    weight = weight,
    height = height,
    bmi = bmi,
    age = age
  )

  colnames(df_tp) <- c(
    paste0("weight_", tp),
    paste0("height_", tp),
    paste0("bmi_", tp),
    paste0("age_", tp)
  )

  final_df <- cbind(final_df, df_tp)
}

# ----------------------------- COMPLETE AND PARTIAL -------------------------------

# Complete case: all values (except IID and sex) must be present
cols_to_check <- setdiff(names(final_df), c("IID", "sex"))
child_complete <- final_df[complete.cases(final_df[, cols_to_check]), ]

# Partial: keep all rows with IID, replace NAs with "."
child_partial <- final_df
child_partial[is.na(child_partial)] <- "."

# ----------------------------- EXPORT ----------------------------------

output_dir <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child"
write.table(child_complete,
            file = file.path(output_dir, "child_anthro_all_timepoints_complete.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

write.table(child_partial,
            file = file.path(output_dir, "child_anthro_all_timepoints_partial.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

cat("Complete and partial child datasets saved to:\n", output_dir, "\n")
