# --------------------------------------------------------------------------------------
# Script Name : clean_child_anthropometrics.R
# Purpose     : Clean and summarise child anthropometric data with unit checks
# Date created: 15-05-2025
# Last updated: 15-05-2025
# Author      : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# ----------------------------- SETUP ----------------------------------

child_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-05-07_Mikko/child_anthropometrics_raw.gz"
con <- gzfile(child_file, "rt")
child <- read.delim(con, stringsAsFactors = FALSE)
close(con)

verbose <- TRUE  # Toggle for output summaries

# ----------------------------- FUNCTIONS -------------------------------

clean_numeric <- function(x) {
  x <- gsub(",", ".", x)
  x[x %in% c("", ".", "NA")] <- NA
  suppressWarnings(as.numeric(x))
}

convert_weight <- function(w) {
  w_clean <- clean_numeric(w)
  w_clean[!is.finite(w_clean) | w_clean < 2 | w_clean > 200] <- NA
  return(w_clean)
}

convert_height <- function(h) {
  h_clean <- clean_numeric(h)
  guess_unit <- median(h_clean[h_clean > 0], na.rm = TRUE)

  if (!is.na(guess_unit) && guess_unit > 20) {
    h_clean <- h_clean / 100
    if (verbose) cat("↪ Converted height from cm to meters\n")
  } else {
    if (verbose) cat("↪ Height appears to be in meters\n")
  }

  h_clean[!is.finite(h_clean) | h_clean < 0.3 | h_clean > 2.2] <- NA
  return(h_clean)
}

convert_age <- function(a) {
  a_clean <- clean_numeric(a)
  # Remove anything below 0 or above 25 years
  a_clean[!is.finite(a_clean) | a_clean < 0 | a_clean > 25] <- NA
  return(a_clean)
}

summarize_var <- function(x, varname) {
  x_clean <- x[!is.na(x)]
  summary <- data.frame(
    Variable = varname,
    N = length(x_clean),
    Missing = sum(is.na(x)),
    Mean = round(mean(x_clean), 2),
    SD = round(sd(x_clean), 2),
    Min = round(min(x_clean), 2),
    Max = round(max(x_clean), 2)
  )
  print(summary, row.names = FALSE)
}

check_variable_units <- function(x, varname, range, unit_label) {
  x_num <- suppressWarnings(as.numeric(gsub(",", ".", x)))
  valid <- x_num >= range[1] & x_num <= range[2]
  percent_valid <- round(sum(valid, na.rm = TRUE) / sum(!is.na(x_num)) * 100, 1)
  if (is.nan(percent_valid)) percent_valid <- 0
  status <- if (percent_valid > 90) "Likely valid" else "Possibly incorrect"

  cat(sprintf("\nVariable: %s\n", varname))
  cat(sprintf("  Range checked: [%s–%s] (%s)\n", range[1], range[2], unit_label))
  cat(sprintf("  Valid entries: %d of %d (%.1f%%) → %s\n",
              sum(valid, na.rm = TRUE), sum(!is.na(x_num)), percent_valid, status))
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

# ----------------------------- CLEANING AND SUMMARIZING -------------------------------

child$IID <- rownames(child)
final_df <- data.frame(IID = child$IID, sex = child$sex)

for (tp in names(timepoints)) {
  vars <- timepoints[[tp]]

  w_raw <- child[[vars$weight]]
  h_raw <- child[[vars$height]]
  a_raw <- if (!is.null(vars$age)) child[[vars$age]] else NA

  if (verbose) {
    cat("\n--- Checking Timepoint:", tp, "---\n")
    check_variable_units(w_raw, vars$weight, c(2, 200), "kg")
    check_variable_units(h_raw, vars$height, c(0.3, 2.2), "meters")
    if (!is.null(vars$age)) check_variable_units(a_raw, vars$age, c(0, 25), "years")
  }

  weight <- convert_weight(w_raw)
  height <- convert_height(h_raw)
  age <- convert_age(a_raw)

  bmi <- ifelse(!is.na(weight) & !is.na(height) & is.finite(weight) & is.finite(height),
                weight / (height^2), NA)

  if (verbose) {
    summarize_var(weight, paste0("weight_", tp))
    summarize_var(height, paste0("height_", tp))
    summarize_var(age, paste0("age_", tp))
    summarize_var(bmi, paste0("bmi_", tp))
  }

  df_tp <- data.frame(weight, height, bmi, age)
  colnames(df_tp) <- c(
    paste0("weight_", tp),
    paste0("height_", tp),
    paste0("bmi_", tp),
    paste0("age_", tp)
  )

  final_df <- cbind(final_df, df_tp)
}

# ----------------------------- FINAL CLEAN-UP & EXPORT -------------------------------

# Replace any Inf, -Inf, NaN with NA
final_df[!is.finite(as.matrix(final_df))] <- NA

# Export with "." for NA
output_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_anthro_all_timepoints_cleaned.txt"
write.table(final_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

cat("\n Final cleaned dataset saved to:\n", output_file, "\n"

