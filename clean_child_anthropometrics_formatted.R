# --------------------------------------------------------------------------------------
# Script Name : clean_child_anthropometrics_formatted.R
# Purpose     : Clean and summarise child anthropometric data with unit checks using formatted data
# Date created: 16-05-2025
# Last updated: 16-05-2025
# Author      : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# This includes all child data including <8 years already partially cleaned/formatted

# ----------------------------- SETUP ----------------------------------

child_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-05-07_Mikko/child_anthropometrics.gz"
con <- gzfile(child_file, "rt")
child <- read.delim(con, stringsAsFactors = FALSE)
close(con)

verbose <- TRUE

# ----------------------------- FUNCTIONS -------------------------------

clean_numeric <- function(x) {
  x <- gsub(",", ".", x)
  x[x %in% c("", ".", "NA")] <- NA
  suppressWarnings(as.numeric(x))
}

convert_weight <- function(w) {
  w_clean <- clean_numeric(w)
  w_clean[!is.finite(w_clean) | w_clean < 3 | w_clean > 150] <- NA
  return(w_clean)
}

convert_height <- function(h) {
  h_clean <- clean_numeric(h)
  guess_unit <- median(h_clean[h_clean > 0], na.rm = TRUE)
  if (!is.na(guess_unit) && guess_unit > 20) {
    h_clean <- h_clean / 100
    if (verbose) cat("↪ Converted height from cm to metres\n")
  }
  h_clean[!is.finite(h_clean) | h_clean < 0.45 | h_clean > 2.1] <- NA
  return(h_clean)
}

convert_age <- function(a) {
  a_clean <- clean_numeric(a)
  a_clean <- a_clean / 365.25  # Convert from days to years
  if (verbose) cat("↪ Converted age from days to years\n")
  a_clean[!is.finite(a_clean) | a_clean < 0 | a_clean > 25] <- NA
  return(a_clean)
}

summarise_var <- function(x, varname) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) {
    cat(sprintf("Variable: %s — all values are missing\n", varname))
    return(NULL)
  }
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

# ----------------------------- TIMEPOINT DEFINITIONS -------------------------------

timepoints <- list(
  "3m"    = list(weight = "weight_3m",    height = "length_3m",   age = "age_3m"),
  "6m"    = list(weight = "weight_6m",    height = "length_6m",   age = "age_6m"),
  "8m"    = list(weight = "weight_8m",    height = "length_8m",   age = "age_8m"),
  "1y"    = list(weight = "weight_1y",    height = "length_1y",   age = "age_1y"),
  "16m"   = list(weight = "weight_16m",   height = "length_16m",  age = "age_16m"),
  "2y"    = list(weight = "weight_2y",    height = "length_2y",   age = "age_2y"),
  "3y"    = list(weight = "weight_3y",    height = "length_3y",   age = "age_3y"),
  "5y"    = list(weight = "weight_5y",    height = "length_5y",   age = "age_5y"),
  "7y"    = list(weight = "weight_7y",    height = "length_7y",   age = "age_7y"),
  "8y"    = list(weight = "weight_8y",    height = "length_8y",   age = "age_8y")
)

# ----------------------------- CLEANING LOOP -------------------------------

# Confirm child_sentrix_id exists
if (!"child_sentrix_id" %in% names(child)) stop("Missing 'child_sentrix_id' in dataset")

# Create output starting from IDs
final_df <- data.frame(IID = child$child_sentrix_id, sex = child$sex, stringsAsFactors = FALSE)

for (tp in names(timepoints)) {
  vars <- timepoints[[tp]]

  w_raw <- child[[vars$weight]]
  h_raw <- child[[vars$height]]
  a_raw <- child[[vars$age]]

  if (verbose) cat("\n--- Timepoint:", tp, "---\n")

  weight <- convert_weight(w_raw)
  height <- convert_height(h_raw)
  age <- convert_age(a_raw)

  bmi <- ifelse(!is.na(weight) & !is.na(height), weight / (height^2), NA)
  bmi[bmi < 8 | bmi > 60] <- NA

  if (verbose) {
    summarise_var(weight, paste0("weight_", tp))
    summarise_var(height, paste0("height_", tp))
    summarise_var(age, paste0("age_", tp))
    summarise_var(bmi, paste0("bmi_", tp))
  }

  df_tp <- data.frame(weight, height, bmi, age)
  colnames(df_tp) <- paste0(c("weight_", "height_", "bmi_", "age_"), tp)

  final_df <- cbind(final_df, df_tp)
}

# ----------------------------- EXPORT -------------------------------

output_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/child/child_anthro_cleaned_new.txt"
write.table(final_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

cat("\n Cleaning complete. Data saved to:\n", output_file, "\n")
