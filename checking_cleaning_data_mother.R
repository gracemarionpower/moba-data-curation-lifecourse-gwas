# --------------------------------------------------------------------------------------
# Script Name : checking_cleaning_data_mother.R
# Purpose     : To identify the best data sources and clean maternal data before curating time-varying datasets
# Date created: 13-05-2025
# Last updated: 14-05-2025
# Author      : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# NB. Maternal data during pregnancy is excluded as it does not reflect typical BMI values 
# and is subject to multiple confounding factors

# Load data
parent_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-05-07_Mikko/parent.gz"

con <- gzfile(parent_file, "rt")
parent <- read.delim(con, stringsAsFactors = FALSE)
close(con)

# Select variables
mother_data <- parent[, c(
  "mother_sentrix_id",
  "mother_height_self",
  "mother_weight_beginning_self",
  "mother_height",
  "mother_age_15w"
)]
colnames(mother_data)[1] <- "IID"

# Convert all variables to numeric (do this FIRST)
mother_data$mother_height_self <- as.numeric(gsub(",", ".", mother_data$mother_height_self))
mother_data$mother_height <- as.numeric(gsub(",", ".", mother_data$mother_height))
mother_data$mother_weight_beginning_self <- as.numeric(gsub(",", ".", mother_data$mother_weight_beginning_self))
mother_data$mother_age_15w <- as.numeric(gsub(",", ".", mother_data$mother_age_15w))

# Set implausible height values (< 30 cm or > 500 cm) to NA BEFORE converting to meters
mother_data$mother_height_self[
  mother_data$mother_height_self < 30 | mother_data$mother_height_self > 500
] <- NA
mother_data$mother_height[
  mother_data$mother_height < 30 | mother_data$mother_height > 500
] <- NA

# Convert height from cm to meters
mother_data$mother_height_self <- mother_data$mother_height_self / 100
mother_data$mother_height <- mother_data$mother_height / 100

# Subtract 15 weeks (in years) from age
mother_data$age_minus_15w <- mother_data$mother_age_15w - (15 / 52.1775)

# Summary function
summarize_var <- function(x) {
  x_clean <- x[!is.na(x)]
  data.frame(
    Sample_Size = length(x_clean),
    Missing = sum(is.na(x)),
    Mean = round(mean(x_clean), 2),
    Median = round(median(x_clean), 2),
    SD = round(sd(x_clean), 2),
    Min = round(min(x_clean), 2),
    Max = round(max(x_clean), 2)
  )
}

# Summary table for height
height_self_summary <- summarize_var(mother_data$mother_height_self)
height_summary <- summarize_var(mother_data$mother_height)

summary_table_heights <- rbind(
  cbind(Variable = "mother_height_self", height_self_summary),
  cbind(Variable = "mother_height", height_summary)
)
print(summary_table_heights, row.names = FALSE)

# Prepregnancy BMI data
final_mother_data <- mother_data[, c("IID", "mother_weight_beginning_self", "mother_height", "age_minus_15w")]
colnames(final_mother_data) <- c("IID", "weight_prepreg", "height_prepreg", "age_prepreg")
final_mother_data$bmi_prepreg <- final_mother_data$weight_prepreg / (final_mother_data$height_prepreg^2)
final_mother_data <- final_mother_data[, c("IID", "weight_prepreg", "height_prepreg", "bmi_prepreg", "age_prepreg")]

# Complete and partial-case datasets
mother_data_complete <- final_mother_data[complete.cases(final_mother_data), ]
mother_data_partial <- final_mother_data[!is.na(final_mother_data$IID), ]
mother_data_partial[is.na(mother_data_partial)] <- "."

# Save files
write.table(mother_data_complete,
            file = "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_anthro_prepreg_complete.txt",
            sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
write.table(mother_data_partial,
            file = "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_anthro_prepreg_partial.txt",
            sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

# Summary for complete-case
summarize_dataframe <- function(df, vars = NULL) {
  if (is.null(vars)) vars <- names(df)
  summaries <- lapply(vars, function(var) {
    x <- as.numeric(df[[var]])
    x_clean <- x[!is.na(x)]
    data.frame(
      Variable = var,
      Sample_Size = length(x_clean),
      Missing = sum(is.na(df[[var]])),
      Mean = round(mean(x_clean), 2),
      Median = round(median(x_clean), 2),
      SD = round(sd(x_clean), 2),
      Min = round(min(x_clean), 2),
      Max = round(max(x_clean), 2)
    )
  })
  do.call(rbind, summaries)
}

summary_table_complete <- summarize_dataframe(
  mother_data_complete,
  vars = c("weight_prepreg", "height_prepreg", "bmi_prepreg", "age_prepreg")
)
cat("\nSummary of complete-case dataset:\n")
print(summary_table_complete, row.names = FALSE)


# Postnatal maternal data: BMI and age at 14m, 3y, 5y, 8y


# Define timepoints with variable names and child ages (years)
timepoints <- list(
  "14m" = list(weight = "weight_mother_14m", height = "height_mother_14m", child_age = 14 / 12),
  "3y" = list(weight = "mother_weight_3y", height = "mother_height_3y", child_age = 3),
  "5y" = list(weight = "mother_weight_5y", height = "mother_height_5y", child_age = 5),
  "8y" = list(weight = "mother_weight_8y", height = "mother_height_8y", child_age = 8)
)

# Loop through timepoints
for (tp in names(timepoints)) {
  wt_var <- timepoints[[tp]]$weight
  ht_var <- timepoints[[tp]]$height
  child_age <- timepoints[[tp]]$child_age

  if (!all(c(wt_var, ht_var) %in% colnames(parent))) {
    warning(paste("Skipping", tp, "- variables not found"))
    next
  }

  temp_data <- parent[, c("mother_sentrix_id", wt_var, ht_var, "mother_age_15w")]
  colnames(temp_data) <- c("IID", "weight", "height", "mother_age_15w")

  temp_data$weight <- as.numeric(gsub(",", ".", temp_data$weight))
  temp_data$height <- as.numeric(gsub(",", ".", temp_data$height))

  temp_data$height[temp_data$height < 30 | temp_data$height > 500] <- NA
  temp_data$height <- temp_data$height / 100

  # Calculate maternal age at the timepoint
  temp_data[[paste0("age_", tp)]] <- temp_data$mother_age_15w - (15 / 52.1775) + child_age
  temp_data[[paste0("bmi_", tp)]] <- temp_data$weight / (temp_data$height^2)

  output_data <- temp_data[, c("IID", "weight", "height", paste0("bmi_", tp), paste0("age_", tp))]
  colnames(output_data) <- c(paste0(c("IID", "weight_", "height_", "bmi_", "age_"), tp))

  output_complete <- output_data[complete.cases(output_data), ]
  output_partial <- output_data[!is.na(output_data[[1]]), ]
  output_partial[is.na(output_partial)] <- "."

  write.table(output_complete,
              file = paste0("/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_anthro_", tp, "_complete.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
  
  write.table(output_partial,
              file = paste0("/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_anthro_", tp, "_partial.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
}

