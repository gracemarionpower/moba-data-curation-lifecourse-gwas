
# --------------------------------------------------------------------------------------
# Script Name : checking_cleaning_data.R
# Purpose     : To identify the best data sources and clean parental data before curating time-varying datasets
# Date created: 13-05-2025
# Last updated: 13-05-2025
# Author      : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# NB. Maternal data during pregnancy is excluded as it does not reflect typical BMI values 
# and is subject to multiple confounding factors

# ---------------------------
# Load data
# ---------------------------

parent_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-05-07_Mikko/parent.gz"
output_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_data_check.txt"

con <- gzfile(parent_file, "rt")
parent <- read.delim(con, stringsAsFactors = FALSE)
close(con)

# ---------------------------
# Select variables
# ---------------------------

mother_data <- parent[, c(
  "mother_sentrix_id",
  "mother_height_self",
  "mother_weight_beginning_self",
  "mother_height",
  "mother_age_15w"
)]
colnames(mother_data)[1] <- "IID"

# ---------------------------
# Convert all variables to numeric (do this FIRST)
# ---------------------------

mother_data$mother_height_self <- as.numeric(gsub(",", ".", mother_data$mother_height_self))
mother_data$mother_height <- as.numeric(gsub(",", ".", mother_data$mother_height))
mother_data$mother_weight_beginning_self <- as.numeric(gsub(",", ".", mother_data$mother_weight_beginning_self))
mother_data$mother_age_15w <- as.numeric(gsub(",", ".", mother_data$mother_age_15w))

# ---------------------------
# Set implausible values (< 30 cm or > 500 cm) to NA BEFORE converting to meters
# ---------------------------

mother_data$mother_height_self[
  mother_data$mother_height_self < 30 | mother_data$mother_height_self > 500
] <- NA

mother_data$mother_height[
  mother_data$mother_height < 30 | mother_data$mother_height > 500
] <- NA

# Convert from cm to meters
mother_data$mother_height_self <- mother_data$mother_height_self / 100
mother_data$mother_height <- mother_data$mother_height / 100

# Subtract 15 weeks (in years)
mother_data$age_minus_15w <- mother_data$mother_age_15w - (15 / 52.1775)

# ---------------------------
# Summary table for height variables
# ---------------------------

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

height_self_summary <- summarize_var(mother_data$mother_height_self)
height_summary <- summarize_var(mother_data$mother_height)

summary_table <- rbind(
  cbind(Variable = "mother_height_self", height_self_summary),
  cbind(Variable = "mother_height", height_summary)
)

print(summary_table, row.names = FALSE)

# Create final cleaned dataset

final_mother_data <- mother_data[, c("IID", "mother_weight_beginning_self", "mother_height")]

# Rename columns to final names
colnames(final_mother_data) <- c("IID", "weight_prepreg", "height_prepreg")

# Create final cleaned dataset
final_mother_data <- mother_data[, c("IID", "mother_weight_beginning_self", "mother_height")]
colnames(final_mother_data) <- c("IID", "weight_prepreg", "height_prepreg")

# Calculate BMI
final_mother_data$bmi_prepreg <- final_mother_data$weight_prepreg / (final_mother_data$height_prepreg^2)

# Preview
head(final_mother_data)

# Save cleaned file
write.table(final_mother_data,
            file = "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_cleaned_anthro_prepreg.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


