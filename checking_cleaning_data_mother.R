# --------------------------------------------------------------------------------------
# Script Name : checking_cleaning_data_mother.R
# Purpose     : Clean maternal data and generate timepoint-specific datasets
# Date created: 13-05-2025
# Last updated: 15-05-2025
# Author      : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# NB: Maternal data during pregnancy is excluded as physiological changes during this period may introduce systematic differences

# ----------------------------- SETUP ----------------------------------

# Load data
parent_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-05-07_Mikko/parent.gz"
parent <- read.delim(gzfile(parent_file, "rt"), stringsAsFactors = FALSE)

# ----------------------------- FUNCTIONS -------------------------------

clean_numeric <- function(x) as.numeric(gsub(",", ".", x))

convert_height <- function(h_cm) {
  h <- clean_numeric(h_cm)
  h[h < 30 | h > 500] <- NA
  return(h / 100)
}

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

summarise_dataframe <- function(df, vars = NULL) {
  if (is.null(vars)) vars <- names(df)
  do.call(rbind, lapply(vars, function(var) {
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
  }))
}

# ----------------------------- VARIABLE DEFINITIONS -------------------------------

# Core identifiers (always needed)
core_vars <- c("mother_sentrix_id", "mother_age_15w")

# Prepregnancy-specific variables
prepreg_vars <- c("mother_weight_beginning_self", "mother_height", "mother_height_self")

# Postnatal timepoints 
timepoints <- list(
  "3y" = list(weight = "mother_weight_3y", height = "mother_height_3y", age_offset = 3+0.75), # adding 9 months of pregnancy == 0.75
  "5y" = list(weight = "mother_weight_5y", height = "mother_height_5y", age_offset = 5+0.75), # adding 9 months of pregnancy == 0.75
  "8y" = list(weight = "mother_weight_8y", height = "mother_height_8y", age_offset = 8+0.75) # adding 9 months of pregnancy == 0.75
)

# Combine all variables needed
all_vars <- unique(c(core_vars, prepreg_vars, unlist(lapply(timepoints, function(x) c(x$weight, x$height)))))
mother_data <- parent[, all_vars]
colnames(mother_data)[colnames(mother_data) == "mother_sentrix_id"] <- "IID"

# ----------------------------- CLEAN PREPREG AGE DATA -------------------------------

mother_data$mother_age_15w <- clean_numeric(mother_data$mother_age_15w)
mother_data$age_prepreg <- mother_data$mother_age_15w - (15 / 52.1775)  # adjust for 15 weeks

# ----------------------------- CLEAN & SUMMARISE HEIGHT VARIANTS -------------------------------

mother_data$mother_height_self <- convert_height(mother_data$mother_height_self)
mother_data$mother_height      <- convert_height(mother_data$mother_height)

height_self_summary <- summarize_var(mother_data$mother_height_self)
height_summary      <- summarize_var(mother_data$mother_height)

summary_table_heights <- rbind(
  cbind(Variable = "mother_height_self", height_self_summary),
  cbind(Variable = "mother_height", height_summary)
)
cat("\nHeight variable summary:\n")
print(summary_table_heights, row.names = FALSE)

# ----------------------------- PREPREGNANCY BMI DATA -------------------------------

mother_data$mother_weight_beginning_self <- clean_numeric(mother_data$mother_weight_beginning_self)

prepreg <- data.frame(
  IID = mother_data$IID,
  weight_prepreg = mother_data$mother_weight_beginning_self,
  height_prepreg = mother_data$mother_height,
  age_prepreg = mother_data$age_prepreg
)
prepreg$bmi_prepreg <- prepreg$weight_prepreg / (prepreg$height_prepreg^2)
prepreg <- prepreg[, c("IID", "weight_prepreg", "height_prepreg", "bmi_prepreg", "age_prepreg")]

# Save files
prepreg_complete <- prepreg[complete.cases(prepreg), ]
prepreg_partial  <- prepreg[!is.na(prepreg$IID), ]
prepreg_partial[is.na(prepreg_partial)] <- "."

# Add sex and rename columns (except IID)
prepreg_complete$sex <- 2
prepreg_partial$sex <- 2
colnames(prepreg_complete)[colnames(prepreg_complete) != "IID"] <- paste0("mum_", colnames(prepreg_complete)[colnames(prepreg_complete) != "IID"])
colnames(prepreg_partial)[colnames(prepreg_partial) != "IID"] <- paste0("mum_", colnames(prepreg_partial)[colnames(prepreg_partial) != "IID"])

write.table(prepreg_complete, "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_anthro_prepreg_complete.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
write.table(prepreg_partial, "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_anthro_prepreg_partial.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

cat("\nSummary of complete-case dataset (Prepregnancy):\n")
print(summarize_dataframe(prepreg_complete, c("mum_weight_prepreg", "mum_height_prepreg", "mum_bmi_prepreg", "mum_age_prepreg")), row.names = FALSE)

# ----------------------------- POSTNATAL TIMEPOINT DATA -------------------------------

for (tp in names(timepoints)) {
  vars <- timepoints[[tp]]
  weight <- clean_numeric(mother_data[[vars$weight]])
  height <- convert_height(mother_data[[vars$height]])
  age <- mother_data$age_prepreg + vars$age_offset
  bmi <- weight / (height^2)

  temp <- data.frame(
    IID = mother_data$IID,
    weight = weight,
    height = height,
    bmi = bmi,
    age = age
  )
  colnames(temp) <- c("IID", 
                      paste0("weight_", tp),
                      paste0("height_", tp),
                      paste0("bmi_", tp),
                      paste0("age_", tp))

  temp$sex <- 2
  colnames(temp)[colnames(temp) != "IID"] <- paste0("mum_", colnames(temp)[colnames(temp) != "IID"])

  temp_complete <- temp[complete.cases(temp), ]
  temp_partial  <- temp[!is.na(temp[[1]]), ]
  temp_partial[is.na(temp_partial)] <- "."

  write.table(temp_complete, paste0("/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_anthro_", tp, "_complete.txt"), sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
  write.table(temp_partial,  paste0("/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_anthro_", tp, "_partial.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
}
