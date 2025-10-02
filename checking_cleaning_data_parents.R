# --------------------------------------------------------------------------------------
# Script Name : checking_cleaning_data_parents.R
# Purpose     : Clean maternal data and generate timepoint-specific datasets
# Date created: 13-05-2025
# Last updated: 02-10-2025
# Author      : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# NB: Maternal data during pregnancy is excluded as physiological changes during this period may introduce systematic differences

# ----------------------------- SETUP ----------------------------------

# Load data
parent_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-05-07_Mikko/parent.gz"
parent <- read.delim(parent_file, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")

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

core_vars <- c("mother_sentrix_id", "mother_age_15w")
prepreg_vars <- c("mother_weight_beginning_self", "mother_height", "mother_height_self")

# NOTE: `mother_weight_14m` and `mother_height_14m` are treated as 14 YEARS (not months)
# `mother_weight_18m` is assumed to occur at 1.5 years (18 months) and we use height from 3y timepoint

timepoints <- list(
  "3y"   = list(weight = "mother_weight_3y", height = "mother_height_3y", age_offset = 3 + 0.75),
  "5y"   = list(weight = "mother_weight_5y", height = "mother_height_5y", age_offset = 5 + 0.75),
  "8y"   = list(weight = "mother_weight_8y", height = "mother_height_8y", age_offset = 8 + 0.75),
  "14y"  = list(weight = "weight_mother_14m", height = "height_mother_14m", age_offset = 14 + 0.75),  # Note: "14m" interpreted as 14 years
  "1.5y" = list(weight = "mother_weight_18m", height = "mother_height_3y", age_offset = 1.5 + 0.75)   # 1.5 years (18 months) with 3y height
)

all_vars <- unique(c(core_vars, prepreg_vars, unlist(lapply(timepoints, function(x) c(x$weight, x$height)))))
mother_data <- parent[, all_vars]
colnames(mother_data)[colnames(mother_data) == "mother_sentrix_id"] <- "IID"

# ----------------------------- CLEAN PREPREG AGE DATA -------------------------------

mother_data$mother_age_15w <- clean_numeric(mother_data$mother_age_15w)
mother_data$age_prepreg <- mother_data$mother_age_15w - (15 / 52.1775)

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

prepreg_complete <- prepreg[complete.cases(prepreg), ]
prepreg_partial  <- prepreg[!is.na(prepreg$IID), ]
prepreg_partial[is.na(prepreg_partial)] <- "."

prepreg_complete$sex <- 2
prepreg_partial$sex <- 2
colnames(prepreg_complete)[colnames(prepreg_complete) != "IID"] <- paste0("mum_", colnames(prepreg_complete)[colnames(prepreg_complete) != "IID"])
colnames(prepreg_partial)[colnames(prepreg_partial) != "IID"] <- paste0("mum_", colnames(prepreg_partial)[colnames(prepreg_partial) != "IID"])

write.table(prepreg_complete, "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/parents/mother_anthro_prepreg_complete.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
write.table(prepreg_partial, "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/parents/mother_anthro_prepreg_partial.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

cat("\nSummary of complete-case dataset (Prepregnancy):\n")
print(summarise_dataframe(prepreg_complete, c("mum_weight_prepreg", "mum_height_prepreg", "mum_bmi_prepreg", "mum_age_prepreg")), row.names = FALSE)

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

  write.table(temp_complete, paste0("/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/parents/mother_anthro_", tp, "_complete.txt"), sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
  write.table(temp_partial,  paste0("/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/parents/mother_anthro_", tp, "_partial.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
}

# ============================= Fathers (age required; 45f BMI uses hf height) =============================

father_id_col <- "father_sentrix_id"
if (!father_id_col %in% names(parent)) stop("Missing column: father_sentrix_id")

# Timepoints with age available
dad_timepoints <- list(
  "hf" = list(height = "height_hf",      weight = "weight_hf",      age = "age_answering_q_hf",  height_proxy = NA_character_),
  "45f"= list(height = NA_character_,     weight = "weight_now_45f", age = "age_answering_q_45f", height_proxy = "height_hf")  # use hf height
)

make_and_write_dad_tp <- function(tp, spec, df_parent) {
  # Require age column
  if (is.na(spec$age) || !(spec$age %in% names(df_parent))) {
    message("Skipping ", tp, " — age column missing."); return(invisible(NULL))
  }

  a <- clean_numeric(df_parent[[spec$age]])
  # Height: direct if present, else proxy (for 45f -> height_hf)
  has_height_col   <- !is.na(spec$height)      && (spec$height      %in% names(df_parent))
  has_height_proxy <- !is.na(spec$height_proxy) && (spec$height_proxy %in% names(df_parent))

  h <- if (has_height_col) {
    convert_height(df_parent[[spec$height]])
  } else if (has_height_proxy) {
    convert_height(df_parent[[spec$height_proxy]])
  } else {
    NA_real_
  }
  height_source <- if (has_height_col) spec$height else if (has_height_proxy) spec$height_proxy else NA_character_

  # Weight
  w <- if (!is.na(spec$weight) && spec$weight %in% names(df_parent)) clean_numeric(df_parent[[spec$weight]]) else NA_real_

  # BMI (height is in meters from convert_height)
  bmi <- ifelse(!is.na(h) & h > 0 & !is.na(w), w / (h^2), NA_real_)

  # Build frame (keep rows with AGE only; completeness defined below)
  tmp <- data.frame(
    IID    = df_parent[[father_id_col]],
    weight = w,
    height = h,
    bmi    = bmi,
    age    = a,
    stringsAsFactors = FALSE
  )
  names(tmp) <- c("IID",
                  paste0("weight_", tp),
                  paste0("height_", tp),
                  paste0("bmi_", tp),
                  paste0("age_", tp))
  tmp$sex <- 1
  names(tmp)[names(tmp) != "IID"] <- paste0("dad_", names(tmp)[names(tmp) != "IID"])

  # Add a provenance flag for 45f height
  if (tp == "45f") {
    tmp$dad_height_45f_source <- ifelse(!is.na(height_source), height_source, ".")
  }

  # Enforce: must have AGE to be kept at all
  tmp <- tmp[!is.na(tmp[[paste0("dad_age_", tp)]]), , drop = FALSE]

  # Define complete rows:
  #  - hf: require age + height + weight
  #  - 45f: require age + weight + (proxy) height so BMI is computed
  comp_vars <- c(paste0("dad_age_", tp), paste0("dad_weight_", tp))
  if (tp %in% c("hf","45f")) comp_vars <- c(comp_vars, paste0("dad_height_", tp))

  tmp_complete <- tmp[complete.cases(tmp[, comp_vars, drop = FALSE]), ]
  tmp_partial  <- tmp[!is.na(tmp$IID), ]
  tmp_partial[is.na(tmp_partial)] <- "."

  out_base <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation"
  write.table(tmp_complete, file.path(out_base, paste0("father_anthro_", tp, "_complete.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
  write.table(tmp_partial,  file.path(out_base, paste0("father_anthro_", tp, "_partial.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE, na = ".")

  # quick console QC
  cat("\nFather ", tp, " — kept (age present): ", nrow(tmp),
      " | complete: ", nrow(tmp_complete),
      if (tp == "45f") paste0(" | height source used: ", unique(na.omit(tmp$dad_height_45f_source))), "\n", sep = "")
}

# Run
for (tp in names(dad_timepoints)) {
  make_and_write_dad_tp(tp, dad_timepoints[[tp]], parent)
}
                                  
