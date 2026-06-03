# --------------------------------------------------------------------------------------
# Script Name : checking_cleaning_data_parents_HDGB_under8FID.R
# Purpose     : Clean parent anthropometric data and generate Lifecourse-GWAS inputs
#               using HDGB raw phenotypes and under-8 genotype FID/IID structure
# Date created: 03-06-2026
# Author      : Grace M. Power
# --------------------------------------------------------------------------------------

library(data.table)

# ----------------------------- SETUP ----------------------------------

parent_file <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_26-03-23_hdgb/parent.gz"

fid_file <- "/home/grace.power/work/gpower/analysis/LifecourseGWAS_pipeline_run_under8s/genotype_input_dir/scratch/tophits.fam"

root <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation"
parents_dir <- file.path(root, "parents/hdgb_under8FID")
dir.create(parents_dir, recursive = TRUE, showWarnings = FALSE)

parent <- read.delim(
  parent_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

# ----------------------------- FID/IID MAP ----------------------------

fidmap <- fread(fid_file, header = FALSE)
setnames(fidmap, c("FID", "IID", "father", "mother", "sex_geno", "pheno"))
fidmap <- fidmap[, .(FID, IID)]

iid2fid <- setNames(fidmap$FID, fidmap$IID)

# ----------------------------- FUNCTIONS -------------------------------

clean_numeric <- function(x) {
  x <- gsub(",", ".", x)
  x[x %in% c("", ".", "NA")] <- NA
  suppressWarnings(as.numeric(x))
}

convert_height <- function(h_cm) {
  h <- clean_numeric(h_cm)
  h[h < 30 | h > 500] <- NA
  h / 100
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
    x <- suppressWarnings(as.numeric(df[[var]]))
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

add_fid <- function(df) {
  df$FID <- iid2fid[df$IID]
  df <- df[!is.na(df$FID), ]
  df <- df[, c("FID", "IID", setdiff(names(df), c("FID", "IID")))]
  df
}

write_tsv <- function(df, filename) {
  write.table(
    df,
    file = file.path(parents_dir, filename),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = "."
  )
}

# ----------------------------- VARIABLE DEFINITIONS --------------------

core_vars <- c("mother_sentrix_id", "mother_age_15w")
prepreg_vars <- c("mother_weight_beginning_self", "mother_height", "mother_height_self")

timepoints <- list(
  "3y"   = list(weight = "mother_weight_3y",  height = "mother_height_3y",  age_offset = 3 + 0.75),
  "5y"   = list(weight = "mother_weight_5y",  height = "mother_height_5y",  age_offset = 5 + 0.75),
  "8y"   = list(weight = "mother_weight_8y",  height = "mother_height_8y",  age_offset = 8 + 0.75),
  "14y"  = list(weight = "weight_mother_14m", height = "height_mother_14m", age_offset = 14 + 0.75),
  "1.5y" = list(weight = "mother_weight_18m", height = "mother_height_3y",  age_offset = 1.5 + 0.75)
)

all_vars <- unique(c(
  core_vars,
  prepreg_vars,
  unlist(lapply(timepoints, function(x) c(x$weight, x$height)))
))

missing_vars <- setdiff(all_vars, names(parent))
if (length(missing_vars) > 0) {
  stop("Missing expected mother columns: ", paste(missing_vars, collapse = ", "))
}

mother_data <- parent[, all_vars]
colnames(mother_data)[colnames(mother_data) == "mother_sentrix_id"] <- "IID"

# ----------------------------- MOTHER: PREPREG -------------------------

mother_data$mother_age_15w <- clean_numeric(mother_data$mother_age_15w)
mother_data$age_prepreg <- mother_data$mother_age_15w - (15 / 52.1775)

mother_data$mother_height_self <- convert_height(mother_data$mother_height_self)
mother_data$mother_height <- convert_height(mother_data$mother_height)

cat("\nHeight variable summary:\n")
print(rbind(
  cbind(Variable = "mother_height_self", summarize_var(mother_data$mother_height_self)),
  cbind(Variable = "mother_height", summarize_var(mother_data$mother_height))
), row.names = FALSE)

mother_data$mother_weight_beginning_self <- clean_numeric(mother_data$mother_weight_beginning_self)

prepreg <- data.frame(
  IID = mother_data$IID,
  weight_prepreg = mother_data$mother_weight_beginning_self,
  height_prepreg = mother_data$mother_height,
  age_prepreg = mother_data$age_prepreg,
  stringsAsFactors = FALSE
)

prepreg$bmi_prepreg <- prepreg$weight_prepreg / (prepreg$height_prepreg^2)
prepreg <- prepreg[, c("IID", "weight_prepreg", "height_prepreg", "bmi_prepreg", "age_prepreg")]

prepreg_complete <- prepreg[complete.cases(prepreg), ]
prepreg_partial <- prepreg[!is.na(prepreg$IID), ]

prepreg_complete$sex <- 2
prepreg_partial$sex <- 2

colnames(prepreg_complete)[colnames(prepreg_complete) != "IID"] <-
  paste0("mum_", colnames(prepreg_complete)[colnames(prepreg_complete) != "IID"])

colnames(prepreg_partial)[colnames(prepreg_partial) != "IID"] <-
  paste0("mum_", colnames(prepreg_partial)[colnames(prepreg_partial) != "IID"])

prepreg_complete <- add_fid(prepreg_complete)
prepreg_partial <- add_fid(prepreg_partial)
prepreg_partial[is.na(prepreg_partial)] <- "."

write_tsv(prepreg_complete, "mother_anthro_prepreg_complete.txt")
write_tsv(prepreg_partial, "mother_anthro_prepreg_partial.txt")

cat("\nSummary of complete-case dataset, prepregnancy:\n")
print(summarise_dataframe(
  prepreg_complete,
  c("mum_weight_prepreg", "mum_height_prepreg", "mum_bmi_prepreg", "mum_age_prepreg")
), row.names = FALSE)

# ----------------------------- MOTHER: POSTNATAL -----------------------

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
    age = age,
    stringsAsFactors = FALSE
  )

  colnames(temp) <- c(
    "IID",
    paste0("weight_", tp),
    paste0("height_", tp),
    paste0("bmi_", tp),
    paste0("age_", tp)
  )

  temp$sex <- 2

  colnames(temp)[colnames(temp) != "IID"] <-
    paste0("mum_", colnames(temp)[colnames(temp) != "IID"])

  temp_complete <- temp[complete.cases(temp), ]
  temp_partial <- temp[!is.na(temp$IID), ]

  temp_complete <- add_fid(temp_complete)
  temp_partial <- add_fid(temp_partial)
  temp_partial[is.na(temp_partial)] <- "."

  write_tsv(temp_complete, paste0("mother_anthro_", tp, "_complete.txt"))
  write_tsv(temp_partial, paste0("mother_anthro_", tp, "_partial.txt"))
}

# ----------------------------- FATHERS --------------------------------

father_id_col <- "father_sentrix_id"
if (!father_id_col %in% names(parent)) stop("Missing column: father_sentrix_id")

dad_timepoints <- list(
  "hf" = list(
    height = "height_hf",
    weight = "weight_hf",
    age = "age_answering_q_hf",
    height_proxy = NA_character_
  ),
  "45f" = list(
    height = NA_character_,
    weight = "weight_now_45f",
    age = "age_answering_q_45f",
    height_proxy = "height_hf"
  )
)

make_and_write_dad_tp <- function(tp, spec, df_parent) {
  if (is.na(spec$age) || !(spec$age %in% names(df_parent))) {
    message("Skipping ", tp, " — age column missing.")
    return(invisible(NULL))
  }

  a <- clean_numeric(df_parent[[spec$age]])

  has_height_col <- !is.na(spec$height) && spec$height %in% names(df_parent)
  has_height_proxy <- !is.na(spec$height_proxy) && spec$height_proxy %in% names(df_parent)

  h <- if (has_height_col) {
    convert_height(df_parent[[spec$height]])
  } else if (has_height_proxy) {
    convert_height(df_parent[[spec$height_proxy]])
  } else {
    NA_real_
  }

  height_source <- if (has_height_col) {
    spec$height
  } else if (has_height_proxy) {
    spec$height_proxy
  } else {
    NA_character_
  }

  w <- if (!is.na(spec$weight) && spec$weight %in% names(df_parent)) {
    clean_numeric(df_parent[[spec$weight]])
  } else {
    NA_real_
  }

  bmi <- ifelse(!is.na(h) & h > 0 & !is.na(w), w / (h^2), NA_real_)

  tmp <- data.frame(
    IID = df_parent[[father_id_col]],
    weight = w,
    height = h,
    bmi = bmi,
    age = a,
    stringsAsFactors = FALSE
  )

  names(tmp) <- c(
    "IID",
    paste0("weight_", tp),
    paste0("height_", tp),
    paste0("bmi_", tp),
    paste0("age_", tp)
  )

  tmp$sex <- 1

  names(tmp)[names(tmp) != "IID"] <-
    paste0("dad_", names(tmp)[names(tmp) != "IID"])

  if (tp == "45f") {
    tmp$dad_height_45f_source <- ifelse(!is.na(height_source), height_source, ".")
  }

  tmp <- tmp[!is.na(tmp[[paste0("dad_age_", tp)]]), , drop = FALSE]

  comp_vars <- c(paste0("dad_age_", tp), paste0("dad_weight_", tp))
  if (tp %in% c("hf", "45f")) {
    comp_vars <- c(comp_vars, paste0("dad_height_", tp))
  }

  tmp_complete <- tmp[complete.cases(tmp[, comp_vars, drop = FALSE]), ]
  tmp_partial <- tmp[!is.na(tmp$IID), ]

  tmp_complete <- add_fid(tmp_complete)
  tmp_partial <- add_fid(tmp_partial)
  tmp_partial[is.na(tmp_partial)] <- "."

  write_tsv(tmp_complete, paste0("father_anthro_", tp, "_complete.txt"))
  write_tsv(tmp_partial, paste0("father_anthro_", tp, "_partial.txt"))

  cat("\nFather ", tp,
      " — kept age present: ", nrow(tmp),
      " | complete after FID merge: ", nrow(tmp_complete),
      if (tp == "45f") paste0(" | height source used: ", unique(na.omit(tmp$dad_height_45f_source))),
      "\n", sep = "")
}

for (tp in names(dad_timepoints)) {
  make_and_write_dad_tp(tp, dad_timepoints[[tp]], parent)
}

# ----------------------------- LONG FILES ------------------------------

make_long <- function(measure = c("bmi", "height")) {
  measure <- match.arg(measure)

  files <- list.files(parents_dir, pattern = "_complete\\.txt$", full.names = TRUE)
  out <- list()
  k <- 0L

  for (f in files) {
    df <- tryCatch(
      read.table(
        f,
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(df) || !all(c("FID", "IID") %in% names(df))) next

    mcols <- grep(paste0("(^|_)", measure, "_"), names(df), value = TRUE)
    acols <- grep("(^|_)age_", names(df), value = TRUE)

    if (!length(mcols) || !length(acols)) next

    tp_measure <- sub(paste0(".*?_?", measure, "_"), "", mcols)
    tp_age <- sub(".*?_?age_", "", acols)
    tps <- intersect(tp_measure, tp_age)

    if (!length(tps)) next

    for (tp in tps) {
      mc <- mcols[tp_measure == tp][1]
      ac <- acols[tp_age == tp][1]

      tmp <- df[, c("FID", "IID", mc, ac)]
      names(tmp) <- c("FID", "IID", "value", "age")

      tmp$value[tmp$value %in% c(".", "")] <- NA
      tmp$age[tmp$age %in% c(".", "")] <- NA

      suppressWarnings({
        tmp$value <- as.numeric(tmp$value)
        tmp$age <- as.numeric(tmp$age)
      })

      tmp <- tmp[!is.na(tmp$value) & !is.na(tmp$age), , drop = FALSE]

      if (nrow(tmp)) {
        k <- k + 1L
        out[[k]] <- tmp
      }
    }
  }

  if (!length(out)) {
    return(data.frame(FID = character(), IID = character(), value = double(), age = double()))
  }

  ans <- do.call(rbind, out)
  ans <- ans[order(ans$FID, ans$IID, ans$age, ans$value), ]
  unique(ans)
}

bmi_long <- make_long("bmi")
height_long <- make_long("height")

write_tsv(bmi_long, "bmi.txt")
write_tsv(height_long, "height.txt")

# ----------------------------- FINAL CHECKS ----------------------------

cat("\nFinal outputs:\n")
cat("BMI rows:", nrow(bmi_long), " | unique IIDs:", length(unique(bmi_long$IID)), "\n")
cat("Height rows:", nrow(height_long), " | unique IIDs:", length(unique(height_long$IID)), "\n")

cat("\nFiles written to:\n", parents_dir, "\n")
