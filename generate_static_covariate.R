# In bash, run: R

# --------------------------------------------------------------------------------------
# Script Name : generate_static_covariates.R
# Purpose     : To generate a static_covariates.txt file for the Lifecourse GWAS pipeline
# Date created: 13-05-2025
# Last updated: 06-10-2025
# Author      : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# Define input file paths
input_file  <- "/home/grace.power/archive/moba/pheno/v12/pheno_anthropometrics_25-09-04_HDGB_compatible/mfr.gz"
psam_file   <- "/home/grace.power/archive/moba/geno/HDGB-MoBaGenetics/2025.01.30_beta/moba_genotypes_2025.01.30_common.psam"
output_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/static_covariates.txt"

# Read compressed mfr file
con <- gzfile(input_file, "rt")
mfr <- read.delim(con, stringsAsFactors = FALSE)
close(con)

# Read the .psam file to get FID-IID mapping
psam <- read.table(psam_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "")
colnames(psam)[1] <- "FID"  # Rename '#FID' to 'FID'

# Create list of rows from mfr
rows <- list()
for (i in 1:nrow(mfr)) {
  child_id   <- mfr$child_sentrix_id[i]
  mother_id  <- mfr$mother_sentrix_id[i]
  father_id  <- mfr$father_sentrix_id[i]
  child_sex  <- mfr$sex[i]
  child_yob  <- mfr$birth_year[i]
  mother_yob <- mfr$mother_birth_year[i]
  father_yob <- mfr$father_birth_year[i]

  # Add child
  rows[[length(rows) + 1]] <- c(IID = child_id, sex = child_sex, yob = child_yob)
  # Add mother
  if (!is.na(mother_id)) {
    rows[[length(rows) + 1]] <- c(IID = mother_id, sex = 2, yob = mother_yob)
  }
  # Add father
  if (!is.na(father_id)) {
    rows[[length(rows) + 1]] <- c(IID = father_id, sex = 1, yob = father_yob)
  }
}

# Convert to data frame
static_covariate <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)

# Convert types
static_covariate$sex <- as.numeric(static_covariate$sex)
static_covariate$yob <- as.numeric(static_covariate$yob)

# Merge FID using psam file
static_covariate <- merge(psam[, c("FID", "IID")], static_covariate, by = "IID")

# Reorder columns
static_covariate <- static_covariate[, c("FID", "IID", "sex", "yob")]

# Remove duplicates
static_covariate <- unique(static_covariate)

# Write to file
write.table(static_covariate,
            file = output_file,
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("✅ static_covariates.txt written with", nrow(static_covariate), "unique rows\n")




