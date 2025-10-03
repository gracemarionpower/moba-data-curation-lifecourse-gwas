# In bash, run: R

WIP!

# --------------------------------------------------------------------------------------
# Script Name : generate_time-varying.R
# Purpose     : To generate bmi.txt and height.txt files for the Lifecourse GWAS pipeline
# Date created: 13-05-2025
# Last updated: 13-05-2025
# Author      : Grace M. Power
# Collaborators: Marc Vaudel and Stefan Johansson, University of Bergen
# --------------------------------------------------------------------------------------

# Setup file paths

input_file <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/mother_anthro_prepreg_complete.txt"
psam_file <- "/home/grace.power/archive/moba/geno/HDGB-MoBaGenetics/2025.01.30_beta/moba_genotypes_2025.01.30_common.psam"

output_height <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/height.txt"
output_bmi <- "/home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/bmi.txt"

# Load data

mother_data <- read.delim(input_file, stringsAsFactors = FALSE)
psam <- read.table(psam_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "")
colnames(psam)[1] <- "FID"

# Merge to add FID from psam

merged_data <- merge(psam[, c("FID", "IID")], mother_data, by = "IID")

# Create height file

height_data <- merged_data[, c("FID", "IID", "height_prepreg", "age_prepreg")]
colnames(height_data) <- c("FID", "IID", "value", "age")

write.table(height_data, file = output_height, sep = "\t", row.names = FALSE, quote = FALSE)

# Create BMI file

bmi_data <- merged_data[, c("FID", "IID", "bmi_prepreg", "age_prepreg")]
colnames(bmi_data) <- c("FID", "IID", "value", "age")

write.table(bmi_data, file = output_bmi, sep = "\t", row.names = FALSE, quote = FALSE)


