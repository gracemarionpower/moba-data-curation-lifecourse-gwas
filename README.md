# data-curation-lifecourse-gwas

This repository documents the data curation steps carried out on the HUNT Cloud to prepare MoBa BMI and height data for the Lifecourse GWAS processing pipeline. For requirements, see: https://github.com/MRCIEU/Lifecourse-GWAS/wiki.

The pipeline performs the following tasks:

1. Generates a static covariate file containing FID, IID, sex, yob for children, mothers, and fathers. Each individual is treated independently.
2. Generates two separate long-format time-varying trait files (one for BMI and one for height) with columns: FID, IID, value, age. Individuals may have multiple repeated measurements, each recorded independently and will be treated in the pipeline as seperate individuals.

NB. These datasets are curated to run in GWAS with pgen files and the FID is derived by merging with the psam file 

Collaborators include: Marc Vaudel and Stefan Johansson, University of Bergen
