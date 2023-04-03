###### Script to compute CSF biomarker PC scores
# Author: Alexander Neumann a.neumann@erasmusmc.nl

# Load packages
library(psych)    #For descriptives and PCA functions
library(missMDA)  #PCA-based imputation

#### Load dataset
load("...")

# The script assumes that the biomarker levels are included in
# a data.frame called biomarker.data, which contains following columns:
# ID, TAU, PTAU, CSFNFL, ABETA, CSFNG, YKL40
# Note, the CSF biomarkers should have their original distribution, they will
# be later transformed with INT

# Define sets of CSF biomarker variable names
csf_phenotypes <- c("TAU","PTAU","CSFNFL","ABETA","CSFNG","YKL40")

# Extract CSF phenotypes and reduce to observations with at least 4 biomarkers
# Note: pTau and Tau are considered 1 biomarker in this case
biomarker.data <- biomarker.data[rowMeans(!is.na(biomarker.data[c("TAU","CSFNFL","ABETA","CSFNG","YKL40")])) >= 4/5, c("ID", csf_phenotypes)]

# INT transformation of the CSF biomarkers
biomarker.data[csf_phenotypes] <- lapply(biomarker.data[csf_phenotypes], function(x) {
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
})

# Check, whether successful
multi.hist(biomarker.data[csf_phenotypes])

### Perform PCA based imputation
# Set seed to make analysis replicable
set.seed(20200617)
# Impute based on 5 PCs
biomarker_int.imputePCA <- imputePCA(biomarker.data[csf_phenotypes], ncp = 5)
# Extract imputed data
biomarker_int_imputed.data <- as.data.frame(biomarker_int.imputePCA$completeObs)
# Add back the participant ID
biomarker_int_imputed.data$ID <- biomarker.data$ID

### compute PC scores
# Load fitted PCA from EMIF + ADNI analyses
load("results/pca.Rdata")
# Compute PC scores based on weights in discovery
pca_biomarker.data <- data.frame(predict(pca.fit, data = biomarker_int_imputed.data[csf_phenotypes]))
# ADD ID
pca_biomarker.data$ID <- biomarker.data$ID


### Descriptives and sensitivity checks
pca_biomarker_single.data <- merge(pca_biomarker.data, biomarker_int_imputed.data, by = "ID")
describe(pca_biomarker_single.data)
cor(pca_biomarker_single.data[c(csf_phenotypes,paste0("RC",1:5))])