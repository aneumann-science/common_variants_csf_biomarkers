library(data.table)
library(dplyr)
library(lavaan)
library(NormPsy)
library(ggplot2)

# Load phenotype
load("EMIF_AD_CSF_PCA_GWAS/phenotypes/pca_cov_emif.Rdata")
load("EMIF_AD_CSF_PCA_GWAS/phenotypes/pca_cov_adni.Rdata")

# Merge EMIF and ADNI 
pca_cov_emif.data$MCI <- pca_cov_emif.data$MCI_Baseline_Diagnosis
pca_cov_emif.data$AD <- pca_cov_emif.data$AD_Baseline_Diagnosis
pca_cov_emif.data$GENDER <- pca_cov_emif.data$Gender
pca_cov_emif.data$AGE <- pca_cov_emif.data$ageatcsfpet
emif_adni.data <- bind_rows(pca_cov_emif.data, pca_cov_adni.data)

# Create unified ordered diagnosis variable
# NC = 0, MCI = 1, AD = 2
emif_adni.data$diagnosis <- emif_adni.data$MCI
emif_adni.data$diagnosis[emif_adni.data$AD == 1] <- 2
emif_adni.data$diagnosis <- ordered(emif_adni.data$diagnosis)

# Sex interactions
emif_adni.data$SEX <- emif_adni.data$GENDER
emif_adni.data$sexXrc1 <- emif_adni.data$SEX*emif_adni.data$RC1
emif_adni.data$sexXrc2 <- emif_adni.data$SEX*emif_adni.data$RC2
emif_adni.data$sexXrc3 <- emif_adni.data$SEX*emif_adni.data$RC3
emif_adni.data$sexXrc4 <- emif_adni.data$SEX*emif_adni.data$RC4
emif_adni.data$sexXrc5 <- emif_adni.data$SEX*emif_adni.data$RC5

# Sex interaction model
biomarker_sex_int.model <- '
  diagnosis ~ sexXrc1 + sexXrc2 + sexXrc3 + sexXrc4 + sexXrc5 + rc1_diagnosis*RC1 + rc2_diagnosis*RC2 + rc3_diagnosis*RC3 + rc4_diagnosis*RC4 + rc5_diagnosis*RC5 + SEX + AGE + study 
'
biomarker_sex_int.fit <- sem(biomarker_sex_int.model, data = emif_adni.data, estimator = "WLSMV")
summary(biomarker_sex_int.fit)
biomarker_sex_int.results <- parameterestimates(biomarker_sex_int.fit)
write.csv(biomarker_sex_int.results, file = "results/biomarker_sex_int.csv")

# Sensitvity logistic regression
summary(glm(AD ~ sexXrc1 + sexXrc2 + sexXrc3 + sexXrc4 + sexXrc5 + RC1 + RC2 + RC3 + RC4 + RC5 + SEX + AGE + study, data = emif_adni.data[emif_adni.data$MCI != 1, ]))
summary(glm(MCI ~ sexXrc1 + sexXrc2 + sexXrc3 + sexXrc4 + sexXrc5 + RC1 + RC2 + RC3 + RC4 + RC5 + SEX + AGE + study, data = emif_adni.data[emif_adni.data$AD != 1, ]))
