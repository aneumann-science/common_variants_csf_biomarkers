library(data.table)
library(dplyr)
library(lavaan)
library(NormPsy)
library(poolr)
library(BinaryDosage)

# Load phenotype
load("EMIF_AD_CSF_PCA_GWAS/phenotypes/pca_cov_emif.Rdata")
load("EMIF_AD_CSF_PCA_GWAS/phenotypes/pca_cov_adni.Rdata")

# Load covariates and AD
covariates_emif.data <- fread("EMIF_AD_CSF_PCA_GWAS/genotypes/EMIF/EMIF_covar_file.txt")
covariates_adni.data <- fread("EMIF_AD_CSF_PCA_GWAS/genotypes/ADNI/ADNI_covar_file.txt")

# Genotypes
genotype_emif.data <- fread("data/EMIF.raw")
genotype_adni.data <- fread("data/ADNI.raw")

# Merge everything
emif_phenotype_covariates.data <- merge(pca_cov_emif.data[c("DNA_ID","RC1","RC2","RC3","RC4","RC5")], covariates_emif.data, by.x = "DNA_ID", by.y = "#IID")
emif.data <- merge(emif_phenotype_covariates.data, genotype_emif.data[,c(2,7:28)], by.x = "DNA_ID", by.y = "IID")
names(emif.data)[names(emif.data) == 'DNA_ID'] <- 'ID'

adni_phenotype_covariates.data <- merge(pca_cov_adni.data[c("ID","RC1","RC2","RC3","RC4","RC5")], covariates_adni.data, by.x = "ID", by.y = "#IID")
adni.data <- merge(adni_phenotype_covariates.data, genotype_adni.data[,c(2,7:28)], by.x = "ID", by.y = "IID")

# Add study identifier
emif_adni.data <- bind_rows(emif.data, adni.data)

# ADD EMIF ID and fix NA for ADNI_genotyping
emif_adni.data$EMIF <- 0 
emif_adni.data$EMIF[1:672] <- 1
emif_adni.data$ADNI_genotyping[emif_adni.data$EMIF == 1] <- 0

# Create unified ordered diagnosis variable
# NC = 0, MCI = 1, AD = 2
emif_adni.data$diagnosis <- emif_adni.data$MCI
emif_adni.data$diagnosis[emif_adni.data$AD == 1] <- 2
emif_adni.data$diagnosis <- ordered(emif_adni.data$diagnosis)

snp_names <- names(genotype_emif.data[,c(7:28)])

# Get MAF and imputation quality
eaf_rsq.list <- lapply(snp_names, function(snp_name) {
  emif_adni.data$snp <- emif_adni.data[,snp_name]
  eac <- colSums(emif_adni.data$snp)
  eaf <- mean(emif_adni.data$snp, na.rm = TRUE)/2
  p2 <- emif_adni.data$snp - 1
  p2[p2<0] <- 0
  rsq <- getrsq(emif_adni.data$snp, p2=p2)
  data.frame(eaf,rsq)
})

eaf_rsq.data <- do.call(rbind, eaf_rsq.list)
write.csv(eaf_rsq.data, file = "results/eaf_rsq.csv")

### Check correlation between YKL-40 hits
cor(emif_adni.data$`1_203150756_1:203150756_T_C_T`, emif_adni.data$`1_203156080_1:203156080_T_C_T`)^2

# For testing purposes
emif_adni.data$snp <- emif_adni.data$`1_203150756_1:203150756_T_C_T`
emif_adni.data$snp <- emif_adni.data$`19_45411941_19:45411941_T_C_T`
emif_adni.data$snp <- emif_adni.data$`7_12270770_7:12270770_T_A_T`

###### SEM
### Main effects mediation model
mediation.model <- '
  # SNP -> PC
  RC1 ~ snp_RC1*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC2 ~ snp_RC2*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC3 ~ snp_RC3*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC4 ~ snp_RC4*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC5 ~ snp_RC5*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF

  # PC -> Latent AD  
  diagnosis ~ rc1_diagnosis*RC1 + rc2_diagnosis*RC2 + rc3_diagnosis*RC3 + rc4_diagnosis*RC4 + rc5_diagnosis*RC5 + snp_diagnosis*snp  + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  
  # Mediation effect
  rc1_mediation := snp_RC1*rc1_diagnosis
  rc2_mediation := snp_RC2*rc2_diagnosis
  rc3_mediation := snp_RC3*rc3_diagnosis
  rc4_mediation := snp_RC4*rc4_diagnosis
  rc5_mediation := snp_RC5*rc5_diagnosis
  
  # All mediation effects
  biomarker_mediation := rc1_mediation+rc2_mediation+rc3_mediation+rc4_mediation+rc5_mediation
  
  # Total effect
  total := biomarker_mediation + snp_diagnosis
  
  # Proportions
  rc1_mediation_prop := rc1_mediation/total
  rc2_mediation_prop := rc2_mediation/total
  rc3_mediation_prop := rc3_mediation/total
  rc4_mediation_prop := rc4_mediation/total
  rc5_mediation_prop := rc5_mediation/total
  direct_prop := snp_diagnosis/total
  biomarker_mediation_prop := biomarker_mediation/total
'

# Loop through all SNPs
mediation.list <- lapply(snp_names, function(snp_name) {
  # Assign the current SNP name as "snp"
  emif_adni.data$snp <- emif_adni.data[,snp_name]
  
  # Fit model
  mediation.fit <- sem(mediation.model, data = emif_adni.data, estimator = "WLSMV")
  
  # Format results
  results.data <- parameterestimates(mediation.fit)
  rc1_mediation.data <- results.data[results.data$label == "rc1_mediation",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc1_mediation.data) <- paste0(names(rc1_mediation.data),"_rc1_mediation")
  rc2_mediation.data <- results.data[results.data$label == "rc2_mediation",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc2_mediation.data) <- paste0(names(rc2_mediation.data),"_rc2_mediation")
  rc3_mediation.data <- results.data[results.data$label == "rc3_mediation",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc3_mediation.data) <- paste0(names(rc3_mediation.data),"_rc3_mediation")
  rc4_mediation.data <- results.data[results.data$label == "rc4_mediation",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc4_mediation.data) <- paste0(names(rc4_mediation.data),"_rc4_mediation")
  rc5_mediation.data <- results.data[results.data$label == "rc5_mediation",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc5_mediation.data) <- paste0(names(rc5_mediation.data),"_rc5_mediation")
  biomarker_mediation.data <- results.data[results.data$label == "biomarker_mediation",c("est","se","pvalue","ci.lower","ci.upper")]
  names(biomarker_mediation.data) <- paste0(names(biomarker_mediation.data),"_biomarker_mediation")
  direct.data <- results.data[results.data$label == "snp_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(direct.data) <- paste0(names(direct.data),"_direct")
  total.data <- results.data[results.data$label == "total",c("est","se","pvalue","ci.lower","ci.upper")]
  names(total.data) <- paste0(names(total.data),"_total")
  
  # Combine all results
  cbind(snp_name,total.data,direct.data,biomarker_mediation.data,rc1_mediation.data,rc2_mediation.data,rc3_mediation.data,rc4_mediation.data,rc5_mediation.data)
})

# Merge all results
mediation.data <- do.call(rbind, mediation.list)

# Format and export
mediation_formatted_rounded <- round(mediation.data[c("est_rc1_mediation","est_rc2_mediation","est_rc3_mediation","est_rc4_mediation","est_rc5_mediation","est_biomarker_mediation","est_direct","est_total")],2)
mediation_formatted.data <- cbind(mediation.data$snp_name,mediation_formatted_rounded)
mediation_formatted.data$prop <- round(mediation.data$est_biomarker_mediation/mediation.data$est_total, 2)
write.csv(mediation.data, file = "results/mediation.csv", row.names = F, quote = F)
write.csv(mediation_formatted.data, file = "results/mediation_formatted.csv", row.names = F, quote = F)

mediation_formatted_rounded_p <- round(mediation.data[c("pvalue_rc1_mediation","pvalue_rc2_mediation","pvalue_rc3_mediation","pvalue_rc4_mediation","pvalue_rc5_mediation","pvalue_biomarker_mediation","pvalue_direct","pvalue_total")],3)
mediation_formatted_p.data <- cbind(mediation.data$snp_name,mediation_formatted_rounded_p)
mediation_formatted_p.data[mediation_formatted_p.data < (0.05/6)] <- "sig"
write.csv(mediation_formatted_p.data, file = "results/mediation_formateed_p.csv", row.names = F, quote = F)

### Sex interaction mediaion models
# Make male = 0 and female 1
emif_adni.data$SEX <- emif_adni.data$SEX - 1

# Sex interaction mdediation model
mediation_int.model <- '
  # SNP -> PC
  RC1 ~ intRC1*sexXsnp + snp_RC1*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC2 ~ intRC2*sexXsnp + snp_RC2*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC3 ~ intRC3*sexXsnp + snp_RC3*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC4 ~ intRC4*sexXsnp + snp_RC4*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC5 ~ intRC5*sexXsnp + snp_RC5*snp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
 
  # PC -> Latent AD   
  diagnosis ~ rc1_diagnosis*RC1 + rc2_diagnosis*RC2 + rc3_diagnosis*RC3 + rc4_diagnosis*RC4 + rc5_diagnosis*RC5 + snp_diagnosis*snp + snp_diagnosis_int*sexXsnp + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  
  # Mediation effect males
  rc1_mediation_male := snp_RC1*rc1_diagnosis
  rc2_mediation_male := snp_RC2*rc2_diagnosis
  rc3_mediation_male := snp_RC3*rc3_diagnosis
  rc4_mediation_male := snp_RC4*rc4_diagnosis
  rc5_mediation_male := snp_RC5*rc5_diagnosis
  
  # SNP -> PC females
  snp_RC1_female_total := snp_RC1 + intRC1
  snp_RC2_female_total := snp_RC2 + intRC2
  snp_RC3_female_total := snp_RC3 + intRC3
  snp_RC4_female_total := snp_RC4 + intRC4
  snp_RC5_female_total := snp_RC5 + intRC5
  # Direct females
  snp_diagnosis_female_total := snp_diagnosis + snp_diagnosis_int
  
  # Mediation effect females
  rc1_mediation_female := snp_RC1_female_total*rc1_diagnosis
  rc2_mediation_female := snp_RC2_female_total*rc2_diagnosis
  rc3_mediation_female := snp_RC3_female_total*rc3_diagnosis
  rc4_mediation_female := snp_RC4_female_total*rc4_diagnosis
  rc5_mediation_female := snp_RC5_female_total*rc5_diagnosis
  
  # All mediation effects male
  biomarker_mediation_male := rc1_mediation_male+rc2_mediation_male+rc3_mediation_male+rc4_mediation_male+rc5_mediation_male
  # Total male
  total_male := biomarker_mediation_female + snp_diagnosis

  # All mediation effects female  
  biomarker_mediation_female := rc1_mediation_female+rc2_mediation_female+rc3_mediation_female+rc4_mediation_female+rc5_mediation_female
  # Totalmale
  total_female := biomarker_mediation_female + snp_diagnosis_female_total
  
  # Sex differences
  rc1_mediation_dif := rc1_mediation_female - rc1_mediation_male
  rc2_mediation_dif := rc2_mediation_female - rc2_mediation_male
  rc3_mediation_dif := rc3_mediation_female - rc3_mediation_male
  rc4_mediation_dif := rc4_mediation_female - rc4_mediation_male
  rc5_mediation_dif := rc5_mediation_female - rc5_mediation_male
  mediation_dif := total_female - total_male
  total_dif := total_female - total_male
'

# For testing
emif_adni.data$sexXsnp <- emif_adni.data$SEX*emif_adni.data$snp
mediation_int.fit <- sem(mediation_int.model, data = emif_adni.data, estimator = "WLSMV")

# Loop through all SNPs
mediation_int.list <- lapply(snp_names, function(snp_name) {
  # Assign the current SNP name as "snp"
  emif_adni.data$snp <- emif_adni.data[,snp_name]
  # Create Sex*SNP interaction term
  emif_adni.data$sexXsnp <- emif_adni.data$SEX*emif_adni.data$snp
  
  # Fit interaction model
  mediation_int.fit <- sem(mediation_int.model, data = emif_adni.data, estimator = "WLSMV")
  
  ### Format results male
  results.data <- parameterestimates(mediation_int.fit)
  # SNP -> RC
  rc1_snp_male.data <- results.data[results.data$label == "snp_RC1",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc1_snp_male.data) <- paste0(names(rc1_snp_male.data),"_rc1_snp")
  rc2_snp_male.data <- results.data[results.data$label == "snp_RC2",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc2_snp_male.data) <- paste0(names(rc2_snp_male.data),"_rc2_snp")
  rc3_snp_male.data <- results.data[results.data$label == "snp_RC3",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc3_snp_male.data) <- paste0(names(rc3_snp_male.data),"_rc3_snp")
  rc4_snp_male.data <- results.data[results.data$label == "snp_RC4",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc4_snp_male.data) <- paste0(names(rc4_snp_male.data),"_rc4_snp")
  rc5_snp_male.data <- results.data[results.data$label == "snp_RC5",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc5_snp_male.data) <- paste0(names(rc5_snp_male.data),"_rc5_snp")
  # RC -> diagnosis
  rc1_diagnosis_male.data <- results.data[results.data$label == "rc1_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc1_diagnosis_male.data) <- paste0(names(rc1_diagnosis_male.data),"_rc1_diagnosis")
  rc2_diagnosis_male.data <- results.data[results.data$label == "rc2_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc2_diagnosis_male.data) <- paste0(names(rc2_diagnosis_male.data),"_rc2_diagnosis")
  rc3_diagnosis_male.data <- results.data[results.data$label == "rc3_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc3_diagnosis_male.data) <- paste0(names(rc3_diagnosis_male.data),"_rc3_diagnosis")
  rc4_diagnosis_male.data <- results.data[results.data$label == "rc4_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc4_diagnosis_male.data) <- paste0(names(rc4_diagnosis_male.data),"_rc4_diagnosis")
  rc5_diagnosis_male.data <- results.data[results.data$label == "rc5_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc5_diagnosis_male.data) <- paste0(names(rc5_diagnosis_male.data),"_rc5_diagnosis")
  # Mediation
  rc1_mediation_male.data <- results.data[results.data$label == "rc1_mediation_male",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc1_mediation_male.data) <- paste0(names(rc1_mediation_male.data),"_rc1_mediation")
  rc2_mediation_male.data <- results.data[results.data$label == "rc2_mediation_male",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc2_mediation_male.data) <- paste0(names(rc2_mediation_male.data),"_rc2_mediation")
  rc3_mediation_male.data <- results.data[results.data$label == "rc3_mediation_male",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc3_mediation_male.data) <- paste0(names(rc3_mediation_male.data),"_rc3_mediation")
  rc4_mediation_male.data <- results.data[results.data$label == "rc4_mediation_male",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc4_mediation_male.data) <- paste0(names(rc4_mediation_male.data),"_rc4_mediation")
  rc5_mediation_male.data <- results.data[results.data$label == "rc5_mediation_male",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc5_mediation_male.data) <- paste0(names(rc5_mediation_male.data),"_rc5_mediation")
  # Mediation difference
  rc1_mediation_dif.data <- results.data[results.data$label == "rc1_mediation_dif",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc1_mediation_dif.data) <- paste0(names(rc1_mediation_dif.data),"_rc1_mediation_dif")
  rc2_mediation_dif.data <- results.data[results.data$label == "rc2_mediation_dif",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc2_mediation_dif.data) <- paste0(names(rc2_mediation_dif.data),"_rc2_mediation_dif")
  rc3_mediation_dif.data <- results.data[results.data$label == "rc3_mediation_dif",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc3_mediation_dif.data) <- paste0(names(rc3_mediation_dif.data),"_rc3_mediation_dif")
  rc4_mediation_dif.data <- results.data[results.data$label == "rc4_mediation_dif",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc4_mediation_dif.data) <- paste0(names(rc4_mediation_dif.data),"_rc4_mediation_dif")
  rc5_mediation_dif.data <- results.data[results.data$label == "rc5_mediation_dif",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc5_mediation_dif.data) <- paste0(names(rc5_mediation_dif.data),"_rc5_mediation_dif")
  mediation_dif.data <- results.data[results.data$label == "mediation_dif",c("est","se","pvalue","ci.lower","ci.upper")]
  names(mediation_dif.data) <- paste0(names(mediation_dif.data),"_mediation_dif")
  total_dif.data <- results.data[results.data$label == "total_dif",c("est","se","pvalue","ci.lower","ci.upper")]
  names(total_dif.data) <- paste0(names(total_dif.data),"_total_dif")
  # Additional data
  biomarker_mediation_male.data <- results.data[results.data$label == "biomarker_mediation_male",c("est","se","pvalue","ci.lower","ci.upper")]
  names(biomarker_mediation_male.data) <- paste0(names(biomarker_mediation_male.data),"_biomarker_mediation")
  mediation_dif.data <- results.data[results.data$label == "mediation_dif",c("est","se","pvalue","ci.lower","ci.upper")]
  names(mediation_dif.data) <- paste0(names(mediation_dif.data),"_mediation_dif")
  direct_male.data <- results.data[results.data$label == "snp_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(direct_male.data) <- paste0(names(direct_male.data),"_direct")
  direct_dif.data <- results.data[results.data$label == "snp_diagnosis_int",c("est","se","pvalue","ci.lower","ci.upper")]
  names(direct_dif.data) <- paste0(names(direct_dif.data),"_direct_dif")
  total_male.data <- results.data[results.data$label == "total_male",c("est","se","pvalue","ci.lower","ci.upper")]
  names(total_male.data) <- paste0(names(total_male.data),"_total")
  
  ### Format results female
  # SNP -> RC
  rc1_snp_female.data <- results.data[results.data$label == "snp_RC1_female_total",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc1_snp_female.data) <- paste0(names(rc1_snp_female.data),"_rc1_snp")
  rc2_snp_female.data <- results.data[results.data$label == "snp_RC2_female_total",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc2_snp_female.data) <- paste0(names(rc2_snp_female.data),"_rc2_snp")
  rc3_snp_female.data <- results.data[results.data$label == "snp_RC3_female_total",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc3_snp_female.data) <- paste0(names(rc3_snp_female.data),"_rc3_snp")
  rc4_snp_female.data <- results.data[results.data$label == "snp_RC4_female_total",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc4_snp_female.data) <- paste0(names(rc4_snp_female.data),"_rc4_snp")
  rc5_snp_female.data <- results.data[results.data$label == "snp_RC5_female_total",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc5_snp_female.data) <- paste0(names(rc5_snp_female.data),"_rc5_snp")
  
  # SNP effect difference
  rc1_snp_dif.data <- results.data[results.data$label == "intRC1",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc1_snp_dif.data) <- paste0(names(rc1_snp_dif.data),"_rc1_snp_dif")
  rc2_snp_dif.data <- results.data[results.data$label == "intRC2",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc2_snp_dif.data) <- paste0(names(rc2_snp_dif.data),"_rc2_snp_dif")
  rc3_snp_dif.data <- results.data[results.data$label == "intRC3",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc3_snp_dif.data) <- paste0(names(rc3_snp_dif.data),"_rc3_snp_dif")
  rc4_snp_dif.data <- results.data[results.data$label == "intRC4",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc4_snp_dif.data) <- paste0(names(rc4_snp_dif.data),"_rc4_snp_dif")
  rc5_snp_dif.data <- results.data[results.data$label == "intRC5",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc5_snp_dif.data) <- paste0(names(rc5_snp_dif.data),"_rc5_snp_dif")
  
  # RC -> diagnosis
  rc1_diagnosis_female.data <- results.data[results.data$label == "rc1_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc1_diagnosis_female.data) <- paste0(names(rc1_diagnosis_female.data),"_rc1_diagnosis")
  rc2_diagnosis_female.data <- results.data[results.data$label == "rc2_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc2_diagnosis_female.data) <- paste0(names(rc2_diagnosis_female.data),"_rc2_diagnosis")
  rc3_diagnosis_female.data <- results.data[results.data$label == "rc3_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc3_diagnosis_female.data) <- paste0(names(rc3_diagnosis_female.data),"_rc3_diagnosis")
  rc4_diagnosis_female.data <- results.data[results.data$label == "rc4_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc4_diagnosis_female.data) <- paste0(names(rc4_diagnosis_female.data),"_rc4_diagnosis")
  rc5_diagnosis_female.data <- results.data[results.data$label == "rc5_diagnosis",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc5_diagnosis_female.data) <- paste0(names(rc5_diagnosis_female.data),"_rc5_diagnosis")
  
  rc1_mediation_female.data <- results.data[results.data$label == "rc1_mediation_female",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc1_mediation_female.data) <- paste0(names(rc1_mediation_female.data),"_rc1_mediation")
  rc2_mediation_female.data <- results.data[results.data$label == "rc2_mediation_female",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc2_mediation_female.data) <- paste0(names(rc2_mediation_female.data),"_rc2_mediation")
  rc3_mediation_female.data <- results.data[results.data$label == "rc3_mediation_female",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc3_mediation_female.data) <- paste0(names(rc3_mediation_female.data),"_rc3_mediation")
  rc4_mediation_female.data <- results.data[results.data$label == "rc4_mediation_female",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc4_mediation_female.data) <- paste0(names(rc4_mediation_female.data),"_rc4_mediation")
  rc5_mediation_female.data <- results.data[results.data$label == "rc5_mediation_female",c("est","se","pvalue","ci.lower","ci.upper")]
  names(rc5_mediation_female.data) <- paste0(names(rc5_mediation_female.data),"_rc5_mediation")
  biomarker_mediation_female.data <- results.data[results.data$label == "biomarker_mediation_female",c("est","se","pvalue","ci.lower","ci.upper")]
  names(biomarker_mediation_female.data) <- paste0(names(biomarker_mediation_female.data),"_biomarker_mediation")
  direct_female.data <- results.data[results.data$label == "snp_diagnosis_female_total",c("est","se","pvalue","ci.lower","ci.upper")]
  names(direct_female.data) <- paste0(names(direct_female.data),"_direct")
  total_dif.data <- results.data[results.data$label == "total_dif",c("est","se","pvalue","ci.lower","ci.upper")]
  names(total_dif.data) <- paste0(names(total_dif.data),"_total_dif")
  total_female.data <- results.data[results.data$label == "total_female",c("est","se","pvalue","ci.lower","ci.upper")]
  names(total_female.data) <- paste0(names(total_female.data),"_total")
  
  # Combine all results
  results_male.data <- cbind(snp_name,rc1_snp_male.data,rc1_snp_dif.data,rc1_diagnosis_male.data,rc1_mediation_male.data,rc1_mediation_dif.data,rc2_snp_male.data,rc2_snp_dif.data,rc2_diagnosis_male.data,rc2_mediation_male.data,rc2_mediation_dif.data,rc3_snp_male.data,rc3_snp_dif.data,rc3_diagnosis_male.data,rc3_mediation_male.data,rc3_mediation_dif.data,rc4_snp_male.data,rc4_snp_dif.data,rc4_diagnosis_male.data,rc4_mediation_male.data,rc4_mediation_dif.data,rc5_snp_male.data,rc5_snp_dif.data,rc5_diagnosis_male.data,rc5_mediation_male.data,rc5_mediation_dif.data,biomarker_mediation_male.data,mediation_dif.data,direct_male.data,direct_dif.data,total_male.data,total_dif.data)
  results_female.data <- cbind(snp_name,rc1_snp_female.data,rc1_snp_dif.data,rc1_diagnosis_female.data,rc1_mediation_female.data,rc1_mediation_dif.data,rc2_snp_female.data,rc2_snp_dif.data,rc2_diagnosis_female.data,rc2_mediation_female.data,rc2_mediation_dif.data,rc3_snp_female.data,rc3_snp_dif.data,rc3_diagnosis_female.data,rc3_mediation_female.data,rc3_mediation_dif.data,rc4_snp_female.data,rc4_snp_dif.data,rc4_diagnosis_female.data,rc4_mediation_female.data,rc4_mediation_dif.data,rc5_snp_female.data,rc5_snp_dif.data,rc5_diagnosis_female.data,rc5_mediation_female.data,rc5_mediation_dif.data,biomarker_mediation_female.data,mediation_dif.data,direct_female.data,direct_dif.data,total_female.data,total_dif.data)
  names(results_female.data) <- names(results_male.data)
  results_formatted.data <- rbind(results_male.data,results_female.data)
  results_formatted.data$sex <- c("male","female")
  
  return(results_formatted.data)
})

# Merge results from all SNPs
mediation_int.data <- do.call(rbind, mediation_int.list)

# Format and export results
mediation_int_formatted_rounded <- round(mediation_int.data[c("est_rc1_mediation","est_rc2_mediation","est_rc3_mediation","est_rc4_mediation","est_rc5_mediation","est_biomarker_mediation","est_direct","est_total")],2)
mediation_int_formatted.data <- cbind(mediation_int.data$snp_name,mediation_int_formatted_rounded)
mediation_int_formatted.data$prop <- round(mediation_int.data$est_biomarker_mediation/mediation_int.data$est_total, 2)
write.csv(mediation_int.data, file = "results/mediation.csv", row.names = F, quote = F)
write.csv(mediation_int_formatted.data, file = "results/mediation_int_formatted.csv", row.names = F, quote = F)

mediation_int_formatted_rounded_p <- round(mediation_int.data[c("pvalue_rc1_mediation","pvalue_rc2_mediation","pvalue_rc3_mediation","pvalue_rc4_mediation","pvalue_rc5_mediation","pvalue_biomarker_mediation","pvalue_direct","pvalue_total")],3)
mediation_int_formatted_p.data <- cbind(mediation_int.data$snp_name,mediation_int_formatted_rounded_p)
mediation_int_formatted_p.data[mediation_int_formatted_p.data < 0.05] <- "sig"
write.csv(mediation_int_formatted_p.data, file = "results/mediation_int_formatted_p.csv", row.names = F, quote = F)

mediation_int_formatted_rounded_p_dif <- round(mediation_int.data[c("pvalue_rc1_mediation_dif","pvalue_rc2_mediation_dif","pvalue_rc3_mediation_dif","pvalue_rc4_mediation_dif","pvalue_rc5_mediation_dif","pvalue_mediation_dif","pvalue_direct_dif","pvalue_total_dif")],3)
mediation_int_formatted_p_dif.data <- cbind(mediation.data$snp_name,mediation_int_formatted_rounded_p_dif)
mediation_int_formatted_p_dif.data[mediation_int_formatted_p_dif.data < 0.05] <- "sig"
write.csv(mediation_int_formatted_p_dif.data, file = "results/mediation_int_formatted_p_dif.csv", row.names = F, quote = F)


### Sensitivity check, stratified analysis
# SEM
mediation_stratified.model <- '
  RC1 ~ snp_RC1*snp + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC2 ~ snp_RC2*snp + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC3 ~ snp_RC3*snp + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC4 ~ snp_RC4*snp + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  RC5 ~ snp_RC5*snp + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  
  diagnosis ~ rc1_diagnosis*RC1 + rc2_diagnosis*RC2 + rc3_diagnosis*RC3 + rc4_diagnosis*RC4 + rc5_diagnosis*RC5 + snp_diagnosis*snp  + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  snp ~ PC1 + PC2 + PC3 + PC4 + PC5 + ADNI_genotyping + EMIF
  
  rc1_mediation := snp_RC1*rc1_diagnosis
  rc2_mediation := snp_RC2*rc2_diagnosis
  rc3_mediation := snp_RC3*rc3_diagnosis
  rc4_mediation := snp_RC4*rc4_diagnosis
  rc5_mediation := snp_RC5*rc5_diagnosis
  
  biomarker_mediation := rc1_mediation+rc2_mediation+rc3_mediation+rc4_mediation+rc5_mediation
  total := biomarker_mediation + snp_diagnosis
'


emif_adni.data$snp <- emif_adni.data$`7_12270770_7:12270770_T_A_T`
mediation_male.fit <- sem(mediation_stratified.model, data = emif_adni.data[emif_adni.data$SEX == 0, ], estimator = "WLSMV")
mediation_female.fit <- sem(mediation_stratified.model, data = emif_adni.data[emif_adni.data$SEX == 1, ], estimator = "WLSMV")
summary(mediation_male.fit)
summary(mediation_female.fit)

