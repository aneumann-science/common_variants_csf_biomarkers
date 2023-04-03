###### ADNI
library(data.table)
library(psych)

### Load ADNI datasets and PCA
load("phenotypes/adni/adnimerge.rdata")
load("phenotypes/adni/blennowcsfnfl.rdata")
load("phenotypes/adni/blennowcsfng.rdata")
ykl40.data <- read.csv("phenotypes/adni/csfproteomics/adni_csfproteomics.csv")

### Prepare ADNI merge file
# Reduce to baseline measurements only
adnimerge <- adnimerge[adnimerge$VISCODE == "bl", ]

# Replace detection limit with numeric values 
adnimerge$TAU[adnimerge$TAU == "<80"] <- 80
adnimerge$TAU[adnimerge$TAU == ">1300"] <- 1300
adnimerge$TAU <- as.numeric(adnimerge$TAU)
hist(adnimerge$TAU)

adnimerge$PTAU[adnimerge$PTAU == "<8"] <- 8
adnimerge$PTAU[adnimerge$PTAU == ">120"] <- 120
adnimerge$PTAU <- as.numeric(adnimerge$PTAU)
hist(adnimerge$PTAU)

adnimerge$ABETA[adnimerge$ABETA == "<200"] <- 200
adnimerge$ABETA[adnimerge$ABETA == ">1700"] <- 1700
adnimerge$ABETA <- as.numeric(adnimerge$ABETA)
hist(adnimerge$ABETA)

# Convert gender to a dummy variable
adnimerge$GENDER <- model.matrix(~ PTGENDER, data = adnimerge)[,-1]
adnimerge$GENDER <- 1-adnimerge$GENDER

# Convert diagnosis to dummy variables, merge EMCI and LMCI, treat SMC as CN
options(na.action='na.pass')
adnimerge[c("SMC","EMCI","LMCI","AD")] <- model.matrix(~ DX.bl, data = adnimerge)[,-1]
adnimerge$MCI <- adnimerge$LMCI
adnimerge$MCI[adnimerge$EMCI == 1] <- 1

# Select relevant variables
adnimerge <- adnimerge[c("RID","PTID","TAU","PTAU","ABETA","MMSE","AGE","GENDER","MCI","AD")]

### Prepare NFL
# Check comments
table(blennowcsfnfl$COMMENTS)
# Check duplicates
table(duplicated(blennowcsfnfl$RID))
# Average duplicate values
blennowcsfnfl <- aggregate(CSFNFL ~ RID, blennowcsfnfl, mean)

### Prepare neurogranin
# Check comments
table(blennowcsfng$GOTCOMMENT)

# Assign lowest value to observations below detection value
# First check min value without 0
min(blennowcsfng$CSFNG[blennowcsfng$CSFNG!= 0])
blennowcsfng$CSFNG[blennowcsfng$GOTCOMMENT== "Very low level, below detection"] <- 5.774937
# Check duplicates
table(duplicated(blennowcsfng$CSFNG))
# Average duplicate values
blennowcsfng <- aggregate(CSFNG ~ RID, blennowcsfng, mean)

### Prepare YKL-40 data
# Only baseline
ykl40.data <- ykl40.data[ykl40.data$VISCODE2 == "bl", ]
# YKL-40 names
ykl40_names <- c("CH3L1_HUMAN.ILGQQVPYATK.Y5.","CH3L1_HUMAN.ILGQQVPYATK.Y9.","CH3L1_HUMAN.VTIDSSYDIAK.Y7.","CH3L1_HUMAN.VTIDSSYDIAK.Y8.")

# Standardize the measures and average
ykl40.data[ykl40_names] <- lapply(ykl40.data[ykl40_names], scale)
ykl40.data$YKL40 <- rowMeans(ykl40.data[ykl40_names], na.rm = T)

# Check correlations to see whether everything is ok
cor(ykl40.data[c(ykl40_names,"YKL40")], use = "complete.obs")

# Reduce to relevant variables
ykl40.data <- ykl40.data[c("RID","YKL40")]

# Merge all data
adnimerge_nfl <- merge(adnimerge, blennowcsfnfl, by = "RID", all = T)
adnimerge_nfl_ng <- merge(adnimerge_nfl, blennowcsfng, by = "RID", all = T)
adni.data <- merge(adnimerge_nfl_ng, ykl40.data, by = "RID", all = T)

# Check which ID have a wrong lower case letter in the VCF file
adni.fam <- fread("../adni_emif_genotypes/adni/adni_rarevariants.fam")
wrong_id <- adni.fam[grep("s", adni.fam[,"V2"]), "V2"]

# Check which row numbers in the adnimerge have the non-matching ID
wrong_id_row <- grep(paste(wrong_id,collapse="|"), adni.data$PTID, ignore.case=T)

# Replace S with s, so that it matches with the VCF files
adni.data$PTID[wrong_id_row] <- gsub("S", "s", adni.data[wrong_id_row, "PTID"])

### CSF biomarkers
# Define sets of variables, which will be useful to refer to later
csf_phenotypes <- c("TAU","PTAU","CSFNFL","ABETA","CSFNG","YKL40")

# Extract CSF phenotypes and reduce to observations with at least three biomarkers
csf_adni.data <- adni.data[rowMeans(!is.na(adni.data[c("TAU","CSFNFL","ABETA","CSFNG","YKL40")])) >= 4/5, c("PTID", csf_phenotypes)]


# INT transformation of the CSF biomarkers
csf_adni.data[csf_phenotypes] <- lapply(csf_adni.data[csf_phenotypes], function(x) {
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
})
multi.hist(csf_adni.data[csf_phenotypes])

###### EMIF
# Load EMIF phenotype data
emif_ad.data <- read.csv("phenotypes/emif/emif_ad_full.csv")

# Only baseline
emif_ad.data <- emif_ad.data[emif_ad.data$Fu == 0, ]

### CSF biomarkers
# Define sets of variables, which will be useful to refer to later
csf_phenotypes <- c("Ptau_ASSAY_Zscore","Ttau_ASSAY_Zscore","Central_CSF_NFL","Central_CSF_AB42","Central_CSF_Neurogranin","Central_CSF_YKL40")

# Extract CSF phenotypes and reduce to observations with at least four biomarkers
csf_emif.data <- emif_ad.data[rowMeans(!is.na(emif_ad.data[c("Ptau_ASSAY_Zscore","Central_CSF_NFL","Central_CSF_AB42","Central_CSF_Neurogranin","Central_CSF_YKL40")])) >= 4/5, c("SubjectId", csf_phenotypes)]

# Check data and distribution
str(csf_emif.data)
describe(csf_emif.data)
multi.hist(csf_emif.data[csf_phenotypes])

# INT transformation of the CSF biomarkers
csf_emif.data[csf_phenotypes] <- lapply(csf_emif.data[csf_phenotypes], function(x) {
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
})
multi.hist(csf_emif.data[csf_phenotypes])


###### PCA
library(missMDA)
library(psych)
library(dplyr)

# Add study identifier
csf_adni.data$study <- "adni"
csf_emif.data$study <- "emif"

# Reduce to relevant variables
adni_biomarkers <- c("TAU","PTAU","CSFNFL","ABETA","YKL40","CSFNG")
emif_biomarkers <- c("Ttau_ASSAY_Zscore","Ptau_ASSAY_Zscore","Central_CSF_NFL","Central_CSF_AB42","Central_CSF_YKL40","Central_CSF_Neurogranin")

adni_biomarkers.data <- csf_adni.data[c("study","PTID",adni_biomarkers)]
emif_biomarkers.data <- csf_emif.data[c("study","SubjectId",emif_biomarkers)]

# Harmonize CSF biomarker names
names(adni_biomarkers.data) <- c("study","ID",adni_biomarkers)
names(emif_biomarkers.data) <- c("study","ID",adni_biomarkers)

# Make direction in EMIF for Tau more conventional
emif_biomarkers.data$TAU <- -1*emif_biomarkers.data$TAU
emif_biomarkers.data$PTAU <- -1*emif_biomarkers.data$PTAU

# Check correlations
cor(adni_biomarkers.data[adni_biomarkers], use = "pairwise.complete.obs")
cor(emif_biomarkers.data[adni_biomarkers], use = "pairwise.complete.obs")

# Merge both studies
adni_emif.data <- as.data.frame(rbind(adni_biomarkers.data, emif_biomarkers.data))
adni_emif.data$study <- as.factor(adni_emif.data$study)

# Descriptives
multi.hist(adni_emif.data[adni_biomarkers])
describe(adni_emif.data)

# Use cross-validation to determine optimum number of componencts
set.seed(20200617)
nb <- estim_ncpPCA(adni_emif.data[adni_biomarkers], method.cv = "loo")
nb$ncp
png(file="figures/pca/elbow_plot_full.png")
plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")
dev.off()

# Perform PCA based imputation
set.seed(20200617)
adni_emif_int.imputePCA <- imputePCA(adni_emif.data[adni_biomarkers], ncp = 5)
adni_emif_int_imputed.data <- as.data.frame(adni_emif_int.imputePCA$completeObs)
adni_emif_int_imputed.data$study <- adni_emif.data$study
adni_emif_int_imputed.data$ID <- adni_emif.data$ID

# PCA
pca.fit <- pca(adni_emif_int_imputed.data[adni_biomarkers], nfactors = 5)
# Save loadings
write.csv(unclass(pca.fit$loadings), "results/pca_biomarkers_int_adni_emif_full.csv")
save(pca.fit, file = "results/pca.Rdata")
# compute PC scores
pca_adni_emif.data <- data.frame(predict(pca.fit, data = adni_emif_int_imputed.data[adni_biomarkers]))
pca_adni_emif.data$study <- adni_emif.data$study
pca_adni_emif.data$ID <- adni_emif.data$ID
# Check PC distribution
multi.hist(pca_adni_emif.data[c("RC1","RC3","RC2","RC4","RC5")])

# check data
describe(pca_adni_emif.data)

# Save PC data
save(pca_adni_emif.data, file = "phenotypes/joint/pca_adni_emif_all.Rdata")

### Heatmap
library(corrplot)

loadings.mat <- t(unclass(pca.fit$loadings))
rownames(loadings.mat) <- c("Tau pathology/degeneration","Injury/inflammation","Aβ pathology","Non-AD inflammation","Non-AD synaptic functioning")
colnames(loadings.mat) <- c("Tau","pTau","NFL","Amyloid","YKL-40","Ng")

col3 <- colorRampPalette(c("red", "white", "blue")) 

png(file="figures/pca/heatmap.png", width = 4000, height = 3000, pointsize = 80)
corrplot(loadings.mat, tl.srt = 0, col = col3(20), method = "color")
dev.off()

loadings_cor.mat <- cbind(loadings.mat,c(-0.14,-0.20,0.28,-0.02,0.01))
colnames(loadings_cor.mat) <- c("Tau","pTau","NFL","Amyloid","YKL-40","Ng","MMSE cor")

png(file="figures/pca/heatmap_cor.png", width = 4000, height = 3000, pointsize = 80)
corrplot(loadings_cor.mat, tl.srt = 0, col = col3(20), method = "color")
dev.off()

### Female
# Female ids
adni_female <- adni.data[adni.data$GENDER == 1,"PTID"]
emif_female <- emif_ad.data[emif_ad.data$Gender == 1,"SubjectId"]

# PCA
adni_emif_int_imputed_female.data <- adni_emif_int_imputed.data[adni_emif_int_imputed.data$ID %in% c(adni_female,emif_female), ]
pca_female.fit <- pca(adni_emif_int_imputed_female.data[adni_biomarkers], nfactors = 5)
# Save loadings
write.csv(unclass(pca_female.fit$loadings), "results/pca_biomarkers_int_adni_emif_full_female.csv")

### Male
# Male ids
adni_male <- adni.data[adni.data$GENDER == 0,"PTID"]
emif_male <- emif_ad.data[emif_ad.data$Gender == 0,"SubjectId"]

# PCA
adni_emif_int_imputed_male.data <- adni_emif_int_imputed.data[adni_emif_int_imputed.data$ID %in% c(adni_male,emif_male), ]
pca_male.fit <- pca(adni_emif_int_imputed_male.data[adni_biomarkers], nfactors = 5)
# Save loadings
write.csv(unclass(pca_male.fit$loadings), "results/pca_biomarkers_int_adni_emif_full_male.csv")

### Check overlap with genetic data in EMIF
gwas_id <- fread("phenotypes/joint/gwas_list.txt", header = F, col.names = "gwas_id")
gwas_id_adni <- gwas_id$gwas_id[1:1106] 
gwas_id_emif <- gwas_id$gwas_id[1107:2001] 

# Convert to SubjectID
gwas_id_emif <- merge(gwas_id_emif, emif_ad.data[c("SubjectId","DNA_ID")], by.x = "x", by.y = "DNA_ID")
names(gwas_id_emif) <- c("DNA_ID", "SubjectId")

# Check overlap
pca_adni_emif_gwas.data <- pca_adni_emif.data[pca_adni_emif.data$ID %in% gwas_id_adni | pca_adni_emif.data$ID %in% gwas_id_emif$SubjectId, ]
dim(pca_adni_emif_gwas.data[pca_adni_emif_gwas.data$study == "adni", ])
dim(pca_adni_emif_gwas.data[pca_adni_emif_gwas.data$study == "emif", ])

### Merge with covariate information
## ADNI
pca_adni.data <- pca_adni_emif.data[pca_adni_emif.data$study == "adni", ]
cov_adni.data <- adni.data[rowMeans(!is.na(adni.data[c("TAU","CSFNFL","ABETA","CSFNG","YKL40")])) >= 4/5, c("PTID","AGE","GENDER","AD","MCI")] 
pca_cov_adni.data <- merge(pca_adni.data, cov_adni.data, by.x = "ID", by.y = "PTID")

# Sanity checks
describe(pca_cov_adni.data)
round(cor(pca_cov_adni.data[c("RC1","RC3","RC2","RC4","RC5","MCI","AD","AGE","GENDER")]),2)

# Save
save(pca_cov_adni.data, file = "phenotypes/GWAS/pca_cov_adni.Rdata")

## EMIF
pca_emif.data <- pca_adni_emif.data[pca_adni_emif.data$study == "emif", ] 
cov_emif.data <- emif_ad.data[rowMeans(!is.na(emif_ad.data[c("Ptau_ASSAY_Zscore","Central_CSF_NFL","Central_CSF_AB42","Central_CSF_Neurogranin","Central_CSF_YKL40")])) >= 4/5, ] 

# Convert diagnosis to dummy variables, treat SMC as CN
options(na.action='na.pass')
cov_emif.data[c("AD_Baseline_Diagnosis","MCI_Baseline_Diagnosis","NL_Baseline_Diagnosis","SCI_Baseline_Diagnosis")] <- model.matrix(~ Diagnosis, data = cov_emif.data)[,-1]

# Compute age at CSF measurement
cov_emif.data$ageatcsfpet <- rowSums(cbind(cov_emif.data$Age, cov_emif.data$DateCSFCollection), na.rm = T)
cov_emif.data <- cov_emif.data[rowMeans(!is.na(cov_emif.data[c("Ptau_ASSAY_Zscore","Central_CSF_NFL","Central_CSF_AB42","Central_CSF_Neurogranin","Central_CSF_YKL40")])) >= 4/5, c("SubjectId","DNA_ID","ageatcsfpet","Gender","AD_Baseline_Diagnosis","MCI_Baseline_Diagnosis")]
pca_cov_emif.data <- merge(pca_emif.data, cov_emif.data, by.x = "ID",  by.y = "SubjectId")

# Sanity checks
describe(pca_cov_emif.data)
round(cor(pca_cov_emif.data[c("RC1","RC3","RC2","RC4","RC5","MCI_Baseline_Diagnosis","AD_Baseline_Diagnosis","ageatcsfpet","Gender")]),2)

# Save
save(pca_cov_emif.data, file = "phenotypes/GWAS/pca_cov_emif.Rdata")

### Save file for imputed/INT transformed single biomarker
# ADNI
pca_cov_single_adni.data <- merge(pca_cov_adni.data, adni_emif_int_imputed.data[c("ID",adni_biomarkers)], by = "ID")
save(pca_cov_single_adni.data, file = "phenotypes/GWAS/pca_cov_single_adni.Rdata")

# EMIF
pca_cov_single_emif.data <- merge(pca_cov_emif.data, adni_emif_int_imputed.data[c("ID",adni_biomarkers)], by = "ID")
save(pca_cov_single_emif.data, file = "phenotypes/GWAS/pca_cov_single_emif.Rdata")

                                      
