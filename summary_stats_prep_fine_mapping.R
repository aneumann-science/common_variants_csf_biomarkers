library(data.table)

# Map data
map.data <- fread("../HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz")
map.data <- map.data[,c("#CHROM","POS","ID")]; gc()
lift_over.data <- fread("../1000G_p3v5.TOPMed_Imputed.add_UKB.Liftover_hg19to38.txt.gz")

lift_map.data <- merge(lift_over.data, map.data, by.x = c("chr37","pos37"), by.y = c("#CHROM","POS"))
map.data <- NULL; gc()
lift_over.data <- NULL; gc()

lift_map.data$chr37 <- as.numeric(lift_map.data$chr37)
lift_map.data$chr38 <- as.numeric(lift_map.data$chr38)

# PC1 
pc1.data <- fread("EMIF_AD_CSF_PCA_GWAS/GWAS_results/all_meta-analysis/PC1_meta_all1.meta.gz")
pc1.data$N <- 973
pc1.data <- pc1.data[,c("CHROM","POS","A1","A2","BETA","P","N")]
pc1.data$CHROM <- as.numeric(pc1.data$CHROM)

pc1.data <- merge(lift_map.data, pc1.data, by.x = c("chr37", "pos37"),  by.y = c("CHROM", "POS"))
pc1.data <- pc1.data[,.(ID,chr38,pos38,N,P,A1,A2,BETA)]
names(pc1.data) <- c("rsID","CHR","BP","N","P","A1","A2","BETA")
fwrite(pc1.data, file = "summary_stats/pc1.txt", sep = "\t")

pc1.data <- NULL; gc()

# PC2
pc2.data <- fread("EMIF_AD_CSF_PCA_GWAS/GWAS_results/all_meta-analysis/PC2_meta_all1.meta.gz")
pc2.data$N <- 973
pc2.data <- pc2.data[,c("CHROM","POS","A1","A2","BETA","P","N")]
pc2.data$CHROM <- as.numeric(pc2.data$CHROM)

pc2.data <- merge(lift_map.data, pc2.data, by.x = c("chr37", "pos37"),  by.y = c("CHROM", "POS"))
pc2.data <- pc2.data[,.(ID,chr38,pos38,N,P,A1,A2,BETA)]
names(pc2.data) <- c("rsID","CHR","BP","N","P","A1","A2","BETA")
fwrite(pc2.data, file = "summary_stats/pc2.txt", sep = "\t")

pc2.data <- NULL; gc()

# PC3
pc3.data <- fread("EMIF_AD_CSF_PCA_GWAS/GWAS_results/all_meta-analysis/PC3_meta_all1.meta.gz")
pc3.data$N <- 973
pc3.data <- pc3.data[,c("CHROM","POS","A1","A2","BETA","P","N")]
pc3.data$CHROM <- as.numeric(pc3.data$CHROM)

pc3.data <- merge(lift_map.data, pc3.data, by.x = c("chr37", "pos37"),  by.y = c("CHROM", "POS"))
pc3.data <- pc3.data[,.(ID,chr38,pos38,N,P,A1,A2,BETA)]
names(pc3.data) <- c("rsID","CHR","BP","N","P","A1","A2","BETA")
fwrite(pc3.data, file = "summary_stats/pc3.txt", sep = "\t")

pc3.data <- NULL; gc()

# PC4
pc4.data <- fread("EMIF_AD_CSF_PCA_GWAS/GWAS_results/all_meta-analysis/PC4_meta_all1.meta.gz")
pc4.data$N <- 973
pc4.data <- pc4.data[,c("CHROM","POS","A1","A2","BETA","P","N")]
pc4.data$CHROM <- as.numeric(pc4.data$CHROM)

pc4.data <- merge(lift_map.data, pc4.data, by.x = c("chr37", "pos37"),  by.y = c("CHROM", "POS"))
pc4.data <- pc4.data[,.(ID,chr38,pos38,N,P,A1,A2,BETA)]
names(pc4.data) <- c("rsID","CHR","BP","N","P","A1","A2","BETA")
fwrite(pc4.data, file = "summary_stats/pc4.txt", sep = "\t")

pc4.data <- NULL; gc()

# PC5
pc5.data <- fread("EMIF_AD_CSF_PCA_GWAS/GWAS_results/all_meta-analysis/PC5_meta_all1.meta.gz")
pc5.data$N <- 973
pc5.data <- pc5.data[,c("CHROM","POS","A1","A2","BETA","P","N")]
pc5.data$CHROM <- as.numeric(pc5.data$CHROM)

pc5.data <- merge(lift_map.data, pc5.data, by.x = c("chr37", "pos37"),  by.y = c("CHROM", "POS"))
pc5.data <- pc5.data[,.(ID,chr38,pos38,N,P,A1,A2,BETA)]
names(pc5.data) <- c("rsID","CHR","BP","N","P","A1","A2","BETA")
fwrite(pc5.data, file = "summary_stats/pc5.txt", sep = "\t")

pc5.data <- NULL; gc()

# PC5 replication
pc5_replication.data <- fread("EMIF_AD_CSF_PCA_GWAS/replication/GWAS_results/RC5_GWAS.RC5.glm.linear.gz")
pc5_replication.data$N <- 786
pc5_replication.data <- pc5_replication.data[,c("#CHROM","POS","REF","ALT","BETA","P","N")]
pc5_replication.data$`#CHROM` <- as.numeric(pc5_replication.data$`#CHROM`)

pc5_replication.data <- merge(lift_map.data, pc5_replication.data, by.x = c("chr38", "pos38"),  by.y = c("#CHROM", "POS"))
pc5_replication.data <- pc5_replication.data[,.(ID,chr38,pos38,N,P,REF,ALT,BETA)]
names(pc5_replication.data) <- c("rsID","CHR","BP","N","P","A1","A2","BETA")
fwrite(pc5_replication.data, file = "summary_stats/pc5_replication.txt", sep = "\t")

pc5_replication.data <- NULL; gc()

### Wightman AD GWAS
wightman.data <- fread("summary_stats/PGCALZ2sumstatsExcluding23andMe.txt.gz")
wightman.data <- wightman.data[,c("chr","PosGRCh37","testedAllele","otherAllele","z","p","N")]

wightman.data <- merge(lift_map.data, wightman.data, by.x = c("chr37", "pos37"),  by.y = c("chr", "PosGRCh37"))
wightman.data <- wightman.data[,.(ID,chr38,pos38,N,p,testedAllele,otherAllele,z)]
names(wightman.data) <- c("SNP","CHR","BP","N","P","A1","A2","Z")
fwrite(wightman.data, file = "summary_stats/wightman.txt", sep = "\t")

wightman.data <- NULL; gc()

# Educational attainment
ea.data <- fread("summary_stats/EA4_additive_excl_23andMe.txt.gz")
ea.data$N <- 765283
ea.data <- ea.data[,c("rsID","Chr","BP","Effect_allele","Other_allele","Beta","P_unadj","N")]

ea.data <- merge(lift_map.data, ea.data, by.x = c("chr37", "pos37"),  by.y = c("Chr", "BP"))
ea.data <- ea.data[,.(rsID,chr38,pos38,N,P_unadj,Effect_allele,Other_allele,Beta)]
names(ea.data) <- c("rsID","CHR","BP","N","P","A1","A2","BETA")
fwrite(ea.data, file = "summary_stats/ea.txt", sep = "\t")

ea.data <- NULL; gc()

# Cognitive ability
cognitve.data <- fread("summary_stats/Davies2018_OPEN_DATASET_summary_results.txt")
cognitve.data$N <- 282014
cognitve.data <- cognitve.data[,c("Markername","Chromosome","Position","Effect_allele","Other_allele","Zscore","P-value","N")]

cognitve.data <- merge(lift_map.data, cognitve.data, by.x = c("chr37", "pos37"),  by.y = c("Chromosome", "Position"))
cognitve.data <- cognitve.data[,c("Markername","chr38","pos38","N","P-value","Effect_allele","Other_allele","Zscore")]
names(cognitve.data) <- c("rsID","CHR","BP","N","P","A1","A2","Z")
fwrite(cognitve.data, file = "summary_stats/cognitive.txt", sep = "\t")

cognitive.data <- NULL; gc()

# MDD
mdd.data <- fread("summary_stats/PGC_UKB_depression_genome-wide.txt")
mdd.data$N <- 807553
mdd.data <- mdd.data[,c("MarkerName","A1","A2","LogOR","P","N")]

mdd.data <- merge(lift_map.data, mdd.data, by.x = "ID",  by.y = "MarkerName")
mdd.data <- mdd.data[,.(ID,chr38,pos38,N,P,A1,A2,LogOR)]
names(mdd.data) <- c("rsID","CHR","BP","N","P","A1","A2","BETA")
fwrite(mdd.data, file = "summary_stats/mdd.txt", sep = "\t")

mdd.data <- NULL; gc()