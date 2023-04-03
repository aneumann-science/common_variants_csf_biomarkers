# GWAS power analysis from https://github.com/kaustubhad/gwas-power
source("power_calc_functions.R")

# Power for genome-wide significance and a sample size of 973, dependong on
# maf and beta
pow <- power_beta_maf(beta = (1:100)/100, maf = (1:50)/100, n = 973, pval=5E-8)
write.csv(pow, file = "results/power_curve.csv")

pdf("figures/power_heatmap.pdf", width=10, height=9)
power_plot(pow, "Beta", "MAF")
dev.off()