### PC5
# Convert summary statistics to focus format
focus munge summary_stats/pc5_fusion.txt --output summary_stats/pc5_fusion
# Fine-mapping analysis for GRIN2D locus
focus finemap summary_stats/pc5_fusion.sumstats.gz /home/fahri/TWAS_GRCh38/LDREF_GRCh38_1KG_Phased/hwe1e6.1000G.EURn404.GRCh38_fk.chr19 weights/DatabaseForFOCUS_FULL_heritable.ROSMAP_DLPFC_n560_eQTL.31January2023.db --chr 19 --locations weights/GRIN2D_locus_10k.bed --out results/GRIN2D_locus_10k --p-threshold 1E-4 --prior-prob gencode38 --plot

### AD (EADB)
focus finemap /home/fahri/EADB_TWAS_Manuscript/FineMapping_and_ConditionalAnalyses/FOCUS_FineMapping/EADBv3_MetaProxy_SummaryForFOCUS_rsIDversion_19August2020.tsv.gz /home/fahri/TWAS_GRCh38/LDREF_GRCh38_1KG_Phased/hwe1e6.1000G.EURn404.GRCh38_fk.chr19 weights/DatabaseForFOCUS_FULL_heritable.ROSMAP_DLPFC_n560_eQTL.31January2023.db --chr 19 --locations weights/GRIN2D_locus_10k.bed --out results/GRIN2D_locus_10k_eadb --p-threshold 1E-4 --prior-prob gencode38 --plot

### AD (Wightman)
# Convert summary statistics to focus format
focus munge summary_stats/wightman.txt --output summary_stats/wightman
# Fine-mapping analysis for GRIN2D locus
focus finemap summary_stats/wightman.sumstats.gz /home/fahri/TWAS_GRCh38/LDREF_GRCh38_1KG_Phased/hwe1e6.1000G.EURn404.GRCh38_fk.chr19 weights/DatabaseForFOCUS_FULL_heritable.ROSMAP_DLPFC_n560_eQTL.31January2023.db --chr 19 --locations weights/GRIN2D_locus_10k.bed --out results/GRIN2D_locus_10k_wightman --p-threshold 1E-2 --prior-prob gencode38 --plot

### Educational attainment
# Convert summary statistics to focus format
focus munge summary_stats/ea.txt --output summary_stats/ea
# Fine-mapping analysis for GRIN2D locus
focus finemap summary_stats/ea.sumstats.gz /home/fahri/TWAS_GRCh38/LDREF_GRCh38_1KG_Phased/hwe1e6.1000G.EURn404.GRCh38_fk.chr19 weights/DatabaseForFOCUS_FULL_heritable.ROSMAP_DLPFC_n560_eQTL.31January2023.db --chr 19 --locations weights/GRIN2D_locus_10k.bed --out results/GRIN2D_locus_10k_ea --p-threshold 1E-2 --prior-prob gencode38 --plot

### Cognitive ability
# Convert summary statistics to focus format
focus munge summary_stats/cognitive.txt --output summary_stats/cognitive
# Fine-mapping analysis for GRIN2D locus
focus finemap summary_stats/cognitive.sumstats.gz /home/fahri/TWAS_GRCh38/LDREF_GRCh38_1KG_Phased/hwe1e6.1000G.EURn404.GRCh38_fk.chr19 weights/DatabaseForFOCUS_FULL_heritable.ROSMAP_DLPFC_n560_eQTL.31January2023.db --chr 19 --locations weights/GRIN2D_locus_10k.bed --out results/GRIN2D_locus_10k_cognitive --p-threshold 1E-2 --prior-prob gencode38 --plot

### Major depressive disorder
# Convert summary statistics to focus format
focus munge summary_stats/mdd.txt --output summary_stats/mdd
# Fine-mapping analysis for GRIN2D locus
focus finemap summary_stats/mdd.sumstats.gz /home/fahri/TWAS_GRCh38/LDREF_GRCh38_1KG_Phased/hwe1e6.1000G.EURn404.GRCh38_fk.chr19 weights/DatabaseForFOCUS_FULL_heritable.ROSMAP_DLPFC_n560_eQTL.31January2023.db --chr 19 --locations weights/GRIN2D_locus_10k.bed --out results/GRIN2D_locus_10k_mdd --p-threshold 0.5 --prior-prob gencode38 --plot