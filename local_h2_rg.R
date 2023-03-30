# Load LAVA for local heritability and genetic correlation analyses
library(LAVA)

### Read in summary statistics and related info
input = process.input(input.info.file="local_rg/input.info.txt",           # Input files
                      sample.overlap.file="local_rg/sample.overlap.txt",   # Sample overlap
                      ref.prefix="local_rg/g1000_eur",
                      phenos=c("pc1","pc2","pc3","pc4","pc5"))             

### Read in locus info file
loci = read.loci("local_rg/hits.loci")

### APOE
APOE = process.locus(loci[2,], input)
# Local heritability
APOE.univ <- run.univ(APOE)
write.csv(APOE.univ, file = "results/local_h2/discovery/APOE_h2.csv")
# Local genetic correlation
APOE.biv <- run.bivar(APOE, target = "pc2")
write.csv(APOE.univ, file = "results/local_rg/discovery/APOE_rg.csv")

### TMEM106B
TMEM106B = process.locus(loci[3,], input)
# Local heritability
TMEM106B.univ <- run.univ(TMEM106B)
write.csv(TMEM106B.univ, file = "results/local_h2/discovery/TMEM106B_h2.csv")
# Local genetic correlation
TMEM106B.biv <- run.bivar(TMEM106B, target = "pc3")
write.csv(TMEM106B.biv, file = "results/local_rg/discovery/TMEM106B_rg.csv")

### rs145791381
rs145791381 = process.locus(loci[4,], input)
# Local heritability
rs145791381.univ <- run.univ(rs145791381)
write.csv(rs145791381.univ, file = "results/local_h2/discovery/rs145791381_h2.csv")
# Local genetic correlation
rs145791381.biv <- run.bivar(rs145791381, target = "pc4")
write.csv(rs145791381.biv, file = "results/local_rg/discovery/rs145791381_rg.csv")

### CHI3L1
CHI3L1 = process.locus(loci[1,], input)
# Local heritability
CHI3L1.univ <- run.univ(CHI3L1)
write.csv(CHI3L1.univ, file = "results/local_h2/discovery/CHI3L1_h2.csv")
# Local genetic correlation
CHI3L1.biv <- run.bivar(CHI3L1, target = "pc4")
write.csv(CHI3L1.biv, file = "results/local_rg/discovery/CHI3L1_rg.csv")

### GRIN2D
GRIN2D = process.locus(loci[5,], input)
# Local heritability
GRIN2D.univ <- run.univ(GRIN2D)
write.csv(GRIN2D.univ, file = "results/local_h2/discovery/GRIN2D_h2.csv")
# Local genetic correlation
GRIN2D.biv <- run.bivar(GRIN2D, target = "pc5")
write.csv(GRIN2D.biv, file = "results/local_rg/discovery/GRIN2D_rg.csv")


