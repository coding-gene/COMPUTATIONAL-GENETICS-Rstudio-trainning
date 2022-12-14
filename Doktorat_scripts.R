# POSTAVLJANJE RADNOG DIREKTORIJA
setwd("C:\\Users\\Ivan\\OneDrive\\Desktop\\R_WORKING_SPACE\\RSTUDIO IMPORT")
getwd() # provjera radnog direktorija

# UCITAVANJE PAKETA
library("adegenet")

# UVOZ PODATAKA
FRAX <- read.genetix("Genetix.gtx")

# KONVERTIRANJE IZ genind U genpop FORMAT
data(FRAX)
frax_genpop <- genind2genpop(FRAX)

# BASIC DATA ANALYSIS
frax_genpop <- summary(FRAX)
names(frax_genpop)

par(mfrow=c(2,2))

plot(frax_genpop$n.by.pop, frax_genpop$pop.n.all, xlab="Colonies sample size",
     ylab="Number of alleles",main="Alleles numbers and sample sizes",
     type="n")

text(frax_genpop$n.by.pop,frax_genpop$pop.n.all,lab=names(frax_genpop$n.by.pop))

barplot(frax_genpop$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")

barplot(frax_genpop$Hexp-frax_genpop$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")

barplot(frax_genpop$n.by.pop, main="Sample sizes per population", 
        ylab="Number of genotypes",las=3)

# HARDY-WEINBERG ZAKON RAVNOTEZE
library(pegas)
data(FRAX)
frax.hwt <- hw.test(FRAX, B=0)
frax.hwt

# F-STATISTICS
library("hierfstat") # Po populacijama
fstat(FRAX)

library(pegas) # Po lokusima
Fst(as.loci(FRAX))

Gtest <- gstat.randtest(FRAX,nsim=99)
Gtest
plot(Gtest)

matFst <- pairwise.fst(FRAX[1:50,])
matFst

# INBREEDING


# PCA ANALYSIS
data(FRAX)
sum(is.na(FRAX$tab))

X <- tab(FRAX, freq = TRUE, NA.method = "mean")
class(X)

dim(X)

pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

s.label(pca1$li)
title("PCA of FRAX datasetnnaxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(FRAX))
title("PCA of Fraxinus angustifolia CSO clones")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li,pop(FRAX),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA of Fraxinus angustifolia CSO clones")
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)

col <- funky(15)
s.class(pca1$li, pop(FRAX),xax=1,yax=3, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)

colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA of Fraxinus angustifolia datasetnnaxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

colorplot(pca1$li[c(1,3)], pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA of Fraxinus angustifolia datasetnnaxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

# BAYESIAN
getwd()



