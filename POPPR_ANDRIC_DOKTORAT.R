# Uvoz library POPPR
library("poppr")

# Provjera radnog direktorija i eventualno postavljanje odgovarajuceg
getwd()

# Uvoz podataka
fraxinus <- read.genalex("POPPR.csv")


# The genotype accumulation curve
gac <- genotype_curve(fraxinus, sample = 1000, quiet = TRUE)

# Allele frequencies, missing data, and ploidy
data("fraxinus")
(fraxinuslt <- locus_table(fraxinus)) # frekvencija alela

info_table(fraxinus, type = "missing", plot = TRUE) # missing data plot

tail(genind2df(fraxinus, sep = "/")) # transforms the data from genind object into a dataframe and adds '/' as separators between alleles at each locus

fraxinus.ploidy <- info_table(fraxinus, type = "ploidy", plot = TRUE, low = "black", high = "orange")

tail(fraxinus.ploidy)

# Genotypic richness, diversity, and evenness

data("fraxinus")

monpop_diversity <- poppr(fraxinus)

library("vegan")
mon.tab <- mlg.table(fraxinus, plot = FALSE)
min_sample <- min(rowSums(mon.tab))
rarecurve(mon.tab, sample = min_sample, xlab = "Sample Size", ylab = "Expected MLGs")
title("Rarefaction of Fruit Rot and Blossom Blight")

N      <- monpop_diversity$N      # number of samples
lambda <- monpop_diversity$lambda # Simpson's index
(N/(N - 1)) * lambda              # Corrected Simpson's index

mon.tab <- mlg.table(fraxinus)

# The index of association

library("poppr")
library("magrittr")

MX <- popsub(fraxinus, "Population I")
ia(MX, sample = 999)

SA <- popsub(fraxinus, "Population II")
ia(SA, sample = 999)

# Population structure: GST, genetic distance, and clustering
library("mmod")
Gst_Hedrick(fraxinus) # GST an example with Felis catus data

library("poppr") # Genetic Distance
library("ape") # To visualize the tree using the "nj" function
library("magrittr")
data(fraxinus)
set.seed(10)
nt_samples <- sample(nInd(fraxinus), 92)
frax92       <- fraxinus[nt_samples]
(fraxdist    <- provesti.dist(frax92))

theTree <- fraxdist %>% # neighbor-joining tree
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
plot(theTree)
add.scale.bar(length = 0.05) # add a scale bar showing 5% difference

set.seed(999)
aboot(frax92, dist = provesti.dist, sample = 200, tree = "nj", cutoff = 50, quiet = TRUE)

# K-means hierarchical clustering
library("poppr") # Cazma
CA <- popsub(fraxinus, "Population I")
CAclust <- find.clusters(CA)

CAclust

library("poppr") # Nova Gradiska
NG <- popsub(fraxinus, "Population II")
NGclust <- find.clusters(NG)

NGclust






