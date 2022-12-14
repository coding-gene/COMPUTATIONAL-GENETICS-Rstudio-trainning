# POSTAVLJANJE RADNOG DIREKTORIJA
setwd("E:\\R_WORKING_SPACE\\RSTUDIO IMPORT")
getwd() # provjera radnog direktorija

# UVOZ PODATAKA
library("poppr")
PAPUK <- read.genalex("GEN.csv")
PAPUK

# The genotype accumulation curve
gac <- genotype_curve(PAPUK, sample = 1000, quiet = TRUE)

# Allele frequencies, missing data, and ploidy
data("PAPUK")
(PAPUKlt <- locus_table(PAPUK))

info_table(PAPUK, type = "missing", plot = TRUE)

tail(genind2df(PAPUK, sep = "/"))

PAPUK.ploidy <- info_table(PAPUK, type = "ploidy", plot = TRUE, low = "black", high = "orange")