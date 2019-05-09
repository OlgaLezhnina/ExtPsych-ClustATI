# set wd
setwd("~/R")
# get packages
# for dealing with the data
library(tidyverse)
# for protocol6
library(psych)
library(ltm)
library(MASS)
library(lavaan)
library(mirt)
library(eRm)
library(mokken)
library(MBESS)
# for additional factor analysis
library(nFactors)
# for additional multivariate normality tests
library(QuantPsyc)
library(MVN)
# for cluster analysis
library(ClustOfVar)

## Get all data and subsets
# import the data
mturk_data <- read_csv("~/R/mturk_data.csv")
# add rfs mean
rfs.data <- dplyr::select(mturk_data, starts_with("rfs"))
mturk_data$rfs_mean <- rowMeans(rfs.data)
all.data <- mturk_data
# select the ATI scale without the mean column
ati.data <- dplyr::select(all.data, starts_with("ati"), -"ati_mean")
# select items of all scales, without means
scales.data <- dplyr::select(all.data, -c("id", ends_with("date"), starts_with("education"), "lastpage", "startlanguage", "age", "sex", ends_with("mean"), starts_with("bfi10_"), starts_with("ip")))
# select means of all scales for cluster analysis
means.data <- dplyr::select(mturk_data, ends_with("mean"), starts_with("bfi10_"))

# explore missing data
sum(is.na(scales.data))
# no missing data in any scales

## Start six steps of analysis for ATI scale according to Dima (2018) https://doi.org/10.1080/21642850.2018.1472602

# define a function to partition an item set into mokken scales - lowerbound from .05 to .80
moscales.for.lowerbounds <- function(x, lowerbounds = seq(from = 0.05, to = 0.80, by = 0.05)) {
  ret.value <- NULL
  for (lowerbound in lowerbounds)
  {
    tmp <- aisp(x, lowerbound = lowerbound)
    if (is.null(ret.value)) {
      ret.value <- data.frame("Item" = rownames(tmp), "Scales." = tmp[, 1])
    }
    else {
      ret.value <- cbind(ret.value, "Scales." = tmp[, 1])
    }
    names(ret.value)[ncol(ret.value)] <- sprintf("%.2f", lowerbound)
  }
  rownames(ret.value) <- NULL
  ret.value
}

# Step 1
# get descriptives
descrmyitems <- as.data.frame(round(psych::describe(ati.data), 2))
# barplots for items
ati.items <- names(ati.data)
par(mfrow = c(3, 3))
for (n in ati.items)
{
  distr <- table(ati.data[, n])
  barplot(distr,
    main = n,
    ylim = c(0, 80),
    col = gray.colors(20),
    ylab = "Number of respondents",
    xlab = "Response (1 to 6)"
  )
}
# corelation matrix
bluesqs <- cor(ati.data, method = "pearson")
psych::cor.plot(bluesqs,
  numbers = TRUE, main = "correlations between items",
  cex = 0.5, cex.axis = 0.7
)
# multivariate outliers
par(mfrow = c(1, 1))
d2mydata <- outlier(ati.data, cex = .6, bad = 3, ylim = c(0, 50))
hist(d2mydata)
# number of outliers
sum((1 - pchisq(d2mydata, ncol(ati.data))) < .001)
# which exactly
1 - pchisq(d2mydata, ncol(ati.data)) < .001 # 16 and 142
# max D2 value of the outliers
round(max(d2mydata), 2)

## My own checking multivariate normality, not from Dima (2018) https://doi.org/10.1080/21642850.2018.1472602

# check multivariate normality with QuantPsyc
mvnorm.ati1 <- QuantPsyc::mult.norm(ati.data)
mvnorm.ati1$mult.test
# outliers
mvnorm.ati1$Dsq > mvnorm.ati1$CriticalDsq
sum(mvnorm.ati1$Dsq > mvnorm.ati1$CriticalDsq) # 5 outliers
# check with MVN
mvnorm.ati2 <- MVN::mvn(ati.data, multivariatePlot = "qq")
mvnorm.ati2
# Return to Dima (2018) https://doi.org/10.1080/21642850.2018.1472602

## Step 2
ati.dataN <- as.data.frame(ati.data)
# Guttman errors outliers
xPlus <- rowSums(ati.dataN)
gPlus <- mokken::check.errors(ati.dataN)$Gplus
hist(gPlus)
# get the outliers
Q3 <- summary(gPlus)[[5]]
IQR <- Q3 - summary(gPlus)[[2]]
outlier <- gPlus > Q3 + 1.5 * IQR
# which exactly is the cut off value
Q3 + 1.5 * IQR
# see outliers
outlier
sum(outlier)

# homoheneity
Hvalues <- coefH(ati.dataN)
# homogeneity for items
Hvalues$Hi
# homogeneity for the complete set
Hvalues$H
# examine aisp for increasing c levels (run the function defined above)
motable.mydata <- moscales.for.lowerbounds(ati.dataN)
# save it as a data frame
aispmydata <- as.data.frame(motable.mydata)

## from now on, ATI and ATI8 will be explored
ati.data8 <- as.data.frame(ati.data[, -3])

# homogeneity for ATI8
Hvalues8 <- coefH(ati.data8)
# homogeneity for items
Hvalues8$Hi
# homogeneity for the complete set
Hvalues8$H

# check conditional association (local independence) for ATI
CA.def.atiscale <- check.ca(ati.dataN, TRUE)
CA.def.atiscale$InScale[[1]]

# check monotonicity with default minsize for ATI
monotonicity.def.atiscale <- check.monotonicity(ati.dataN, minvi = .03)
monotonicity.def.atiscale$results
summary(monotonicity.def.atiscale)
plot(monotonicity.def.atiscale)

# check monotonicity with default minsize for ATI8
monotonicity.def.atiscale8 <- check.monotonicity(ati.data8, minvi = .03)
summary(monotonicity.def.atiscale8)
plot(monotonicity.def.atiscale8)

# Investigate the assumption of non-intersecting item step response functions (ISRFs)
# MIIO for ATI
miio.atiscale <- check.iio(ati.dataN)
miio.atiscale
miio.atiscale$items.removed
summary(miio.atiscale)$item.summary

# MIIO for ATI8
miio.atiscale8 <- check.iio(ati.data8)
miio.atiscale8$items.removed
summary(miio.atiscale8)$item.summary

# Step 3
# Rating Scale model for ATI
fit.atiscale <- eRm::RSM(ati.dataN)
summary(fit.atiscale)

# Rating Scale model for ATI8
fit.ati8scale <- RSM(ati.data8)
summary(fit.ati8scale)

# person parameters for ATI
ppr <- person.parameter(fit.atiscale)
# person parameters for ATI8
ppr8 <- person.parameter(fit.ati8scale)

# Information criteria for ATI and ATI8
IC(ppr)
IC(ppr8)

# item-pair residuals for testing local dependencies
fit.ltm.atiscale <- ltm::grm(ati.dataN, constrained = TRUE)
fit.ltm.ati8scale <- grm(ati.data8, constrained = TRUE)

# model summary
summary(fit.ltm.atiscale)
summary(fit.ltm.ati8scale)

# Fit on the Two-Way Margins
margins(fit.ltm.atiscale)
margins(fit.ltm.ati8scale)

# plot Item Characteristic Curves, all nine plots together
par(mfrow = c(3, 3))
plotICC(fit.atiscale, legend = T, legpos = T)

## My own way of plotting ICCs, not from Dima (2018) https://doi.org/10.1080/21642850.2018.1472602

# use package mirt to get all plots together
mod <- mirt(ati.dataN, 1)
plot(mod, type = "trace")
plot(mod, type = "infotrace")

# and separate item plots, e.g.:
plot(mod, type = "trace", which.items = 1, facet_items = FALSE)
plot(mod, type = "trace", which.items = 3, facet_items = FALSE)

## Returnt to Dima(2018) https://doi.org/10.1080/21642850.2018.1472602

# plot person-item map
par(mfrow = c(1, 1))
plotPImap(fit.atiscale, sorted = TRUE)

# plot item difficulty & infit statistics
plotPWmap(fit.atiscale)

# separation reliability for ATI and ATI8
round(as.numeric(SepRel(ppr))[1], 2)
round(as.numeric(SepRel(ppr8))[1], 2)

# person separation for ATI and ATI8
round(sqrt(SepRel(ppr)$sep.rel / (1 - SepRel(ppr)$sep.rel)), 2)
round(sqrt(SepRel(ppr8)$sep.rel / (1 - SepRel(ppr8)$sep.rel)), 2)

# item fit (infit and outfit) for ATI
itemfit.ati <- itemfit(ppr)
itemfit.ati.print <- as.data.frame(print(itemfit.ati))

# Personfit
personfit.atiscale <- personfit(ppr)
PersonFitTBL <- as.data.frame(print(personfit.atiscale))
# misfitting persons
misoutfits <- nrow(PersonFitTBL[(PersonFitTBL$p.outfitZ > 1.96 | PersonFitTBL$p.outfitZ < -1.96)
& (PersonFitTBL$p.outfitMSQ < 0.6 | PersonFitTBL$p.outfitMSQ > 1.4), ])
misinfits <- nrow(PersonFitTBL[(PersonFitTBL$p.infitZ > 1.96 | PersonFitTBL$p.infitZ < -1.96)
& (PersonFitTBL$p.infitMSQ < 0.6 | PersonFitTBL$p.infitMSQ > 1.4), ])
# number and percentage of misfitting persons
misoutfits
round(misoutfits * 100 / nrow(PersonFitTBL), 2)
misinfits
round(misinfits * 100 / nrow(PersonFitTBL), 2)

# Step 4
## Add my own KMO and Barlett analysis with psych package, as it was missing in Dima (2018) https://doi.org/10.1080/21642850.2018.1472602
# KMO and Barlett
print(KMO(ati.data))
print(cortest.bartlett(ati.data))
# extracting factors with nFactors package
ev <- eigen(cor(ati.data))
ap <- parallel(subject = nrow(ati.data), var = ncol(ati.data), rep = 1000, cent = .05)
nS <- nScree(x = ev$values, aparallel = ap$eigen$qevpea)
plotnScree(nS)

## Return to Dima (2018) https://doi.org/10.1080/21642850.2018.1472602
# very simple structure analysis
vss(ati.data, 5)

# factor solution
fa(ati.data, nfactors = 1, fm = "minres", n.iter = 10)

# ICLUST (groups items) for ATI
summary(iclust(ati.dataN, title = "ICLUST for ATI"))

# CFA for ATI
CFA.ati <- "
Factor1 =~ ati01 + ati02 + ati03R + ati04 + ati05 + ati06R + ati07 + ati08R + ati09"
fitCFA.ati <- lavaan::cfa(CFA.ati, data = ati.data)
summary(fitCFA.ati, standardized = TRUE, fit.measures = TRUE)

## My own addition, see modification indices
modificationIndices(fitCFA.ati, sort. = TRUE, power = TRUE)

## Return to Dima (2018) https://doi.org/10.1080/21642850.2018.1472602
# Step5
# CTT item properties for ATI
psych::alpha(ati.data)
# for ATI8
psych::alpha(ati.data8)

## A recommended splithalf test is used, as there was a warning for guttman(ati.data):
# "Guttman has been deprecated.  The use of the splitHalf function is recommended"
splitHalf(ati.data) # for ATI
splitHalf(ati.data8) # for ATI8

# omega with CI for ATI
omega.ati <- MBESS::ci.reliability(
  data = ati.data, type = "omega", conf.level = 0.95,
  interval.type = "perc", B = 100
)
omega.ati
# omega with CI for ATI8
omega.ati8 <- ci.reliability(
  data = ati.data8, type = "omega", conf.level = 0.95,
  interval.type = "perc", B = 100
)
omega.ati8

## STEP6
# Item frequencies etc for ATI
mydata.scores <- scoreItems(ati.items, ati.data)
print(mydata.scores, short = FALSE)
# Item frequencies etc for ati8
ati.items8 <- ati.items[-3]
mydata.scores8 <- scoreItems(ati.items8, ati.data8)
print(mydata.scores8, short = FALSE)

# descriptives for ATI
psych::describe(ati.data)
psych::alpha(ati.data)

# descriptives for ATI8
psych::describe(ati.data8)
psych::alpha(ati.data8)

# histogram for ATI
ati.mean <- rowMeans(ati.data)
hist(ati.mean)
# histogram for ATI8
ati.mean8 <- rowMeans(ati.data8)
hist(ati.mean8)

# FINISH Dima (2918) https://doi.org/10.1080/21642850.2018.1472602

## Hierarchical clustering on variables with ClustOfVar package
# transform data (scales and means) into matrix
meansM.data <- data.matrix(means.data)
scalesM.data <- data.matrix(scales.data)
# dendogram of means
tree.means <- hclustvar(meansM.data)
plot(tree.means)
# dendogram of all items
tree.scales <- hclustvar(scalesM.data)
plot(tree.scales)
# stability of partitions, means
stab.means <- stability(tree.means)
# stability of partitions, all items
stab.scales <- stability(tree.scales)
# get partition into K clusters for means
cutree.means2 <- cutreevar(tree.means, k = 2, matsim = TRUE)
par(mfrow = c(2, 1))
plot(cutree.means2) # plots of clusters, correlation with the synthetic variable
cutree.means2$E # gain in cohesion
# for 8 clusters as shown by stability plot
cutree.means8 <- cutreevar(tree.means, k = 8, matsim = TRUE)
par(mfrow = c(4, 2))
plot(cutree.means8)
cutree.means8$E
# for 9 clusters as shown by stability plot
cutree.means9 <- cutreevar(tree.means, k = 9, matsim = TRUE)
par(mfrow = c(3, 3))
plot(cutree.means9)
cutree.means9$E
# get partition into 2 clusters for all items
cutree.scales2 <- cutreevar(tree.scales, k = 2, matsim = TRUE)
par(mfrow = c(1, 1))
plot(cutree.scales2)
cutree.scales2$E
# for next maximum of stability plot, 26 clusters
cutree.scales26 <- cutreevar(tree.scales, k = 26, matsim = TRUE)
par(mfrow = c(2, 2))
plot(cutree.scales26)
cutree.scales26$E
# DONE
