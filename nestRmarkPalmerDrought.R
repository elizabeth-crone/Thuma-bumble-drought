rm(list = ls(all= T)) #clears environment 
library(RMark)
library(car)
library(gplots)

setwd("")
MarkPath = "C:/Program Files/"
dat<-read.csv("all years capture history.csv", colClasses = c("factor", "factor", "factor", "integer", "integer", "integer", "integer", "integer", "character")) #import capture histories and keep the leading zeroes

table(dat$ch)
summary(dat)
bees = dat[, c(1,9)]
head(bees)
bees.processed <- process.data(bees, model = "POPAN", groups = c("year"), begin.time = 0, time.intervals = rep(1,4))

bees.processed$group.covariates
bees.ddl <- make.design.data(bees.processed)
bees.ddl

# we are following an example online
# https://www.montana.edu/rotella/documents/502/lab05Rmark.html  (thanks, Jay Rotella!)
# palmer drought severity https://www.ncei.noaa.gov/access/monitoring/weekly-palmers/time-series/1903
# we are using Palmer Drought Index (not hydrological) for Jun-Aug +/- 1 month from survey start
# we used coastal MA, not west or central
precip.dat = c(-0.7467, 2.8143, -0.7807, 0.6504, -2.0800, 0.5623, 1.3987) #PDI averaged over Jun-Aug each year

#create precipitation variable associated with p (detection probability)
bees.ddl$p$precip = precip.dat[1] 
bees.ddl$p$precip[bees.ddl$p$year == "2019"] = precip.dat[2]  
bees.ddl$p$precip[bees.ddl$p$year == "2020"] = precip.dat[3]  
bees.ddl$p$precip[bees.ddl$p$year == "2021"] = precip.dat[4] 
bees.ddl$p$precip[bees.ddl$p$year == "2022"] = precip.dat[5]
bees.ddl$p$precip[bees.ddl$p$year == "2023"] = precip.dat[6]  
bees.ddl$p$precip[bees.ddl$p$year == "2024"] = precip.dat[7]

#create precipitation variable associated with Phi (survival)
bees.ddl$Phi$precip = precip.dat[1] 
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2019"] = precip.dat[2]  
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2020"] = precip.dat[3]  
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2021"] = precip.dat[4] 
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2022"] = precip.dat[5] 
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2023"] = precip.dat[6]  
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2024"] = precip.dat[7] 

# this next bit of code sets the detection probability for the fifth survey to 0 for 2018 and 2019
# there were only four surveys in this year
bees.ddl$p$fix = NA 
bees.ddl$p$fix[bees.ddl$p$year == "2018" & bees.ddl$p$time == 4] = 0
bees.ddl$p$fix[bees.ddl$p$year == "2019" & bees.ddl$p$time == 4] = 0

bees.ddl

################################################
# simple models ############################
################################################

# the next lines of code define possible predictors of survival (phi) and detection (p) probability
Phi.dot = list(formula = ~1)
Phi.year = list(formula = ~0+year)
Phi.precip = list(formula = ~precip)
Phi.1 = list(formula = ~1, fixed = 1)

p.dot = list(formula = ~1)
p.year = list(formula = ~year)
p.precip = list(formula = ~precip)

# "N" in this model is the number of missed nests in each year, not the actual population size
# it is not sensible to have this be the same every year
N.year = list(formula = ~0+ year)

# testing whether nests "entered" the study during the surveys
pent.0 = list(formula = ~1, fixed = 0) # no nests entered
pent.varies = list(formula = ~ time*year) # sometimes nests entered

# this code creates all possible combinations of possible predictors of survival, abundance and detection probability
bees.model.list = create.model.list("POPAN")
bees.model.list

# set 'silent = F' to print out all the models as they are run
bees.results = mark.wrapper(bees.model.list, data = bees.processed, ddl = bees.ddl, silent = T) 


bees.results # the top model has constant survival, constant detection rates and no new nests appearing during the study

# getting estimates for each year
bees.results$Phi.year.p.dot.pent.0.N.year$results$real

# interpreting the top model
bees.results$Phi.dot.p.dot.pent.0.N.year$results$real
plogis(1.88) # average nest survival = 0.87
plogis(-0.46) # average detection probability = 0.39
popan.derived(bees.processed, bees.results$Phi.dot.p.dot.pent.0.N.year)$NGross # this is the total population size
# 29, 36, 10, 10, 2, 24, 24 (total population each year, 2018-2024)
popan.derived(bees.processed, bees.results$Phi.dot.p.dot.pent.0.N.year)$NGross.vcv # variance-covariance of total population size

# post-hoc test of population size vs. precipitation
Ns = popan.derived(bees.processed, bees.results$Phi.dot.p.dot.pent.0.N.year)$NGross 
Ns.SE = sqrt(diag(popan.derived(bees.processed, bees.results$Phi.dot.p.dot.pent.0.N.year)$NGross.vcv))

# calculating observed vs. actual number of nests
N.obs = c(23,28,9.2,8.8,2,20.2,19.8)
plotCI(N.obs, Ns, uiw = 1.96*Ns.SE, xlim = c(0,30), ylim = c(0,45))
abline(a = 0, b = 1)
m.obs = glm(N.obs ~ 0+ Ns)
summary(m.obs)
cor(Ns, N.obs)
library(msm)
deltamethod(~x1+x2+x3+x4+x5+x6+x7, 
            mean = popan.derived(bees.processed, bees.results$Phi.dot.p.dot.pent.0.N.year)$NGross,
            cov = popan.derived(bees.processed, bees.results$Phi.dot.p.dot.pent.0.N.year)$NGross.vcv,
            ses = T)

# pdf = 2.83 x 3.N.obs# pdf = 2.83 x 3.70
plot(precip.dat, Ns)
m.N = lm(Ns ~ precip.dat, weights = 1/(1+Ns.SE)) #weights the one precise estimate from 2022
summary(m.N)
Anova(m.N) #significant effect of precipitation on population size/abundance
# F value = 10.99, df=1, p=0.0211

# this relationship is weaker without weighting the one very precise estimate from 2022
m.N.no_wts = lm(Ns ~ precip.dat)
summary(m.N.no_wts)
Anova(m.N.no_wts) #now pvalue = 0.0611

PDI = seq(-2.1, 2.9, 0.1)
preds = predict(m.N, newdata = data.frame(precip.dat = PDI), se.fit = T)

plot(PDI, preds$fit, type = "l", lwd = 4, xlab = "", ylab = "", ylim = c(0,45))
points(PDI, preds$fit+1.96*preds$se.fit, type = "l")
points(PDI, preds$fit-1.96*preds$se.fit, type = "l")
points(precip.dat, Ns, pch = 19, cex = sqrt(3/(1+Ns.SE)), col = "goldenrod")
mtext(side = 1, line = 2, "Drought (Palmer Index)")
mtext(side = 2, line = 2, "# nests found")
#nest abundance increases with increasing precipitation/higher Palmer Index value

# NB; survival is lowest in 2020 - the year with a wet spring and really dry Jul-Aug - possibly for discussion
bees.results$Phi.year.p.dot.pent.0.N.year$results$real

####################################
#  analyses comparing gris and imp #
####################################


rm(list = ls(all= T)) #clears environment - this is a handy way to get rid of all the exploratory models
library(RMark)
library(car)
library(gplots)

setwd("")
MarkPath = "C:/Program Files/"

dat<-read.csv("all years capture history.csv", colClasses = c("factor", "factor", "factor", "integer", "integer", "integer", "integer", "integer", "character")) #import capture histories and keep the leading zeroes
# need to add columns for habitat type and bee species

table(dat$ch)
summary(dat)
#bees = dat[dat$year != "2024", 2:5]
bees = dat[dat$Species %in% c("Griseocollis", "Impatiens"), c(1,2,9)]
head(bees)
bees.processed <- process.data(bees, model = "POPAN", groups = c("year", "Species"), begin.time = 0, time.intervals = rep(1,4))

bees.processed$group.covariates
bees.ddl <- make.design.data(bees.processed)
bees.ddl

# this seems like there should be a better way to do it, but I'm following an example online
# https://www.montana.edu/rotella/documents/502/lab05Rmark.html  (thanks, Jay Rotella!)
# palmer drought severity https://www.ncei.noaa.gov/access/monitoring/weekly-palmers/time-series/1903
# we are using Palmer Drought Index (not hydrological) for Jun-Aug +/1 1 month from survey start
# we used coastal MA, not west or central
precip.dat = c(-0.7467, 2.8143, -0.7807, 0.6504, -2.0800, 0.5623, 1.3987) #PDI averaged over Jun-Aug each year

#create precipitation variable associated with p (detection probability)
bees.ddl$p$precip = precip.dat[1] 
bees.ddl$p$precip[bees.ddl$p$year == "2019"] = precip.dat[2]  
bees.ddl$p$precip[bees.ddl$p$year == "2020"] = precip.dat[3]  
bees.ddl$p$precip[bees.ddl$p$year == "2021"] = precip.dat[4] 
bees.ddl$p$precip[bees.ddl$p$year == "2022"] = precip.dat[5]
bees.ddl$p$precip[bees.ddl$p$year == "2023"] = precip.dat[6]  
bees.ddl$p$precip[bees.ddl$p$year == "2024"] = precip.dat[7]

#create precipitation variable associated with Phi (survival)
bees.ddl$Phi$precip = precip.dat[1] 
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2019"] = precip.dat[2]  
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2020"] = precip.dat[3]  
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2021"] = precip.dat[4] 
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2022"] = precip.dat[5] 
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2023"] = precip.dat[6]  
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2024"] = precip.dat[7] 


bees.ddl$N$precip = precip.dat[1] 
bees.ddl$N$precip[bees.ddl$N$year == "2019"] = precip.dat[2]  
bees.ddl$N$precip[bees.ddl$N$year == "2020"] = precip.dat[3]  
bees.ddl$N$precip[bees.ddl$N$year == "2021"] = precip.dat[4] 
bees.ddl$N$precip[bees.ddl$N$year == "2022"] = precip.dat[5] 
bees.ddl$N$precip[bees.ddl$N$year == "2023"] = precip.dat[6]  
bees.ddl$N$precip[bees.ddl$N$year == "2024"] = precip.dat[7] 


# this next bit of code sets the capture probability for the fifth survey to 0 for 2018 and 2019
# there were only four surveys in this year
bees.ddl$p$fix = NA 
bees.ddl$p$fix[bees.ddl$p$year == "2018" & bees.ddl$p$time == 4] = 0
bees.ddl$p$fix[bees.ddl$p$year == "2019" & bees.ddl$p$time == 4] = 0

bees.ddl

################################################
# simple models with species####################
################################################

# the next lines of code define possible predictors of survival and capture probability
Phi.dot = list(formula = ~1)
Phi.year = list(formula = ~0+year)
Phi.SP = list(formula = ~0+Species)
Phi.precip = list(formula = ~precip)
Phi.yearSP = list(formula = ~0+year + Species)
Phi.precipSP = list(formula = ~0+precip + Species)
Phi.precipxSP = list(formula = ~0+precip*Species)

p.dot = list(formula = ~1)
p.year = list(formula = ~Species)


# "N" in this model is the number of missed nests in each year, not the actual population size
# it is not sensible to have this be the same every year
N.year = list(formula = ~0+ group)

# testing assume no nests "entered" the study during the surveys
pent.0 = list(formula = ~1, fixed = 0) # no nests entered

# this code creates all possible combinations of possible predictors of survival, abundance and capture probability
bees.model.list = create.model.list("POPAN")
bees.model.list

# set 'silent = F' to print out all the models as they are run
bees.results = mark.wrapper(bees.model.list, data = bees.processed, ddl = bees.ddl, silent = T) 


bees.results # the top model (for the data I used) has constant survival, constant detection rates and no new nests appearing during the study
#best model also has species as predictor but not precipitation or year
#so there are differences between species but those differences do not vary with year or precipitation

# interpreting the top model
bees.results$Phi.SP.p.dot.pent.0.N.year$results$beta
bees.results$Phi.SP.p.dot.pent.0.N.year$results$real
plogis(0.7933) # average nest survival for Gris = 0.69
plogis(17.2) # average nest survival for Imp = 1
plogis(-0.346) # detection probability (p) = 0.42
popan.derived(bees.processed, bees.results$Phi.SP.p.dot.pent.0.N.year)$NGross # this is the total population size
#13 groups, N estimates listed in groups order (run just bees.results$Phi.SP.p.dot.pent.0.N.year to see full model output)
#1 gris 2018: 8 
#2 gris 2019: 17
#3 gris 2020: 2
#4 gris 2021: 5
#5 gris 2023: 12
#6 gris 2024: 18
#7 imp 2018: 16
#8 imp 2019: 15
#9 imp 2020: 4
#10 imp 2021: 5
#11 imp 2022: 2
#12 imp 2023: 11
#13 imp 2024: 7
#gris counts 2018-2024:  8, 17, 2, 5, 0, 12, 18
#imp  counts 2018-2024: 16, 15, 4, 5, 2, 11, 7
popan.derived(bees.processed, bees.results$Phi.SP.p.dot.pent.0.N.year)$NGross.vcv # variance-covariance of total population size

# post-hoc test of population size vs. precipitation
Ns = popan.derived(bees.processed, bees.results$Phi.dot.p.dot.pent.0.N.year)$NGross 
Ns.SE = sqrt(diag(popan.derived(bees.processed, bees.results$Phi.dot.p.dot.pent.0.N.year)$NGross.vcv))
SPs = bees.ddl$N$Species
precip = bees.ddl$N$precip
#hack to add no Gris in 2022
Ns = c(Ns, 0)
Ns.SE = c(Ns.SE, Ns.SE[11]) # setting SE for Gris to be the same as SE for Imp in 2022
SPs = c(SPs, SPs[1])
precip = c(precip, precip[11])

# comparing the estimated number of nests in each year
use.precip = precip[7:13]
m.N.SP2 = glm(cbind(round(Ns[c(1:6,14)]), round(Ns[7:13]))~use.precip, family = binomial)
summary(m.N.SP2)
Anova(m.N.SP2) #so precipitation has no effect on abundance of each species relative to each other, chis = 0.7055, df=1, p=0.4009

PDI = seq(-2.1, 2.9, 0.1)
preds = predict(m.N.SP2, se.fit = T, newdata = data.frame(use.precip = PDI))
plot(PDI, 100*plogis(preds$fit), xlab = "", ylab = "", ylim = c(0,100), type = "l", lwd = 4, col = "gray")
points(PDI, 100*plogis(preds$fit+1.96*preds$se.fit), type = "l", col = "gray")
points(PDI, 100*plogis(preds$fit-1.96*preds$se.fit), type = "l", col = "gray")
points(precip.dat, 100*(Ns[c(1:6,14)])/(Ns[c(1:6,14)]+Ns[7:13]), cex = sqrt((Ns[c(1:6,14)] + Ns[7:13])/10), pch = 19, lwd = 2, col = "purple1")
mtext(side = 1, line = 2, "Drought (Palmer Index)")
mtext(side = 2, line = 2, "% B. griseocollis")

m.N.SP1 = glm(Ns ~ SPs*precip, weights = 1/(Ns.SE + 1))
summary(m.N.SP1)
Anova(m.N.SP1) #so precipitation has no effect on abundance of each species relative to each other, chis = 0.7055, df=1, p=0.4009
m.N.SP1a = glm(Ns ~ SPs+precip, weights = 1/(Ns.SE + 1))
summary(m.N.SP1a)

#species generally found in different numbers from each other, but no trend seen between years or based on drought index
PDI = seq(-2.1, 2.9, 0.1)
predsI = predict(m.N.SP1, newdata = data.frame(precip = PDI, SPs = "Impatiens"), se.fit = T)
plot(PDI, predsI$fit, type = "l", lwd = 4, xlab = "Drought (Palmer Index)", ylab = "# nests found", ylim = c(0,20), col = "goldenrod")
points(PDI, predsI$fit+1.96*predsI$se.fit, type = "l", col = "goldenrod")
points(PDI, predsI$fit-1.96*predsI$se.fit, type = "l", col = "goldenrod")
points(precip[7:13], Ns[7:13], pch = 19, col = "goldenrod")

predsG = predict(m.N.SP1, newdata = data.frame(precip = PDI, SPs = "Griseocollis"), se.fit = T)
points(PDI, predsG$fit, type = "l", lwd = 4, col = "purple4")
points(PDI, predsG$fit+1.96*predsG$se.fit, type = "l", col = "purple4")
points(PDI, predsG$fit-1.96*predsG$se.fit, type = "l", col = "purple4")
points(precip[c(1:6,14)], Ns[c(1:6,14)], pch = 19, col = "purple4")

####################################
#  analyses comparing forests and meadows #
####################################

rm(list = ls(all= T)) #clears environment 
library(RMark)
library(car)
library(gplots)

setwd("")
MarkPath = "C:/Program Files/"
dat<-read.csv("all years capture history.csv", colClasses = c("factor", "factor", "factor", "integer", "integer", "integer", "integer", "integer", "character")) #import capture histories and keep the leading zeroes

table(dat$ch)
summary(dat)
bees = dat[dat$Species =="Impatiens", c(1,3,9)]
head(bees)
summary(bees)
bees$Habitat.Type[bees$Habitat.Type == "Hayfield"] = "Meadow"
bees.processed <- process.data(bees, model = "POPAN", groups = c("year", "Habitat.Type"), begin.time = 0, time.intervals = rep(1,4))

bees.processed$group.covariates
bees.ddl <- make.design.data(bees.processed)
bees.ddl

# we are following an example online
# https://www.montana.edu/rotella/documents/502/lab05Rmark.html  (thanks, Jay Rotella!)
# palmer drought severity https://www.ncei.noaa.gov/access/monitoring/weekly-palmers/time-series/1903
# we are using Palmer Drought Index (not hydrological) for Jun-Aug +/- 1 month from survey start
# we used coastal MA, not west or central
precip.dat = c(-0.7467, 2.8143, -0.7807, 0.6504, -2.0800, 0.5623, 1.3987) #PDI averaged over Jun-Aug each year

#create precipitation variable associated with p (detection probability)
bees.ddl$p$precip = precip.dat[1] 
bees.ddl$p$precip[bees.ddl$p$year == "2019"] = precip.dat[2]  
bees.ddl$p$precip[bees.ddl$p$year == "2020"] = precip.dat[3]  
bees.ddl$p$precip[bees.ddl$p$year == "2021"] = precip.dat[4] 
bees.ddl$p$precip[bees.ddl$p$year == "2022"] = precip.dat[5]
bees.ddl$p$precip[bees.ddl$p$year == "2023"] = precip.dat[6]  
bees.ddl$p$precip[bees.ddl$p$year == "2024"] = precip.dat[7]

#create precipitation variable associated with Phi (survival)
bees.ddl$Phi$precip = precip.dat[1] 
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2019"] = precip.dat[2]  
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2020"] = precip.dat[3]  
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2021"] = precip.dat[4] 
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2022"] = precip.dat[5] 
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2023"] = precip.dat[6]  
bees.ddl$Phi$precip[bees.ddl$Phi$year == "2024"] = precip.dat[7] 

# this next bit of code sets the detection probability for the fifth survey to 0 for 2018 and 2019
# there were only four surveys in this year
bees.ddl$p$fix = NA 
bees.ddl$p$fix[bees.ddl$p$year == "2018" & bees.ddl$p$time == 4] = 0
bees.ddl$p$fix[bees.ddl$p$year == "2019" & bees.ddl$p$time == 4] = 0

bees.ddl

################################################
# simple models with habitat####################
################################################

# the next lines of code define possible predictors of survival and detection probability
Phi.dot = list(formula = ~1)
Phi.year = list(formula = ~0+year)
Phi.SP = list(formula = ~0+Habitat.Type)
Phi.precip = list(formula = ~precip)
Phi.yearSP = list(formula = ~0+year + Habitat.Type)
Phi.precipSP = list(formula = ~0+precip + Habitat.Type)
Phi.precipxSP = list(formula = ~0+precip*Habitat.Type)

p.dot = list(formula = ~1)
p.hab = list(formula = ~0+Habitat.Type)


# "N" in this model is the number of missed nests in each year, not the actual population size
# it is not sensible to have this be the same every year
N.sep = list(formula = ~0+ group)

# testing assume no nests "entered" the study during the surveys
pent.0 = list(formula = ~1, fixed = 0) # no nests entered

# this code creates all possible combinations of possible predictors of survival, abundance and detection probability
bees.model.list = create.model.list("POPAN")
bees.model.list

# set 'silent = F' to print out all the models as they are run
bees.results = mark.wrapper(bees.model.list, data = bees.processed, ddl = bees.ddl, silent = T) 


bees.results # the top model has constant survival, constant detection rates and no new nests appearing during the study
#has habitat as predictor but not year or PDSI

# interpreting the top model
bees.results$Phi.dot.p.hab.pent.0.N.sep$results$beta
bees.results$Phi.dot.p.hab.pent.0.N.sep$results$real

plogis(5.9) # average nest survival (for Imp) = 0.997
plogis(-0.7453803) # detection probability in Forests = 0.32
plogis(-0.0360526) # detection probability in Meadows = 0.49
popan.derived(bees.processed, bees.results$Phi.dot.p.hab.pent.0.N.sep)$NGross # this is the total population size
#forest nest counts 2018-2024: 10, 3, 1, 2, 1, 2, 2
#meadow nest counts 2018-2024: 6, 11, 3, 3, 1, 9, 5
popan.derived(bees.processed, bees.results$Phi.dot.p.hab.pent.0.N.sep)$NGross.vcv # variance-covariance of total population size

# post-hoc test of population size vs. drought
Ns = round(popan.derived(bees.processed, bees.results$Phi.dot.p.hab.pent.0.N.sep)$NGross)
precip = precip.dat
habitat = bees.ddl$N$Habitat.Type

# comparing the estimated number of nests in each year
m.N.hab = glm(cbind(round(Ns[c(1:7)]), round(Ns[8:14]))~precip.dat, family = binomial) # modeling proportion of nests in forests
summary(m.N.hab) 
plogis(-0.32) 
Anova(m.N.hab) #significant effect of precipitation/drought index on nesting chis=4.79, df=1, p=0.0286

# how different would it be if we used raw data
N.obs.M = c(6,11,3,3,1,7,5)
N.obs.F = c(9,3,1,2,1,2,2)
m.N.hab.obs = glm(cbind(N.obs.F, N.obs.M)~precip.dat, family = binomial) # modeling proportion of nests in forests
summary(m.N.hab.obs) 
Anova(m.N.hab.obs) #significant effect of precipitation/drought index on nesting chis=4.79, df=1, p=0.0286


PDI = seq(-2.1, 2.9, 0.01)
preds = predict(m.N.hab, se.fit = T, newdata = data.frame(precip.dat = PDI))
plot(PDI, 100*plogis(preds$fit), xlab = "", ylab = "", ylim = c(0,100), type = "l", lwd = 4)
points(PDI, 100*plogis(preds$fit+1.96*preds$se.fit), type = "l")
points(PDI, 100*plogis(preds$fit-1.96*preds$se.fit), type = "l")
points(precip.dat, 100*(Ns[c(1:7)])/(Ns[8:14]+Ns[c(1:7)]), pch = 19, cex = sqrt(Ns[8:14]+Ns[c(1:7)])/2, col = "darkolivegreen4")
mtext(side = 1, line = 2, "drought (Palmer Index)")
mtext(side = 2, line = 2, "% in forest")
#see more impatiens nesting in forest when drought index is lower/there is more drought


#####################3
# post-hoc test - did the proportion of rare species decrease after 2020 drought?
######################

num.rare = c(2, 2, 2, 0,0,0,0)
num.nests = c(23,28,8,9,2,20,20) # excludes one unidentified nest in 2020
after.drought = c(0,0,0,1,1,1,1)

m.rare = glm(cbind(num.rare, num.nests-num.rare)~after.drought, family = binomial)
summary(m.rare)
Anova(m.rare)

#####################3
# figure of other mark-recapture estimates
######################

phis = c(0.868, 0.689, 1.00)
phis.low = c(0.776, 0.561, 0.999)
phis.upp = c(0.926, 0.793, 1.000)

ps = c(0.388, 0.322, 0.491)
ps.low = c(0.326, 0.218, 0.401)
ps.upp = c(0.453, 0.447, 0.582)


library(plotrix)
xvals = c(1, 1.7, 2.3)
plotCI(xvals, phis, ui = phis.upp, li = phis.low, xlim = c(0.7, 2.5), pch = 21, cex = 1.25, pt.bg = "steelblue", xaxt = "n", xlab = "", ylab = "")
axis(side = 1, at = xvals, labels = c("all", "gris.", "imp."), cex = 0.8)
mtext(side = 1, line = 2, "Bombus species")
mtext(side = 2, line = 3, "Survival during survey period")
mtext(side = 2, line = 2, "(probability per week)")

plotCI(xvals, ps, ui = ps.upp, li = ps.low, xlim = c(0.7, 2.5), pch = 21, cex = 1.25, pt.bg = "steelblue", xaxt = "n", xlab = "", ylab = "")
axis(side = 1, at = xvals, labels = c("all", "forest", "open"), cex = 0.8)
axis(side = 1, at = xvals, labels = c("", "", "open"), cex = 0.8)
mtext(side = 1, line = 2, "Habitat type")
mtext(side = 2, line = 3, "Detection probability")
mtext(side = 2, line = 2, "(per survey)")
