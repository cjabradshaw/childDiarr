# Analysing DHS data for child health-climate relationships
# Hira Fatima & Corey Bradshaw
# September 2024

#rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(dismo)
library(gbm)
library(spatstat.random)
library(ggpubr)
library(usdm)
library(truncnorm)

# functions
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# reading the preprocessed DHS data for all countries - imputed data
u5dat <- try(read.csv("DHSclusterLevelDiarrData.csv"), silent = TRUE)
head(u5dat)
str(u5dat)

# check data distributions & formats
# continuous variables

hist(u5dat$diarp) # RESPONSE: probability of child diarrhoea
hist(u5dat$diarv)

hist(u5dat$hhmemmn) # household members
hist(u5dat$hhmemsd)

hist(u5dat$wbmimn) # woman bmi
hist(u5dat$wbmisd)

hist(u5dat$wagemn) # woman age
hist(u5dat$wagesd)

hist(u5dat$wedumn) # woman edu
hist(u5dat$wedusd)

hist(u5dat$stuntmn) # child stunting
hist(u5dat$stuntsd)

hist(u5dat$uweighmn) # child underweight
hist(u5dat$uweighsd)

hist(u5dat$wastmn) # child wasting
hist(u5dat$wastsd)

# binomial probabilities
hist(u5dat$csectp) # delivery by c section (pr(yes))
hist(u5dat$csectv)

hist(u5dat$pncp) # postnatal checks (pr(yes))
hist(u5dat$pncv)

hist(u5dat$hhsexp) # postnatal checks (pr(female))
hist(u5dat$hhsexv)

hist(u5dat$dplacep) # displaced? (pr(yes))
hist(u5dat$dplacev)

hist(u5dat$acchcp) # access to healthcare? (pr(yes))
hist(u5dat$acchcv)

hist(u5dat$regionp) # urban or regional (pr(rural))
hist(u5dat$regionv)

  # change to binomial
  u5dat$regbin <- ifelse(u5dat$regionp < 0.5, 0, 1)
  u5dat$regbin <- as.factor(u5dat$regbin) # factor
  table(u5dat$regbin)
  str(u5dat$regbin)

hist(u5dat$femp) # child gender (pr(female))
hist(u5dat$femv)
  
hist(u5dat$drugipp) # taking drugs for intestinal parasites? (pr(yes))
hist(u5dat$drugipv)

hist(u5dat$dwimprp) # improved drinking water? (pr(yes))
hist(u5dat$dwimprv)

hist(u5dat$sfimpp) # improved sanitation facilities? (pr(yes))
hist(u5dat$sfimpv)

# ordinal categorical
str(u5dat$wlthmd) # wealth index (1 -> 5)
table(u5dat$wlthmd)

  # change to ordinal
  u5dat$wlthord <- factor(u5dat$wlthmd, ordered=T)
  table(u5dat$wlthord)
  str(u5dat$wlthord)

str(u5dat$birthsmd) # birth size (1 -> 3)
table(u5dat$birthsmd)
  
  # change to ordinal
  u5dat$birthsord <- factor(u5dat$birthsmd, ordered=T)
  table(u5dat$birthsord)
  str(u5dat$birthsord)
  
## climate variables
u5dat$MATmn <- u5dat$wc2.1_30s_bio_1_mean_summary # mean annual temperature
hist(u5dat$MATmn)
u5dat$MATsd <- u5dat$wc2.1_30s_bio_1_sd_summary
hist(u5dat$MATsd)

u5dat$TARmn <- u5dat$wc2.1_30s_bio_7_mean_summary # temperature annual range
hist(u5dat$TARmn)
u5dat$TARsd <- u5dat$wc2.1_30s_bio_7_sd_summary
hist(u5dat$TARsd)

u5dat$PRCPmn <- u5dat$wc2.1_30s_bio_12_mean_summary # annual precipitation
hist(u5dat$PRCPmn)
u5dat$PRCPsd <- u5dat$wc2.1_30s_bio_12_sd_summary
hist(u5dat$PRCPsd)

u5dat$PRCPSNSmn <- u5dat$wc2.1_30s_bio_15_mean_summary # precipitation seasonality
hist(u5dat$PRCPSNSmn)
u5dat$PRCPSNSsd <- u5dat$wc2.1_30s_bio_15_sd_summary
hist(u5dat$PRCPSNSsd)

u5dat$PRCPDQmn <- u5dat$wc2.1_30s_bio_17_mean_summary # precipitation driest quarter
hist(u5dat$PRCPDQmn)
u5dat$PRCPDQsd <- u5dat$wc2.1_30s_bio_17_sd_summary
hist(u5dat$PRCPDQsd)

# create country variable
u5dat$cntry <- substr(u5dat$clust, start=1, stop=2)
table(u5dat$cntry)


############################
## boosted regression trees
############################

## PHASE 1 - socio-economic
# household members, woman education, household head gender, regional/urban, wealth, displacement status,
# access to improved drinking water, access to improved sanitation, access to healthcare
# hhmemmn (hhmemsd); wedumn (wedusd); hhsexp (hhsexv); regbin; wlthmd; dplacep (dplacev);
# dwimprp (dwimprv); sfimpp (sfimpv); acchcp (acchcv)

PHASE1vars <- u5dat %>%
  dplyr::select(cntry, clust, diarp, diarv, hhmemmn, hhmemsd, wedumn, wedusd, hhsexp, hhsexv, regbin,
                wlthmd, dplacep, dplacev, dwimprp, dwimprv, sfimpp, sfimpv, acchcp, acchcv)
head(PHASE1vars)
head(PHASE1vars[,c(5,7,9,11,12,13,15,17,19)])
dim(PHASE1vars)
PHASE1vars.naomit <- na.omit(PHASE1vars)
dim(PHASE1vars.naomit)

# variance inflation factor
usdm::vif(PHASE1vars.naomit[,c(5,7,9,13,15,17,19)])

# deterministic BRT with means only
PHASE1.brt.dtrm.mn <- gbm.step(PHASE1vars.naomit, gbm.x = attr(PHASE1vars.naomit, "names")[c(5,7,9,11,12,13,15,17,19)],
                         gbm.y = attr(PHASE1vars.naomit, "names")[3],
                         family="gaussian", max.trees=1000000, tolerance = 0.01,
                         learning.rate = 0.05, bag.fraction=0.75, tree.complexity = 2)
summary(PHASE1.brt.dtrm.mn)
gbm.plot(PHASE1.brt.dtrm.mn)
gbm.plot.fits(PHASE1.brt.dtrm.mn)

PHASE1.brt.dtrm.mn.CV.cor <- 100 * PHASE1.brt.dtrm.mn$cv.statistics$correlation.mean
PHASE1.brt.dtrm.mn.CV.cor.se <- 100 * PHASE1.brt.dtrm.mn$cv.statistics$correlation.se
print(c(PHASE1.brt.dtrm.mn.CV.cor, PHASE1.brt.dtrm.mn.CV.cor))

plot(PHASE1vars$sfimpp, PHASE1vars$diarp, pch=19, xlab="improved sanitation pr", ylab="diarrhoea pr")
plot(PHASE1vars$wedumn, PHASE1vars$diarp, pch=19, xlab="woman education", ylab="diarrhoea pr")
plot(PHASE1vars$hhmemmn, PHASE1vars$diarp, pch=19, xlab="household size", ylab="diarrhoea pr")
plot(PHASE1vars$hhsexp, PHASE1vars$diarp, pch=19, xlab="household head gender", ylab="diarrhoea pr")

## iterate resampled boosted regression trees
iter <- 100
eq.sp.points <- 100

# create storage arrays
out.colnames <- c("hhmem", "wedu", "hhsex", "dplace", "dwimpr", "sfimp", "acchc")
val.arr <- pred.arr <- array(data = NA, dim=c(eq.sp.points, length=7, iter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                             out.colnames, paste("i",1:iter, sep="")))
reg.mat <- matrix(data=NA, nrow=iter, ncol=2)
wlth.mat <- matrix(data=NA, nrow=iter, ncol=5)  

# create storage vectors
CV.cor.vec <- CV.cor.se.vec <- hhmem.ri <- hhsex.ri <- wedu.ri <- reg.ri <- wlth.ri <- dplace.ri <- dwimpr.ri <-
  sfimp.ri <- acchc.ri <- rep(NA, iter)

for (i in 1:iter) {
  
  # country-resampling loop
  cntry.vec <- attr(table(PHASE1vars.naomit$cntry), "names")
  cntry.samp <- round(min(table(PHASE1vars.naomit$cntry)) / 2, 0)
  dat.cresamp <- PHASE1vars.naomit[1,]
  
  for (c in 1:length(cntry.vec)) {
    cntry.dat <- PHASE1vars.naomit[PHASE1vars.naomit$cntry == cntry.vec[c],]
    cntry.resamp <- cntry.dat[sort(sample(1:dim(cntry.dat)[1], cntry.samp, replace=F)), ]
    dat.cresamp <- rbind(dat.cresamp, cntry.resamp)
  } # end c
  dat.cresamp <- dat.cresamp[-1,]
  
  # resample
  # response
  diar.alpha <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$alpha
  diar.beta <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$beta
  diar.stch <- rbeta(length(diar.alpha), diar.alpha, diar.beta)
  diar.stoch <- ifelse(is.na(diar.stch==T), 0, diar.stch) # diarrhoea probability
  
  # continuous
  hhmem.stoch <- truncnorm::rtruncnorm(length(dat.cresamp$hhmemmn), a=0, b=Inf, mean=dat.cresamp$hhmemmn, sd=dat.cresamp$hhmemsd)
  wedu.stoch <- truncnorm::rtruncnorm(length(dat.cresamp$wedumn), a=0, b=Inf, mean=dat.cresamp$wedumn, sd=dat.cresamp$wedusd)
  
  # binomial probabilities
  hhsex.alpha <- estBetaParams(dat.cresamp$hhsexp, dat.cresamp$hhsexv^2)$alpha
  hhsex.beta <- estBetaParams(dat.cresamp$hhsexp, dat.cresamp$hhsexv^2)$beta
  hhsex.stch <- rbeta(length(hhsex.alpha), hhsex.alpha, hhsex.beta)
  hhsex.stoch <- ifelse(is.na(hhsex.stch==T), 0, hhsex.stch) # household head gender (female) probability

  dplace.alpha <- estBetaParams(dat.cresamp$dplacep, dat.cresamp$dplacev^2)$alpha
  dplace.beta <- estBetaParams(dat.cresamp$dplacep, dat.cresamp$dplacev^2)$beta
  dplace.stch <- rbeta(length(dplace.alpha), dplace.alpha, dplace.beta)
  dplace.stoch <- ifelse(is.na(dplace.stch==T), 0, dplace.stch) # displacement probability

  dwimpr.alpha <- estBetaParams(dat.cresamp$dwimprp, dat.cresamp$dwimprv^2)$alpha
  dwimpr.beta <- estBetaParams(dat.cresamp$dwimprp, dat.cresamp$dwimprv^2)$beta
  dwimpr.stch <- rbeta(length(dwimpr.alpha), dwimpr.alpha, dwimpr.beta)
  dwimpr.stoch <- ifelse(is.na(dwimpr.stch==T), 0, dwimpr.stch) # access to improved water
  
  sfimp.alpha <- estBetaParams(dat.cresamp$sfimpp, dat.cresamp$sfimpv^2)$alpha
  sfimp.beta <- estBetaParams(dat.cresamp$sfimpp, dat.cresamp$sfimpv^2)$beta
  sfimp.stch <- rbeta(length(sfimp.alpha), sfimp.alpha, sfimp.beta)
  sfimp.stoch <- ifelse(is.na(sfimp.stch==T), 0, sfimp.stch) # access to improved sanitation

  acchc.alpha <- estBetaParams(dat.cresamp$acchcp, dat.cresamp$acchcv^2)$alpha
  acchc.beta <- estBetaParams(dat.cresamp$acchcp, dat.cresamp$acchcv^2)$beta
  acchc.stch <- rbeta(length(acchc.alpha), acchc.alpha, acchc.beta)
  acchc.stoch <- ifelse(is.na(acchc.stch==T), 0, acchc.stch) # access to healthcare
  acchc.stoch
  
  # binomial categories
  reg <- dat.cresamp$regbin
  
  # ordinal factors
  wlth.stch <- spatstat.random::rpoistrunc(length(dat.cresamp$wlthmd), lambda=dat.cresamp$wlthmd, minimum=1, method="harding")
  wlth.stch1 <- ifelse(wlth.stch > 5, 5, wlth.stch)
  wlth.stoch <- factor(wlth.stch1, ordered=T)
  
  ## collect resampled variables into single dataframe
  dat.resamp <- data.frame(diar.stoch,hhmem.stoch,wedu.stoch,hhsex.stoch,dplace.stoch,
                           dwimpr.stoch,sfimp.stoch,acchc.stoch,reg,wlth.stoch)
  colnames(dat.resamp) <- c("diar", "hhmem", "wedu", "hhsex", "dplace", "dwimpr", "sfimp", "acchc",
                            "reg", "wlth")
  # scale
  dat.resamp.cont.sc <- as.data.frame(scale(dat.resamp[,1:8], center=T, scale=T))
  dat.resamp.sc <- data.frame(dat.resamp.cont.sc, "reg"=dat.resamp$reg, "wlth"=dat.resamp$wlth)
  #head(dat.resamp.sc)
  #str(dat.resamp.sc)
  
  # boosted regression tree
  brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:10)],
                           gbm.y = attr(dat.resamp.sc, "names")[1],
                           family="gaussian", max.trees=1000000, tolerance = 0.001,
                           learning.rate = 0.001, bag.fraction=0.5, tree.complexity = 2,
                           silent=T, tolerance.method = "auto")

  if (i == 1 & is.null(brt.fit.rsmp)==F) {
    brt.fit.rsmp.old <- brt.fit.rsmp
  }
  
  if (is.null(brt.fit.rsmp) == T) {
    brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:10)],
                             gbm.y = attr(dat.resamp.sc, "names")[1],
                             family="gaussian", max.trees=1000000, tolerance = 0.0001,
                             learning.rate = 0.0001, bag.fraction=0.7, tree.complexity = 2,
                             silent=T, tolerance.method = "auto", step.size=10)
  }
  if (is.null(brt.fit.rsmp) == T) {
    brt.fit.rsmp <- brt.fit.rsmp.old
  }
  summ.fit <- summary(brt.fit.rsmp)
  
  if (is.null(brt.fit.rsmp) == F) {
    brt.fit.rsmp.old <- brt.fit.rsmp
  }
  
  # variable relative importance
  hhmem.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[1])]
  wedu.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[2])]
  hhsex.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[3])]
  dplace.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[4])]
  dwimpr.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[5])]
  sfimp.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[6])]
  acchc.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[7])]
  reg.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == "reg")]
  wlth.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == "wlth")]
  
  # goodness of fit
  CV.cor.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.mean
  CV.cor.se.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=7)
  for (p in 1:7) {
    RESP.val[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit.rsmp$var.names[1:7]
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit.rsmp$var.names[1:7]
  
  val.arr[, , i] <- as.matrix(RESP.val.dat)
  pred.arr[, , i] <- as.matrix(RESP.pred.dat)
  
  # categorical predictions
  RESP.reg <- plot.gbm(brt.fit.rsmp, i.var=8, continuous.resolution = eq.sp.points, return.grid=T)
  colnames(RESP.reg)[2] <- "pred"
  reg.mat[i,] <- RESP.reg[,2]
  
  RESP.wlth <- plot.gbm(brt.fit.rsmp, i.var=9, continuous.resolution = eq.sp.points, return.grid=T)
  colnames(RESP.wlth)[2] <- "pred"
  wlth.mat[i,] <- RESP.wlth[,2]
  
  print(i)
  
} # end i

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 7
pred.update <- pred.arr[,,1:iter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:iter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

PHASE1.boot.reg.mean <- apply(reg.mat, MARGIN=1, mean, na.rm=T)
PHASE1.boot.reg.sd <- apply(reg.mat, MARGIN=1, sd, na.rm=T)
PHASE1.reg.pred.med <- apply(reg.mat, MARGIN=1, median, na.rm=T)
PHASE1.reg.pred.lo <- apply(reg.mat, MARGIN=1, quantile, probs=0.025, na.rm=T)
PHASE1.reg.pred.up <- apply(reg.mat, MARGIN=1, quantile, probs=0.975, na.rm=T)

PHASE1.boot.wlth.mean <- apply(wlth.mat, MARGIN=1, mean, na.rm=T)
PHASE1.boot.wlth.sd <- apply(wlth.mat, MARGIN=1, sd, na.rm=T)
PHASE1.wlth.pred.med <- apply(wlth.mat, MARGIN=1, median, na.rm=T)
PHASE1.wlth.pred.lo <- apply(wlth.mat, MARGIN=1, quantile, probs=0.025, na.rm=T)
PHASE1.wlth.pred.up <- apply(wlth.mat, MARGIN=1, quantile, probs=0.975, na.rm=T)

PHASE1.pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
PHASE1.pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
PHASE1.pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
PHASE1.val.med <- apply(val.arr[,,1:iter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
PHASE1.CV.cor.update <- CV.cor.vec[1:iter]
PHASE1.CV.cor.se.update <- CV.cor.se.vec[1:iter]

hhmem.ri.update <- hhmem.ri[1:iter]
wedu.ri.update <- wedu.ri[1:iter]
hhsex.ri.update <- hhsex.ri[1:iter]
dplace.ri.update <- dplace.ri[1:iter]
dwimpr.ri.update <- dwimpr.ri[1:iter]
sfimp.ri.update <- sfimp.ri[1:iter]
acchc.ri.update <- acchc.ri[1:iter]
reg.ri.update <- reg.ri[1:iter]
wlth.ri.update <- wlth.ri[1:iter]

for (k in 1:kappa.n) {
  PHASE1.CV.cor.mean <- mean(PHASE1.CV.cor.update, na.rm=T); PHASE1.CV.cor.sd <- sd(PHASE1.CV.cor.update, na.rm=T)
  PHASE1.CV.cor.se.mean <- mean(PHASE1.CV.cor.se.update, na.rm=T); PHASE1.CV.cor.se.sd <- sd(PHASE1.CV.cor.se.update, na.rm=T)
  
  hhmem.mean <- mean(hhmem.ri.update, na.rm=T); hhmem.sd <- sd(hhmem.ri.update, na.rm=T)
  wedu.mean <- mean(wedu.ri.update, na.rm=T); wedu.sd <- sd(wedu.ri.update, na.rm=T)
  hhsex.mean <- mean(hhsex.ri.update, na.rm=T); hhsex.sd <- sd(hhsex.ri.update, na.rm=T)
  dplace.mean <- mean(dplace.ri.update, na.rm=T); dplace.sd <- sd(dplace.ri.update, na.rm=T)
  dwimpr.mean <- mean(dwimpr.ri.update, na.rm=T); dwimpr.sd <- sd(dwimpr.ri.update, na.rm=T)
  sfimp.mean <- mean(sfimp.ri.update, na.rm=T); sfimp.sd <- sd(sfimp.ri.update, na.rm=T)
  acchc.mean <- mean(acchc.ri.update, na.rm=T); acchc.sd <- sd(acchc.ri.update, na.rm=T)
  reg.mean <- mean(reg.ri.update, na.rm=T); reg.sd <- sd(reg.ri.update, na.rm=T)
  wlth.mean <- mean(wlth.ri.update, na.rm=T); wlth.sd <- sd(wlth.ri.update, na.rm=T)
  
  for (u in 1:iter) {
    PHASE1.CV.cor.update[u] <- ifelse((PHASE1.CV.cor.update[u] < (PHASE1.CV.cor.mean-kappa*PHASE1.CV.cor.sd) | PHASE1.CV.cor.update[u] >
                                  (PHASE1.CV.cor.mean+kappa*PHASE1.CV.cor.sd)), NA, PHASE1.CV.cor.update[u])
    PHASE1.CV.cor.se.update[u] <- ifelse((PHASE1.CV.cor.se.update[u] < (PHASE1.CV.cor.se.mean-kappa*PHASE1.CV.cor.se.sd) | PHASE1.CV.cor.se.update[u]
                                   > (PHASE1.CV.cor.se.mean+kappa*PHASE1.CV.cor.se.sd)), NA, PHASE1.CV.cor.se.update[u])
    
    hhmem.ri.update[u] <- ifelse((hhmem.ri.update[u] < (hhmem.mean-kappa*hhmem.sd) | hhmem.ri.update[u] > (hhmem.mean+kappa*hhmem.sd)), NA, hhmem.ri.update[u])
    wedu.ri.update[u] <- ifelse((wedu.ri.update[u] < (wedu.mean-kappa*wedu.sd) | wedu.ri.update[u] > (wedu.mean+kappa*wedu.sd)), NA, wedu.ri.update[u])
    hhsex.ri.update[u] <- ifelse((hhsex.ri.update[u] < (hhsex.mean-kappa*hhsex.sd) | hhsex.ri.update[u] > (hhsex.mean+kappa*hhsex.sd)), NA, hhsex.ri.update[u])
    dplace.ri.update[u] <- ifelse((dplace.ri.update[u] < (dplace.mean-kappa*dplace.sd) | dplace.ri.update[u] > (dplace.mean+kappa*dplace.sd)), NA, dplace.ri.update[u])
    dwimpr.ri.update[u] <- ifelse((dwimpr.ri.update[u] < (dwimpr.mean-kappa*dwimpr.sd) | dwimpr.ri.update[u] > (dwimpr.mean+kappa*dwimpr.sd)), NA, dwimpr.ri.update[u])
    sfimp.ri.update[u] <- ifelse((sfimp.ri.update[u] < (sfimp.mean-kappa*sfimp.sd) | sfimp.ri.update[u] > (sfimp.mean+kappa*sfimp.sd)), NA, sfimp.ri.update[u])
    acchc.ri.update[u] <- ifelse((acchc.ri.update[u] < (acchc.mean-kappa*acchc.sd) | acchc.ri.update[u] > (acchc.mean+kappa*acchc.sd)), NA, acchc.ri.update[u])
    reg.ri.update[u] <- ifelse((reg.ri.update[u] < (reg.mean-kappa*reg.sd) | reg.ri.update[u] > (reg.mean+kappa*reg.sd)), NA, reg.ri.update[u])
    wlth.ri.update[u] <- ifelse((wlth.ri.update[u] < (wlth.mean-kappa*wlth.sd) | wlth.ri.update[u] > (wlth.mean+kappa*wlth.sd)), NA, wlth.ri.update[u])
    
  } # end u
  
  print(k)
  
} # end k

# summaries
PHASE1.CV.cor.med <- median(PHASE1.CV.cor.update, na.rm=TRUE)
PHASE1.CV.cor.lo <- quantile(PHASE1.CV.cor.update, probs=0.025, na.rm=TRUE)
PHASE1.CV.cor.up <- quantile(PHASE1.CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(PHASE1.CV.cor.lo,PHASE1.CV.cor.med,PHASE1.CV.cor.up))

hhmem.ri.lo <- quantile(hhmem.ri.update, probs=0.025, na.rm=TRUE)
hhmem.ri.med <- median(hhmem.ri.update, na.rm=TRUE)
hhmem.ri.up <- quantile(hhmem.ri.update, probs=0.975, na.rm=TRUE)

wedu.ri.lo <- quantile(wedu.ri.update, probs=0.025, na.rm=TRUE)
wedu.ri.med <- median(wedu.ri.update, na.rm=TRUE)
wedu.ri.up <- quantile(wedu.ri.update, probs=0.975, na.rm=TRUE)

hhsex.ri.lo <- quantile(hhsex.ri.update, probs=0.025, na.rm=TRUE)
hhsex.ri.med <- median(hhsex.ri.update, na.rm=TRUE)
hhsex.ri.up <- quantile(hhsex.ri.update, probs=0.975, na.rm=TRUE)

dplace.ri.lo <- quantile(dplace.ri.update, probs=0.025, na.rm=TRUE)
dplace.ri.med <- median(dplace.ri.update, na.rm=TRUE)
dplace.ri.up <- quantile(dplace.ri.update, probs=0.975, na.rm=TRUE)

dwimpr.ri.lo <- quantile(dwimpr.ri.update, probs=0.025, na.rm=TRUE)
dwimpr.ri.med <- median(dwimpr.ri.update, na.rm=TRUE)
dwimpr.ri.up <- quantile(dwimpr.ri.update, probs=0.975, na.rm=TRUE)

sfimp.ri.lo <- quantile(sfimp.ri.update, probs=0.025, na.rm=TRUE)
sfimp.ri.med <- median(sfimp.ri.update, na.rm=TRUE)
sfimp.ri.up <- quantile(sfimp.ri.update, probs=0.975, na.rm=TRUE)

acchc.ri.lo <- quantile(acchc.ri.update, probs=0.025, na.rm=TRUE)
acchc.ri.med <- median(acchc.ri.update, na.rm=TRUE)
acchc.ri.up <- quantile(acchc.ri.update, probs=0.975, na.rm=TRUE)

reg.ri.lo <- quantile(reg.ri.update, probs=0.025, na.rm=TRUE)
reg.ri.med <- median(reg.ri.update, na.rm=TRUE)
reg.ri.up <- quantile(reg.ri.update, probs=0.975, na.rm=TRUE)

wlth.ri.lo <- quantile(wlth.ri.update, probs=0.025, na.rm=TRUE)
wlth.ri.med <- median(wlth.ri.update, na.rm=TRUE)
wlth.ri.up <- quantile(wlth.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(hhmem.ri.lo,wedu.ri.lo,hhsex.ri.lo,dplace.ri.lo,dwimpr.ri.lo,sfimp.ri.lo,acchc.ri.lo,reg.ri.lo,wlth.ri.lo)
ri.med <- c(hhmem.ri.med,wedu.ri.med,hhsex.ri.med,dplace.ri.med,dwimpr.ri.med,sfimp.ri.med,acchc.ri.med,reg.ri.med,wlth.ri.med)
ri.up <- c(hhmem.ri.up,wedu.ri.up,hhsex.ri.up,dplace.ri.up,dwimpr.ri.up,sfimp.ri.up,acchc.ri.up,reg.ri.up,wlth.ri.up)
ri.out <- as.data.frame(cbind(ri.med,ri.up,ri.lo))
rownames(ri.out) <- c(out.colnames, "reg", "wlth")
PHASE1.ri.sort <- ri.out[order(ri.out[,1], decreasing=T), ]
PHASE1.ri.sort

PHASE1.ri.plt <- ggplot(PHASE1.ri.sort) +
  geom_bar(aes(x=reorder(row.names(PHASE1.ri.sort), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar( aes(x=row.names(PHASE1.ri.sort), ymin=ri.lo, ymax=ri.up),
                 width=0.4, colour="black", alpha=0.9, size=0.3)
PHASE1.ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")


## PHASE 2 - female traits
# woman bmi, woman age, delivery by c-section, postnatal checks
# wbmimn (wbmisd); wagemn (wagesd); csectp (csectv); pncp (pncv)

PHASE2vars <- u5dat %>%
  dplyr::select(cntry, clust, diarp, diarv, wbmimn, wbmisd, wagemn, wagesd, csectp, csectv, pncp, pncv)
head(PHASE2vars)
head(PHASE2vars[,c(5,7,9,11)])
dim(PHASE2vars)
PHASE2vars.naomit <- na.omit(PHASE2vars)
dim(PHASE2vars.naomit)

# variance inflation factor
usdm::vif(PHASE2vars.naomit[,c(5,7,9,11)])

# deterministic BRT with means only
PHASE2.brt.dtrm.mn <- gbm.step(PHASE2vars.naomit, gbm.x = attr(PHASE2vars.naomit, "names")[c(5,7,9,11)],
                               gbm.y = attr(PHASE2vars.naomit, "names")[3],
                               family="gaussian", max.trees=1000000, tolerance = 0.01,
                               learning.rate = 0.05, bag.fraction=0.75, tree.complexity = 2)
summary(PHASE2.brt.dtrm.mn)
gbm.plot(PHASE2.brt.dtrm.mn)
gbm.plot.fits(PHASE2.brt.dtrm.mn)

PHASE2.brt.dtrm.mn.CV.cor <- 100 * PHASE2.brt.dtrm.mn$cv.statistics$correlation.mean
PHASE2.brt.dtrm.mn.CV.cor.se <- 100 * PHASE2.brt.dtrm.mn$cv.statistics$correlation.se
print(c(PHASE2.brt.dtrm.mn.CV.cor, PHASE2.brt.dtrm.mn.CV.cor))

plot(PHASE2vars$wbmimn, PHASE2vars$diarp, pch=19, xlab="woman bmi", ylab="diarrhoea pr")
  
  
  ## iterate resampled boosted regression trees
  iter <- 100
  eq.sp.points <- 100
  
  # create storage arrays
  out.colnames <- c("wbmi", "wage", "csect", "pnc")
  val.arr <- pred.arr <- array(data = NA, dim=c(eq.sp.points, length=4, iter),
                               dimnames=list(paste("x",1:eq.sp.points,sep=""),
                                             out.colnames, paste("i",1:iter, sep="")))
  
  # create storage vectors
  CV.cor.vec <- CV.cor.se.vec <- wbmi.ri <- wage.ri <- csect.ri <- pnc.ri <- rep(NA, iter)
  
  for (i in 1:iter) {
    
    # country-resampling loop
    cntry.vec <- attr(table(PHASE2vars.naomit$cntry), "names")
    cntry.samp <- round(min(table(PHASE2vars.naomit$cntry)) / 2, 0)
    dat.cresamp <- PHASE2vars.naomit[1,]
    
    for (c in 1:length(cntry.vec)) {
      cntry.dat <- PHASE2vars.naomit[PHASE2vars.naomit$cntry == cntry.vec[c],]
      cntry.resamp <- cntry.dat[sort(sample(1:dim(cntry.dat)[1], cntry.samp, replace=F)), ]
      dat.cresamp <- rbind(dat.cresamp, cntry.resamp)
    } # end c
    dat.cresamp <- dat.cresamp[-1,]
    
    # resample
    # response
    diar.alpha <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$alpha
    diar.beta <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$beta
    diar.stch <- rbeta(length(diar.alpha), diar.alpha, diar.beta)
    diar.stoch <- ifelse(is.na(diar.stch==T), 0, diar.stch) # diarrhoea probability
    
    # continuous
    wbmi.stoch <- truncnorm::rtruncnorm(length(dat.cresamp$wbmimn), a=0, b=Inf, mean=dat.cresamp$wbmimn, sd=dat.cresamp$wbmisd) # woman bmi
    wage.stoch <- truncnorm::rtruncnorm(length(dat.cresamp$wagemn), a=0, b=Inf, mean=dat.cresamp$wagemn, sd=dat.cresamp$wagesd) # woman age
    
    # binomial probabilities
    csect.alpha <- estBetaParams(dat.cresamp$csectp, dat.cresamp$csectv^2)$alpha
    csect.beta <- estBetaParams(dat.cresamp$csectp, dat.cresamp$csectv^2)$beta
    csect.stch <- rbeta(length(csect.alpha), csect.alpha, csect.beta)
    csect.stoch <- ifelse(is.na(csect.stch==T), 0, csect.stch) # C-section probability
    
    pnc.alpha <- estBetaParams(dat.cresamp$pncp, dat.cresamp$pncv^2)$alpha
    pnc.beta <- estBetaParams(dat.cresamp$pncp, dat.cresamp$pncv^2)$beta
    pnc.stch <- rbeta(length(pnc.alpha), pnc.alpha, pnc.beta)
    pnc.stoch <- ifelse(is.na(pnc.stch==T), 0, pnc.stch) # postnatal checks probability
    
    ## collect resampled variables into single dataframe
    dat.resamp <- data.frame(diar.stoch,wbmi.stoch,wage.stoch,csect.stoch,pnc.stoch)
    colnames(dat.resamp) <- c("diar", "wbmi", "wage", "csect", "pnc")
    
    # scale
    dat.resamp.sc <- as.data.frame(scale(dat.resamp, center=T, scale=T))
    head(dat.resamp.sc)
    
    
    # boosted regression tree
    brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:5)],
                             gbm.y = attr(dat.resamp.sc, "names")[1],
                             family="gaussian", max.trees=1000000, tolerance = 0.00001,
                             learning.rate = 0.00025, bag.fraction=0.5, tree.complexity = 2,
                             silent=T, tolerance.method = "auto")
    
    if (i == 1 & is.null(brt.fit.rsmp)==F) {
      brt.fit.rsmp.old <- brt.fit.rsmp
    }
    
    
    if (is.null(brt.fit.rsmp) == T) {
      brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:5)],
                               gbm.y = attr(dat.resamp.sc, "names")[1],
                               family="gaussian", max.trees=1000000, tolerance = 0.000001,
                               learning.rate = 0.00001, bag.fraction=0.5, tree.complexity = 2,
                               silent=T, tolerance.method = "auto")
    }
    if (is.null(brt.fit.rsmp) == T) {
      brt.fit.rsmp <- brt.fit.rsmp.old
    }
    summ.fit <- summary(brt.fit.rsmp)
    
    if (is.null(brt.fit.rsmp) == F) {
      brt.fit.rsmp.old <- brt.fit.rsmp
    }
    
    # variable relative importance
    wbmi.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[1])]
    wage.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[2])]
    csect.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[3])]
    pnc.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[4])]
    
    # goodness of fit
    CV.cor.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.mean
    CV.cor.se.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.se
    
    # response curves
    RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=4)
    for (p in 1:4) {
      RESP.val[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,1]
      RESP.pred[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,2]
    } # end p
    RESP.val.dat <- as.data.frame(RESP.val)
    colnames(RESP.val.dat) <- brt.fit.rsmp$var.names
    RESP.pred.dat <- as.data.frame(RESP.pred)
    colnames(RESP.pred.dat) <- brt.fit.rsmp$var.names
    
    val.arr[, , i] <- as.matrix(RESP.val.dat)
    pred.arr[, , i] <- as.matrix(RESP.pred.dat)
    
    print(i)
    
  } # end i
  
  # kappa method to reduce effects of outliers on bootstrap estimates
  kappa <- 2
  kappa.n <- 4
  pred.update <- pred.arr[,,1:iter]
  
  for (k in 1:kappa.n) {
    boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
    boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
    
    for (z in 1:iter) {
      pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                    (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
    } # end z
    print(k)
  } # end k
  
  PHASE2.pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
  PHASE2.pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
  PHASE2.pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
  PHASE2.val.med <- apply(val.arr[,,1:iter], MARGIN=c(1,2), median, na.rm=T)
  
  # kappa method for output vectors
  PHASE2.CV.cor.update <- CV.cor.vec[1:iter]
  PHASE2.CV.cor.se.update <- CV.cor.se.vec[1:iter]
  
  wbmi.ri.update <- hhmem.ri[1:iter]
  wage.ri.update <- wage.ri[1:iter]
  csect.ri.update <- csect.ri[1:iter]
  pnc.ri.update <- pnc.ri[1:iter]
  
  for (k in 1:kappa.n) {
    PHASE2.CV.cor.mean <- mean(PHASE2.CV.cor.update, na.rm=T); PHASE2.CV.cor.sd <- sd(PHASE2.CV.cor.update, na.rm=T)
    PHASE2.CV.cor.se.mean <- mean(PHASE2.CV.cor.se.update, na.rm=T); PHASE2.CV.cor.se.sd <- sd(PHASE2.CV.cor.se.update, na.rm=T)
    
    wbmi.mean <- mean(wbmi.ri.update, na.rm=T); wbmi.sd <- sd(wbmi.ri.update, na.rm=T)
    wage.mean <- mean(wage.ri.update, na.rm=T); wage.sd <- sd(wage.ri.update, na.rm=T)
    csect.mean <- mean(csect.ri.update, na.rm=T); csect.sd <- sd(csect.ri.update, na.rm=T)
    pnc.mean <- mean(pnc.ri.update, na.rm=T); pnc.sd <- sd(pnc.ri.update, na.rm=T)
  
    for (u in 1:iter) {
      PHASE2.CV.cor.update[u] <- ifelse((PHASE2.CV.cor.update[u] < (PHASE2.CV.cor.mean-kappa*PHASE2.CV.cor.sd) | PHASE2.CV.cor.update[u] >
                                           (PHASE2.CV.cor.mean+kappa*PHASE2.CV.cor.sd)), NA, PHASE2.CV.cor.update[u])
      PHASE2.CV.cor.se.update[u] <- ifelse((PHASE2.CV.cor.se.update[u] < (PHASE2.CV.cor.se.mean-kappa*PHASE2.CV.cor.se.sd) | PHASE2.CV.cor.se.update[u]
                                            > (PHASE2.CV.cor.se.mean+kappa*PHASE2.CV.cor.se.sd)), NA, PHASE2.CV.cor.se.update[u])
      
      wbmi.ri.update[u] <- ifelse((wbmi.ri.update[u] < (wbmi.mean-kappa*wbmi.sd) | wbmi.ri.update[u] > (wbmi.mean+kappa*wbmi.sd)), NA, wbmi.ri.update[u])
      wage.ri.update[u] <- ifelse((wage.ri.update[u] < (wage.mean-kappa*wedu.sd) | wage.ri.update[u] > (wage.mean+kappa*wedu.sd)), NA, wage.ri.update[u])
      csect.ri.update[u] <- ifelse((csect.ri.update[u] < (csect.mean-kappa*csect.sd) | csect.ri.update[u] > (csect.mean+kappa*csect.sd)), NA, csect.ri.update[u])
      pnc.ri.update[u] <- ifelse((pnc.ri.update[u] < (pnc.mean-kappa*pnc.sd) | pnc.ri.update[u] > (pnc.mean+kappa*pnc.sd)), NA, pnc.ri.update[u])
      
    } # end u
    
    print(k)
    
  } # end k
  
  # summaries
  PHASE2.CV.cor.med <- median(PHASE2.CV.cor.update, na.rm=TRUE)
  PHASE2.CV.cor.lo <- quantile(PHASE2.CV.cor.update, probs=0.025, na.rm=TRUE)
  PHASE2.CV.cor.up <- quantile(PHASE2.CV.cor.update, probs=0.975, na.rm=TRUE)
  print(c(PHASE2.CV.cor.lo,PHASE2.CV.cor.med,PHASE2.CV.cor.up))
  
  wbmi.ri.lo <- quantile(wbmi.ri.update, probs=0.025, na.rm=TRUE)
  wbmi.ri.med <- median(wbmi.ri.update, na.rm=TRUE)
  wbmi.ri.up <- quantile(wbmi.ri.update, probs=0.975, na.rm=TRUE)
  
  wage.ri.lo <- quantile(wage.ri.update, probs=0.025, na.rm=TRUE)
  wage.ri.med <- median(wage.ri.update, na.rm=TRUE)
  wage.ri.up <- quantile(wage.ri.update, probs=0.975, na.rm=TRUE)
  
  csect.ri.lo <- quantile(csect.ri.update, probs=0.025, na.rm=TRUE)
  csect.ri.med <- median(csect.ri.update, na.rm=TRUE)
  csect.ri.up <- quantile(csect.ri.update, probs=0.975, na.rm=TRUE)
  
  pnc.ri.lo <- quantile(pnc.ri.update, probs=0.025, na.rm=TRUE)
  pnc.ri.med <- median(pnc.ri.update, na.rm=TRUE)
  pnc.ri.up <- quantile(pnc.ri.update, probs=0.975, na.rm=TRUE)
  
  ri.lo <- c(wbmi.ri.lo,wage.ri.lo,csect.ri.lo,pnc.ri.lo)
  ri.med <- c(wbmi.ri.med,wage.ri.med,csect.ri.med,pnc.ri.med)
  ri.up <- c(wbmi.ri.up,wage.ri.up,csect.ri.up,pnc.ri.up)
  ri.out <- as.data.frame(cbind(ri.med,ri.up,ri.lo))
  rownames(ri.out) <- out.colnames
  PHASE2.ri.sort <- ri.out[order(ri.out[,1], decreasing=T), ]
  PHASE2.ri.sort
  
  PHASE2.ri.plt <- ggplot(PHASE2.ri.sort) +
    geom_bar(aes(x=reorder(row.names(PHASE2.ri.sort), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
    geom_errorbar( aes(x=row.names(PHASE2.ri.sort), ymin=ri.lo, ymax=ri.up),
                   width=0.4, colour="black", alpha=0.9, size=0.3)
  PHASE2.ri.plt + coord_flip() +
    xlab("relative influence") + ylab("")
  


## PHASE 3 - child traits
# gender, birth size, stunting, underweight, wasting, drugs for intestinal parasites
# femp (femv); birthsord; stuntmn (stundsd); uweighmn (uweighsd); wastmn (wastsd); drugipp (drugipv)

PHASE3vars <- u5dat %>%
  dplyr::select(cntry, clust, diarp, diarv, femp, femv, birthsord, stuntmn, stuntsd, uweighmn, uweighsd,
                wastmn, wastsd, drugipp, drugipv)
head(PHASE3vars)
head(PHASE3vars[,c(5,7,8,10,12,14)])
dim(PHASE3vars)
PHASE3vars.naomit <- na.omit(PHASE3vars)
dim(PHASE3vars.naomit)

# variance inflation factor
usdm::vif(PHASE3vars.naomit[,c(5,8,10,12,14)])

# deterministic BRT with means only
PHASE3.brt.dtrm.mn <- gbm.step(PHASE3vars.naomit, gbm.x = attr(PHASE3vars.naomit, "names")[c(5,7,8,10,12,14)],
                               gbm.y = attr(PHASE3vars.naomit, "names")[3],
                               family="gaussian", max.trees=1000000, tolerance = 0.01,
                               learning.rate = 0.01, bag.fraction=0.75, tree.complexity = 2)
summary(PHASE3.brt.dtrm.mn)
gbm.plot(PHASE3.brt.dtrm.mn)
gbm.plot.fits(PHASE3.brt.dtrm.mn)

PHASE3.brt.dtrm.mn.CV.cor <- 100 * PHASE3.brt.dtrm.mn$cv.statistics$correlation.mean
PHASE3.brt.dtrm.mn.CV.cor.se <- 100 * PHASE3.brt.dtrm.mn$cv.statistics$correlation.se
print(c(PHASE3.brt.dtrm.mn.CV.cor, PHASE3.brt.dtrm.mn.CV.cor))

plot(PHASE3vars$uweighmn, PHASE3vars$diarp, pch=19, xlab="underweight pr", ylab="diarrhoea pr")
plot(PHASE3vars$stuntmn, PHASE3vars$diarp, pch=19, xlab="stunting pr", ylab="diarrhoea pr")
plot(PHASE3vars$wastmn, PHASE3vars$diarp, pch=19, xlab="wasting pr", ylab="diarrhoea pr")


## iterate resampled boosted regression trees
iter <- 100
eq.sp.points <- 100

# create storage arrays
out.colnames <- c("fem", "stunt", "uweigh", "wast", "drugip")
val.arr <- pred.arr <- array(data = NA, dim=c(eq.sp.points, length=5, iter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                                           out.colnames, paste("i",1:iter, sep="")))
brthsz.mat <- matrix(data=NA, nrow=iter, ncol=3)

# create storage vectors
CV.cor.vec <- CV.cor.se.vec <- fem.ri <- brthsz.ri <- stunt.ri <- uweigh.ri <- wast.ri <- drugip.ri <- rep(NA, iter)

for (i in 1:iter) {
  
  # country-resampling loop
  cntry.vec <- attr(table(PHASE3vars.naomit$cntry), "names")
  cntry.samp <- round(min(table(PHASE3vars.naomit$cntry)) / 2, 0)
  dat.cresamp <- PHASE3vars.naomit[1,]
  
  for (c in 1:length(cntry.vec)) {
    cntry.dat <- PHASE3vars.naomit[PHASE3vars.naomit$cntry == cntry.vec[c],]
    cntry.resamp <- cntry.dat[sort(sample(1:dim(cntry.dat)[1], cntry.samp, replace=T)), ]
    dat.cresamp <- rbind(dat.cresamp, cntry.resamp)
  } # end c
  dat.cresamp <- dat.cresamp[-1,]
  
  # resample
  # response
  diar.alpha <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$alpha
  diar.beta <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$beta
  diar.stch <- rbeta(length(diar.alpha), diar.alpha, diar.beta)
  diar.stoch <- ifelse(is.na(diar.stch==T), 0, diar.stch) # diarrhoea probability
  
  # continuous
  stunt.stoch <- rnorm(length(dat.cresamp$stuntmn), mean=dat.cresamp$stuntmn, sd=dat.cresamp$stuntsd)
  uweigh.stoch <- rnorm(length(dat.cresamp$uweighmn), mean=dat.cresamp$uweighmn, sd=dat.cresamp$uweighsd)
  wast.stoch <- rnorm(length(dat.cresamp$wastmn), mean=dat.cresamp$wastmn, sd=dat.cresamp$wastsd)
  
  # binomial probabilities
  fem.alpha <- estBetaParams(dat.cresamp$femp, dat.cresamp$femv^2)$alpha
  fem.beta <- estBetaParams(dat.cresamp$femp, dat.cresamp$femv^2)$beta
  fem.stch <- rbeta(length(fem.alpha), fem.alpha, fem.beta)
  fem.stoch <- ifelse(is.na(fem.stch==T), 0, fem.stch) # probability of child gender being female
  
  drugip.alpha <- estBetaParams(dat.cresamp$drugipp, dat.cresamp$drugipv^2)$alpha
  drugip.beta <- estBetaParams(dat.cresamp$drugipp, dat.cresamp$drugipv^2)$beta
  drugip.stch <- rbeta(length(drugip.alpha), drugip.alpha, drugip.beta)
  drugip.stoch <- ifelse(is.na(drugip.stch==T), 0, drugip.stch) # drugs for intestinal parasites
  
  # ordinal factors
  brthsz.stch <- spatstat.random::rpoistrunc(length(dat.cresamp$birthsord), lambda=as.integer(dat.cresamp$birthsord), minimum=1, method="harding")
  brthsz.stch1 <- ifelse(brthsz.stch > 3, 3, brthsz.stch)
  brthsz.stoch <- factor(brthsz.stch1, ordered=T)

  ## collect resampled variables into single dataframe
  dat.resamp <- data.frame(diar.stoch,stunt.stoch,uweigh.stoch,wast.stoch,fem.stoch,drugip.stoch,brthsz.stoch)
  colnames(dat.resamp) <- c("diar", "stunt", "uweigh", "wast", "fem", "drugip", "brthsz")
  
  # scale
  dat.resamp.cont.sc <- as.data.frame(scale(dat.resamp[,1:6], center=T, scale=T))
  dat.resamp.sc <- data.frame(dat.resamp.cont.sc, "brthsz"=dat.resamp$brthsz)
  #head(dat.resamp.sc)
  #dim(dat.resamp.sc)
  #str(dat.resamp.sc)
  
  # boosted regression tree
  brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:7)],
                           gbm.y = attr(dat.resamp.sc, "names")[1],
                           family="gaussian", max.trees=1000000, tolerance = 0.0001,
                           learning.rate = 0.000001, bag.fraction=0.7, tree.complexity = 2,
                           silent=T, tolerance.method = "auto", step.size=40)
  
  if (i == 1 & is.null(brt.fit.rsmp)==F) {
    brt.fit.rsmp.old <- brt.fit.rsmp
  }
  
  if (is.null(brt.fit.rsmp) == T) {
    brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:7)],
                             gbm.y = attr(dat.resamp.sc, "names")[1],
                             family="gaussian", max.trees=1000000, tolerance = 0.0001,
                             learning.rate = 0.0000001, bag.fraction=0.7, tree.complexity = 2,
                             silent=T, tolerance.method = "auto", step.size=10)
  }
  if (is.null(brt.fit.rsmp) == T) {
    brt.fit.rsmp <- brt.fit.rsmp.old
  }
  summ.fit <- summary(brt.fit.rsmp)
  
  if (is.null(brt.fit.rsmp) == F) {
    brt.fit.rsmp.old <- brt.fit.rsmp
  }
  
  # variable relative importance
  stunt.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[2])]
  uweigh.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[3])]
  wast.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[4])]
  fem.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[1])]
  drugip.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[5])]
  brthsz.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == "brthsz")]
  
  # goodness of fit
  CV.cor.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.mean
  CV.cor.se.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.se
  
  # response curves (continuous & probabilities)
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=5)
  for (p in 1:5) {
    RESP.val[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit.rsmp$var.names[1:5]
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit.rsmp$var.names[1:5]
  #RESP.pred.dat

  val.arr[, , i] <- as.matrix(RESP.val.dat)
  pred.arr[, , i] <- as.matrix(RESP.pred.dat)
  
  # categorical predictions
  RESP.brthsz <- plot.gbm(brt.fit.rsmp, i.var=6, continuous.resolution = eq.sp.points, return.grid=T)
  colnames(RESP.brthsz)[2] <- "pred"

  brthsz.mat[i,] <- RESP.brthsz[,2]
  
  print(i)
  
} # end i

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:iter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:iter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

boot.brthsz.mean <- apply(brthsz.mat, MARGIN=1, mean, na.rm=T)
boot.brthsz.sd <- apply(brthsz.mat, MARGIN=1, sd, na.rm=T)
PHASE3.brthsz.pred.med <- apply(brthsz.mat, MARGIN=1, median, na.rm=T)
PHASE3.brthsz.pred.lo <- apply(brthsz.mat, MARGIN=1, quantile, probs=0.025, na.rm=T)
PHASE3.brthsz.pred.up <- apply(brthsz.mat, MARGIN=1, quantile, probs=0.975, na.rm=T)

PHASE3.pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
PHASE3.pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
PHASE3.pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
PHASE3.val.med <- apply(val.arr[,,1:iter], MARGIN=c(1,2), median, na.rm=T)


# kappa method for output vectors
PHASE3.CV.cor.update <- CV.cor.vec[1:iter]
PHASE3.CV.cor.se.update <- CV.cor.se.vec[1:iter]

stunt.ri.update <- stunt.ri[1:iter]
uweigh.ri.update <- uweigh.ri[1:iter]
wast.ri.update <- wast.ri[1:iter]
fem.ri.update <- fem.ri[1:iter]
drugip.ri.update <- drugip.ri[1:iter]
brthsz.ri.update <- brthsz.ri[1:iter]

for (k in 1:kappa.n) {
  PHASE3.CV.cor.mean <- mean(PHASE3.CV.cor.update, na.rm=T); PHASE3.CV.cor.sd <- sd(PHASE3.CV.cor.update, na.rm=T)
  PHASE3.CV.cor.se.mean <- mean(PHASE3.CV.cor.se.update, na.rm=T); PHASE3.CV.cor.se.sd <- sd(PHASE3.CV.cor.se.update, na.rm=T)
  
  stunt.mean <- mean(stunt.ri.update, na.rm=T); stunt.sd <- sd(stunt.ri.update, na.rm=T)
  uweigh.mean <- mean(uweigh.ri.update, na.rm=T); uweigh.sd <- sd(uweigh.ri.update, na.rm=T)
  wast.mean <- mean(wast.ri.update, na.rm=T); wast.sd <- sd(wast.ri.update, na.rm=T)
  fem.mean <- mean(fem.ri.update, na.rm=T); fem.sd <- sd(fem.ri.update, na.rm=T)
  drugip.mean <- mean(drugip.ri.update, na.rm=T); drugip.sd <- sd(drugip.ri.update, na.rm=T)
  brthsz.mean <- mean(brthsz.ri.update, na.rm=T); brthsz.sd <- sd(brthsz.ri.update, na.rm=T)
  
  for (u in 1:iter) {
    PHASE3.CV.cor.update[u] <- ifelse((PHASE3.CV.cor.update[u] < (PHASE3.CV.cor.mean-kappa*PHASE3.CV.cor.sd) | PHASE3.CV.cor.update[u] >
                                         (PHASE3.CV.cor.mean+kappa*PHASE3.CV.cor.sd)), NA, PHASE3.CV.cor.update[u])
    PHASE3.CV.cor.se.update[u] <- ifelse((PHASE3.CV.cor.se.update[u] < (PHASE3.CV.cor.se.mean-kappa*PHASE3.CV.cor.se.sd) | PHASE3.CV.cor.se.update[u]
                                          > (PHASE3.CV.cor.se.mean+kappa*PHASE3.CV.cor.se.sd)), NA, PHASE3.CV.cor.se.update[u])
    
    stunt.ri.update[u] <- ifelse((stunt.ri.update[u] < (stunt.mean-kappa*stunt.sd) | stunt.ri.update[u] > (stunt.mean+kappa*stunt.sd)), NA, stunt.ri.update[u])
    uweigh.ri.update[u] <- ifelse((uweigh.ri.update[u] < (uweigh.mean-kappa*uweigh.sd) | uweigh.ri.update[u] > (uweigh.mean+kappa*uweigh.sd)), NA, uweigh.ri.update[u])
    wast.ri.update[u] <- ifelse((wast.ri.update[u] < (wast.mean-kappa*wast.sd) | wast.ri.update[u] > (wast.mean+kappa*wast.sd)), NA, wast.ri.update[u])
    fem.ri.update[u] <- ifelse((fem.ri.update[u] < (fem.mean-kappa*fem.sd) | fem.ri.update[u] > (fem.mean+kappa*fem.sd)), NA, fem.ri.update[u])
    drugip.ri.update[u] <- ifelse((drugip.ri.update[u] < (drugip.mean-kappa*drugip.sd) | drugip.ri.update[u] > (drugip.mean+kappa*drugip.sd)), NA, drugip.ri.update[u])
    brthsz.ri.update[u] <- ifelse((brthsz.ri.update[u] < (brthsz.mean-kappa*brthsz.sd) | brthsz.ri.update[u] > (brthsz.mean+kappa*brthsz.sd)), NA, brthsz.ri.update[u])
    
  } # end u
  
  print(k)
  
} # end k

# summaries
PHASE3.CV.cor.med <- median(PHASE3.CV.cor.update, na.rm=TRUE)
PHASE3.CV.cor.lo <- quantile(PHASE3.CV.cor.update, probs=0.025, na.rm=TRUE)
PHASE3.CV.cor.up <- quantile(PHASE3.CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(PHASE3.CV.cor.lo,PHASE3.CV.cor.med,PHASE3.CV.cor.up))

stunt.ri.lo <- quantile(stunt.ri.update, probs=0.025, na.rm=TRUE)
stunt.ri.med <- median(stunt.ri.update, na.rm=TRUE)
stunt.ri.up <- quantile(stunt.ri.update, probs=0.975, na.rm=TRUE)

uweigh.ri.lo <- quantile(uweigh.ri.update, probs=0.025, na.rm=TRUE)
uweigh.ri.med <- median(uweigh.ri.update, na.rm=TRUE)
uweigh.ri.up <- quantile(uweigh.ri.update, probs=0.975, na.rm=TRUE)

wast.ri.lo <- quantile(wast.ri.update, probs=0.025, na.rm=TRUE)
wast.ri.med <- median(wast.ri.update, na.rm=TRUE)
wast.ri.up <- quantile(wast.ri.update, probs=0.975, na.rm=TRUE)

fem.ri.lo <- quantile(fem.ri.update, probs=0.025, na.rm=TRUE)
fem.ri.med <- median(fem.ri.update, na.rm=TRUE)
fem.ri.up <- quantile(fem.ri.update, probs=0.975, na.rm=TRUE)

drugip.ri.lo <- quantile(drugip.ri.update, probs=0.025, na.rm=TRUE)
drugip.ri.med <- median(drugip.ri.update, na.rm=TRUE)
drugip.ri.up <- quantile(drugip.ri.update, probs=0.975, na.rm=TRUE)

brthsz.ri.lo <- quantile(brthsz.ri.update, probs=0.025, na.rm=TRUE)
brthsz.ri.med <- median(brthsz.ri.update, na.rm=TRUE)
brthsz.ri.up <- quantile(brthsz.ri.update, probs=0.975, na.rm=TRUE)


ri.lo <- c(stunt.ri.lo,uweigh.ri.lo,wast.ri.lo,fem.ri.lo,drugip.ri.lo,brthsz.ri.lo)
ri.med <- c(stunt.ri.med,uweigh.ri.med,wast.ri.med,fem.ri.med,drugip.ri.med,brthsz.ri.med)
ri.up <- c(stunt.ri.up,uweigh.ri.up,wast.ri.up,fem.ri.up,drugip.ri.up,brthsz.ri.up)
ri.out <- as.data.frame(cbind(ri.med,ri.up,ri.lo))
rownames(ri.out) <- c(out.colnames,"brthsz")
PHASE3.ri.sort <- ri.out[order(ri.out[,1], decreasing=T), ]
PHASE3.ri.sort

PHASE3.ri.plt <- ggplot(PHASE3.ri.sort) +
  geom_bar(aes(x=reorder(row.names(PHASE3.ri.sort), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar( aes(x=row.names(PHASE3.ri.sort), ymin=ri.lo, ymax=ri.up),
                 width=0.4, colour="black", alpha=0.9, size=0.3)
PHASE3.ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")


## PHASE 4 - climate variables
# mean annual temperature, temperature annual range, annual precipitation, precipitation seasonality,
# precipitation driest quarter
# MATmn (MATsd); TARmn (TARsd); PRCPmn (PRCPsd); PRCPSNSmn (PRCPSNSsd); PRCPDQmn (PRCPDQsd)

PHASE4vars <- u5dat %>%
  dplyr::select(cntry, clust, diarp, diarv, MATmn, MATsd, TARmn, TARsd, PRCPmn, PRCPsd,
                PRCPSNSmn, PRCPSNSsd, PRCPDQmn, PRCPDQsd)
head(PHASE4vars)
head(PHASE4vars[,c(5,7,9,11,13)])
dim(PHASE4vars)
PHASE4vars.naomit <- na.omit(PHASE4vars)
dim(PHASE4vars.naomit)

# variance inflation factor
usdm::vif(PHASE4vars.naomit[,c(5,7,9,11,13)])

# deterministic BRT with means only
PHASE4.brt.dtrm.mn <- gbm.step(PHASE4vars.naomit, gbm.x = attr(PHASE4vars.naomit, "names")[c(5,7,9,11,13)],
                               gbm.y = attr(PHASE4vars.naomit, "names")[3],
                               family="gaussian", max.trees=1000000, tolerance = 0.01,
                               learning.rate = 0.01, bag.fraction=0.75, tree.complexity = 2)
summary(PHASE4.brt.dtrm.mn)
gbm.plot(PHASE4.brt.dtrm.mn)
gbm.plot.fits(PHASE4.brt.dtrm.mn)

PHASE4.brt.dtrm.mn.CV.cor <- 100 * PHASE4.brt.dtrm.mn$cv.statistics$correlation.mean
PHASE4.brt.dtrm.mn.CV.cor.se <- 100 * PHASE4.brt.dtrm.mn$cv.statistics$correlation.se
print(c(PHASE4.brt.dtrm.mn.CV.cor, PHASE4.brt.dtrm.mn.CV.cor))

plot(PHASE4vars$TARmn, PHASE4vars$diarp, pch=19, xlab="temperature annual range", ylab="diarrhoea pr")
plot(PHASE4vars$PRCPmn, PHASE4vars$diarp, pch=19, xlab="annual precipitation", ylab="diarrhoea pr")
plot(PHASE4vars$MATmn, PHASE4vars$diarp, pch=19, xlab="mean annual temperature", ylab="diarrhoea pr")


## iterate resampled boosted regression trees
iter <- 100
eq.sp.points <- 100

# create storage arrays
out.colnames <- c("MAT", "TAR", "PRCP", "PRCPSNS", "PRCPDQ")
val.arr <- pred.arr <- array(data = NA, dim=c(eq.sp.points, length=5, iter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                                           out.colnames, paste("i",1:iter, sep="")))

# create storage vectors
CV.cor.vec <- CV.cor.se.vec <- MAT.ri <- TAR.ri <- PRCP.ri <- PRCPSNS.ri <- PRCPDQ.ri <- rep(NA, iter)

for (i in 1:iter) {
  
  # country-resampling loop
  cntry.vec <- attr(table(PHASE4vars.naomit$cntry), "names")
  cntry.samp <- round(min(table(PHASE4vars.naomit$cntry)) / 2, 0)
  dat.cresamp <- PHASE4vars.naomit[1,]
  
  for (c in 1:length(cntry.vec)) {
    cntry.dat <- PHASE4vars.naomit[PHASE4vars.naomit$cntry == cntry.vec[c],]
    cntry.resamp <- cntry.dat[sort(sample(1:dim(cntry.dat)[1], cntry.samp, replace=F)), ]
    dat.cresamp <- rbind(dat.cresamp, cntry.resamp)
  } # end c
  dat.cresamp <- dat.cresamp[-1,]
  
  # resample
  # response
  diar.alpha <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$alpha
  diar.beta <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$beta
  diar.stch <- rbeta(length(diar.alpha), diar.alpha, diar.beta)
  diar.stoch <- ifelse(is.na(diar.stch==T), 0, diar.stch) # diarrhoea probability
  
  # continuous
  MAT.stoch <- rnorm(length(dat.cresamp$MATmn), mean=dat.cresamp$MATmn, sd=dat.cresamp$MATsd) # mean annual temperature
  TAR.stoch <- rnorm(length(dat.cresamp$TARmn), mean=dat.cresamp$TARmn, sd=dat.cresamp$TARsd) # temperature annual range
  PRCP.stoch <- rnorm(length(dat.cresamp$PRCPmn), mean=dat.cresamp$PRCPmn, sd=dat.cresamp$PRCPsd) # annual precipitation
  PRCPSNS.stoch <- rnorm(length(dat.cresamp$PRCPSNSmn), mean=dat.cresamp$PRCPSNSmn, sd=dat.cresamp$PRCPSNSsd) # precipitation seasonality
  PRCPDQ.stoch <- rnorm(length(dat.cresamp$PRCPDQmn), mean=dat.cresamp$PRCPDQmn, sd=dat.cresamp$PRCPDQsd) # precipitation driest quarter
  
  ## collect resampled variables into single dataframe
  dat.resamp <- data.frame(diar.stoch,MAT.stoch,TAR.stoch,PRCP.stoch,PRCPSNS.stoch,PRCPDQ.stoch)
  colnames(dat.resamp) <- c("diar", "MAT", "TAR", "PRCP", "PRCPSNS", "PRCPDQ")
  
  # scale
  dat.resamp.sc <- as.data.frame(scale(dat.resamp, center=T, scale=T))
  head(dat.resamp.sc)
  
  # boosted regression tree
  brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:6)],
                           gbm.y = attr(dat.resamp.sc, "names")[1],
                           family="gaussian", max.trees=1000000, tolerance = 0.00001,
                           learning.rate = 0.001, bag.fraction=0.5, tree.complexity = 2,
                           silent=T, tolerance.method = "auto")
  
  if (i == 1 & is.null(brt.fit.rsmp)==F) {
    brt.fit.rsmp.old <- brt.fit.rsmp
  }
  
  
  if (is.null(brt.fit.rsmp) == T) {
    brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:6)],
                             gbm.y = attr(dat.resamp.sc, "names")[1],
                             family="gaussian", max.trees=1000000, tolerance = 0.000001,
                             learning.rate = 0.0001, bag.fraction=0.5, tree.complexity = 2,
                             silent=T, tolerance.method = "auto", step.size = 20)
  }
  if (is.null(brt.fit.rsmp) == T) {
    brt.fit.rsmp <- brt.fit.rsmp.old
  }
  summ.fit <- summary(brt.fit.rsmp)
  
  if (is.null(brt.fit.rsmp) == F) {
    brt.fit.rsmp.old <- brt.fit.rsmp
  }
  
  # variable relative importance
  MAT.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[1])]
  TAR.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[2])]
  PRCP.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[3])]
  PRCPSNS.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[4])]
  PRCPDQ.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[5])]
  
  # goodness of fit
  CV.cor.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.mean
  CV.cor.se.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=5)
  for (p in 1:5) {
    RESP.val[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit.rsmp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit.rsmp$var.names
  
  val.arr[, , i] <- as.matrix(RESP.val.dat)
  pred.arr[, , i] <- as.matrix(RESP.pred.dat)
  
  print(i)
  
} # end i

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:iter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:iter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

PHASE4.pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
PHASE4.pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
PHASE4.pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
PHASE4.val.med <- apply(val.arr[,,1:iter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
PHASE4.CV.cor.update <- CV.cor.vec[1:iter]
PHASE4.CV.cor.se.update <- CV.cor.se.vec[1:iter]

MAT.ri.update <- MAT.ri[1:iter]
TAR.ri.update <- TAR.ri[1:iter]
PRCP.ri.update <- PRCP.ri[1:iter]
PRCPSNS.ri.update <- PRCPSNS.ri[1:iter]
PRCPDQ.ri.update <- PRCPDQ.ri[1:iter]

for (k in 1:kappa.n) {
  PHASE4.CV.cor.mean <- mean(PHASE4.CV.cor.update, na.rm=T); PHASE4.CV.cor.sd <- sd(PHASE4.CV.cor.update, na.rm=T)
  PHASE4.CV.cor.se.mean <- mean(PHASE4.CV.cor.se.update, na.rm=T); PHASE4.CV.cor.se.sd <- sd(PHASE4.CV.cor.se.update, na.rm=T)
  
  MAT.mean <- mean(MAT.ri.update, na.rm=T); MAT.sd <- sd(MAT.ri.update, na.rm=T)
  TAR.mean <- mean(wedu.ri.update, na.rm=T); TAR.sd <- sd(TAR.ri.update, na.rm=T)
  PRCP.mean <- mean(PRCP.ri.update, na.rm=T); PRCP.sd <- sd(PRCP.ri.update, na.rm=T)
  PRCPSNS.mean <- mean(PRCPSNS.ri.update, na.rm=T); PRCPSNS.sd <- sd(PRCPSNS.ri.update, na.rm=T)
  PRCPDQ.mean <- mean(PRCPDQ.ri.update, na.rm=T); PRCPDQ.sd <- sd(PRCPDQ.ri.update, na.rm=T)
  
  for (u in 1:iter) {
    PHASE4.CV.cor.update[u] <- ifelse((PHASE4.CV.cor.update[u] < (PHASE4.CV.cor.mean-kappa*PHASE4.CV.cor.sd) | PHASE4.CV.cor.update[u] >
                                         (PHASE4.CV.cor.mean+kappa*PHASE4.CV.cor.sd)), NA, PHASE4.CV.cor.update[u])
    PHASE4.CV.cor.se.update[u] <- ifelse((PHASE4.CV.cor.se.update[u] < (PHASE4.CV.cor.se.mean-kappa*PHASE4.CV.cor.se.sd) | PHASE4.CV.cor.se.update[u]
                                          > (PHASE4.CV.cor.se.mean+kappa*PHASE4.CV.cor.se.sd)), NA, PHASE4.CV.cor.se.update[u])
    
    MAT.ri.update[u] <- ifelse((MAT.ri.update[u] < (MAT.mean-kappa*MAT.sd) | MAT.ri.update[u] > (MAT.mean+kappa*MAT.sd)), NA, MAT.ri.update[u])
    TAR.ri.update[u] <- ifelse((TAR.ri.update[u] < (TAR.mean-kappa*TAR.sd) | TAR.ri.update[u] > (TAR.mean+kappa*TAR.sd)), NA, TAR.ri.update[u])
    PRCP.ri.update[u] <- ifelse((PRCP.ri.update[u] < (PRCP.mean-kappa*PRCP.sd) | PRCP.ri.update[u] > (PRCP.mean+kappa*PRCP.sd)), NA, PRCP.ri.update[u])
    PRCPSNS.ri.update[u] <- ifelse((PRCPSNS.ri.update[u] < (PRCPSNS.mean-kappa*PRCPSNS.sd) | PRCPSNS.ri.update[u] > (PRCPSNS.mean+kappa*PRCPSNS.sd)), NA, PRCPSNS.ri.update[u])
    PRCPDQ.ri.update[u] <- ifelse((PRCPDQ.ri.update[u] < (PRCPDQ.mean-kappa*PRCPDQ.sd) | PRCPDQ.ri.update[u] > (PRCPDQ.mean+kappa*PRCPDQ.sd)), NA, PRCPDQ.ri.update[u])
    
  } # end u
  
  print(k)
  
} # end k

# summaries
PHASE4.CV.cor.med <- median(PHASE4.CV.cor.update, na.rm=TRUE)
PHASE4.CV.cor.lo <- quantile(PHASE4.CV.cor.update, probs=0.025, na.rm=TRUE)
PHASE4.CV.cor.up <- quantile(PHASE4.CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(PHASE4.CV.cor.lo,PHASE4.CV.cor.med,PHASE4.CV.cor.up))

MAT.ri.lo <- quantile(MAT.ri.update, probs=0.025, na.rm=TRUE)
MAT.ri.med <- median(MAT.ri.update, na.rm=TRUE)
MAT.ri.up <- quantile(MAT.ri.update, probs=0.975, na.rm=TRUE)

TAR.ri.lo <- quantile(TAR.ri.update, probs=0.025, na.rm=TRUE)
TAR.ri.med <- median(TAR.ri.update, na.rm=TRUE)
TAR.ri.up <- quantile(TAR.ri.update, probs=0.975, na.rm=TRUE)

PRCP.ri.lo <- quantile(PRCP.ri.update, probs=0.025, na.rm=TRUE)
PRCP.ri.med <- median(PRCP.ri.update, na.rm=TRUE)
PRCP.ri.up <- quantile(PRCP.ri.update, probs=0.975, na.rm=TRUE)

PRCPSNS.ri.lo <- quantile(PRCPSNS.ri.update, probs=0.025, na.rm=TRUE)
PRCPSNS.ri.med <- median(PRCPSNS.ri.update, na.rm=TRUE)
PRCPSNS.ri.up <- quantile(PRCPSNS.ri.update, probs=0.975, na.rm=TRUE)

PRCPDQ.ri.lo <- quantile(PRCPDQ.ri.update, probs=0.025, na.rm=TRUE)
PRCPDQ.ri.med <- median(PRCPDQ.ri.update, na.rm=TRUE)
PRCPDQ.ri.up <- quantile(PRCPDQ.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(MAT.ri.lo,TAR.ri.lo,PRCP.ri.lo,PRCPSNS.ri.lo,PRCPDQ.ri.lo)
ri.med <- c(MAT.ri.med,TAR.ri.med,PRCP.ri.med,PRCPSNS.ri.med,PRCPDQ.ri.med)
ri.up <- c(MAT.ri.up,TAR.ri.up,PRCP.ri.up,PRCPSNS.ri.up,PRCPDQ.ri.up)
ri.out <- as.data.frame(cbind(ri.med,ri.up,ri.lo))
rownames(ri.out) <- out.colnames
PHASE4.ri.sort <- ri.out[order(ri.out[,1], decreasing=T), ]
PHASE4.ri.sort

PHASE4.ri.plt <- ggplot(PHASE4.ri.sort) +
  geom_bar(aes(x=reorder(row.names(PHASE4.ri.sort), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar( aes(x=row.names(PHASE4.ri.sort), ymin=ri.lo, ymax=ri.up),
                 width=0.4, colour="black", alpha=0.9, size=0.3)
PHASE4.ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")


#################################
## compare RI plots among phases
#################################
P1 <- PHASE1.ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")
P2 <- PHASE2.ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")
P3 <- PHASE3.ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")
P4 <- PHASE4.ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")

ggarrange(P1, P2, P3, P4, 
          labels = c("SOCIO-ECONOMIC", "FEMALE TRAITS", "CHILD TRAITS", "CLIMATE"),
          label.x=0.4, label.y=0.2, hjust=-0.7, vjust=3,
          ncol = 2, nrow = 2)

# summary CVs by phase
print(c(PHASE1.CV.cor.lo,PHASE1.CV.cor.med,PHASE1.CV.cor.up)) # socio-economic
print(c(PHASE2.CV.cor.lo,PHASE2.CV.cor.med,PHASE2.CV.cor.up)) # mother traits
print(c(PHASE3.CV.cor.lo,PHASE3.CV.cor.med,PHASE3.CV.cor.up)) # child traits
print(c(PHASE4.CV.cor.lo,PHASE4.CV.cor.med,PHASE4.CV.cor.up)) # climate


## FINAL COMBINED PHASE5
## variables to keep: wedu, hhmem, pnc, fem, drugip, TAR, PRCP
PHASE5vars <- u5dat %>%
  dplyr::select(cntry, clust, diarp, diarv, wedumn, wedusd, hhmemmn, hhmemsd, pncp, pncv,
                femp, femv, drugipp, drugipv, TARmn, TARsd, PRCPmn, PRCPsd)
head(PHASE5vars)
head(PHASE5vars[,c(5,7,9,11,13,15,17)])
dim(PHASE5vars)
PHASE5vars.naomit <- na.omit(PHASE5vars)
dim(PHASE5vars.naomit)

# variance inflation factor
usdm::vif(PHASE5vars.naomit[,c(5,7,9,11,13,15,17)])

# deterministic BRT with means only
PHASE5.brt.dtrm.mn <- gbm.step(PHASE5vars.naomit, gbm.x = attr(PHASE5vars.naomit, "names")[c(5,7,9,11,13,15,17)],
                               gbm.y = attr(PHASE5vars.naomit, "names")[3],
                               family="gaussian", max.trees=1000000, tolerance = 0.01,
                               learning.rate = 0.001, bag.fraction=0.75, tree.complexity = 2, step.size=30)
summary(PHASE5.brt.dtrm.mn)
gbm.plot(PHASE5.brt.dtrm.mn)
gbm.plot.fits(PHASE5.brt.dtrm.mn)

PHASE5.brt.dtrm.mn.CV.cor <- 100 * PHASE5.brt.dtrm.mn$cv.statistics$correlation.mean
PHASE5.brt.dtrm.mn.CV.cor.se <- 100 * PHASE5.brt.dtrm.mn$cv.statistics$correlation.se
print(c(PHASE5.brt.dtrm.mn.CV.cor, PHASE5.brt.dtrm.mn.CV.cor))

plot(PHASE5vars$wedumn, PHASE5vars$diarp, pch=19, xlab="woman's education", ylab="diarrhoea pr")
plot(PHASE5vars$TARmn, PHASE5vars$diarp, pch=19, xlab="temperature annual range", ylab="diarrhoea pr")
plot(PHASE5vars$PRCPmn, PHASE5vars$diarp, pch=19, xlab="annual precipitation", ylab="diarrhoea pr")


## iterate resampled boosted regression trees
iter <- 1000
eq.sp.points <- 100

# create storage arrays
out.colnames <- c("wedu", "hhmem", "pnc", "fem", "drugip", "TAR", "PRCP")
val.arr <- pred.arr <- array(data = NA, dim=c(eq.sp.points, length=7, iter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                                           out.colnames, paste("i",1:iter, sep="")))

# create storage vectors
CV.cor.vec <- CV.cor.se.vec <- wedu.ri <- hhmem.ri <- pnc.ri <- fem.ri <- drugip.ri <-
    TAR.ri <- PRCP.ri <- rep(NA, iter)

for (i in 1:iter) {
  
  # country-resampling loop
  cntry.vec <- attr(table(PHASE5vars.naomit$cntry), "names")
  cntry.samp <- round(min(table(PHASE5vars.naomit$cntry)) / 2, 0)
  dat.cresamp <- PHASE5vars.naomit[1,]
  
  for (c in 1:length(cntry.vec)) {
    cntry.dat <- PHASE5vars.naomit[PHASE5vars.naomit$cntry == cntry.vec[c],]
    cntry.resamp <- cntry.dat[sort(sample(1:dim(cntry.dat)[1], cntry.samp, replace=F)), ]
    dat.cresamp <- rbind(dat.cresamp, cntry.resamp)
  } # end c
  dat.cresamp <- dat.cresamp[-1,]
  
  # resample
  # response
  diar.alpha <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$alpha
  diar.beta <- estBetaParams(dat.cresamp$diarp, dat.cresamp$diarv^2)$beta
  diar.stch <- rbeta(length(diar.alpha), diar.alpha, diar.beta)
  diar.stoch <- ifelse(is.na(diar.stch==T), 0, diar.stch) # diarrhoea probability
  
  # continuous
  wedu.stoch <- truncnorm::rtruncnorm(length(dat.cresamp$wedumn), a=0, b=Inf, mean=dat.cresamp$wedumn, sd=dat.cresamp$wedusd) # woman's education
  hhmem.stoch <- truncnorm::rtruncnorm(length(dat.cresamp$hhmemmn), a=0, b=Inf, mean=dat.cresamp$hhmemmn, sd=dat.cresamp$hhmemsd) # household members
  TAR.stoch <- rnorm(length(dat.cresamp$TARmn), mean=dat.cresamp$TARmn, sd=dat.cresamp$TARsd) # temperature annual range
  PRCP.stoch <- rnorm(length(dat.cresamp$PRCPmn), mean=dat.cresamp$PRCPmn, sd=dat.cresamp$PRCPsd) # annual precipitation
  
  # binomial probabilities
  fem.alpha <- estBetaParams(dat.cresamp$femp, dat.cresamp$femv^2)$alpha
  fem.beta <- estBetaParams(dat.cresamp$femp, dat.cresamp$femv^2)$beta
  fem.stch <- rbeta(length(fem.alpha), fem.alpha, fem.beta)
  fem.stoch <- ifelse(is.na(fem.stch==T), 0, fem.stch) # probability of child gender being female
  
  drugip.alpha <- estBetaParams(dat.cresamp$drugipp, dat.cresamp$drugipv^2)$alpha
  drugip.beta <- estBetaParams(dat.cresamp$drugipp, dat.cresamp$drugipv^2)$beta
  drugip.stch <- rbeta(length(drugip.alpha), drugip.alpha, drugip.beta)
  drugip.stoch <- ifelse(is.na(drugip.stch==T), 0, drugip.stch) # drugs for intestinal parasites
  
  pnc.alpha <- estBetaParams(dat.cresamp$pncp, dat.cresamp$pncv^2)$alpha
  pnc.beta <- estBetaParams(dat.cresamp$pncp, dat.cresamp$pncv^2)$beta
  pnc.stch <- rbeta(length(pnc.alpha), pnc.alpha, pnc.beta)
  pnc.stoch <- ifelse(is.na(pnc.stch==T), 0, pnc.stch) # postnatal checks probability
  
  ## collect resampled variables into single dataframe
  dat.resamp <- data.frame(diar.stoch,wedu.stoch,hhmem.stoch,pnc.stoch,fem.stoch,
                           drugip.stoch,TAR.stoch,PRCP.stoch)
  colnames(dat.resamp) <- c("diar", "wedu", "hhmem", "pnc", "fem", "drugip", "TAR", "PRCP")
  
  # scale
  dat.resamp.sc <- na.omit(as.data.frame(scale(dat.resamp, center=T, scale=T)))
  head(dat.resamp.sc)
  
  
  # boosted regression tree
  brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:8)],
                           gbm.y = attr(dat.resamp.sc, "names")[1],
                           family="gaussian", max.trees=1000000, tolerance = 0.00001,
                           learning.rate = 0.001, bag.fraction=0.5, tree.complexity = 2,
                           silent=T, tolerance.method = "auto", step.size=40)
  
  if (i == 1 & is.null(brt.fit.rsmp)==F) {
    brt.fit.rsmp.old <- brt.fit.rsmp
  }
  
  
  if (is.null(brt.fit.rsmp) == T) {
    brt.fit.rsmp <- gbm.step(dat.resamp.sc, gbm.x = attr(dat.resamp.sc, "names")[c(2:8)],
                             gbm.y = attr(dat.resamp.sc, "names")[1],
                             family="gaussian", max.trees=1000000, tolerance = 0.000001,
                             learning.rate = 0.0001, bag.fraction=0.5, tree.complexity = 2,
                             silent=T, tolerance.method = "auto", step.size=20)
  }
  if (is.null(brt.fit.rsmp) == T) {
    brt.fit.rsmp <- brt.fit.rsmp.old
  }
  summ.fit <- summary(brt.fit.rsmp)
  
  if (is.null(brt.fit.rsmp) == F) {
    brt.fit.rsmp.old <- brt.fit.rsmp
  }
  
  # variable relative importance
  wedu.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[1])]
  hhmem.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[2])]
  pnc.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[3])]
  fem.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[4])]
  drugip.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[5])]
  TAR.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[6])]
  PRCP.ri[i] <- summ.fit$rel.inf[which(summ.fit$var == out.colnames[7])]
  
  # goodness of fit
  CV.cor.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.mean
  CV.cor.se.vec[i] <- 100*brt.fit.rsmp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=7)
  for (p in 1:7) {
    RESP.val[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit.rsmp, i.var=p, continuous.resolution = eq.sp.points, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit.rsmp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit.rsmp$var.names
  
  val.arr[, , i] <- as.matrix(RESP.val.dat)
  pred.arr[, , i] <- as.matrix(RESP.pred.dat)
  
  print(i)
  
} # end i

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 7
pred.update <- pred.arr[,,1:iter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:iter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

PHASE5.pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
PHASE5.pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
PHASE5.pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
PHASE5.val.med <- apply(val.arr[,,1:iter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
PHASE5.CV.cor.update <- CV.cor.vec[1:iter]
PHASE5.CV.cor.se.update <- CV.cor.se.vec[1:iter]

wedu.ri.update <- wedu.ri[1:iter]
hhmem.ri.update <- hhmem.ri[1:iter]
pnc.ri.update <- pnc.ri[1:iter]
fem.ri.update <- fem.ri[1:iter]
drugip.ri.update <- drugip.ri[1:iter]
TAR.ri.update <- TAR.ri[1:iter]
PRCP.ri.update <- PRCP.ri[1:iter]

for (k in 1:kappa.n) {
  PHASE5.CV.cor.mean <- mean(PHASE5.CV.cor.update, na.rm=T); PHASE5.CV.cor.sd <- sd(PHASE5.CV.cor.update, na.rm=T)
  PHASE5.CV.cor.se.mean <- mean(PHASE5.CV.cor.se.update, na.rm=T); PHASE5.CV.cor.se.sd <- sd(PHASE5.CV.cor.se.update, na.rm=T)
  
  wedu.mean <- mean(wedu.ri.update, na.rm=T); wedu.sd <- sd(wedu.ri.update, na.rm=T)
  hhmem.mean <- mean(hhmem.ri.update, na.rm=T); hhmem.sd <- sd(hhmem.ri.update, na.rm=T)
  pnc.mean <- mean(pnc.ri.update, na.rm=T); pnc.sd <- sd(pnc.ri.update, na.rm=T)
  fem.mean <- mean(fem.ri.update, na.rm=T); fem.sd <- sd(fem.ri.update, na.rm=T)
  drugip.mean <- mean(drugip.ri.update, na.rm=T); drugip.sd <- sd(drugip.ri.update, na.rm=T)
  TAR.mean <- mean(TAR.ri.update, na.rm=T); TAR.sd <- sd(TAR.ri.update, na.rm=T)
  PRCP.mean <- mean(PRCP.ri.update, na.rm=T); PRCP.sd <- sd(PRCP.ri.update, na.rm=T)
  
  for (u in 1:iter) {
    PHASE5.CV.cor.update[u] <- ifelse((PHASE5.CV.cor.update[u] < (PHASE5.CV.cor.mean-kappa*PHASE5.CV.cor.sd) | PHASE5.CV.cor.update[u] >
                                         (PHASE5.CV.cor.mean+kappa*PHASE5.CV.cor.sd)), NA, PHASE5.CV.cor.update[u])
    PHASE5.CV.cor.se.update[u] <- ifelse((PHASE5.CV.cor.se.update[u] < (PHASE5.CV.cor.se.mean-kappa*PHASE5.CV.cor.se.sd) | PHASE5.CV.cor.se.update[u]
                                          > (PHASE5.CV.cor.se.mean+kappa*PHASE5.CV.cor.se.sd)), NA, PHASE5.CV.cor.se.update[u])
    
    wedu.ri.update[u] <- ifelse((wedu.ri.update[u] < (wedu.mean-kappa*wedu.sd) | wedu.ri.update[u] > (wedu.mean+kappa*wedu.sd)), NA, wedu.ri.update[u])
    hhmem.ri.update[u] <- ifelse((hhmem.ri.update[u] < (hhmem.mean-kappa*wedu.sd) | hhmem.ri.update[u] > (hhmem.mean+kappa*wedu.sd)), NA, hhmem.ri.update[u])
    pnc.ri.update[u] <- ifelse((pnc.ri.update[u] < (pnc.mean-kappa*pnc.sd) | pnc.ri.update[u] > (pnc.mean+kappa*pnc.sd)), NA, pnc.ri.update[u])
    fem.ri.update[u] <- ifelse((fem.ri.update[u] < (fem.mean-kappa*fem.sd) | fem.ri.update[u] > (fem.mean+kappa*fem.sd)), NA, fem.ri.update[u])
    drugip.ri.update[u] <- ifelse((drugip.ri.update[u] < (drugip.mean-kappa*drugip.sd) | drugip.ri.update[u] > (drugip.mean+kappa*drugip.sd)), NA, drugip.ri.update[u])
    TAR.ri.update[u] <- ifelse((TAR.ri.update[u] < (TAR.mean-kappa*TAR.sd) | TAR.ri.update[u] > (TAR.mean+kappa*TAR.sd)), NA, TAR.ri.update[u])
    PRCP.ri.update[u] <- ifelse((PRCP.ri.update[u] < (PRCP.mean-kappa*PRCP.sd) | PRCP.ri.update[u] > (PRCP.mean+kappa*PRCP.sd)), NA, PRCP.ri.update[u])
    
  } # end u
  
  print(k)
  
} # end k

# summaries
PHASE5.CV.cor.med <- median(PHASE5.CV.cor.update, na.rm=TRUE)
PHASE5.CV.cor.lo <- quantile(PHASE5.CV.cor.update, probs=0.025, na.rm=TRUE)
PHASE5.CV.cor.up <- quantile(PHASE5.CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(PHASE5.CV.cor.lo,PHASE5.CV.cor.med,PHASE5.CV.cor.up))

wedu.ri.lo <- quantile(wedu.ri.update, probs=0.025, na.rm=TRUE)
wedu.ri.med <- median(wedu.ri.update, na.rm=TRUE)
wedu.ri.up <- quantile(wedu.ri.update, probs=0.975, na.rm=TRUE)

hhmem.ri.lo <- quantile(hhmem.ri.update, probs=0.025, na.rm=TRUE)
hhmem.ri.med <- median(hhmem.ri.update, na.rm=TRUE)
hhmem.ri.up <- quantile(hhmem.ri.update, probs=0.975, na.rm=TRUE)

pnc.ri.lo <- quantile(pnc.ri.update, probs=0.025, na.rm=TRUE)
pnc.ri.med <- median(pnc.ri.update, na.rm=TRUE)
pnc.ri.up <- quantile(pnc.ri.update, probs=0.975, na.rm=TRUE)

fem.ri.lo <- quantile(fem.ri.update, probs=0.025, na.rm=TRUE)
fem.ri.med <- median(fem.ri.update, na.rm=TRUE)
fem.ri.up <- quantile(fem.ri.update, probs=0.975, na.rm=TRUE)

drugip.ri.lo <- quantile(drugip.ri.update, probs=0.025, na.rm=TRUE)
drugip.ri.med <- median(drugip.ri.update, na.rm=TRUE)
drugip.ri.up <- quantile(drugip.ri.update, probs=0.975, na.rm=TRUE)

TAR.ri.lo <- quantile(TAR.ri.update, probs=0.025, na.rm=TRUE)
TAR.ri.med <- median(TAR.ri.update, na.rm=TRUE)
TAR.ri.up <- quantile(TAR.ri.update, probs=0.975, na.rm=TRUE)

PRCP.ri.lo <- quantile(PRCP.ri.update, probs=0.025, na.rm=TRUE)
PRCP.ri.med <- median(PRCP.ri.update, na.rm=TRUE)
PRCP.ri.up <- quantile(PRCP.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(wedu.ri.lo,hhmem.ri.lo,pnc.ri.lo,fem.ri.lo,drugip.ri.lo,TAR.ri.lo,PRCP.ri.lo)
ri.med <- c(wedu.ri.med,hhmem.ri.med,pnc.ri.med,fem.ri.med,drugip.ri.med,TAR.ri.med,PRCP.ri.med)
ri.up <- c(wedu.ri.up,hhmem.ri.up,pnc.ri.up,fem.ri.up,drugip.ri.up,TAR.ri.up,PRCP.ri.up)
ri.out <- as.data.frame(cbind(ri.med,ri.up,ri.lo))
rownames(ri.out) <- out.colnames
PHASE5.ri.sort <- ri.out[order(ri.out[,1], decreasing=T), ]
PHASE5.ri.sort

PHASE5.ri.plt <- ggplot(PHASE5.ri.sort) +
  geom_bar(aes(x=reorder(row.names(PHASE5.ri.sort), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar( aes(x=row.names(PHASE5.ri.sort), ymin=ri.lo, ymax=ri.up),
                 width=0.4, colour="black", alpha=0.9, size=0.3)
PHASE5.ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(PHASE5.CV.cor.lo,PHASE5.CV.cor.med,PHASE5.CV.cor.up))


## plot predicted relationships
par(mfrow=c(2,2))
plot(RESP.val.dat$TAR, PHASE5.pred.med[,6], type="l", lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="temperature annual range",
     ylim=c(min(PHASE5.pred.lo[,6]), max(PHASE5.pred.up[,6])))
polygon(c(RESP.val.dat$TAR, rev(RESP.val.dat$TAR)), c(PHASE5.pred.up[,6], rev(PHASE5.pred.lo[,6])),
        col="lightblue", density=30)

plot(RESP.val.dat$PRCP, PHASE5.pred.med[,7], type="l", lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="annual precipitation",
     ylim=c(min(PHASE5.pred.lo[,7]), max(PHASE5.pred.up[,7])))
polygon(c(RESP.val.dat$PRCP, rev(RESP.val.dat$PRCP)), c(PHASE5.pred.up[,7], rev(PHASE5.pred.lo[,7])),
        col="lightblue", density=30)

plot(RESP.val.dat$wedu, PHASE5.pred.med[,1], type="l", lwd=2, lty=2, 
     ylab="pr diarrhoea", xlab="woman's education",
     ylim=c(min(PHASE5.pred.lo[,1]), max(PHASE5.pred.up[,1])))
polygon(c(RESP.val.dat$wedu, rev(RESP.val.dat$wedu)), c(PHASE5.pred.up[,1], rev(PHASE5.pred.lo[,1])),
        col="lightblue", density=30)

plot(RESP.val.dat$hhmem, PHASE5.pred.med[,2], type="l",  lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="household members",
     ylim=c(min(PHASE5.pred.lo[,2]), max(PHASE5.pred.up[,2])))
polygon(c(RESP.val.dat$hhmem, rev(RESP.val.dat$hhmem)), c(PHASE5.pred.up[,2], rev(PHASE5.pred.lo[,2])),
        col="lightblue", density=30)
par(mfrow=c(1,1))



plot(RESP.val.dat$fem, PHASE5.pred.med[,4], type="l",  lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="pr female child",
     ylim=c(min(PHASE5.pred.lo[,4]), max(PHASE5.pred.up[,4])))
polygon(c(RESP.val.dat$fem, rev(RESP.val.dat$fem)), c(PHASE5.pred.up[,4], rev(PHASE5.pred.lo[,4])),
        col="lightblue", density=30)


## transform back to original scales
head(PHASE5vars.naomit[,c(5,7,9,11,13,15,17)])
PHASE5vars.sc <- scale(PHASE5vars.naomit[,c(5,7,9,11,13,15,17)], center=T, scale=T)
centres <- as.numeric(attr(PHASE5vars.sc,"scaled:center"))
scales <- as.numeric(attr(PHASE5vars.sc,"scaled:scale"))
mins <- as.numeric(apply(PHASE5vars.naomit[,c(5,7,9,11,13,15,17)], MARGIN=2, min, na.rm=T))
maxs <- as.numeric(apply(PHASE5vars.naomit[,c(5,7,9,11,13,15,17)], MARGIN=2, max, na.rm=T))

TARval.bsc <- (scales[6] * RESP.val.dat$TAR) + centres[6]
PRCPval.bsc <- (scales[7] * RESP.val.dat$PRCP) + centres[7]
weduval.bsc <- (scales[1] * RESP.val.dat$wedu) + centres[1]
hhmemval.bsc <- (scales[2] * RESP.val.dat$hhmem) + centres[2]

diarp.sc <- scale(PHASE5vars.naomit$diarp, center=T, scale=T)
diarp.centre <- as.numeric(attr(diarp.sc,"scaled:center"))
diarp.scale <- as.numeric(attr(diarp.sc,"scaled:scale"))
diarp.min <- min(PHASE5vars.naomit$diarp, na.rm=T)
diarp.max <- max(PHASE5vars.naomit$diarp, na.rm=T)

TARdiarp.med.bsc <- as.numeric((diarp.scale * PHASE5.pred.med[,6]) + diarp.centre)
TARdiarp.lo.bsc <- as.numeric((diarp.scale * PHASE5.pred.lo[,6]) + diarp.centre)
TARdiarp.up.bsc <- as.numeric((diarp.scale * PHASE5.pred.up[,6]) + diarp.centre)

PRCPdiarp.med.bsc <- as.numeric((diarp.scale * PHASE5.pred.med[,7]) + diarp.centre)
PRCPdiarp.lo.bsc <- as.numeric((diarp.scale * PHASE5.pred.lo[,7]) + diarp.centre)
PRCPdiarp.up.bsc <- as.numeric((diarp.scale * PHASE5.pred.up[,7]) + diarp.centre)

wedudiarp.med.bsc <- as.numeric((diarp.scale * PHASE5.pred.med[,1]) + diarp.centre)
wedudiarp.lo.bsc <- as.numeric((diarp.scale * PHASE5.pred.lo[,1]) + diarp.centre)
wedudiarp.up.bsc <- as.numeric((diarp.scale * PHASE5.pred.up[,1]) + diarp.centre)

hhmemdiarp.med.bsc <- as.numeric((diarp.scale * PHASE5.pred.med[,2]) + diarp.centre)
hhmemdiarp.lo.bsc <- as.numeric((diarp.scale * PHASE5.pred.lo[,2]) + diarp.centre)
hhmemdiarp.up.bsc <- as.numeric((diarp.scale * PHASE5.pred.up[,2]) + diarp.centre)

ReScale <- function(x,first,last){(last-first)/(max(x)-min(x))*(x-min(x))+first}

TARval.bscrs <- ReScale(TARval.bsc, mins[6], maxs[6])
PRCPval.bscrs <- ReScale(PRCPval.bsc, mins[7], maxs[7])
weduval.bscrs <- ReScale(weduval.bsc, mins[1], maxs[1])
hhmemval.bscrs <- ReScale(hhmemval.bsc, mins[2], maxs[2])

par(mfrow=c(2,2))
plot(TARval.bscrs, TARdiarp.med.bsc, type="l", lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="temperature annual range (°C)",
     ylim=c(min(TARdiarp.lo.bsc), max(TARdiarp.up.bsc)))
polygon(c(TARval.bscrs, rev(TARval.bscrs)), c(TARdiarp.up.bsc, rev(TARdiarp.lo.bsc)),
        col="lightblue", density=30)

plot(PRCPval.bscrs, PRCPdiarp.med.bsc, type="l", lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="annual precipitation (mm)",
     ylim=c(min(PRCPdiarp.lo.bsc), max(PRCPdiarp.up.bsc)))
polygon(c(PRCPval.bscrs, rev(PRCPval.bscrs)), c(PRCPdiarp.up.bsc, rev(PRCPdiarp.lo.bsc)),
        col="lightblue", density=30)

plot(weduval.bscrs, wedudiarp.med.bsc, type="l", lwd=2, lty=2, 
     ylab="pr diarrhoea", xlab="woman's education (years)",
     ylim=c(min(wedudiarp.lo.bsc), max(wedudiarp.up.bsc)))
polygon(c(weduval.bscrs, rev(weduval.bscrs)), c(wedudiarp.up.bsc, rev(wedudiarp.lo.bsc)),
        col="lightblue", density=30)

plot(hhmemval.bscrs, hhmemdiarp.med.bsc, type="l",  lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="household members",
     ylim=c(min(hhmemdiarp.lo.bsc), max(hhmemdiarp.up.bsc)))
polygon(c(hhmemval.bscrs, rev(hhmemval.bscrs)), c(hhmemdiarp.up.bsc, rev(hhmemdiarp.lo.bsc)),
        col="lightblue", density=30)
par(mfrow=c(1,1))

## all with y axis on same scale
yaxs.min <- min(TARdiarp.lo.bsc,PRCPdiarp.lo.bsc,wedudiarp.lo.bsc,hhmemdiarp.lo.bsc)
yaxs.max <- max(TARdiarp.up.bsc,PRCPdiarp.up.bsc,wedudiarp.up.bsc,hhmemdiarp.up.bsc)

par(mfrow=c(2,2))
plot(TARval.bscrs, TARdiarp.med.bsc, type="l", lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="temperature annual range (°C)",
     ylim=c(yaxs.min, yaxs.max))
polygon(c(TARval.bscrs, rev(TARval.bscrs)), c(TARdiarp.up.bsc, rev(TARdiarp.lo.bsc)),
        col="lightblue", density=30)

plot(PRCPval.bscrs, PRCPdiarp.med.bsc, type="l", lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="annual precipitation (mm)",
     ylim=c(yaxs.min, yaxs.max))
polygon(c(PRCPval.bscrs, rev(PRCPval.bscrs)), c(PRCPdiarp.up.bsc, rev(PRCPdiarp.lo.bsc)),
        col="lightblue", density=30)

plot(weduval.bscrs, wedudiarp.med.bsc, type="l", lwd=2, lty=2, 
     ylab="pr diarrhoea", xlab="woman's education (years)",
     ylim=c(yaxs.min, yaxs.max))
polygon(c(weduval.bscrs, rev(weduval.bscrs)), c(wedudiarp.up.bsc, rev(wedudiarp.lo.bsc)),
        col="lightblue", density=30)

plot(hhmemval.bscrs, hhmemdiarp.med.bsc, type="l",  lwd=2, lty=2,
     ylab="pr diarrhoea", xlab="household members",
     ylim=c(yaxs.min, yaxs.max))
polygon(c(hhmemval.bscrs, rev(hhmemval.bscrs)), c(hhmemdiarp.up.bsc, rev(hhmemdiarp.lo.bsc)),
        col="lightblue", density=30)
par(mfrow=c(1,1))


save.image('DHSresampledBRT.RData')
