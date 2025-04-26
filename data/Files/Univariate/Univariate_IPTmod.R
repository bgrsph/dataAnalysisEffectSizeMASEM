
#IPT Trials Univariate Meta (Study 2) 
#   + additional analyses (sensitivity, publication bias, plots)
#   + moderator analyses by treatment type (CBT v. IPT)

#Load metafor
library(metafor)
library(dplyr)

#Load full dataset--select unimeta.cvs when window opens
unimeta1<-read.csv(file.choose())

#Create independent matrices dataset for the CME analyses

#Flag with 0 the less representative matrices first
unimeta1$indep <- 1
unimeta1$indep[(unimeta1$studyid==78) & (unimeta1$TCid==1)] <- 0
unimeta1$indep[(unimeta1$studyid==463) & (unimeta1$TCid==1)] <- 0
unimeta1$indep[(unimeta1$studyid==844) & (unimeta1$TCid==2)] <- 0
unimeta1$indep[(unimeta1$studyid==10065) & (unimeta1$TCid==2)] <- 0
unimeta1$indep[(unimeta1$studyid==10114) & (unimeta1$TCid==2)] <- 0

#Filter in indep matrices
unimeta2 <- filter(unimeta1, indep==1) 



#######IPT Path a: tx-cm

#Fit random-effects models
#   to the full data unimeta1 Using default estimator REML
#   to the full data unimeta1, with Trim & Fill to estimate the number of missing NULL studies
#       and with estimated NULL datapoints added to the observed data
#Compute total no. of participants for each analysis

#MA = negative cognition

IPT_MA <- rma(yi, vi, data=unimeta1, subset=(IPT==1 & VarPair=="0_1_3_MA"), slab = paste(author, comparison, sep = ", "))
summary(IPT_MA)
N_IPT_MA <- filter(unimeta1, IPT==1 & VarPair=="0_1_3_MA")
sum(N_IPT_MA$n_r1, na.rm=TRUE)

trimfill(IPT_MA) 
funnel(trimfill(IPT_MA))

#ME = social engagement

IPT_ME <- rma(yi, vi, data=unimeta1, subset=(IPT==1 &  VarPair=="0_1_3_ME"), slab = paste(author, comparison, sep = ", "))
summary(IPT_ME)
N_IPT_ME <- filter(unimeta1, IPT==1 & VarPair=="0_1_3_ME")
sum(N_IPT_ME$n_r1, na.rm=TRUE)

trimfill(IPT_ME) 
funnel(trimfill(IPT_ME))

#ML2 = family functioning

IPT_ML2 <- rma(yi, vi, data=unimeta1, subset=(IPT==1 &  VarPair=="0_1_3_ML2"), slab = paste(author, comparison, sep = ", "))
summary(IPT_ML2)
N_IPT_ML2 <- filter(unimeta1, IPT==1 & VarPair=="0_1_3_ML2")
sum(N_IPT_ML2$n_r1, na.rm=TRUE)

trimfill(IPT_ML2) 
funnel(trimfill(IPT_ML2))



#######IPT Path c: tx-outcome

#Fit random-effects models
#   to the full data unimeta1 Using default estimator REML
#   to the full data unimeta1, with Trim & Fill to estimate the number of missing NULL studies
#       and with estimated NULL datapoints added to the observed data
#Compute total no. of participants for each analysis

#MA = negative cognition

IPT_MAO <- rma(yi, vi, data=unimeta1, subset=(IPT==1 & CMA_O==1), slab = paste(author, comparison, sep = ", "))
summary(IPT_MAO)
N_IPT_MAO <- filter(unimeta1, IPT==1 & CMA_O==1)
sum(N_IPT_MAO$n_r1, na.rm=TRUE)

trimfill(IPT_MAO) 
funnel(trimfill(IPT_MAO))

#ME = social engagement

IPT_MEO <- rma(yi, vi, data=unimeta1, subset=(IPT==1 & CME_O==1), slab = paste(author, comparison, sep = ", "))
summary(IPT_MEO)
N_IPT_MEO <- filter(unimeta1, IPT==1 & CME_O==1)
sum(N_IPT_MEO$n_r1, na.rm=TRUE)

trimfill(IPT_MEO) 
funnel(trimfill(IPT_MEO))

#ML2 = family functioning DOES NOT MATCH UP in meta result, N matches

IPT_ML2O <- rma(yi, vi, data=unimeta1, subset=(IPT==1 & CML2_O==1), slab = paste(author, comparison, sep = ", "))
summary(IPT_ML2O) 
N_IPT_ML2O <- filter(unimeta1, IPT==1 & CML2_O==1)
sum(N_IPT_ML2O$n_r1, na.rm=TRUE)

trimfill(IPT_ML2O) 
funnel(trimfill(IPT_ML2O))



#######CBT Path b: cm-outcome

#Compute sampling variances
#   for full data unimeta1
#   for independent comparisons subset unimeta2 for sensitivity analyses
dat1 <- escalc(measure="COR", ri=rga, ni=n_r1, vtype = "LS",
               data=unimeta1, append = TRUE, replace = FALSE)
dat2 <- escalc(measure="COR", ri=rga, ni=n_r1, vtype = "LS",
               data=unimeta2, append = TRUE, replace = FALSE)

#Fit random-effects models
#   to the full data unimeta1 Using default estimator REML
#   to the full data unimeta1, with Trim & Fill to estimate the number of missing NULL studies
#       and with estimated NULL datapoints added to the observed data
#Compute total no. of participants for each analysis

#MA = negative cognition

MAO_b_IPT <- rma(yi, vi, data=dat1, subset=(IPT==1 & VarPair=="3_MA_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MAO_b_IPT)
N_IPT_MAO_b <- filter(dat1, IPT==1 & VarPair=="3_MA_3_O")
sum(N_IPT_MAO_b$n_r1, na.rm=TRUE)

#ME = social engagement

MEO_b_IPT <- rma(yi, vi, data=dat1, subset=(IPT==1 & VarPair=="3_ME_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MEO_b_IPT)
N_IPT_MEO_b <- filter(dat1, IPT==1 & VarPair=="3_ME_3_O")
sum(N_IPT_MEO_b$n_r1, na.rm=TRUE)

#ML2 = family functioning

ML2O_b_IPT <- rma(yi, vi, data=dat1, subset=(IPT==1 & VarPair=="3_ML2_3_O"), slab = paste(author, comparison, sep = ", "))
summary(ML2O_b_IPT)
N_IPT_ML2O_b <- filter(dat1, IPT==1 & VarPair=="3_ML2_3_O")
sum(N_IPT_ML2O_b$n_r1, na.rm=TRUE)



#######CBT v. IPT Path a: tx-cm

#Fit mixed effects moderator meta-analysis
#   to the full data unimeta1 Using default estimator REML for CMA and CML2
#   to the independent comparisons subset unimeta2 for CME
#Compute total no. of participants for each analysis

#MA = negative cognition

Mod_MA <- rma(yi, vi, data=unimeta1, subset=(VarPair=="0_1_3_MA"), slab = paste(author, comparison, sep = ", "), mods = ~ IPT)
summary(Mod_MA) 
N_mod_MA <- filter(unimeta1, VarPair=="0_1_3_MA")
sum(N_mod_MA$n_r1, na.rm=TRUE)

#ME = social engagement

Mod_ME <- rma(yi, vi, data=unimeta2, subset=(VarPair=="0_1_3_ME"), slab = paste(author, comparison, sep = ", "), mods = ~ IPT)
summary(Mod_ME) 
N_mod_ME <- filter(unimeta2, VarPair=="0_1_3_ME")
sum(N_mod_ME$n_r1, na.rm=TRUE)

#ML2 = family functioning

Mod_ML2 <- rma(yi, vi, data=unimeta1, subset=(VarPair=="0_1_3_ML2"), slab = paste(author, comparison, sep = ", "), mods = ~ IPT)
summary(Mod_ML2) 
N_mod_ML2 <- filter(unimeta1, VarPair=="0_1_3_ML2")
sum(N_mod_ML2$n_r1, na.rm=TRUE)



#######CBT v. IPT Path b: cm-outcome

#Compute sampling variances
#   for full data unimeta1
#   for independent comparisons subset unimeta2 for sensitivity analyses
dat1 <- escalc(measure="COR", ri=rga, ni=n_r1, vtype = "LS",
               data=unimeta1, append = TRUE)
dat2 <- escalc(measure="COR", ri=rga, ni=n_r1, vtype = "LS",
               data=unimeta2, append = TRUE)

#Fit mixed effects moderator meta-analysis
#   to the full data unimeta1 Using default estimator REML for CMA and CML2
#   to the independent comparisons subset unimeta2 for CME
#Compute total no. of participants for each analysis

#MA = negative cognition

Mod_MAO_b <- rma(yi, vi, data=dat1, subset=(VarPair=="3_MA_3_O"), slab = paste(author, comparison, sep = ", "), mods = ~ IPT)
summary(Mod_MAO_b)
N_mod_MAO_b <- filter(dat1, VarPair=="3_MA_3_O")
sum(N_mod_MAO_b$n_r1, na.rm=TRUE)

#ME = social engagement

Mod_MEO_b <- rma(yi, vi, data=dat2, subset=(VarPair=="3_ME_3_O"), slab = paste(author, comparison, sep = ", "), mods = ~ IPT)
summary(Mod_MEO_b)
N_mod_MEO_b <- filter(dat2, VarPair=="3_ME_3_O")
sum(N_mod_MEO_b$n_r1, na.rm=TRUE)

#ML2 = family functioning

Mod_ML2O_b <- rma(yi, vi, data=dat1, subset=(VarPair=="3_ML2_3_O"), slab = paste(author, comparison, sep = ", "), mods = ~ IPT)
summary(Mod_ML2O_b)
N_mod_ML2O_b <- filter(dat1, VarPair=="3_ML2_3_O")
sum(N_mod_ML2O_b$n_r1, na.rm=TRUE)

