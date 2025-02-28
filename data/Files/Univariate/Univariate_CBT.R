
#CBT Trials Univariate Meta-analysis (Study 1) 
#   + additional analyses (sensitivity, publication bias, plots)


#Load metafor
library(metafor)
library(dplyr)

#Load full dataset--select unimeta.cvs when window opens
unimeta1<-read.csv(file.choose())

#Create subset of only independent comparisons
#   Flag with 1 all comparisons
#   Then change to 0 the dependent and less representative comparisons
unimeta1$indep <- 1
unimeta1$indep[(unimeta1$studyid==78) & (unimeta1$TCid==1)] <- 0
unimeta1$indep[(unimeta1$studyid==463) & (unimeta1$TCid==1)] <- 0
unimeta1$indep[(unimeta1$studyid==844) & (unimeta1$TCid==2)] <- 0
unimeta1$indep[(unimeta1$studyid==10065) & (unimeta1$TCid==2)] <- 0
unimeta1$indep[(unimeta1$studyid==10114) & (unimeta1$TCid==2)] <- 0

#Filter in indep comparisons (remove dependent comparisons flagged above)
unimeta2 <- filter(unimeta1, indep==1) 



#######CBT Path a: tx-cm

#Fit random-effects models
#   to the full data unimeta1 Using default estimator REML
#   to the full data unimeta1, with Trim & Fill to estimate the number of missing NULL studies
#       and with estimated NULL datapoints added to the observed data
#   to independent comparisons subset unimeta2 for sensitivity analyses
#Compute total no. of participants for each analysis

#MA = negative cognition

CBT_MA <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & VarPair=="0_1_3_MA"), slab = paste(author, comparison, sep = ", "))
summary(CBT_MA) 
N_CBT_MA <- filter(unimeta1, IPT==0 & VarPair=="0_1_3_MA")
sum(N_CBT_MA$n_r1, na.rm=TRUE) 

trimfill(CBT_MA) 
funnel(trimfill(CBT_MA))

CBT_MAi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & VarPair=="0_1_3_MA"))
summary(CBT_MAi)
N_CBT_MAi <- filter(unimeta2, IPT==0 & VarPair=="0_1_3_MA")
sum(N_CBT_MAi$n_r1, na.rm=TRUE)

#ME = social engagement

CBT_ME <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & VarPair=="0_1_3_ME"), slab = paste(author, comparison, sep = ", "))
summary(CBT_ME)
N_CBT_ME <- filter(unimeta1, IPT==0 & VarPair=="0_1_3_ME")
sum(N_CBT_ME$n_r1, na.rm=TRUE)

trimfill(CBT_ME) 
funnel(trimfill(CBT_ME))

CBT_MEi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & VarPair=="0_1_3_ME"))
summary(CBT_MEi)
N_CBT_MEi <- filter(unimeta2, IPT==0 & VarPair=="0_1_3_ME")
sum(N_CBT_MEi$n_r1, na.rm=TRUE)

#ML2 = family functioning

CBT_ML2 <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & VarPair=="0_1_3_ML2"), slab = paste(author, comparison, sep = ", "))
summary(CBT_ML2)
N_CBT_ML2 <- filter(unimeta1, IPT==0 & VarPair=="0_1_3_ML2")
sum(N_CBT_ML2$n_r1, na.rm=TRUE) 

trimfill(CBT_ML2) 
funnel(trimfill(CBT_ML2))

CBT_ML2i <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & VarPair=="0_1_3_ML2"))
summary(CBT_ML2i)
N_CBT_ML2i <- filter(unimeta2, IPT==0 & VarPair=="0_1_3_ML2")
sum(N_CBT_ML2i$n_r1, na.rm=TRUE) 


#MC = problem solving

CBT_MC <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & VarPair=="0_1_3_MC"), slab = paste(author, comparison, sep = ", "))
summary(CBT_MC) 
N_CBT_MC <- filter(unimeta1, IPT==0 & VarPair=="0_1_3_MC")
sum(N_CBT_MC$n_r1, na.rm=TRUE) 

trimfill(CBT_MC) 
funnel(trimfill(CBT_MC))

CBT_MCi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & VarPair=="0_1_3_MC"))
summary(CBT_MCi) 
N_CBT_MCi <- filter(unimeta2, IPT==0 & VarPair=="0_1_3_MC")
sum(N_CBT_MCi$n_r1, na.rm=TRUE) 

#MB = pleasant activities

CBT_MB <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & VarPair=="0_1_3_MB"), slab = paste(author, comparison, sep = ", "))
summary(CBT_MB)
N_CBT_MB <- filter(unimeta1, IPT==0 & VarPair=="0_1_3_MB")
sum(N_CBT_MB$n_r1, na.rm=TRUE)

trimfill(CBT_MB) 
funnel(trimfill(CBT_MB))

CBT_MBi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & VarPair=="0_1_3_MB"))
summary(CBT_MBi)
N_CBT_MBi <- filter(unimeta2, IPT==0 & VarPair=="0_1_3_MB")
sum(N_CBT_MBi$n_r1, na.rm=TRUE)

#MD = reframing

CBT_MD <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & VarPair=="0_1_3_MD"), slab = paste(author, comparison, sep = ", "))
summary(CBT_MD) 
N_CBT_MD <- filter(unimeta1, IPT==0 & VarPair=="0_1_3_MD")
sum(N_CBT_MD$n_r1, na.rm=TRUE) 

trimfill(CBT_MD) 
funnel(trimfill(CBT_MD))

CBT_MDi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & VarPair=="0_1_3_MD"))
summary(CBT_MDi) 
N_CBT_MDi <- filter(unimeta2, IPT==0 & VarPair=="0_1_3_MD")
sum(N_CBT_MDi$n_r1, na.rm=TRUE)

#ML29 = avoidance

CBT_ML29 <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & VarPair=="0_1_3_ML29"), slab = paste(author, comparison, sep = ", "))
summary(CBT_ML29) 
N_CBT_ML29 <- filter(unimeta1, IPT==0 & VarPair=="0_1_3_ML29")
sum(N_CBT_ML29$n_r1, na.rm=TRUE)

trimfill(CBT_ML29) 
funnel(trimfill(CBT_ML29))

#all ML29 comparisons are independent so no need to run sensitivity analysis



#######CBT Path c: tx-outcome

#Fit random-effects models
#   to the full data unimeta1 Using default estimator REML
#   to the full data unimeta1, with Trim & Fill to estimate the number of missing NULL studies
#       and with estimated NULL datapoints added to the observed data
#   to independent comparisons subset unimeta2 for sensitivity analyses
#Compute total no. of participants for each analysis

#MA = negative cognition

CBT_MAO <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & CMA_O==1), slab = paste(author, comparison, sep = ", "))
summary(CBT_MAO) 
N_CBT_MAO <- filter(unimeta1, IPT==0 & CMA_O==1)
sum(N_CBT_MAO$n_r1, na.rm=TRUE) 

trimfill(CBT_MAO) 
funnel(trimfill(CBT_MAO))

CBT_MAOi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & CMA_O==1))
summary(CBT_MAOi)
N_CBT_MAOi <- filter(unimeta2, IPT==0 & CMA_O==1)
sum(N_CBT_MAOi$n_r1, na.rm=TRUE)

#ME = social engagement

CBT_MEO <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & CME_O==1), slab = paste(author, comparison, sep = ", "))
summary(CBT_MEO)
N_CBT_MEO <- filter(unimeta1, IPT==0 & CME_O==1)
sum(N_CBT_MEO$n_r1, na.rm=TRUE)

trimfill(CBT_MEO) 
funnel(trimfill(CBT_MEO))

CBT_MEOi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & CME_O==1))
summary(CBT_MEOi)
N_CBT_MEOi <- filter(unimeta2, IPT==0 & CME_O==1)
sum(N_CBT_MEOi$n_r1, na.rm=TRUE)

#ML2 = family functioning

CBT_ML2O <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & CML2_O==1), slab = paste(author, comparison, sep = ", "))
summary(CBT_ML2O)
N_CBT_ML2O <- filter(unimeta1, IPT==0 & CML2_O==1)
sum(N_CBT_ML2O$n_r1, na.rm=TRUE) 

trimfill(CBT_ML2O) 
funnel(trimfill(CBT_ML2O))

CBT_ML2Oi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & CML2_O==1))
summary(CBT_ML2Oi)
N_CBT_ML2Oi <- filter(unimeta2, IPT==0 & CML2_O==1)
sum(N_CBT_ML2Oi$n_r1, na.rm=TRUE) 

#MC = problem solving

CBT_MCO <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & CMC_O==1), slab = paste(author, comparison, sep = ", "))
summary(CBT_MCO) 
N_CBT_MCO <- filter(unimeta1, IPT==0 & CMC_O==1)
sum(N_CBT_MCO$n_r1, na.rm=TRUE) 

trimfill(CBT_MCO) 
funnel(trimfill(CBT_MCO))

CBT_MCOi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & CMC_O==1))
summary(CBT_MCOi) 
N_CBT_MCOi <- filter(unimeta2, IPT==0 & CMC_O==1)
sum(N_CBT_MCOi$n_r1, na.rm=TRUE) 

#MB = pleasant activities

CBT_MBO <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & CMB_O==1), slab = paste(author, comparison, sep = ", "))
summary(CBT_MBO)
N_CBT_MBO <- filter(unimeta1, IPT==0 & CMB_O==1)
sum(N_CBT_MBO$n_r1, na.rm=TRUE)

trimfill(CBT_MBO) 
funnel(trimfill(CBT_MBO))

CBT_MBOi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & CMB_O==1))
summary(CBT_MBOi)
N_CBT_MBOi <- filter(unimeta2, IPT==0 & CMB_O==1)
sum(N_CBT_MBOi$n_r1, na.rm=TRUE)

#MD = reframing

CBT_MDO <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & CMD_O==1), slab = paste(author, comparison, sep = ", "))
summary(CBT_MDO) 
N_CBT_MDO <- filter(unimeta1, IPT==0 & CMD_O==1)
sum(N_CBT_MDO$n_r1, na.rm=TRUE) 

trimfill(CBT_MDO) 
funnel(trimfill(CBT_MDO))

CBT_MDOi <- rma(yi, vi, data=unimeta2, subset=(IPT==0 & CMD_O==1))
summary(CBT_MDOi) 
N_CBT_MDOi <- filter(unimeta2, IPT==0 & CMD_O==1)
sum(N_CBT_MDOi$n_r1, na.rm=TRUE)

#ML29 = avoidance

CBT_ML29O <- rma(yi, vi, data=unimeta1, subset=(IPT==0 & CML29_O==1), slab = paste(author, comparison, sep = ", "))
summary(CBT_ML29O) 
N_CBT_ML29O <- filter(unimeta1, IPT==0 & CML29_O==1)
sum(N_CBT_ML29O$n_r1, na.rm=TRUE)

trimfill(CBT_ML29O) 
funnel(trimfill(CBT_ML29O))

#all ML29 comparisons are independent so no need to run sensitivity analysis



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
#   to independent comparisons subset unimeta2 for sensitivity analyses
#Compute total no. of participants for each analysis

#MA = negative cognition

MAO_b <- rma(yi, vi, data=dat1, subset=(IPT==0 & VarPair=="3_MA_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MAO_b) 
N_CBT_MAO_b <- filter(dat1, IPT==0 & VarPair=="3_MA_3_O")
sum(N_CBT_MAO_b$n_r1, na.rm=TRUE)

MAO_bi <- rma(yi, vi, data=dat2, subset=(IPT==0 & VarPair=="3_MA_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MAO_bi)
N_CBT_MAO_bi <- filter(dat2, IPT==0 & VarPair=="3_MA_3_O")
sum(N_CBT_MAO_bi$n_r1, na.rm=TRUE) 

#ME = social engagement

MEO_b <- rma(yi, vi, data=dat1, subset=(IPT==0 & VarPair=="3_ME_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MEO_b) 
N_CBT_MEO_b <- filter(dat1, IPT==0 & VarPair=="3_ME_3_O")
sum(N_CBT_MEO_b$n_r1, na.rm=TRUE)

MEO_bi <- rma(yi, vi, data=dat2, subset=(IPT==0 & VarPair=="3_ME_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MEO_bi)
N_CBT_MEO_bi <- filter(dat2, IPT==0 & VarPair=="3_ME_3_O")
sum(N_CBT_MEO_bi$n_r1, na.rm=TRUE)

#ML2 = family functioning

ML2O_b <- rma(yi, vi, data=dat1, subset=(IPT==0 & VarPair=="3_ML2_3_O"), slab = paste(author, comparison, sep = ", "))
summary(ML2O_b) 
N_CBT_ML2O_b <- filter(dat1, IPT==0 & VarPair=="3_ML2_3_O")
sum(N_CBT_ML2O_b$n_r1, na.rm=TRUE)

ML2O_bi <- rma(yi, vi, data=dat2, subset=(IPT==0 & VarPair=="3_ML2_3_O"), slab = paste(author, comparison, sep = ", "))
summary(ML2O_bi) 
N_CBT_ML2O_bi <- filter(dat2, IPT==0 & VarPair=="3_ML2_3_O")
sum(N_CBT_ML2O_bi$n_r1, na.rm=TRUE)

#MC = problem solving

MCO_b <- rma(yi, vi, data=dat1, subset=(IPT==0 & VarPair=="3_MC_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MCO_b) 
N_CBT_MCO_b <- filter(dat1, IPT==0 & VarPair=="3_MC_3_O")
sum(N_CBT_MCO_b$n_r1, na.rm=TRUE)

#all available MC comparisons are independent so no need to run sensitivity analysis

#MB = pleasant activities

MBO_b <- rma(yi, vi, data=dat1, subset=(IPT==0 & VarPair=="3_MB_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MBO_b)
N_CBT_MBO_b <- filter(dat1, IPT==0 & VarPair=="3_MB_3_O")
sum(N_CBT_MBO_b$n_r1, na.rm=TRUE)

MBO_bi <- rma(yi, vi, data=dat2, subset=(IPT==0 & VarPair=="3_MB_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MBO_bi)
N_CBT_MBO_bi <- filter(dat2, IPT==0 & VarPair=="3_MB_3_O")
sum(N_CBT_MBO_bi$n_r1, na.rm=TRUE)

#MD = reframing

MDO_b <- rma(yi, vi, data=dat1, subset=(IPT==0 & VarPair=="3_MD_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MDO_b)
N_CBT_MDO_b <- filter(dat1, IPT==0 & VarPair=="3_MD_3_O")
sum(N_CBT_MDO_b$n_r1, na.rm=TRUE)

MDO_bi <- rma(yi, vi, data=dat2, subset=(IPT==0 & VarPair=="3_MD_3_O"), slab = paste(author, comparison, sep = ", "))
summary(MDO_bi)
N_CBT_MDO_bi <- filter(dat2, IPT==0 & VarPair=="3_MD_3_O")
sum(N_CBT_MDO_bi$n_r1, na.rm=TRUE)

#ML29 = avoidance

ML29O_b <- rma(yi, vi, data=dat1, subset=(IPT==0 & VarPair=="3_ML29_3_O"), slab = paste(author, comparison, sep = ", "))
summary(ML29O_b)
N_CBT_ML29O_b <- filter(dat1, IPT==0 & VarPair=="3_ML29_3_O")
sum(N_CBT_ML29O_b$n_r1, na.rm=TRUE)

#all ML29 comparisons are independent so no need to run sensitivity analysis


