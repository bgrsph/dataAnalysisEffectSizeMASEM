
#MASEM Models for Testing Moderated Mediation (Study 2)
#   This syntax can be used to test whether treatment type (CBT v. IPT) 
#   moderates pathways in the the mediation model for any candidate mediator.
#   The mediation model is first fitted to CBT trials and then to IPT trials,
#   then paths a and b are tested for equivalence between the CBT and IPT subgroups.
#   Change file names and/or the working directory where the relevant data files are saved.

#Candidate mediators tested with each treatment, full dataset or subset, as reported in Study 2.
#   MA_CBTfull.txt = negative cognition full dataset of CBT matrices (k=28)
#   MA_IPTfull.txt = negative cognition full dataset of IPT matrices (k=3)
#   ME_CBTidp.txt = social engagement subset of independent CBT matrices (k=14)
#   ME_IPTfull.txt = social engagement full dataset of IPT matrices (k=5)
#   ML2_CBTfull.txt = family functioning full dataset of CBT matrices (k=17)
#   ML2_IPTfull.txt = family functioning full dataset of IPT matrices (k=6)

#For corresponding study sample size data files, append "_n" to the end of the text file.

#Load metaSEM and semPlot
library(metaSEM)
library(semPlot)
library(OpenMx)

## Load the functions to facilitate subgroup analysis
source("http://www.suzannejak.nl/subgroup.functions.R")

#Set working directory(library)--enter path to directory in parentheses
setwd ("C:/Users/meiying/OneDrive - Florida International University/Research/MY MS/2 Review, Revise, Resubmit/Mediation Meta Youth Dep CMs/PB/RESUBMIT/OSF_depmedmeta/Data")

#Read in datasets of correlation matrices 
vector_CBT <- readStackVec ("MA_CBTfull.txt") #change file for different candidate mediator
vector_CBT
vector_IPT <- readStackVec ("MA_IPTfull.txt") #change file for different candidate mediator
vector_IPT

#Give column and row headings, enter columns first, check vector
#   Tx=treatment, m1=pre-treatment candidate mediator, o1=pre-treatment outcome
#   m3=post-treatment candidate mediator, o3=post-treatment outcome
  vector_CBT <- lapply(vector_CBT,function(x)
  {dimnames(x) <- list(c("Tx", "m1", "o1", "m3", "o3"),
                       c("Tx", "m1", "o1", "m3", "o3"))
   x})
  vector_CBT 
  vector_IPT <- lapply(vector_IPT,function(x)
  {dimnames(x) <- list(c("Tx", "m1", "o1", "m3", "o3"),
                       c("Tx", "m1", "o1", "m3", "o3"))
  x})
  vector_IPT 

#Read in corresponding study sample sizes, n
n_CBT <- read.table("MA_CBTfull_n.txt") #change file for different candidate mediator
n_CBT
n_IPT <- read.table("MA_IPTfull_n.txt") #change file for different candidate mediator
n_IPT  

#Display sample sizes in the pooled correlation matrix
#   number of treatment-control comparisons k for each correlation
pattern.na(vector_CBT, show.na=FALSE)
pattern.na(vector_IPT, show.na=FALSE)
#   number of participants n for each correlation
pattern.n(vector_CBT, n_CBT)
pattern.n(vector_IPT, n_IPT)

#Test if matrix is positive definite
#   True = positive definite, False = not positive definite
#   NA = missing element in matrices
is.pd(vector_CBT)
is.pd(vector_IPT)


#TSSEM Stage 1

#CBT random effects model, use diagonal matrix for random effects
random1_CBT <- tssem1(vector_CBT, n_CBT, method="REM", RE.type="Diag")
options("scipen"=100, "digits"=20)
summary(random1_CBT)

#Check OpenMx status1 and rerun analysis if needed
#   OpenMx status1 = 0 or 1 means the estimation is considered fine
#   Other values indicate estimation problems, rerun analysis
random1_CBT <-rerun(random1_CBT)
summary(random1_CBT)

#Check OpenMx status1 and rerun analysis if needed
#   OpenMx status1 = 0 or 1 means the estimation is considered fine
#   Other values indicate estimation problems, check Tau2 for very small heterogeneity variances
#   Fix any Tau2 that is < 1 x 10e-10 to 0 in a user-defined model

#Create user-defined model
#   for ME_CBTidp = social engagement subset of independent CBT matrices (k=14)
#RE <- Diag(c(0,0,"0.01*Tau2_3","0.01*Tau2_4","0.01*Tau2_5","0.01*Tau2_6","0.01*Tau2_7",0,0,"0.01*Tau2_10"))

#Show user-defined model and run it
#RE
#random1_CBT <- tssem1(vector_CBT, n_CBT, method="REM", RE.type="User",
#                  RE.constraints = RE)
#summary(random1_CBT)
#Check that OpenMx status1 = 0.

#IPT random effects model, use diagonal matrix for random effects
random1_IPT <- tssem1(vector_IPT, n_IPT, method="REM", RE.type="Diag")
options("scipen"=100, "digits"=20)
summary(random1_IPT)

#Check OpenMx status1 and rerun analysis if needed
#   OpenMx status1 = 0 or 1 means the estimation is considered fine
#   Other values indicate estimation problems, rerun analysis
random1_IPT <-rerun(random1_IPT)
summary(random1_IPT)


#To extract estimated average correlation matrix in matrix form
#   Select the fixed effects and convert it into a correlation matrix
vec2symMat(coef(random1_CBT, select="fixed"),diag=FALSE)
vec2symMat(coef(random1_IPT, select="fixed"),diag=FALSE)


#TSSEM Stage 2

#Create A, asymmetric matrix to estimate regression coefficients needed in analytic model
#   Regress post-treatment mediator onto treatment condition, label Txm3 (path a)
#   Regress post-treatment outcome onto treatment condition, label Txo3 (path c')
#   Regress post-treatment outcome onto post-treatment mediator, label m3o3 (path b)
#   Regress post-treatment mediator onto pre-treatment mediator, label m1m3
#   Regress post-treatment outcome onto pre-treatment outcome, label o1o3
#   Starting value 0.1 for correlations to be estimated
#   Fix the rest of the correlations not needed to 0
A <- create.mxMatrix(c(0, 0,  0, 0, 0,
                        0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 
                        "0.1*Txm3", "0.1*m1m3", 0, 0, 0, 
                        "0.1*Txo3", 0, "0.1*o1o3", "0.1*m3o3", 0),
                      type="Full", byrow=TRUE, ncol=5, nrow=5, as.mxMatrix=FALSE)
## Label matrix for inspecting the model.
dimnames(A)[[1]] <- dimnames(A)[[2]] <- c("Tx", "m1", "o1", "m3", "o3") 
A

#Create S, symmetric variance covariance matrix among variables in analytic model
#   Starting value 0.1 for variances and covariances to be estimated
#   Set variances of Tx, m1, and o1  (predictor variables) to 1
#   Fix the rest of the covariances to 0
S <- create.mxMatrix(c(1, 
                        0, 1,
                        0, "0.1*r_m1o1", 1,
                        0, 0, 0, "0.1*var_m3",
                        0, 0, 0, 0, "0.1*var_o3"), 
                      byrow=TRUE, type="Symm", as.mxMatrix=FALSE)
## label matrix variables, useful for inspecting the model.
dimnames(S)[[1]] <- dimnames(S)[[2]] <- c("Tx", "m1", "o1", "m3", "o3")
S

#Stage 2 analysis without equality constraints between CBT and IPT models
#    Freely estimate all parameters for each subgroup separately

#CBT random effects model, define effects to test, request likelihood-based confidence intervals
#   indirect effect, product of path a, Txm3, and path b, m3o3
#   total effect, path c, Txo3
#   direct effect, total effect minus mediation effect
random2_CBT <- tssem2(random1_CBT, Amatrix=A, Smatrix=S, intervals.type="LB", model.name="TSSEM2 CBT free",
                  mx.algebras=list( indirect_CBT=mxAlgebra(Txm3*m3o3, name="indirect_CBT"),
                                    direct_CBT=mxAlgebra(Txo3, name="direct_CBT"), 
                                    total_CBT=mxAlgebra(Txo3 + Txm3*m3o3, name="total_CBT") ))
summary(random2_CBT)

#Check OpenMx status1 and error messages, rerun analysis if needed
#   OpenMx status1 = 0 or 1 means the estimation is considered fine
#   Other values indicate estimation problems, rerun analysis
random2_CBT <- rerun(random2_CBT)
summary(random2_CBT)

#IPT random effects model, define effects to test, request likelihood-based confidence intervals
#   indirect effect, product of path a, Txm3, and path b, m3o3
#   total effect, path c, Txo3
#   direct effect, total effect minus mediation effect
random2_IPT <- tssem2(random1_IPT, Amatrix=A, Smatrix=S,
                  intervals.type="LB", model.name="TSSEM2 IPT free",
                  mx.algebras=list( indirect_IPT=mxAlgebra(Txm3*m3o3, name="indirect_IPT"),
                                    direct_IPT=mxAlgebra(Txo3, name="direct_IPT"), 
                                    total_IPT=mxAlgebra(Txo3 + Txm3*m3o3, name="total_IPT") ))
summary(random2_IPT)

#Check OpenMx status1 and error messages, rerun analysis if needed
#   OpenMx status1 = 0 or 1 means the estimation is considered fine
#   Other values indicate estimation problems, rerun analysis
random2_IPT <- rerun(random2_IPT)
summary(random2_IPT)


#Moderation Analysis: Testing equality of paths a & b across CBT & IPT models

#Testing the equality of regression coefficients

#Specify a model with the same regression paths to be applied to both CBT & IPT groups
#   We only want to test whether the 2 direct paths, a and b, are the same across groups 
#   For the CBT group, we will fit the model specified in the original A matrix
#   For the IPT group, we will fit a model where paths a and b are constrained to be the same as for the CBT group
#      but allow the other paths to be freely estimated
#   Thus we create a modified A matrix using the same labels for paths a and b (Txm3, m3o3) as those in the original A matrix
#      but different labels for the other paths
A_IPT2 <- create.mxMatrix(c(0, 0,  0, 0, 0,
                            0, 0, 0, 0, 0, 
                            0, 0, 0, 0, 0, 
                            "0.1*Txm3", "0.1*m1m3_IPT2", 0, 0, 0, 
                            "0.3*Txo3_IPT2", 0, "0.1*o1o3_IPT2", "0.1*m3o3", 0),
                          type="Full", byrow=TRUE, ncol=5, nrow=5, as.mxMatrix=FALSE)
#   Label matrix for inspecting the model.
dimnames(A_IPT2)[[1]] <- dimnames(A_IPT2)[[2]] <- c("Tx", "m1", "o1", "m3", "o3") 
A_IPT2

#Allow variances and covariances to differ across the CBT and IPT groups
#   For the CBT group, we will use the original S matrix
#   For the IPT group, we will use a modified S matrix with different labels 
#        for all variances and covariances to be estimated
S_IPT <- create.mxMatrix(c(1,
                           0, 1,
                           0, "0.1*r_m1o1_IPT", 1,
                           0, 0, 0, "0.1*var_m3_IPT",
                           0, 0, 0, 0, "0.1*var_o3_IPT"), 
                         byrow=TRUE, type="Symm", as.mxMatrix=FALSE)
#   Label matrix variables, useful for inspecting the model.
dimnames(S_IPT)[[1]] <- dimnames(S_IPT)[[2]] <- c("Tx", "m1", "o1", "m3", "o3")
S_IPT


#Second stage analysis with equality constraints for 2 direct paths a and b
#   Create the models for the two groups, set the argument run=FALSE
random2_CBTc <- tssem2(random1_CBT, Amatrix=A, Smatrix=S,
                      intervals.type="LB", model.name="TSSEM2 CBT constrained", run=FALSE)
random2_IPTc <- tssem2(random1_IPT, Amatrix=A_IPT2, Smatrix=S_IPT,
                      intervals.type="LB", model.name="TSSEM2 IPT constrained", run=FALSE)

#If error message comes up, rerun analysis
random2_CBTc <- rerun(random2_CBTc)
random2_IPTc <- rerun(random2_IPTc)

#Create the multigroup model for the two groups
random2_constrained <- mxModel(model="same_regression_coef", random2_CBTc, random2_IPTc,
                              mxFitFunctionMultigroup(c("TSSEM2 CBT constrained", "TSSEM2 IPT constrained")))

#Fit multigroup model with equality constraints
random2_constrained.fit <- mxRun(random2_constrained, intervals=TRUE)

#If error message appears stating that one of the parameters has been assigned mutltiple starting values
#   Follow their instructions for randomly selecting of these values (below)
random2_constrained <- omxAssignFirstParameters(random2_constrained)

#Then try again to fit multigroup model with equality constraints
random2_constrained.fit <- mxRun(random2_constrained, intervals=TRUE)

#First make a list of the fitted models in the separate groups
submodels.fit <- list(random2_CBT, random2_IPT)
subgroup.summary(submodels.fit,random2_constrained.fit)

#Clear R environment before rerunning analyses with new data files
rm(list = ls())