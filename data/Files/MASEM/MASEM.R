
#MASEM Models for Testing Mediation (Study 1)
#   This syntax can be used to test whether any candidate mediator 
#   mediates treatment outcome for CBT or IPT
#   in the full dataset or subset of independent/complete matrices.
#   Change file names and/or the working directory where the relevant data files are saved.

#Candidate mediators tested with CBT, full dataset or subset, as reported in Study 1.
#   MA_CBTfull.txt = negative cognition full dataset of CBT matrices (k=28)
#   MA_CBTidpcom.txt = negative cognition subset of independent, complete CBT matrices (k=15)
#   ME_CBTidp.txt = social engagement subset of independent CBT matrices (k=14)
#   ME_CBTidpcom.txt = social engagement subset of independent, complete CBT matrices (k=8)
#   MB_CBTfull.txt = pleasant activities full dataset of CBT matrices (k=7)
#   MB_CBTidpcom.txt = pleasant activities set of independent, complete CBT matrices (k=3)

#For corresponding study sample size data files, append "_n" to the end of the text file.

#Load metaSEM and semPlot
library(metaSEM)
library(semPlot)

#Set working directory(library)--enter path to directory in parentheses
setwd ("C:/Users/meiying/OneDrive - Florida International University/Research/MY MS/2 Review, Revise, Resubmit/Mediation Meta Youth Dep CMs/PB/RESUBMIT/OSF_depmedmeta/Data")

#Read in dataset of correlation matrices 
vector <- readStackVec ("MB_CBTidpcom.txt") #change file for different candidate mediator
vector

#Give column and row headings, enter columns first, check vector
#   Tx=treatment, m1=pre-treatment candidate mediator, o1=pre-treatment outcome
#   m3=post-treatment candidate mediator, o3=post-treatment outcome
  vector <- lapply(vector,function(x)
  {dimnames(x) <- list(c("Tx", "m1", "o1", "m3", "o3"),
                       c("Tx", "m1", "o1", "m3", "o3"))
   x})
  vector 

#Read in corresponding study sample sizes, n
n <- read.table("MB_CBTidpcom_n.txt") #change file for different candidate mediator
n

#Display sample sizes in the pooled correlation matrix
#   number of treatment-control comparisons k for each correlation
pattern.na(vector, show.na=FALSE)
#   number of participants n for each correlation
pattern.n(vector, n)

#Test if matrix is positive definite
#   True = positive definite, False = not positive definite
#   NA = missing element in matrices
is.pd(vector)


#TSSEM Stage 1

#Random effects model, use diagonal matrix for random effects
random1 <- tssem1(vector, n, method="REM", RE.type="Diag")
options("scipen"=100, "digits"=20)
summary(random1)

#Check OpenMx status1 and rerun analysis if needed
#   OpenMx status1 = 0 or 1 means the estimation is considered fine
#   Other values indicate estimation problems, rerun analysis
random1 <-rerun(random1)
summary(random1)

#Check OpenMx status1 and rerun analysis if needed
#   OpenMx status1 = 0 or 1 means the estimation is considered fine
#   Other values indicate estimation problems, check Tau2 for very small heterogeneity variances
#   Fix any Tau2 that is < 1 x 10e-10 to 0 in a user-defined model

#Create user-defined model--select appropriate one

#   for ME_CBTidp = social engagement subset of independent CBT matrices (k=14)
#RE <- Diag(c(0,0,"0.01*Tau2_3","0.01*Tau2_4","0.01*Tau2_5","0.01*Tau2_6","0.01*Tau2_7",0,0,"0.01*Tau2_10"))

#   for MB_CBTfull = pleasant activities full dataset of CBT matrices (k=7)
#RE <- Diag(c(0,0,"0.01*Tau2_3","0.01*Tau2_4", 0, "0.01*Tau2_6", 0,0,"0.01*Tau2_9",0))

#Show user-defined model and run it
#RE
#random1 <- tssem1(vector, n, method="REM", RE.type="User",
#                  RE.constraints = RE)
#summary(random1)
#Check that OpenMx status1 = 0.

#To extract estimated average correlation matrix in matrix form
#   Select the fixed effects and convert it into a correlation matrix
vec2symMat(coef(random1, select="fixed"),diag=FALSE)


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
#   Label matrix for inspecting the model.
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
#   Label matrix variables, useful for inspecting the model.
dimnames(S)[[1]] <- dimnames(S)[[2]] <- c("Tx", "m1", "o1", "m3", "o3")
S

#Random effects model, define effects to test, request likelihood-based confidence intervals
#   indirect effect, product of path a, Txm3, and path b, m3o3
#   total effect, path c, Txo3
#   direct effect, total effect minus mediation effect
random2 <- tssem2(random1, Amatrix=A, Smatrix=S, intervals.type="LB", model.name="TSSEM2 test",
                  mx.algebras=list( indirect=mxAlgebra(Txm3*m3o3, name="indirect"),
                                    direct=mxAlgebra(Txo3, name="direct"), 
                                    total=mxAlgebra(Txo3 + Txm3*m3o3, name="total") ))
summary(random2)

#Check OpenMx status1 and error messages, rerun analysis if needed
#   OpenMx status1 = 0 or 1 means the estimation is considered fine
#   Other values indicate estimation problems, rerun analysis
random2 <- rerun(random2)
summary(random2)

#Plot the model for checking
#   Convert the model to semPlotModel object
#   Add labels, use layout="spring" to avoid overlapping labels
#   Check parameter estimates in the figure
my.plot <- meta2semPlot(random2, manNames=c("Tx", "m1", "o1", "m3", "o3"))
semPaths(my.plot, whatLabels="path", nCharEdges=10, nCharNodes=10,
         layout="spring", color="yellow", edge.label.cex=0.8)
semPaths(my.plot, whatLabels="est", nCharNodes=10, layout="spring",
         color="green", edge.label.cex=1.2)


#Clear R environment before rerunning analyses with new data files
rm(list = ls())