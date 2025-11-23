## Biometrical CLPM - T21 + COVID1-4 + T26 (6 waves) - Full WLS estimator
## Depression Sx + Cigs per Day
## Madhur Singh

rm(list=ls(all=TRUE))
#
library(OpenMx)
# mxOption(NULL, "Default optimizer","CSOLNP")
#

##### Read in the data =================

dataDir <- "Data/"
outDir  <- "Out/"
plotDir <- "Plots/"
logDir  <- "Logs/"

## WLS sum stats from 3_Get_DepSx_CigDay_SumStats_for_WLS.R

load(paste0("Out/MZsumStatDepSmk.RData"))
str(MZsumStatDepSmk)
round(cov2cor(MZsumStatDepSmk@observedStats$cov), 2)

load(paste0("Out/DZsumStatDepSmk.RData"))
str(DZsumStatDepSmk)
round(cov2cor(DZsumStatDepSmk@observedStats$cov), 2)


# Number of time points
Ts <- 6

##### N. vars =========================

## Two phenotypes: 
### DepSxRN (rank-normalized residuals of depression symptom score)
### CigDay3L (ordinal 3-level variable of current cigarettes per day)

nv=2*Ts;  nv
ntv=nv*2; ntv
labs <- paste0(rep(c("DepSxResRN_T","CigDay3L_T"),Ts),
               rep(1:Ts, each=2))
labs
selVars <- c(paste0(labs,"_Tw1"),paste0(labs,"_Tw2"))
selVars

## Ordinal Y variable with 3 levels = 2 Thresholds
varso <- paste0("CigDay3L_T",1:Ts)
varso

ordVars <- selVars[grep("CigDay3L_T",selVars)]
ordVars

nvo  <- length(varso) # N ordinal variables = (Y1, Y2, Y3)
nvo
ntvo <- nvo*2 # total ordinal variables = 2 (Y1, Y2, Y3)
ntvo

nth   <- 2 # number of thresholds

##### start values ====================

svb0 <- rep(0,nv) # intercepts, means  

## parameters SA
vax=(.4/1)    # innov variance A of DepSxResRN_T
vay=(.4/1)    # innov variance A of DrkWkResRN_T
raxy=.2       # correlation between A inno DepSxResRN_T and A innno DrkWkResRN_T
covaxy=sqrt(vax*vay)*raxy
covaxy

# next SE .... same structure
vex=(.4/1)
vey=(.4/1)
rexy=.2      
covexy=sqrt(vex*vey)*rexy
covexy

# next SC .... same structure
vcx=(.2/1)
vcy=(.2/1)
rcxy=.2  # 
covcxy=sqrt(vcx*vcy)*rcxy
covcxy

# phenotypic causal paths
bphxx=.6   #  AR of DepSxResRN_T
bphyy=.6   #  AR of CigDay3L_T
bphxy=.1   #  CigDay3L_T --> DepSxResRN_T causal
bphyx=.1   #  DepSxResRN_T --> CigDay3L_T causal

# pheno on prs parameters (DepSxResRN_T and CigDay3L_T)
# bx=sqrt(.02)   # prsx -> x
# by=sqrt(.02)   # prsy -> y
# vprsx=1
# vprsy=1
# rprsxy=.2 


#### Means matrix ===============

# Create regression model for expected Mean Matrices
B0_  <- mxMatrix( 
  type="Full", nrow=1, ncol=nv, free=TRUE, values=svb0, 
  labels = paste0("b0_",labs),
  name="b0", 
  dimnames = list("b0",labs)
) 
B0_
# Equal means by Twin order
ExpMean <- mxAlgebra(expression=cbind(b0,b0), name='expMean')



#### Thresholds matrix ===============

thLabs <- paste(paste0("th",1:nth),
                rep(varso, each=nth),
                sep = "_")

threG <- mxMatrix( type="Full", nrow=nth, ncol=ntvo, free=F, 
                   values=c(0.5,1), 
                   labels=rep(thLabs,2),
                   dimnames = list(paste0("th",1:nth),
                                   ordVars),
                   name="threG")
threG



#### Phenotypic Beta matrix ===============

## Beta Values
B_vals <- matrix(0, nv, nv)
rownames(B_vals) = colnames(B_vals) = labs

for (i in 1:Ts) {
  ## Lagged first-order paths
  if(i < Ts) {
    ## AR1_x
    B_vals[paste0("DepSxResRN_T",i+1), paste0("DepSxResRN_T",i)] <- bphxx
    ## AR1_y
    B_vals[paste0("CigDay3L_T",i+1), paste0("CigDay3L_T",i)] <- bphyy
    ## DepSxResRN_T_i --> CigDay3L_T_i+1
    B_vals[paste0("CigDay3L_T",i+1), paste0("DepSxResRN_T",i)] <- bphyx
    ## CigDay3L_T_i --> DepSxResRN_T_i+1
    B_vals[paste0("DepSxResRN_T",i+1), paste0("CigDay3L_T",i)] <- bphxy
  }
  ## Cross-sectional causal paths (for later use with DOC)
  ## DepSxResRN_T_i --> CigDay3L_T_i
  B_vals[paste0("CigDay3L_T",i), paste0("DepSxResRN_T",i)] <- 0
  ## CigDay3L_T_i --> DepSxResRN_T_i
  B_vals[paste0("DepSxResRN_T",i), paste0("CigDay3L_T",i)] <- 0
  
} ### END Vals Loop

B_vals


## Beta Free Parameters
B_free <- matrix(FALSE, nv, nv)
rownames(B_free) = colnames(B_free) = labs

for (i in 1:Ts) {
  ## Lagged first-order paths
  if(i < Ts) {
    ## AR1_x
    B_free[paste0("DepSxResRN_T",i+1), paste0("DepSxResRN_T",i)] <- TRUE
    ## AR1_y
    B_free[paste0("CigDay3L_T",i+1), paste0("CigDay3L_T",i)] <- TRUE
    ## DepSxResRN_T_i --> CigDay3L_T_i+1
    B_free[paste0("CigDay3L_T",i+1), paste0("DepSxResRN_T",i)] <- TRUE
    ## CigDay3L_T_i --> DepSxResRN_T_i+1
    B_free[paste0("DepSxResRN_T",i+1), paste0("CigDay3L_T",i)] <- TRUE
  }
  ## Cross-sectional causal paths (for later use with DOC)
  ## DepSxResRN_T_i --> CigDay3L_T_i
  B_free[paste0("CigDay3L_T",i), paste0("DepSxResRN_T",i)] <- FALSE
  ## CigDay3L_T_i --> DepSxResRN_T_i
  B_free[paste0("DepSxResRN_T",i), paste0("CigDay3L_T",i)] <- FALSE
  
} ### END Free Loop
B_free


## Beta Labels
B_labs <- matrix(NA_real_, nv, nv)
rownames(B_labs) = colnames(B_labs) = labs

for (i in 1:Ts) {
  ## Lagged first-order paths
  if(i < Ts) {
    ## AR1_x
    B_labs[paste0("DepSxResRN_T",i+1), paste0("DepSxResRN_T",i)] <- paste0("AR_bx",i+1,"x",i)
    ## AR1_y
    B_labs[paste0("CigDay3L_T",i+1), paste0("CigDay3L_T",i)] <- paste0("AR_by",i+1,"y",i)
    ## DepSxResRN_T_i --> CigDay3L_T_i+1
    B_labs[paste0("CigDay3L_T",i+1), paste0("DepSxResRN_T",i)] <- paste0("CL_by",i+1,"x",i)
    ## CigDay3L_T_i --> DepSxResRN_T_i+1
    B_labs[paste0("DepSxResRN_T",i+1), paste0("CigDay3L_T",i)] <- paste0("CL_bx",i+1,"y",i)
  }
  ## Cross-sectional causal paths
  ## DepSxResRN_T_i --> CigDay3L_T_i
  B_labs[paste0("CigDay3L_T",i), paste0("DepSxResRN_T",i)] <- paste0("Prox_by",i,"x",i)
  ## CigDay3L_T_i --> DepSxResRN_T_i
  B_labs[paste0("DepSxResRN_T",i), paste0("CigDay3L_T",i)] <- paste0("Prox_bx",i,"y",i)
  
} ### END Labs Loop

B_labs


## Combine
BE <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                free = B_free,
                values = B_vals, 
                label = B_labs,
                dimnames = list(labs,labs), 
                name="BE")
BE

# Identity matrix for algebra
I4 = mxMatrix(type="Iden", nrow=nv, ncol=nv, 
              name='I4', dimnames = list(labs,labs))


# Create Matrices for Variance Components

##### covA ==========================

# Values
covA_vals <- matrix(0, nv, nv)
rownames(covA_vals) = colnames(covA_vals) = paste0("A",labs)
covA_vals

for(i in 1:Ts) {
  covA_vals[paste0("ADepSxResRN_T",i),paste0("ADepSxResRN_T",i)] <- vax
  covA_vals[paste0("ACigDay3L_T",i),paste0("ADepSxResRN_T",i)] <- covaxy
  covA_vals[paste0("ADepSxResRN_T",i),paste0("ACigDay3L_T",i)] <- covaxy
  covA_vals[paste0("ACigDay3L_T",i),paste0("ACigDay3L_T",i)] <- vax
}
isSymmetric(covA_vals)
covA_vals


# Labels
covA_labs <- matrix(NA_real_, nv, nv)
rownames(covA_labs) = colnames(covA_labs) = paste0("A",labs)
covA_labs

for(i in 1:Ts) {
  covA_labs[paste0("ADepSxResRN_T",i),paste0("ADepSxResRN_T",i)] <- paste0("VAx",i)
  covA_labs[paste0("ACigDay3L_T",i),paste0("ADepSxResRN_T",i)] <- paste0("covAxy",i)
  covA_labs[paste0("ADepSxResRN_T",i),paste0("ACigDay3L_T",i)] <- paste0("covAxy",i)
  covA_labs[paste0("ACigDay3L_T",i),paste0("ACigDay3L_T",i)] <- paste0("VAy",i)
}
isSymmetric(covA_labs)
covA_labs


# Free
covA_free <- matrix(FALSE, nv, nv)
rownames(covA_free) = colnames(covA_free) = paste0("A",labs)
covA_free

for(i in 1:Ts) {
  covA_free[paste0("ADepSxResRN_T",i),paste0("ADepSxResRN_T",i)] <- TRUE
  covA_free[paste0("ACigDay3L_T",i),paste0("ADepSxResRN_T",i)] <- TRUE
  covA_free[paste0("ADepSxResRN_T",i),paste0("ACigDay3L_T",i)] <- TRUE
  covA_free[paste0("ACigDay3L_T",i),paste0("ACigDay3L_T",i)] <- TRUE
}
isSymmetric(covA_free)
covA_free


## Combine
covPSA <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                    free = covA_free,
                    values = covA_vals, 
                    label = covA_labs,
                    dimnames = list(paste0("A",labs),paste0("A",labs)),
                    name="PSA")
covPSA


##### Beta_A ==========================


## Beta Labels
BA_labs <- matrix(NA_real_, nv, nv)
rownames(BA_labs) = colnames(BA_labs) =  paste0("A",labs)

for (i in 1:(Ts-1)) {
  ## Lagged first-order paths
  ## AR1_x
  BA_labs[paste0("ADepSxResRN_T",i+1), paste0("ADepSxResRN_T",i)] <- paste0("AR_Ax",i+1,"x",i)
  ## AR1_y
  BA_labs[paste0("ACigDay3L_T",i+1), paste0("ACigDay3L_T",i)] <- paste0("AR_Ay",i+1,"y",i)
  ## DepSxResRN_T_i --> CigDay3L_T_i+1
  BA_labs[paste0("ACigDay3L_T",i+1), paste0("ADepSxResRN_T",i)] <- paste0("CL_Ay",i+1,"x",i)
  ## CigDay3L_T_i --> DepSxResRN_T_i+1
  BA_labs[paste0("ADepSxResRN_T",i+1), paste0("ACigDay3L_T",i)] <- paste0("CL_Ax",i+1,"y",i)
} ### END Labs Loop

BA_labs

##combine
BEA <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                 free=F, values=0, # Initial model is phenotypic AR paths only
                 label=BA_labs, 
                 dimnames = list(paste0("A",labs),paste0("A",labs)),
                 name="BEA" )


##### covE ==========================

# Values
covE_vals <- matrix(0, nv, nv)
rownames(covE_vals) = colnames(covE_vals) = paste0("E",labs)
covE_vals

for(i in 1:Ts) {
  covE_vals[paste0("EDepSxResRN_T",i),paste0("EDepSxResRN_T",i)] <- vex
  covE_vals[paste0("ECigDay3L_T",i),paste0("EDepSxResRN_T",i)] <- covexy
  covE_vals[paste0("EDepSxResRN_T",i),paste0("ECigDay3L_T",i)] <- covexy
  covE_vals[paste0("ECigDay3L_T",i),paste0("ECigDay3L_T",i)] <- vex
}
isSymmetric(covE_vals)
covE_vals


# Labels
covE_labs <- matrix(NA_real_, nv, nv)
rownames(covE_labs) = colnames(covE_labs) = paste0("E",labs)
covE_labs

for(i in 1:Ts) {
  covE_labs[paste0("EDepSxResRN_T",i),paste0("EDepSxResRN_T",i)] <- paste0("VEx",i)
  covE_labs[paste0("ECigDay3L_T",i),paste0("EDepSxResRN_T",i)] <- paste0("covExy",i)
  covE_labs[paste0("EDepSxResRN_T",i),paste0("ECigDay3L_T",i)] <- paste0("covExy",i)
  covE_labs[paste0("ECigDay3L_T",i),paste0("ECigDay3L_T",i)] <- paste0("VEy",i)
}
isSymmetric(covE_labs)
covE_labs


# Free
covE_free <- matrix(FALSE, nv, nv)
rownames(covE_free) = colnames(covE_free) = paste0("E",labs)
covE_free

for(i in 1:Ts) {
  covE_free[paste0("EDepSxResRN_T",i),paste0("EDepSxResRN_T",i)] <- TRUE
  covE_free[paste0("ECigDay3L_T",i),paste0("EDepSxResRN_T",i)] <- TRUE
  covE_free[paste0("EDepSxResRN_T",i),paste0("ECigDay3L_T",i)] <- TRUE
  covE_free[paste0("ECigDay3L_T",i),paste0("ECigDay3L_T",i)] <- TRUE
}
isSymmetric(covE_free)
covE_free


## Combine
covPSE <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                    free = covE_free,
                    values = covE_vals, 
                    label = covE_labs,
                    dimnames = list(paste0("E",labs),paste0("E",labs)),
                    name="PSE")
covPSE


##### Beta_E ==========================


## Beta Labels
BE_labs <- matrix(NA_real_, nv, nv)
rownames(BE_labs) = colnames(BE_labs) =  paste0("E",labs)

for (i in 1:(Ts-1)) {
  ## Lagged first-order paths
  ## AR1_x
  BE_labs[paste0("EDepSxResRN_T",i+1), paste0("EDepSxResRN_T",i)] <- paste0("AR_Ex",i+1,"x",i)
  ## AR1_y
  BE_labs[paste0("ECigDay3L_T",i+1), paste0("ECigDay3L_T",i)] <- paste0("AR_Ey",i+1,"y",i)
  ## DepSxResRN_T_i --> CigDay3L_T_i+1
  BE_labs[paste0("ECigDay3L_T",i+1), paste0("EDepSxResRN_T",i)] <- paste0("CL_Ey",i+1,"x",i)
  ## CigDay3L_T_i --> DepSxResRN_T_i+1
  BE_labs[paste0("EDepSxResRN_T",i+1), paste0("ECigDay3L_T",i)] <- paste0("CL_Ex",i+1,"y",i)
} ### END Labs Loop

BE_labs


##combine
BEE <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                 free=F, values=0, # Initial model is phenotypic AR paths only
                 label=BE_labs, 
                 dimnames = list(paste0("E",labs),paste0("E",labs)),
                 name="BEE" )


##### covC ==========================

# Values
covC_vals <- matrix(0, nv, nv)
rownames(covC_vals) = colnames(covC_vals) = paste0("C",labs)
covC_vals

for(i in 1:Ts) {
  covC_vals[paste0("CDepSxResRN_T",i),paste0("CDepSxResRN_T",i)] <- vcx
  covC_vals[paste0("CCigDay3L_T",i),paste0("CDepSxResRN_T",i)] <- covcxy
  covC_vals[paste0("CDepSxResRN_T",i),paste0("CCigDay3L_T",i)] <- covcxy
  covC_vals[paste0("CCigDay3L_T",i),paste0("CCigDay3L_T",i)] <- vcx
}
isSymmetric(covC_vals)
covC_vals


# Labels
covC_labs <- matrix(NA_real_, nv, nv)
rownames(covC_labs) = colnames(covC_labs) = paste0("C",labs)
covC_labs

for(i in 1:Ts) {
  covC_labs[paste0("CDepSxResRN_T",i),paste0("CDepSxResRN_T",i)] <- paste0("VCx",i)
  covC_labs[paste0("CCigDay3L_T",i),paste0("CDepSxResRN_T",i)] <- paste0("covCxy",i)
  covC_labs[paste0("CDepSxResRN_T",i),paste0("CCigDay3L_T",i)] <- paste0("covCxy",i)
  covC_labs[paste0("CCigDay3L_T",i),paste0("CCigDay3L_T",i)] <- paste0("VCy",i)
}
isSymmetric(covC_labs)
covC_labs


# Free
covC_free <- matrix(FALSE, nv, nv)
rownames(covC_free) = colnames(covC_free) = paste0("C",labs)
covC_free

for(i in 1:Ts) {
  covC_free[paste0("CDepSxResRN_T",i),paste0("CDepSxResRN_T",i)] <- TRUE
  covC_free[paste0("CCigDay3L_T",i),paste0("CDepSxResRN_T",i)] <- TRUE
  covC_free[paste0("CDepSxResRN_T",i),paste0("CCigDay3L_T",i)] <- TRUE
  covC_free[paste0("CCigDay3L_T",i),paste0("CCigDay3L_T",i)] <- TRUE
}
isSymmetric(covC_free)
covC_free


## Combine
covPSC <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                    free = covC_free,
                    values = covC_vals, 
                    label = covC_labs,
                    dimnames = list(paste0("C",labs),paste0("C",labs)),
                    name="PSC")
covPSC


##### Beta_C ==========================


## Beta Labels
BC_labs <- matrix(NA_real_, nv, nv)
rownames(BC_labs) = colnames(BC_labs) =  paste0("C",labs)

for (i in 1:(Ts-1)) {
  ## Lagged first-order paths
  ## AR1_x
  BC_labs[paste0("CDepSxResRN_T",i+1), paste0("CDepSxResRN_T",i)] <- paste0("AR_Cx",i+1,"x",i)
  ## AR1_y
  BC_labs[paste0("CCigDay3L_T",i+1), paste0("CCigDay3L_T",i)] <- paste0("AR_Cy",i+1,"y",i)
  ## DepSxResRN_T_i --> CigDay3L_T_i+1
  BC_labs[paste0("CCigDay3L_T",i+1), paste0("CDepSxResRN_T",i)] <- paste0("CL_Cy",i+1,"x",i)
  ## CigDay3L_T_i --> DepSxResRN_T_i+1
  BC_labs[paste0("CDepSxResRN_T",i+1), paste0("CCigDay3L_T",i)] <- paste0("CL_Cx",i+1,"y",i)
} ### END Labs Loop

BC_labs


##combine
BEC <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                 free=F, values=0, # Initial model is phenotypic AR paths only
                 label=BC_labs, 
                 dimnames = list(paste0("C",labs),paste0("C",labs)),
                 name="BEC" )




#### ACE model ==================

### (I-B)^-1
BAi = mxAlgebra(expression=solve(I4-BEA), name='BAi')
BCi = mxAlgebra(expression=solve(I4-BEC), name='BCi')
BEi = mxAlgebra(expression=solve(I4-BEE), name='BEi')

### Expected VA, VC, VE matrices
covA = mxAlgebra(expression= (BAi)%*%PSA%*%t(BAi), name='VA', 
                 dimnames = list(paste0("A",labs), paste0("A",labs)) )
covC = mxAlgebra(expression= (BCi)%*%PSC%*%t(BCi), name='VC', 
                 dimnames = list(paste0("C",labs), paste0("C",labs)) )
covE = mxAlgebra(expression= (BEi)%*%PSE%*%t(BEi), name='VE', 
                 dimnames = list(paste0("E",labs), paste0("E",labs)) )


# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP1      <- mxAlgebra( expression= VA+VC+VE, name="V1", 
                         dimnames = list(labs, labs) )
iBE        <- mxAlgebra(expression = solve(I4-BE), name='iBE') # iBE is not BEi
#
covP       <- mxAlgebra(expression= iBE%*%V1%*%t(iBE), name='V', 
                        dimnames = list(labs, labs) )
covMZ      <- mxAlgebra( expression= iBE%*%(VA+VC)%*%t(iBE), name="cMZ" )
covDZ      <- mxAlgebra( expression= iBE%*%(0.5%x%VA+VC)%*%t(iBE), name="cDZ" )
#
expCovMZ   <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ   <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )


# Create Data Objects for Multiple Groups
dataMZ = mxData( type = "none", observedStats = MZsumStatDepSmk$observedStats, numObs = MZsumStatDepSmk$numObs)
dataDZ = mxData( type = "none", observedStats = DZsumStatDepSmk$observedStats, numObs = DZsumStatDepSmk$numObs)

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars, 
                                  thresholds="threG", threshnames=ordVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars, 
                                  thresholds="threG", threshnames=ordVars )
funWLS    <- mxFitFunctionWLS(type = "WLS", allContinuousMethod = "marginals")

# combine all parameters
pars      <- list(B0_, threG,
                  covPSA, BEA, BAi, covA, 
                  covPSC, BEC, BCi, covC, 
                  covPSE, BEE, BEi, covE, 
                  covP1, covP, I4, BE, iBE )
modelMZ   <- mxModel( pars, ExpMean, covMZ, expCovMZ, ExpMean, dataMZ, expMZ, funWLS, name="MZ" )
modelDZ   <- mxModel( pars, ExpMean, covDZ, expCovDZ, ExpMean, dataDZ, expDZ, funWLS, name="DZ" )	
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )


## Standardized estimates ======================


# Calculate genetic and environmental correlations 
corA      <- mxAlgebra( expression = vec2diag(1/sqrt(diag2vec(VA))) %&% VA, name ="rA", 
                        dimnames = list(paste0("A",labs), paste0("A",labs)) ) 
corC      <- mxAlgebra( expression = vec2diag(1/sqrt(diag2vec(VC))) %&% VC, name ="rC", 
                        dimnames = list(paste0("C",labs), paste0("C",labs)) )
corE      <- mxAlgebra( expression = vec2diag(1/sqrt(diag2vec(VE)))  %&% VE, name ="rE", 
                        dimnames = list(paste0("E",labs), paste0("E",labs)) )

## Standardized covariance paths = gives the same cross-sectional estimates as rA, rC, rE
psiA      <- mxAlgebra( expression = vec2diag(1/sqrt(diag2vec(PSA))) %&% PSA, name ="psiA",
                        dimnames = list(paste0("A",labs), paste0("A",labs)) )
psiC      <- mxAlgebra( expression = vec2diag(1/sqrt(diag2vec(PSC))) %&% PSC, name ="psiC",
                        dimnames = list(paste0("C",labs), paste0("C",labs)) )
psiE      <- mxAlgebra( expression = vec2diag(1/sqrt(diag2vec(PSE))) %&% PSE, name ="psiE",
                        dimnames = list(paste0("E",labs), paste0("E",labs)) )

# Standardized regression paths
stdBeta   <- mxAlgebra( expression= vec2diag(1/sqrt(diag2vec(V))) %*% BE %*% vec2diag(sqrt(diag2vec(V))), name ="stdBeta", 
                        dimnames = list(labs, labs) ) 
stdBetaA  <- mxAlgebra( expression= vec2diag(1/sqrt(diag2vec(VA))) %*% BEA %*% vec2diag(sqrt(diag2vec(VA))), name ="stdBetaA", 
                        dimnames = list(paste0("A",labs), paste0("A",labs)) )
stdBetaC  <- mxAlgebra( expression= vec2diag(1/sqrt(diag2vec(VC))) %*% BEC %*% vec2diag(sqrt(diag2vec(VC))), name ="stdBetaC", 
                        dimnames = list(paste0("C",labs), paste0("C",labs)) )
stdBetaE  <- mxAlgebra( expression= vec2diag(1/sqrt(diag2vec(VE))) %*% BEE %*% vec2diag(sqrt(diag2vec(VE))), name ="stdBetaE", 
                        dimnames = list(paste0("E",labs), paste0("E",labs)) )

## Create Algebra for Unstandardized and Standardized Variance Components
rowUS     <- paste0(rep('US_',nv), labs)
colUS     <- paste0(rep(c('VA_','VC_','VE_','SA_','SC_','SE_'),each=nv), labs)
estUS     <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )

# Combine
calc      <- list( 
  corA, corC, corE,
  psiA, psiC, psiE,
  stdBeta, stdBetaA, stdBetaC, stdBetaE,
  estUS
)


#### RUN THE MODEL =========================
modelcausalA  <- mxModel(name="CLPM_phenoAR", pars, modelMZ, modelDZ, multi, calc )

## Drop all C parameters
modelcausalA <- omxSetParameters(modelcausalA, 
                                 labels = c(paste0("VCx",1:Ts),
                                            paste0("covCxy",1:Ts),
                                            paste0("VCy",1:Ts)), free = F, values = 0)

fitcausalA    <- mxRun( modelcausalA) 
sumcausalA    <- summary( fitcausalA )
sumcausalA

round(mxEval(stdBeta, fitcausalA), 2)
round(  mxSE(stdBeta, fitcausalA), 3)

# Mean Smk --> Dep
SmkDep <- mxAlgebra(expression = cbind(stdBeta[3,2],stdBeta[5,4],stdBeta[7,6],stdBeta[9,8],stdBeta[11,10]), 
                    name ="SmkDep")
mSmkDep   <- mxAlgebra( expression = mean(SmkDep), 
                        name ="mSmkDep" )

# Mean Dep --> Smk
DepSmk   <- mxAlgebra( expression = cbind(stdBeta[4,1],stdBeta[6,3],stdBeta[8,5],stdBeta[10,7],stdBeta[12,9]),
                       name ="DepSmk" )
mDepSmk   <- mxAlgebra( expression = mean(DepSmk),
                        name ="mDepSmk" )

modelcausalA0 <- mxModel(name="CLPM_phenoAR_meanBeta", modelcausalA, 
                         SmkDep, mSmkDep, DepSmk, mDepSmk )

fitcausalA0 <- mxRun( modelcausalA0 )
sumcausalA0 <- summary( fitcausalA0 )
sumcausalA0

round(mxEval(mSmkDep, fitcausalA0), 2) 
round(  mxSE(mSmkDep, fitcausalA0), 3) 

round(mxEval(mDepSmk, fitcausalA0), 2) 
round(  mxSE(mDepSmk, fitcausalA0), 3) 



#### SubModel 1: AR at AE ======================

modelcausalA1 <- mxModel(modelcausalA, name = "CLPM_aeAR")
modelcausalA1 <- omxSetParameters(modelcausalA1, 
                                  labels = c("AR_bx2x1","AR_by2y1",
                                             "AR_bx3x2","AR_by3y2",
                                             "AR_bx4x3","AR_by4y3",
                                             "AR_bx5x4","AR_by5y4",
                                             "AR_bx6x5","AR_by6y5"), free = F, values = 0)
modelcausalA1 <- omxSetParameters(modelcausalA1, 
                                  labels = c("AR_Ax2x1","AR_Ay2y1",
                                             "AR_Ax3x2","AR_Ay3y2",
                                             "AR_Ax4x3","AR_Ay4y3",
                                             "AR_Ax5x4","AR_Ay5y4",
                                             "AR_Ax6x5","AR_Ay6y5"), free = T, values = 0.8)
modelcausalA1 <- omxSetParameters(modelcausalA1, 
                                  labels = c("AR_Ex2x1","AR_Ey2y1",
                                             "AR_Ex3x2","AR_Ey3y2",
                                             "AR_Ex4x3","AR_Ey4y3",
                                             "AR_Ex5x4","AR_Ey5y4",
                                             "AR_Ex6x5","AR_Ey6y5"), free = T, values = 0.3)

fitcausalA1    <- mxRun( modelcausalA1)
sumcausalA1    <- summary( fitcausalA1 )

mxCompare(fitcausalA1, fitcausalA)

sumcausalA1



#### SubModel 2: Drop Phenotypic Causation ======================

modelcausalA2 <- mxModel(fitcausalA1, name = "CLPM_aeAR_noCL")
modelcausalA2 <- omxSetParameters(modelcausalA2,
                                  labels = c("CL_by2x1","CL_bx2y1",
                                             "CL_by3x2","CL_bx3y2",
                                             "CL_by4x3","CL_bx4y3",
                                             "CL_by5x4","CL_bx5y4",
                                             "CL_by6x5","CL_bx6y5"), free = F, values = 0)

fitcausalA2   <- mxRun( modelcausalA2 )

sumcausalA2    <- summary( fitcausalA2 )

mxCompare(fitcausalA1, fitcausalA2)

sumcausalA2


#### SubModel 3: Drop Dep --> Smk Causation ======================

modelcausalA3 <- mxModel(fitcausalA1, name = "CLPM_aeAR_noCLyx")
modelcausalA3 <- omxSetParameters(modelcausalA3,
                                  labels = c("CL_by2x1",
                                             "CL_by3x2",
                                             "CL_by4x3",
                                             "CL_by5x4",
                                             "CL_by6x5"), free = F, values = 0)
fitcausalA3   <- mxRun( modelcausalA3 )

sumcausalA3    <- summary( fitcausalA3 )

mxCompare(fitcausalA1, fitcausalA3)

sumcausalA3


#### SubModel 4: Drop Smk --> Dep Causation ======================

modelcausalA4 <- mxModel(fitcausalA1, name = "CLPM_aeAR_noCLxy")
modelcausalA4 <- omxSetParameters(modelcausalA4,
                                  labels = c("CL_bx2y1",
                                             "CL_bx3y2",
                                             "CL_bx4y3",
                                             "CL_bx5y4",
                                             "CL_bx6y5"), free = F, values = 0)
fitcausalA4   <- mxRun( modelcausalA4 )

sumcausalA4    <- summary( fitcausalA4 )

mxCompare(fitcausalA1, list(fitcausalA2,fitcausalA3,fitcausalA4))

sumcausalA4


### Final Model - Get meta/mean effect size #####

betaCausalA1 <- mxEval(stdBeta, fitcausalA1)
seCausalA1 <- mxSE(stdBeta, fitcausalA1)
colnames(seCausalA1) <- colnames(betaCausalA1)
rownames(seCausalA1) <- rownames(betaCausalA1)

betaCausalA1
seCausalA1

CigDepBeta <- c(betaCausalA1["DepSxResRN_T2","CigDay3L_T1"],
                betaCausalA1["DepSxResRN_T3","CigDay3L_T2"],
                betaCausalA1["DepSxResRN_T4","CigDay3L_T3"],
                betaCausalA1["DepSxResRN_T5","CigDay3L_T4"],
                betaCausalA1["DepSxResRN_T6","CigDay3L_T5"])
CigDepBeta

CigDepSE <- c(seCausalA1["DepSxResRN_T2","CigDay3L_T1"],
              seCausalA1["DepSxResRN_T3","CigDay3L_T2"],
              seCausalA1["DepSxResRN_T4","CigDay3L_T3"],
              seCausalA1["DepSxResRN_T5","CigDay3L_T4"],
              seCausalA1["DepSxResRN_T6","CigDay3L_T5"])
CigDepSE

metaCigDep <- meta::metagen(TE = CigDepBeta, seTE = CigDepSE)
summary(metaCigDep)

DepCigBeta <- c(betaCausalA1["CigDay3L_T2","DepSxResRN_T1"],
                betaCausalA1["CigDay3L_T3","DepSxResRN_T2"],
                betaCausalA1["CigDay3L_T4","DepSxResRN_T3"],
                betaCausalA1["CigDay3L_T5","DepSxResRN_T4"],
                betaCausalA1["CigDay3L_T6","DepSxResRN_T5"])
DepCigBeta

DepCigSE <- c(seCausalA1["CigDay3L_T2","DepSxResRN_T1"],
              seCausalA1["CigDay3L_T3","DepSxResRN_T2"],
              seCausalA1["CigDay3L_T4","DepSxResRN_T3"],
              seCausalA1["CigDay3L_T5","DepSxResRN_T4"],
              seCausalA1["CigDay3L_T6","DepSxResRN_T5"])
DepCigSE

metaDepCig <- meta::metagen(TE = DepCigBeta, seTE = DepCigSE)
summary(metaDepCig)

# Mean Smk --> Dep
SmkDep <- mxAlgebra(expression = cbind(stdBeta[3,2],stdBeta[5,4],stdBeta[7,6],stdBeta[9,8],stdBeta[11,10]), 
                    name ="SmkDep")
mSmkDep   <- mxAlgebra( expression = mean(SmkDep), 
                        name ="mSmkDep" )

# Mean Dep --> Smk
DepSmk   <- mxAlgebra( expression = cbind(stdBeta[4,1],stdBeta[6,3],stdBeta[8,5],stdBeta[10,7],stdBeta[12,9]),
                       name ="DepSmk" )
mDepSmk   <- mxAlgebra( expression = mean(DepSmk),
                        name ="mDepSmk" )

modelcausalA1B <- mxModel(name="CLPM_aeAR_meanBeta", modelcausalA1, 
                          SmkDep, mSmkDep, DepSmk, mDepSmk )

fitcausalA1B <- mxRun( modelcausalA1B )
sumcausalA1B <- summary( fitcausalA1B )
sumcausalA1B

round(mxEval(mSmkDep, fitcausalA1B), 2) 
round(  mxSE(mSmkDep, fitcausalA1B), 3) 

round(mxEval(mDepSmk, fitcausalA1B), 2) 
round(  mxSE(mDepSmk, fitcausalA1B), 3) 

round(mxEval(stdBeta, fitcausalA1B), 2)
round(  mxSE(stdBeta, fitcausalA1B), 3)

round(mxEval(stdBetaA, fitcausalA1B), 2)
round(  mxSE(stdBetaA, fitcausalA1B), 3)

round(mxEval(stdBetaE, fitcausalA1B), 2)
round(  mxSE(stdBetaE, fitcausalA1B), 3)

round(Matrix::nearPD(x = fitcausalA1B$PSA$values, corr = T)$mat, 2)

round(mxEval(psiA, fitcausalA1B), 2)
round(  mxSE(psiA, fitcausalA1B), 3)

round(mxEval(psiE, fitcausalA1B), 2)
round(  mxSE(psiE, fitcausalA1B), 3)

sumcausalA1B
outCausalA1B <- sumcausalA1B$parameters
write.csv(outCausalA1B, paste0(outDir,"CLPM_Biometric_Dep_Smk_TEDS21_Cv14_26.csv"), quote = F)


# END ===========
