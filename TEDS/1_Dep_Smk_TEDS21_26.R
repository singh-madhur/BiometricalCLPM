## Depression + Cigs per Day
## CLPM - T21 + T26

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

load(paste0("Out/MZsumStatDepSmk_T16.RData"))
str(MZsumStatDepSmk_T16)
round(cov2cor(MZsumStatDepSmk_T16@observedStats$cov), 2)

load(paste0("Out/DZsumStatDepSmk_T16.RData"))
str(DZsumStatDepSmk_T16)
round(cov2cor(DZsumStatDepSmk_T16@observedStats$cov), 2)

# Number of time points
Ts <- 2

##### N. vars =========================

## Two phenotypes: 
### DepSxRN (rank-normalized residuals of depression symptom score)
### CigDay3L (ordinal 3-level variable of current cigarettes per day)

nv=2*Ts;  nv
ntv=nv*2; ntv
labs <- paste0(rep(c("DepSxResRN_T","CigDay3L_T"),Ts),
               rep(c(1,6), each=2))
labs
selVars <- c(paste0(labs,"_Tw1"),paste0(labs,"_Tw2"))
selVars

## Ordinal Y variable with 3 levels = 2 Thresholds
varso <- paste0("CigDay3L_T",c(1,6))
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
vax=(.4/1)    # innov variance A of X
vay=(.4/1)    # innov variance A of Y
raxy=.2       # correlation between A inno X and A innno Y
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
bphxx=.6   #  AR of X
bphyy=.6   #  AR of Y
bphxy=.1   #  Y --> X causal
bphyx=.1   #  X --> Y causal

# pheno on prs parameters (X and Y)
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
#
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

for (i in c(1,6)) {
  ## Lagged first-order paths
  if(i < 6) {
    ## AR1_x
    B_vals[paste0("DepSxResRN_T",i+5), paste0("DepSxResRN_T",i)] <- bphxx
    ## AR1_y
    B_vals[paste0("CigDay3L_T",i+5), paste0("CigDay3L_T",i)] <- bphyy
    ## DepSxResRN_T_i --> CigDay3L_T_i+1
    B_vals[paste0("CigDay3L_T",i+5), paste0("DepSxResRN_T",i)] <- bphyx
    ## CigDay3L_T_i --> DepSxResRN_T_i+1
    B_vals[paste0("DepSxResRN_T",i+5), paste0("CigDay3L_T",i)] <- bphxy
  }
  ## Cross-sectional causal paths (for later use with IVs)
  ## DepSxResRN_T_i --> CigDay3L_T_i
  B_vals[paste0("CigDay3L_T",i), paste0("DepSxResRN_T",i)] <- 0
  ## CigDay3L_T_i --> DepSxResRN_T_i
  B_vals[paste0("DepSxResRN_T",i), paste0("CigDay3L_T",i)] <- 0
  
} ### END Vals Loop

B_vals


## Beta Free Parameters
B_free <- matrix(FALSE, nv, nv)
rownames(B_free) = colnames(B_free) = labs

for (i in c(1,6)) {
  ## Lagged first-order paths
  if(i < 6) {
    ## AR1_x
    B_free[paste0("DepSxResRN_T",i+5), paste0("DepSxResRN_T",i)] <- TRUE
    ## AR1_y
    B_free[paste0("CigDay3L_T",i+5), paste0("CigDay3L_T",i)] <- TRUE
    ## DepSxResRN_T_i --> CigDay3L_T_i+1
    B_free[paste0("CigDay3L_T",i+5), paste0("DepSxResRN_T",i)] <- TRUE
    ## CigDay3L_T_i --> DepSxResRN_T_i+1
    B_free[paste0("DepSxResRN_T",i+5), paste0("CigDay3L_T",i)] <- TRUE
  }
  ## Cross-sectional causal paths (for later use with IVs)
  ## DepSxResRN_T_i --> CigDay3L_T_i
  B_free[paste0("CigDay3L_T",i), paste0("DepSxResRN_T",i)] <- FALSE
  ## CigDay3L_T_i --> DepSxResRN_T_i
  B_free[paste0("DepSxResRN_T",i), paste0("CigDay3L_T",i)] <- FALSE
  
} ### END Free Loop
B_free


## Beta Labels
B_labs <- matrix(NA_real_, nv, nv)
rownames(B_labs) = colnames(B_labs) = labs

for (i in c(1,6)) {
  ## Lagged first-order paths
  if(i < 6) {
    ## AR1_x
    B_labs[paste0("DepSxResRN_T",i+5), paste0("DepSxResRN_T",i)] <- paste0("AR_bx",i+5,"x",i)
    ## AR1_y
    B_labs[paste0("CigDay3L_T",i+5), paste0("CigDay3L_T",i)] <- paste0("AR_by",i+5,"y",i)
    ## DepSxResRN_T_i --> CigDay3L_T_i+1
    B_labs[paste0("CigDay3L_T",i+5), paste0("DepSxResRN_T",i)] <- paste0("CL_by",i+5,"x",i)
    ## CigDay3L_T_i --> DepSxResRN_T_i+1
    B_labs[paste0("DepSxResRN_T",i+5), paste0("CigDay3L_T",i)] <- paste0("CL_bx",i+5,"y",i)
  }
  ## Cross-sectional causal paths
  ## DepSxResRN_T_i --> CigDay3L_T_i
  B_labs[paste0("CigDay3L_T",i), paste0("DepSxResRN_T",i)] <- paste0("Inst_by",i,"x",i)
  ## CigDay3L_T_i --> DepSxResRN_T_i
  B_labs[paste0("DepSxResRN_T",i), paste0("CigDay3L_T",i)] <- paste0("Inst_bx",i,"y",i)
  
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

for(i in c(1,6)) {
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

for(i in c(1,6)) {
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

for(i in c(1,6)) {
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
  BA_labs[paste0("ADepSxResRN_T",i+5), paste0("ADepSxResRN_T",i)] <- paste0("AR_Ax",i+5,"x",i)
  ## AR1_y
  BA_labs[paste0("ACigDay3L_T",i+5), paste0("ACigDay3L_T",i)] <- paste0("AR_Ay",i+5,"y",i)
  ## DepSxResRN_T_i --> CigDay3L_T_i+1
  BA_labs[paste0("ACigDay3L_T",i+5), paste0("ADepSxResRN_T",i)] <- paste0("CL_Ay",i+5,"x",i)
  ## CigDay3L_T_i --> DepSxResRN_T_i+1
  BA_labs[paste0("ADepSxResRN_T",i+5), paste0("ACigDay3L_T",i)] <- paste0("CL_Ax",i+5,"y",i)
} ### END Labs Loop

BA_labs

##combine
BEA <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                 free=F, values=0, # Initial model is phenotypic AR paths only
                 label=BA_labs, 
                 dimnames = list(paste0("A",labs),paste0("A",labs)),
                 name="BEA" )
BEA


##### covE ==========================

# Values
covE_vals <- matrix(0, nv, nv)
rownames(covE_vals) = colnames(covE_vals) = paste0("E",labs)
covE_vals

for(i in c(1,6)) {
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

for(i in c(1,6)) {
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

for(i in c(1,6)) {
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
  BE_labs[paste0("EDepSxResRN_T",i+5), paste0("EDepSxResRN_T",i)] <- paste0("AR_Ex",i+5,"x",i)
  ## AR1_y
  BE_labs[paste0("ECigDay3L_T",i+5), paste0("ECigDay3L_T",i)] <- paste0("AR_Ey",i+5,"y",i)
  ## DepSxResRN_T_i --> CigDay3L_T_i+1
  BE_labs[paste0("ECigDay3L_T",i+5), paste0("EDepSxResRN_T",i)] <- paste0("CL_Ey",i+5,"x",i)
  ## CigDay3L_T_i --> DepSxResRN_T_i+1
  BE_labs[paste0("EDepSxResRN_T",i+5), paste0("ECigDay3L_T",i)] <- paste0("CL_Ex",i+5,"y",i)
} ### END Labs Loop

BE_labs


##combine
BEE <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                 free=F, values=0, # Initial model is phenotypic AR paths only
                 label=BE_labs, 
                 dimnames = list(paste0("E",labs),paste0("E",labs)),
                 name="BEE" )
BEE


##### covC ==========================

# Values
covC_vals <- matrix(0, nv, nv)
rownames(covC_vals) = colnames(covC_vals) = paste0("C",labs)
covC_vals

for(i in c(1,6)) {
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

for(i in c(1,6)) {
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

for(i in c(1,6)) {
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
  BC_labs[paste0("CDepSxResRN_T",i+5), paste0("CDepSxResRN_T",i)] <- paste0("AR_Cx",i+5,"x",i)
  ## AR1_y
  BC_labs[paste0("CCigDay3L_T",i+5), paste0("CCigDay3L_T",i)] <- paste0("AR_Cy",i+5,"y",i)
  ## DepSxResRN_T_i --> CigDay3L_T_i+1
  BC_labs[paste0("CCigDay3L_T",i+5), paste0("CDepSxResRN_T",i)] <- paste0("CL_Cy",i+5,"x",i)
  ## CigDay3L_T_i --> DepSxResRN_T_i+1
  BC_labs[paste0("CDepSxResRN_T",i+5), paste0("CCigDay3L_T",i)] <- paste0("CL_Cx",i+5,"y",i)
} ### END Labs Loop

BC_labs


##combine
BEC <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                 free=F, values=0, # Initial model is phenotypic AR paths only
                 label=BC_labs, 
                 dimnames = list(paste0("C",labs),paste0("C",labs)),
                 name="BEC" )

BEC



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
dataMZ = mxData( type = "none", observedStats = MZsumStatDepSmk_T16$observedStats, numObs = MZsumStatDepSmk_T16$numObs)
dataDZ = mxData( type = "none", observedStats = DZsumStatDepSmk_T16$observedStats, numObs = DZsumStatDepSmk_T16$numObs)

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
calc      <- list( corA, corC, corE,
                   psiA, psiC, psiE,
                   stdBeta, stdBetaA, stdBetaC, stdBetaE,
                   estUS
)


#### RUN MODEL =========================
modelcausalA  <- mxModel(name="CLPM_phenoAR", pars, modelMZ, modelDZ, multi, calc )

## Drop all C parameters
modelcausalA <- omxSetParameters(modelcausalA, 
                                 labels = c("VCx1","covCxy1","VCy1",
                                            "VCx6","covCxy6","VCy6"), free = F, values = 0)

fitcausalA    <- mxRun( modelcausalA ) 
sumcausalA    <- summary( fitcausalA )
sumcausalA

round(fitcausalA$V$result, 2)
round(mxSE(V, fitcausalA), 2)

round(cov2cor(fitcausalA$V$result), 3)


round(fitcausalA@algebras$rA$result, 2)
round(fitcausalA@algebras$rE$result, 2)

round(fitcausalA@algebras$stdBeta$result, 2)

outCausalA <- sumcausalA$parameters
write.csv(outCausalA, paste0(outDir,"CLPM_Std_Dep_Smk_TEDS21_26.csv"), quote = F)



#### SubModel 0: Drop Phenotypic Causation/Shared Liabilities ======================

modelcausalA0A <- mxModel(fitcausalA, name = "CLPM_phenoAR_noCL")
modelcausalA0A <- omxSetParameters(modelcausalA0A,
                                   labels = c("CL_by6x1","CL_bx6y1"), free = F, values = 0)
fitcausalA0A   <- mxRun( modelcausalA0A )

modelcausalA0B <- mxModel(fitcausalA, name = "CLPM_phenoAR_noCLyx")
modelcausalA0B <- omxSetParameters(modelcausalA0B,
                                   labels = c("CL_by6x1"), free = F, values = 0)
fitcausalA0B   <- mxRun( modelcausalA0B )

modelcausalA0C <- mxModel(fitcausalA, name = "CLPM_phenoAR_noCLxy")
modelcausalA0C <- omxSetParameters(modelcausalA0C,
                                   labels = c("CL_bx6y1"), free = F, values = 0)
fitcausalA0C   <- mxRun( modelcausalA0C )


mxCompare(fitcausalA, list(fitcausalA0A, fitcausalA0B, fitcausalA0C) )


modelcausalA0D <- mxModel(fitcausalA0B, name = "CLPM_phenoAR_noCovE1")
modelcausalA0D <- omxSetParameters(modelcausalA0D,
                                   labels = c("covExy1"), free = F, values = 0)
fitcausalA0D   <- mxRun( modelcausalA0D )

modelcausalA0E <- mxModel(fitcausalA0B, name = "CLPM_phenoAR_noCovA2")
modelcausalA0E <- omxSetParameters(modelcausalA0E,
                                   labels = c("covAxy6"), free = F, values = 0)
fitcausalA0E   <- mxRun( modelcausalA0E )

modelcausalA0F <- mxModel(fitcausalA0B, name = "CLPM_phenoAR_noCovE2")
modelcausalA0F <- omxSetParameters(modelcausalA0F,
                                   labels = c("covExy6"), free = F, values = 0)
fitcausalA0F   <- mxRun( modelcausalA0F )

modelcausalA0G <- mxModel(fitcausalA0F, name = "CLPM_phenoAR_noCov2")
modelcausalA0G <- omxSetParameters(modelcausalA0G,
                                   labels = c("covAxy6"), free = F, values = 0)
fitcausalA0G   <- mxRun( modelcausalA0G )

modelcausalA0H <- mxModel(fitcausalA0D, name = "CLPM_phenoAR_noCovE1A2")
modelcausalA0H <- omxSetParameters(modelcausalA0H,
                                   labels = c("covAxy6"), free = F, values = 0)
fitcausalA0H   <- mxRun( modelcausalA0H )

modelcausalA0I <- mxModel(fitcausalA0B, name = "CLPM_phenoAR_noCovA1")
modelcausalA0I <- omxSetParameters(modelcausalA0I,
                                   labels = c("covAxy1"), free = F, values = 0)
fitcausalA0I   <- mxRun( modelcausalA0I )

mxCompare(fitcausalA0B, list(fitcausalA0I, fitcausalA0D, 
                             fitcausalA0E, fitcausalA0F, 
                             fitcausalA0G, fitcausalA0H))

sumcausalA0H <- summary(fitcausalA0H)
sumcausalA0H

round(fitcausalA0H@algebras$stdBeta$result, 2)

round(mxSE(stdBeta, fitcausalA0H), 3)

round(fitcausalA0H@algebras$psiA$result, 2)

round(mxSE(psiA, fitcausalA0H), 3)

round(fitcausalA0H@algebras$psiE$result, 2)

round(mxSE(psiE, fitcausalA0H), 3)

outCausalA0H <- sumcausalA0H$parameters
outCausalA0H
write.csv(outCausalA0H, paste0(outDir,"CLPM_Std_Dep_Smk_CovA_TEDS21_26.csv"), quote = F)


#### SubModel 1: AR at AE ======================

modelcausalA1 <- mxModel(fitcausalA, name = "CLPM_aceAR")
modelcausalA1 <- omxSetParameters(modelcausalA1, 
                                  labels = c("AR_bx6x1","AR_by6y1"), free = F, values = 0)
modelcausalA1 <- omxSetParameters(modelcausalA1, 
                                  labels = c("AR_Ax6x1","AR_Ay6y1",
                                             "AR_Ex6x1","AR_Ey6y1"), free = T, values = 0.5)

fitcausalA1    <- mxRun( modelcausalA1)
sumcausalA1    <- summary( fitcausalA1 )
sumcausalA1

mxCompare(fitcausalA1, fitcausalA)

sumcausalA1

round(fitcausalA1@algebras$rA$result, 2)
round(fitcausalA1@algebras$rE$result, 2)
round(fitcausalA1@algebras$stdBeta$result, 2)
round(mxSE(stdBeta, fitcausalA1), 2)
round(fitcausalA1@algebras$stdBetaA$result, 2)
round(mxSE(stdBetaA, fitcausalA1), 2)
round(fitcausalA1@algebras$stdBetaE$result, 2)
round(mxSE(stdBetaE, fitcausalA1), 2)

outCausalA1 <- sumcausalA1$parameters
write.csv(outCausalA1, paste0(outDir,"CLPM_Biometric_Dep_Smk_TEDS21_26.csv"), quote = F)



#### SubModel 2: Cross-Lagged A ======================

modelcausalA2 <- mxModel(fitcausalA1, name = "CLPM_aceAR_CLAxy")
modelcausalA2 <- omxSetParameters(modelcausalA2, 
                                  labels = c("CL_Ay6x1","CL_Ax6y1"), free = T, values = 0.1)

fitcausalA2   <- mxRun( modelcausalA2 )
sumcausalA2   <- summary( fitcausalA2 )
sumcausalA2

round(fitcausalA2@algebras$stdBeta$result, 3)
round(fitcausalA2@algebras$stdBetaA$result, 3)
round(fitcausalA2@algebras$stdBetaE$result, 3)

mxCompare(fitcausalA2, fitcausalA1)


#### SubModel 3: Cross-Lagged E ======================

modelcausalA3 <- mxModel(fitcausalA1, name = "CLPM_aceAR_CLExy")
modelcausalA3 <- omxSetParameters(modelcausalA3, 
                                  labels = c("CL_Ey6x1","CL_Ex6y1"), free = T, values = 0.1)

fitcausalA3   <- mxRun( modelcausalA3 )
sumcausalA3   <- summary( fitcausalA3 )
sumcausalA3

round(fitcausalA3@algebras$stdBeta$result, 3)
round(fitcausalA3@algebras$stdBetaA$result, 3)
round(fitcausalA3@algebras$stdBetaE$result, 3)


#### SubModel 4: Cross-Lagged AE - No Phenotypical Causal Paths ======================

modelcausalA4 <- mxModel(fitcausalA2, name = "CLPM_aceAR_CLAExy")
modelcausalA4 <- omxSetParameters(modelcausalA4, 
                                  labels = c("CL_Ey6x1","CL_Ex6y1"), free = T, values = 0.2)
modelcausalA4 <- omxSetParameters(modelcausalA4, 
                                  labels = c("CL_by6x1","CL_bx6y1"), free = F, values = 0)
fitcausalA4   <- mxRun( modelcausalA4 )
sumcausalA4   <- summary( fitcausalA4 )
sumcausalA4

round(fitcausalA4@algebras$stdBeta$result, 2)
round(fitcausalA4@algebras$stdBetaA$result, 2)
round(mxSE(stdBetaA, model = fitcausalA4), 2)
round(fitcausalA4@algebras$stdBetaE$result, 2)
round(mxSE(stdBetaE, model = fitcausalA4), 2)

#### Compare Models ======================

mxCompare(fitcausalA4, list(fitcausalA1,fitcausalA2,fitcausalA3))

round(fitcausalA1$V$result, 2)
round(cov2cor(fitcausalA1$V$result), 2)
round(fitcausalA1$US$result, 2)

#### SubModel 5: No cross-lagged associations ======================

modelcausalA5 <- mxModel(fitcausalA1, name = "CLPM_aceAR_NoLagged")
modelcausalA5 <- omxSetParameters(modelcausalA5, 
                                  labels = c("CL_by6x1","CL_bx6y1"), free = F, values = 0)
fitcausalA5   <- mxRun( modelcausalA5 )

mxCompare(fitcausalA4, list(fitcausalA1,fitcausalA5))
mxCompare(fitcausalA1, list(fitcausalA5))

sumcausalA5 <- summary(fitcausalA5)
sumcausalA5


#### SubModel 6: DOC Proximal Effects at T1 ======================

modelcausalA6a <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC1yx_rA")
modelcausalA6a <- omxSetParameters(modelcausalA6a, 
                                   labels = c("covExy1"), free = F, values = 0)
modelcausalA6a <- omxSetParameters(modelcausalA6a, 
                                   labels = c("Inst_by1x1"), free = T, values = 0.2)
fitcausalA6a   <- mxRun( modelcausalA6a )
sumcausalA6a   <- summary(fitcausalA6a )
sumcausalA6a


modelcausalA6b <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC1xy_rA")
modelcausalA6b <- omxSetParameters(modelcausalA6b, 
                                   labels = c("covExy1"), free = F, values = 0)
modelcausalA6b <- omxSetParameters(modelcausalA6b, 
                                   labels = c("Inst_bx1y1"), free = T, values = -0.2)
fitcausalA6b   <- mxRun( modelcausalA6b )
sumcausalA6b   <- summary(fitcausalA6b )
sumcausalA6b


modelcausalA6c <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC1yx_rE")
modelcausalA6c <- omxSetParameters(modelcausalA6c, 
                                   labels = c("covAxy1"), free = F, values = 0)
modelcausalA6c <- omxSetParameters(modelcausalA6c, 
                                   labels = c("Inst_by1x1"), free = T, values = 0.2)
fitcausalA6c   <- mxRun( modelcausalA6c )
sumcausalA6c   <- summary(fitcausalA6c )
sumcausalA6c


modelcausalA6d <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC1xy_rE")
modelcausalA6d <- omxSetParameters(modelcausalA6d, 
                                   labels = c("covAxy1"), free = F, values = 0)
modelcausalA6d <- omxSetParameters(modelcausalA6d, 
                                   labels = c("Inst_bx1y1"), free = T, values = 0.1)
fitcausalA6d   <- mxRun( modelcausalA6d )
sumcausalA6d   <- summary(fitcausalA6d )
sumcausalA6d


modelcausalA6e <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC1bidir")
modelcausalA6e <- omxSetParameters(modelcausalA6e, 
                                   labels = c("covAxy1","covExy1"), free = F, values = 0)
modelcausalA6e <- omxSetParameters(modelcausalA6e, 
                                   labels = c("Inst_bx1y1","Inst_by1x1"), free = T, values = -.2)
fitcausalA6e   <- mxRun( modelcausalA6e )
sumcausalA6e   <- summary(fitcausalA6e )
sumcausalA6e

mxCompare(fitcausalA5, list(fitcausalA6e, 
                            fitcausalA6a, fitcausalA6c,
                            fitcausalA6b, fitcausalA6d ))

modelcausalA6f <- mxModel(fitcausalA6a, name = "CLPM_aceAR_rA1")
modelcausalA6f <- omxSetParameters(modelcausalA6f, 
                                   labels = c("Inst_by1x1"), free = F, values = 0)
fitcausalA6f   <- mxRun( modelcausalA6f )
sumcausalA6f   <- summary(fitcausalA6f )
sumcausalA6f


mxCompare(fitcausalA5, list(fitcausalA6a, fitcausalA6f))



#### SubModel 5: DOC Proximal Effects at T2 ======================

modelcausalA5a <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2yx_rA")
modelcausalA5a <- omxSetParameters(modelcausalA5a, 
                                   labels = c("covExy6"), free = F, values = 0)
modelcausalA5a <- omxSetParameters(modelcausalA5a, 
                                   labels = c("Inst_by6x6"), free = T, values = 0.2)
fitcausalA5a   <- mxRun( modelcausalA5a )
sumcausalA5a   <- summary(fitcausalA5a )
sumcausalA5a


modelcausalA5b <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2xy_rA")
modelcausalA5b <- omxSetParameters(modelcausalA5b, 
                                   labels = c("covExy6"), free = F, values = 0)
modelcausalA5b <- omxSetParameters(modelcausalA5b, 
                                   labels = c("Inst_bx6y6"), free = T, values = -0.2)
fitcausalA5b   <- mxRun( modelcausalA5b )
sumcausalA5b   <- summary(fitcausalA5b )
sumcausalA5b


modelcausalA5c <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2yx_rE")
modelcausalA5c <- omxSetParameters(modelcausalA5c, 
                                   labels = c("covAxy6"), free = F, values = 0)
modelcausalA5c <- omxSetParameters(modelcausalA5c, 
                                   labels = c("Inst_by6x6"), free = T, values = 0.2)
fitcausalA5c   <- mxRun( modelcausalA5c )
sumcausalA5c   <- summary(fitcausalA5c )
sumcausalA5c


modelcausalA5d <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2xy_rE")
modelcausalA5d <- omxSetParameters(modelcausalA5d, 
                                   labels = c("covAxy6"), free = F, values = 0)
modelcausalA5d <- omxSetParameters(modelcausalA5d, 
                                   labels = c("Inst_bx6y6"), free = T, values = 0.1)
fitcausalA5d   <- mxRun( modelcausalA5d )
sumcausalA5d   <- summary(fitcausalA5d )
sumcausalA5d


modelcausalA5e <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2bidir")
modelcausalA5e <- omxSetParameters(modelcausalA5e, 
                                   labels = c("covAxy6","covExy6"), free = F, values = 0)
modelcausalA5e <- omxSetParameters(modelcausalA5e, 
                                   labels = c("Inst_bx6y6","Inst_by6x6"), free = T, values = -.2)
fitcausalA5e   <- mxRun( modelcausalA5e )
sumcausalA5e   <- summary(fitcausalA5e )
sumcausalA5e

mxCompare(fitcausalA5, list(fitcausalA5e, 
                            fitcausalA5a, fitcausalA5c,
                            fitcausalA5b, fitcausalA5d ))

modelcausalA5f <- mxModel(modelcausalA5e, name = "CLPM_aceAR_DOC2xy")
modelcausalA5f <- omxSetParameters(modelcausalA5f, 
                                   labels = c("Inst_by6x6"), free = F, values = 0)
fitcausalA5f   <- mxRun( modelcausalA5f )
sumcausalA5f   <- summary(fitcausalA5f )
sumcausalA5f

mxCompare(fitcausalA5d, list(fitcausalA5e,fitcausalA5f))

sumcausalA5d
outCausalA5d <- sumcausalA5d$parameters
write.csv(outCausalA5d, paste0(outDir,"CLPM_Biometric_Dep_Smk_TEDS21_26_DOCSmkDepE.csv"), quote = F)

sumcausalA5e
outCausalA5e <- sumcausalA5e$parameters
write.csv(outCausalA5e, paste0(outDir,"CLPM_Biometric_Dep_Smk_TEDS21_26_DOCbidir.csv"), quote = F)


# Final Models ==========

#### Base Model ========================

round(fitcausalA@algebras$stdBeta$result, 2)
round(mxSE(stdBeta, model = fitcausalA), 3)

round(fitcausalA@algebras$rA$result, 2)
round(mxSE(rA, model = fitcausalA), 3)

round(fitcausalA@algebras$rE$result, 2)
round(mxSE(rE, model = fitcausalA), 3)

round(fitcausalA@algebras$psiA$result, 2)
round(mxSE(psiA, model = fitcausalA), 3)

round(fitcausalA@algebras$psiE$result, 2)
round(mxSE(psiE, model = fitcausalA), 3)

round(fitcausalA$US$result[,1:4], 2)
round(fitcausalA$US$result[,9:12], 2)

round(fitcausalA$US$result[,13:16], 2)
round(fitcausalA$US$result[,21:24], 2)


#### CLPM with AE Simplex ===============

round(fitcausalA1@algebras$stdBeta$result, 2)
round(mxSE(stdBeta, model = fitcausalA1), 3)

round(fitcausalA1@algebras$stdBetaA$result, 2)
round(mxSE(stdBetaA, model = fitcausalA1), 3)

round(fitcausalA1@algebras$stdBetaE$result, 2)
round(mxSE(stdBetaE, model = fitcausalA1), 3)

round(fitcausalA1@algebras$rA$result, 2)
round(mxSE(rA, model = fitcausalA1), 3)

round(fitcausalA1@algebras$rE$result, 2)
round(mxSE(rE, model = fitcausalA1), 3)

round(fitcausalA1$PSA$values, 2)
round(mxSE(PSA, model = fitcausalA1), 3)

round(fitcausalA1@algebras$psiA$result, 2)
round(mxSE(psiA, model = fitcausalA1), 3)

round(fitcausalA1$PSE$values, 2)
round(mxSE(PSE, model = fitcausalA1), 3)

round(fitcausalA1@algebras$psiE$result, 2)
round(mxSE(psiE, model = fitcausalA1), 3)

round(fitcausalA1$US$result[,1:4], 2)
round(fitcausalA1$US$result[,9:12], 2)

round(fitcausalA1$US$result[,13:16], 3)
round(fitcausalA1$US$result[,21:24], 3)
round(mxSE(US, fitcausalA1), 3)


#### CLPM with AE Simplex, No CL association ===============

round(fitcausalA5@algebras$stdBeta$result, 2)
round(mxSE(stdBeta, model = fitcausalA5), 3)

round(fitcausalA5@algebras$stdBetaA$result, 2)
round(mxSE(stdBetaA, model = fitcausalA5), 3)

round(fitcausalA5@algebras$stdBetaE$result, 2)
round(mxSE(stdBetaE, model = fitcausalA5), 3)

round(fitcausalA5@algebras$psiA$result, 2)
round(mxSE(psiA, model = fitcausalA5), 3)

round(fitcausalA5@algebras$psiE$result, 2)
round(mxSE(psiE, model = fitcausalA5), 3)

round(fitcausalA5$US$result[,13:16], 3)
round(fitcausalA5d$US$result[,21:24], 3)

#### CLPM + Bidir DOC model ============================

round(fitcausalA5e@algebras$stdBeta$result, 2)
round(mxSE(stdBeta, model = fitcausalA5e), 3)

round(fitcausalA5e@algebras$stdBetaA$result, 2)
round(mxSE(stdBetaA, model = fitcausalA5e), 3)

round(fitcausalA5e@algebras$stdBetaE$result, 2)
round(mxSE(stdBetaE, model = fitcausalA5e), 3)

round(fitcausalA5e@algebras$psiA$result, 2)
round(mxSE(psiA, model = fitcausalA5e), 3)

round(fitcausalA5e@algebras$psiE$result, 2)
round(mxSE(psiE, model = fitcausalA5e), 3)


#### Best-fitting CLPM + Smk -> Dep DOC model ============

round(fitcausalA5d@algebras$stdBeta$result, 2)
round(mxSE(stdBeta, model = fitcausalA5d), 3)

round(fitcausalA5d@algebras$stdBetaA$result, 2)
round(mxSE(stdBetaA, model = fitcausalA5d), 3)

round(fitcausalA5d@algebras$stdBetaE$result, 2)
round(mxSE(stdBetaE, model = fitcausalA5d), 3)

round(fitcausalA5d@algebras$psiA$result, 2)
round(mxSE(psiA, model = fitcausalA5d), 3)

round(fitcausalA5d@algebras$psiE$result, 2)
round(mxSE(psiE, model = fitcausalA5d), 3)

round(fitcausalA5d$US$result[,1:4], 2)
round(fitcausalA5d$US$result[,9:12], 2)

round(fitcausalA5d$US$result[,13:16], 3)

round(fitcausalA5d$US$result[,21:24], 3)

round(mxSE(US, fitcausalA5d), 3)

# END ===========

