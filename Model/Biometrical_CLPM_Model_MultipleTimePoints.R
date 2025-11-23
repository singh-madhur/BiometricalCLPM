# Biometrical CLPM with Two continuous variables
# MadhurBain Singh, Conor Dolan

# Clear Workspace
rm(list=ls(all=TRUE))

# 
library(MASS)
library(psych)
library(OpenMx)


##### Read in the data =================

## Subset by Zygosity
datmz_
datdz_

# Specify the Number of time points (e.g., 3)
Ts <- 3

##### N. of vars =========================

## Assuming two constructs, X and Y, over Ts time points

nv=2*Ts;  nv
ntv=nv*2; ntv
labs <- paste0(rep(c("X","Y"),Ts),
               rep(1:Ts, each=2))
labs


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


#### Means matrix ===============

# Create regression model for expected Mean Matrices
B0_  <- mxMatrix( 
  type="Full", nrow=1, ncol=nv, free=TRUE, values=svb0, 
  labels = paste0("b0_",labs),
  name="b0", 
  dimnames = list("b0",labs)
) 
B0_

# Equal means by twin order
ExpMean <- mxAlgebra(expression=cbind(b0,b0), name='expMean')


#### Phenotypic Beta matrix ===============

## Beta Values
B_vals <- matrix(0, nv, nv)
rownames(B_vals) = colnames(B_vals) = labs

for (i in 1:Ts) {
  ## Lagged first-order paths
  if(i < Ts) {
    ## AR1_x
    B_vals[paste0("X",i+1), paste0("X",i)] <- bphxx
    ## AR1_y
    B_vals[paste0("Y",i+1), paste0("Y",i)] <- bphyy
    ## X_i --> Y_i+1
    B_vals[paste0("Y",i+1), paste0("X",i)] <- bphyx
    ## Y_i --> X_i+1
    B_vals[paste0("X",i+1), paste0("Y",i)] <- bphxy
  }
  ## Cross-sectional causal paths (for later use with DoC)
  ## X_i --> Y_i
  B_vals[paste0("Y",i), paste0("X",i)] <- 0
  ## Y_i --> X_i
  B_vals[paste0("X",i), paste0("Y",i)] <- 0
  
} ### END Vals Loop

B_vals


## Beta Free Parameters
B_free <- matrix(FALSE, nv, nv)
rownames(B_free) = colnames(B_free) = labs

for (i in 1:Ts) {
  ## Lagged first-order paths
  if(i < Ts) {
    ## AR1_x
    B_free[paste0("X",i+1), paste0("X",i)] <- TRUE
    ## AR1_y
    B_free[paste0("Y",i+1), paste0("Y",i)] <- TRUE
    ## X_i --> Y_i+1
    B_free[paste0("Y",i+1), paste0("X",i)] <- TRUE
    ## Y_i --> X_i+1
    B_free[paste0("X",i+1), paste0("Y",i)] <- TRUE
  }
  ## Cross-sectional causal paths (for later use with DoC)
  ## X_i --> Y_i
  B_free[paste0("Y",i), paste0("X",i)] <- FALSE
  ## Y_i --> X_i
  B_free[paste0("X",i), paste0("Y",i)] <- FALSE
  
} ### END Free Loop
B_free


## Beta Labels
B_labs <- matrix(NA_real_, nv, nv)
rownames(B_labs) = colnames(B_labs) = labs

for (i in 1:Ts) {
  ## Lagged first-order paths
  if(i < Ts) {
    ## AR1_x
    B_labs[paste0("X",i+1), paste0("X",i)] <- paste0("AR_bx",i+1,"x",i)
    ## AR1_y
    B_labs[paste0("Y",i+1), paste0("Y",i)] <- paste0("AR_by",i+1,"y",i)
    ## X_i --> Y_i+1
    B_labs[paste0("Y",i+1), paste0("X",i)] <- paste0("CL_by",i+1,"x",i)
    ## Y_i --> X_i+1
    B_labs[paste0("X",i+1), paste0("Y",i)] <- paste0("CL_bx",i+1,"y",i)
  }
  ## Cross-sectional causal paths
  ## X_i --> Y_i
  B_labs[paste0("Y",i), paste0("X",i)] <- paste0("Inst_by",i,"x",i)
  ## Y_i --> X_i
  B_labs[paste0("X",i), paste0("Y",i)] <- paste0("Inst_bx",i,"y",i)
  
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
  covA_vals[paste0("AX",i),paste0("AX",i)] <- vax
  covA_vals[paste0("AY",i),paste0("AX",i)] <- covaxy
  covA_vals[paste0("AX",i),paste0("AY",i)] <- covaxy
  covA_vals[paste0("AY",i),paste0("AY",i)] <- vax
}
isSymmetric(covA_vals)
covA_vals


# Labels
covA_labs <- matrix(NA_real_, nv, nv)
rownames(covA_labs) = colnames(covA_labs) = paste0("A",labs)
covA_labs

for(i in 1:Ts) {
  covA_labs[paste0("AX",i),paste0("AX",i)] <- paste0("VAx",i)
  covA_labs[paste0("AY",i),paste0("AX",i)] <- paste0("covAxy",i)
  covA_labs[paste0("AX",i),paste0("AY",i)] <- paste0("covAxy",i)
  covA_labs[paste0("AY",i),paste0("AY",i)] <- paste0("VAy",i)
}
isSymmetric(covA_labs)
covA_labs


# Free
covA_free <- matrix(FALSE, nv, nv)
rownames(covA_free) = colnames(covA_free) = paste0("A",labs)
covA_free

for(i in 1:Ts) {
  covA_free[paste0("AX",i),paste0("AX",i)] <- TRUE
  covA_free[paste0("AY",i),paste0("AX",i)] <- TRUE
  covA_free[paste0("AX",i),paste0("AY",i)] <- TRUE
  covA_free[paste0("AY",i),paste0("AY",i)] <- TRUE
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
  BA_labs[paste0("AX",i+1), paste0("AX",i)] <- paste0("AR_Ax",i+1,"x",i)
  ## AR1_y
  BA_labs[paste0("AY",i+1), paste0("AY",i)] <- paste0("AR_Ay",i+1,"y",i)
  ## X_i --> Y_i+1
  BA_labs[paste0("AY",i+1), paste0("AX",i)] <- paste0("CL_Ay",i+1,"x",i)
  ## Y_i --> X_i+1
  BA_labs[paste0("AX",i+1), paste0("AY",i)] <- paste0("CL_Ax",i+1,"y",i)
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
  covE_vals[paste0("EX",i),paste0("EX",i)] <- vex
  covE_vals[paste0("EY",i),paste0("EX",i)] <- covexy
  covE_vals[paste0("EX",i),paste0("EY",i)] <- covexy
  covE_vals[paste0("EY",i),paste0("EY",i)] <- vex
}
isSymmetric(covE_vals)
covE_vals


# Labels
covE_labs <- matrix(NA_real_, nv, nv)
rownames(covE_labs) = colnames(covE_labs) = paste0("E",labs)
covE_labs

for(i in 1:Ts) {
  covE_labs[paste0("EX",i),paste0("EX",i)] <- paste0("VEx",i)
  covE_labs[paste0("EY",i),paste0("EX",i)] <- paste0("covExy",i)
  covE_labs[paste0("EX",i),paste0("EY",i)] <- paste0("covExy",i)
  covE_labs[paste0("EY",i),paste0("EY",i)] <- paste0("VEy",i)
}
isSymmetric(covE_labs)
covE_labs


# Free
covE_free <- matrix(FALSE, nv, nv)
rownames(covE_free) = colnames(covE_free) = paste0("E",labs)
covE_free

for(i in 1:Ts) {
  covE_free[paste0("EX",i),paste0("EX",i)] <- TRUE
  covE_free[paste0("EY",i),paste0("EX",i)] <- TRUE
  covE_free[paste0("EX",i),paste0("EY",i)] <- TRUE
  covE_free[paste0("EY",i),paste0("EY",i)] <- TRUE
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
  BE_labs[paste0("EX",i+1), paste0("EX",i)] <- paste0("AR_Ex",i+1,"x",i)
  ## AR1_y
  BE_labs[paste0("EY",i+1), paste0("EY",i)] <- paste0("AR_Ey",i+1,"y",i)
  ## X_i --> Y_i+1
  BE_labs[paste0("EY",i+1), paste0("EX",i)] <- paste0("CL_Ey",i+1,"x",i)
  ## Y_i --> X_i+1
  BE_labs[paste0("EX",i+1), paste0("EY",i)] <- paste0("CL_Ex",i+1,"y",i)
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
  covC_vals[paste0("CX",i),paste0("CX",i)] <- vcx
  covC_vals[paste0("CY",i),paste0("CX",i)] <- covcxy
  covC_vals[paste0("CX",i),paste0("CY",i)] <- covcxy
  covC_vals[paste0("CY",i),paste0("CY",i)] <- vcx
}
isSymmetric(covC_vals)
covC_vals


# Labels
covC_labs <- matrix(NA_real_, nv, nv)
rownames(covC_labs) = colnames(covC_labs) = paste0("C",labs)
covC_labs

for(i in 1:Ts) {
  covC_labs[paste0("CX",i),paste0("CX",i)] <- paste0("VCx",i)
  covC_labs[paste0("CY",i),paste0("CX",i)] <- paste0("covCxy",i)
  covC_labs[paste0("CX",i),paste0("CY",i)] <- paste0("covCxy",i)
  covC_labs[paste0("CY",i),paste0("CY",i)] <- paste0("VCy",i)
}
isSymmetric(covC_labs)
covC_labs


# Free
covC_free <- matrix(FALSE, nv, nv)
rownames(covC_free) = colnames(covC_free) = paste0("C",labs)
covC_free

for(i in 1:Ts) {
  covC_free[paste0("CX",i),paste0("CX",i)] <- TRUE
  covC_free[paste0("CY",i),paste0("CX",i)] <- TRUE
  covC_free[paste0("CX",i),paste0("CY",i)] <- TRUE
  covC_free[paste0("CY",i),paste0("CY",i)] <- TRUE
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
  BC_labs[paste0("CX",i+1), paste0("CX",i)] <- paste0("AR_Cx",i+1,"x",i)
  ## AR1_y
  BC_labs[paste0("CY",i+1), paste0("CY",i)] <- paste0("AR_Cy",i+1,"y",i)
  ## X_i --> Y_i+1
  BC_labs[paste0("CY",i+1), paste0("CX",i)] <- paste0("CL_Cy",i+1,"x",i)
  ## Y_i --> X_i+1
  BC_labs[paste0("CX",i+1), paste0("CY",i)] <- paste0("CL_Cx",i+1,"y",i)
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
covA = mxAlgebra(expression= (BAi)%*%PSA%*%t(BAi), name='VA')
covC = mxAlgebra(expression= (BCi)%*%PSC%*%t(BCi), name='VC')
covE = mxAlgebra(expression= (BEi)%*%PSE%*%t(BEi), name='VE')


# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP1      <- mxAlgebra( expression= VA+VC+VE, name="V1" )
iBE        <- mxAlgebra(expression = solve(I4-BE), name='iBE') # iBE is not BEi
#
covP       <- mxAlgebra(expression= iBE%*%V1%*%t(iBE), name='V')
covMZ      <- mxAlgebra( expression= iBE%*%(VA+VC)%*%t(iBE), name="cMZ" )
covDZ      <- mxAlgebra( expression= iBE%*%(0.5%x%VA+VC)%*%t(iBE), name="cDZ" )
#
expCovMZ   <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ   <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

# Create Data Objects for Multiple Groups
dataMZ = mxData( observed=datmz_, type="raw")
dataDZ = mxData( observed=datdz_, type="raw")

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars_ )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars_ )
funML     <- mxFitFunctionML()

# combine all parameters
pars      <- list(B0_, covPSA, BEA, BAi, covA, 
                  covPSC, BEC, BCi, covC, 
                  covPSE, BEE, BEi, covE, 
                  covP1, covP, I4, BE, iBE )
modelMZ   <- mxModel( pars, ExpMean, covMZ, expCovMZ, ExpMean, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, ExpMean, covDZ, expCovDZ, ExpMean, dataDZ, expDZ, funML, name="DZ" )	
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Build Model with Confidence Intervals?

modelcausalA  <- mxModel(name="CLPM_phenoAR", pars, modelMZ, modelDZ, multi)

#### RUN MODEL =========================
fitcausalA    <- mxRun( modelcausalA) 
# fitcausalA    <- mxTryHard( modelcausalA) 
sumcausalA    <- summary( fitcausalA )
sumcausalA

# test ID
doID=T
if (doID) {
  testID2=mxCheckIdentification(fitcausalA, details=T)
}  

# reference model
m1_satA=mxRefModels(fitcausalA,run=TRUE) 
mxCompare(m1_satA,fitcausalA) 
summary(fitcausalA, refModels = m1_satA)


#### SubModel 1: AR at ACE ======================
modelcausalA1 <- mxModel(modelcausalA, name = "CLPM_aceAR")
modelcausalA1 <- omxSetParameters(modelcausalA1, 
                                  labels = c("AR_bx2x1","AR_by2y1",
                                             "AR_bx3x2","AR_by3y2"), free = F, values = 0)
modelcausalA1 <- omxSetParameters(modelcausalA1, 
                                  labels = c("AR_Ax2x1","AR_Ay2y1",
                                             "AR_Cx2x1","AR_Cy2y1",
                                             "AR_Ex2x1","AR_Ey2y1",
                                             "AR_Ax3x2","AR_Ay3y2",
                                             "AR_Cx3x2","AR_Cy3y2",
                                             "AR_Ex3x2","AR_Ey3y2"), free = T, values = 0.5)

fitcausalA1    <- mxRun( modelcausalA1)
sumcausalA1    <- summary( fitcausalA1 )
sumcausalA1


#### SubModel 2: Cross-Lagged A ======================

modelcausalA2 <- mxModel(modelcausalA1, name = "CLPM_aceAR_CLAxy")
modelcausalA2 <- omxSetParameters(modelcausalA2, 
                                  labels = c("CL_Ay2x1","CL_Ax2y1",
                                             "CL_Ay3x2","CL_Ax3y2"), free = T, values = 0.1)

fitcausalA2   <- mxRun( modelcausalA2 )
sumcausalA2   <- summary( fitcausalA2 )
sumcausalA2



#### SubModel 3: Cross-Lagged AE ======================

modelcausalA3 <- mxModel(modelcausalA2, name = "CLPM_aceAR_CLAExy")
modelcausalA3 <- omxSetParameters(modelcausalA3, 
                                  labels = c("CL_Ey2x1","CL_Ex2y1",
                                             "CL_Ey3x2","CL_Ex3y2"), free = T, values = 0.1)

fitcausalA3   <- mxRun( modelcausalA3 )
sumcausalA3   <- summary( fitcausalA3 )
sumcausalA3


#### SubModel 4: Cross-Lagged AC ======================

modelcausalA4 <- mxModel(modelcausalA2, name = "CLPM_aceAR_CLACxy")
modelcausalA4 <- omxSetParameters(modelcausalA4, 
                                  labels = c("CL_Cy2x1","CL_Cx2y1",
                                             "CL_Cy3x2","CL_Cx3y2"), free = T, values = 0.1)

fitcausalA4   <- mxRun( modelcausalA4 )
sumcausalA4   <- summary( fitcausalA4 )
sumcausalA4


#### SubModel 5: Cross-Lagged ACE - No Phenotypical Causal Paths ======================
modelcausalA5 <- mxModel(modelcausalA3, name = "CLPM_aceAR_CLACExy")
modelcausalA5 <- omxSetParameters(modelcausalA5, 
                                  labels = c("CL_Cy2x1","CL_Cx2y1",
                                             "CL_Cy3x2","CL_Cx3y2"), free = T, values = 0.1)
modelcausalA5 <- omxSetParameters(modelcausalA5, 
                                  labels = c("CL_by2x1","CL_bx2y1",
                                             "CL_by3x2","CL_bx3y2"), free = F, values = 0)
fitcausalA5   <- mxRun( modelcausalA5 )
sumcausalA5   <- summary( fitcausalA5 )
sumcausalA5


#### SubModel 6: DOC Proximal Effects at T2 ======================

modelcausalA6a <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2yx_rA")
modelcausalA6a <- omxSetParameters(modelcausalA6a, 
                                   labels = c("covExy2"), free = F, values = 0)
modelcausalA6a <- omxSetParameters(modelcausalA6a, 
                                   labels = c("Inst_by2x2"), free = T, values = 0.2)
fitcausalA6a   <- mxRun( modelcausalA6a )
sumcausalA6a   <- summary(fitcausalA6a )
sumcausalA6a


modelcausalA6b <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2xy_rA")
modelcausalA6b <- omxSetParameters(modelcausalA6b, 
                                   labels = c("covExy2"), free = F, values = 0)
modelcausalA6b <- omxSetParameters(modelcausalA6b, 
                                   labels = c("Inst_bx2y2"), free = T, values = -0.2)
fitcausalA6b   <- mxRun( modelcausalA6b )
sumcausalA6b   <- summary(fitcausalA6b )
sumcausalA6b


modelcausalA6c <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2yx_rE")
modelcausalA6c <- omxSetParameters(modelcausalA6c, 
                                   labels = c("covAxy2"), free = F, values = 0)
modelcausalA6c <- omxSetParameters(modelcausalA6c, 
                                   labels = c("Inst_by2x2"), free = T, values = 0.2)
fitcausalA6c   <- mxRun( modelcausalA6c )
sumcausalA6c   <- summary(fitcausalA6c )
sumcausalA6c


modelcausalA6d <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2xy_rE")
modelcausalA6d <- omxSetParameters(modelcausalA6d, 
                                   labels = c("covAxy2"), free = F, values = 0)
modelcausalA6d <- omxSetParameters(modelcausalA6d, 
                                   labels = c("Inst_bx2y2"), free = T, values = 0.1)
fitcausalA6d   <- mxRun( modelcausalA6d )
sumcausalA6d   <- summary(fitcausalA6d )
sumcausalA6d


modelcausalA6e <- mxModel(fitcausalA5, name = "CLPM_aceAR_DOC2bidir")
modelcausalA6e <- omxSetParameters(modelcausalA6e, 
                                   labels = c("covAxy2","covExy2"), free = F, values = 0)
modelcausalA6e <- omxSetParameters(modelcausalA6e, 
                                   labels = c("Inst_bx2y2","Inst_by2x2"), free = T, values = -.2)
fitcausalA6e   <- mxRun( modelcausalA6e )
sumcausalA6e   <- summary(fitcausalA6e )
sumcausalA6e

mxCompare(fitcausalA5, list(fitcausalA6e, 
                            fitcausalA6a, fitcausalA6c,
                            fitcausalA6b, fitcausalA6d ))

# END ===========
