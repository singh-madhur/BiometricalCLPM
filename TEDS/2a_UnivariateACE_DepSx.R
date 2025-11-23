### DepSx univariate twin
### Madhur Singh

## Code adapted from https://hermine-maes.squarespace.com/s/oneACEvca-k9hd.R

library(dplyr)
library(tidyr)
library(ggplot2)

library(OpenMx)
library(psych)
source("miFunctions.R")
mxOption( NULL, "Default optimizer", "CSOLNP" )

## Data ===========================

dataDir <- "Data/"
outDir  <- "Out/"
plotDir <- "Plots/"
logDir  <- "Logs/"

DepSx <- read.table(paste0(dataDir,"TEDS_DepSx_byFam_Jun2024.tsv"), header = T)
str(DepSx)


## Matrices ===========================

# Select Variables for Analysis
vars      <- 'DepSx'                   # list of variables names
nv        <- 1                         # number of variables
ntv       <- nv*2                      # number of total variables
selVars   <- paste(vars,"Tw",c(rep(1,nv),rep(2,nv)),sep="")

# Set Starting Values
svMe      <-  0                        # start value for means
svPa      <- .4                        # start value for variance component for a
svPc      <- .1                        # start value for variance component for c
svPe      <- .5                        # start value for variance component for e


# Create Algebra for expected Mean Matrices
meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean",vars), name="meanG" )

# Create Matrices for Variance Components
covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11", name="VA" ) 
covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPc, label="VC11", name="VC" )
covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VE11", name="VE" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= VA+VC+VE, name="V" )
covMZ     <- mxAlgebra( expression= VA+VC, name="cMZ" )
covDZ     <- mxAlgebra( expression= 0.5%x%VA+ VC, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )


## Start Loop ========================================

timePoints <- 1:6

for (t in 1:6) {
  
  TT <- timePoints[t]
  print(paste("Timepoint",TT))
  
  twinData <- DepSx |> 
    select(zygos, contains(paste0("DepSxResRN_T",TT)))
  
  
  ## Output log file ======================
  
  sink(paste0(logDir,"univar_DepSx_T",TT,".txt"), append=FALSE, split=TRUE)
  print(mxVersion())
  
  print(str(twinData))
  
  colnames(twinData) <- gsub(paste0("ResRN_T",TT,"."), "", colnames(twinData))
  
  twinData <- twinData |> 
    filter(!(is.na(DepSxTw1) & is.na(DepSxTw2)))
  
  print(head(twinData))
  
  ## Select Data for Analysis =======================
  
  mzData    <- subset(twinData, zygos=="MZ", selVars)
  dzData    <- subset(twinData, zygos=="DZ", selVars)
  
  print(dim(mzData))
  print(table(rowSums(is.na(mzData))))
  print(dim(dzData))
  print(table(rowSums(is.na(dzData))))
  
  # Generate Descriptive Statistics
  print(describe(mzData))
  print(describe(dzData))
  print(cov(mzData,use="complete"))
  print(cov(dzData,use="complete"))
  print(corTest(mzData), short = F)
  print(corTest(dzData), short = F)
  
  ## Create Data Objects for Multiple Groups ========================
  
  dataMZ    <- mxData( observed=mzData, type="raw" )
  dataDZ    <- mxData( observed=dzData, type="raw" )
  
  # Create Expectation Objects for Multiple Groups
  expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars )
  expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars )
  funML     <- mxFitFunctionML()
  
  # Create Model Objects for Multiple Groups
  pars      <- list( meanG, covA, covC, covE, covP )
  modelMZ   <- mxModel( pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
  modelDZ   <- mxModel( pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
  multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
  
  # Create Algebra for Unstandardized and Standardized Variance Components
  rowUS     <- vars
  colUS     <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
  estUS     <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )
  
  # Create Confidence Interval Objects
  ciACE     <- mxCI( "US[1,1:6]" )
  
  # Build Model with Confidence Intervals
  modelACE  <- mxModel( "oneACEvc", pars, modelMZ, modelDZ, multi, estUS, ciACE )
  
  ## RUN ACE ==========================
  
  # Run ACE Model
  fitACE    <- mxRun( modelACE, intervals=T )
  refModels <- mxRefModels(fitACE, run = T)
  sumACE    <- summary( fitACE , refModels = refModels)
  print(sumACE)
  
  ## RUN SUBMODELS ===============================
  
  # Run AE model
  modelAE   <- mxModel( fitACE, name="oneAEvc" )
  modelAE   <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
  fitAE     <- mxRun( modelAE, intervals=T )
  sumAE     <- summary( fitAE , refModels = refModels)
  print(sumAE)
  
  # Run CE model
  modelCE   <- mxModel( fitACE, name="oneCEvc" )
  modelCE   <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
  modelCE   <- omxSetParameters( modelCE, labels=c("VE11","VC11"), free=TRUE, values=.6 )
  fitCE     <- mxRun( modelCE, intervals=T )
  sumCE     <- summary( fitCE , refModels = refModels)
  print(sumCE)
  
  # Run E model
  modelE    <- mxModel( fitAE, name="oneEvc" )
  modelE    <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0 )
  fitE      <- mxRun( modelE, intervals=T )
  sumE      <- summary( fitE , refModels = refModels)
  print(sumE)
  
  # Print Comparative Fit Statistics
  print(mxCompare( fitACE, nested <- list(fitAE, fitCE, fitE) ))
  print(round(rbind(fitACE$US$result,fitAE$US$result,fitCE$US$result,fitE$US$result),4))
  
  ## End ====
  sink()
  
  print(paste("Ran time point", TT,"of 6"))
}



