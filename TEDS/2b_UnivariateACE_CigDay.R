### CigDay univariate twin
### Madhur Singh

## Code adapted from https://hermine-maes.squarespace.com/s/oneACEvc-lm27.R

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

CigDay <- read.table(paste0(dataDir,"TEDS_CigDay_byFam_Jun2024.tsv"), header = T)
str(CigDay)


## Matrices ===========================

# Select Variables for Analysis
vars      <- 'CigDay'                  # list of variables names
nv        <- 1                         # number of variables
ntv       <- nv*2                      # number of total variables
nth       <- 2                         # number of thresholds
selVars   <- paste0(vars,"Tw",c(rep(1,nv),rep(2,nv)))

covar     <- c("age", "sex") 
nCov      <- length(covar)
covVars   <- paste0(covar,"Tw",c(rep(1,nCov),rep(2,nCov)))

# Set Starting Values
svBe      <- 0.01                      # start value for regressions
svLTh     <- 1                         # start value for first threshold
svITh     <- 1                         # start value for increments
svTh      <- matrix(rep(c(svLTh,(rep(svITh,nth-1)))),nrow=nth,ncol=nv)     # start value for thresholds
lbTh      <- matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv)     # lower bounds for thresholds

svPa      <- .4                        # start value for variance component for a
svPc      <- .1                        # start value for variance component for c
svPe      <- .5                        # start value for variance component for e



# Create Matrices for Covariates and linear Regression Coefficients
defL   <- mxMatrix( type = "Full", nrow = nCov, ncol = 2, free = FALSE, 
                    labels = paste0("data.",covVars), name = "defL" )
pathBl <- mxMatrix( type = "Full", nrow = 1, ncol = nCov, free = TRUE, values = svBe,
                    labels = paste("b",covar,sep = "_"), name = "pathBl" )


# Create Algebra for expected Mean & Threshold Matrices
meanG     <- mxMatrix( type="Zero", nrow=1, ncol=ntv, name="meanG" )
expMean   <- mxAlgebra( expression= meanG + pathBl%*%defL, name="expMean" )
thinG     <- mxMatrix( type="Full", nrow=nth, ncol=ntv, free=TRUE, values=svTh, lbound=lbTh, labels=labTh("th",vars,nth), name="thinG" )
inc       <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=FALSE, values=1, name="inc" )
threG     <- mxAlgebra( expression= inc %*% thinG, name="threG" )

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

# Constrain Variance of Ordinal Variable Liability
var1     <- mxConstraint( expression=diag2vec(V)==1, name="Var1" )


## Start Loop ========================================

timePoints <- 1:6

for (t in 1:6) {
  
  TT <- timePoints[t]
  print(paste("Timepoint",TT))
  
  twinData <- CigDay |> 
    select(zygos, contains(paste0("CigDay3L_T",TT)), contains(paste0("age_T",TT)), contains("sex"))
  
  
  ## Output log file ======================
  
  sink(paste0(logDir,"univar_CigDay3L_T",TT,".txt"), append=FALSE, split=TRUE)
  print(mxVersion())
  
  print(str(twinData))
  
  colnames(twinData) <- gsub("\\.", "", colnames(twinData))
  colnames(twinData) <- gsub(paste0("_T",TT), "", colnames(twinData))
  colnames(twinData) <- gsub(paste0("3L"), "", colnames(twinData))
  
  ## Remove rows with missing definition variable
  twinData <- twinData |> 
    filter(!(is.na(CigDayTw1) & is.na(CigDayTw2)))
  
  print(dim(twinData))
  
  ## Remove rows with missing pheno on both twins
  twinData[which(is.na(twinData$CigDayTw1)), "ageTw1"] <- -99999
  twinData[which(is.na(twinData$CigDayTw2)), "ageTw2"] <- -99999
  
  ## Remove rows with pheno +nt but age missing
  twinData <- twinData |> 
    filter(!is.na(ageTw1) & !is.na(ageTw2) & !is.na(sexTw1) & !is.na(sexTw2))
  
  print(dim(twinData))
  
  ## Make sex numerical dummy variable
  twinData <- twinData |> 
    mutate(sexTw1 = case_when(sexTw1 == "Female" ~ 0,
                              sexTw1 == "Male" ~ 1),
           sexTw2 = case_when(sexTw2 == "Female" ~ 0,
                              sexTw2 == "Male" ~ 1))
  
  print(head(twinData))
  
  ## Select Data for Analysis =======================
  
  # Select Data for Analysis
  mzData    <- subset(twinData, zygos=="MZ", c(selVars, covVars))
  dzData    <- subset(twinData, zygos=="DZ", c(selVars, covVars))
  mzDataF   <- cbind(mxFactor( x=mzData[,selVars], levels=c(0:nth)), mzData[,covVars] ) 
  dzDataF   <- cbind(mxFactor( x=dzData[,selVars], levels=c(0:nth)), dzData[,covVars] ) 
  
  # Generate Descriptive Statistics
  print(sapply(mzDataF[,1:2],table))
  print(sapply(dzDataF[,1:2],table))
  print(polycor::hetcor(mzDataF[,1:2]))
  print(polycor::hetcor(dzDataF[,1:2]))
  print(dim(mzDataF))
  print(dim(dzDataF))
  
  print(table(rowSums(is.na(mzDataF[,1:2]))))
  print(table(rowSums(is.na(dzDataF[,1:2]))))
  
  ## Create Data Objects for Multiple Groups ========================
  
  dataMZ    <- mxData( observed=mzDataF, type="raw" )
  dataDZ    <- mxData( observed=dzDataF, type="raw" )
  
  # Create Expectation Objects for Multiple Groups
  expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars, thresholds="threG" )
  expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars, thresholds="threG" )
  funML     <- mxFitFunctionML()
  
  # Create Model Objects for Multiple Groups
  pars      <- list( pathBl, meanG, thinG, inc, threG, covA, covC, covE, covP )
  defs      <- list( defL )
  modelMZ   <- mxModel( pars, defs, expMean, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
  modelDZ   <- mxModel( pars, defs, expMean, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
  multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
  
  # Create Algebra for Unstandardized and Standardized Variance Components
  rowUS     <- rep('US',nv)
  colUS     <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
  estUS     <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )
  
  # Create Confidence Interval Objects
  ciACE     <- mxCI( "US[1,1:3]" )
  
  # Build Model with Confidence Intervals
  modelACE  <- mxModel( "oneACEvoa", pars, var1, modelMZ, modelDZ, multi, estUS, ciACE )
  
  ## RUN ACE ==========================
  
  # Run ACE Model
  fitACE    <- mxRun( modelACE, intervals=T )
  sumACE    <- summary( fitACE )
  print(sumACE)
  
  ## RUN SUBMODELS ===============================
  
  # Run AE model
  modelAE   <- mxModel( fitACE, name="oneAEvc" )
  modelAE   <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
  fitAE     <- mxRun( modelAE, intervals=T )
  sumAE     <- summary( fitAE )
  print(sumAE)
  
  # Run CE model
  modelCE   <- mxModel( fitACE, name="oneCEvc" )
  modelCE   <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
  modelCE   <- omxSetParameters( modelCE, labels=c("VE11","VC11"), free=TRUE, values=.6 )
  fitCE     <- mxRun( modelCE, intervals=T )
  sumCE     <- summary( fitCE )
  print(sumCE)
  
  # Run E model
  modelE    <- mxModel( fitAE, name="oneEvc" )
  modelE    <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0 )
  fitE      <- mxRun( modelE )
  sumE      <- summary( fitE )
  print(sumE)
  
  # Print Comparative Fit Statistics
  print(mxCompare( fitACE, nested <- list(fitAE, fitCE, fitE) ))
  print(round(rbind(fitACE$US$result,fitAE$US$result,fitCE$US$result,fitE$US$result),4))
  
  ## End ====
  sink()
  
  print(paste("Ran time point", TT,"of 6"))
}



