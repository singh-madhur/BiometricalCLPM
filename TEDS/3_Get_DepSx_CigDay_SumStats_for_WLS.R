### Get sum stats for Biometrical CLPM with WLS
### DepSx (rank-normalized residuals) and CigDay (3-L ordinal)
### Madhur Singh

library(dplyr)
library(OpenMx)


# Data ===========================

dataDir <- "Data/"
outDir  <- "Out/"
plotDir <- "Plots/"
logDir  <- "Logs/"

DepSx <- read.table(paste0(dataDir,"TEDS_DepSx_byFam_Jun2024.tsv"), header = T)
str(DepSx)

DepSxT <- DepSx |> 
  select(randomfamid:age_T6.Tw2, contains("DepSxResRN"))

str(DepSxT)
dim(DepSxT)

CigDay <- read.table(paste0(dataDir,"TEDS_CigDay_byFam_Jun2024.tsv"), header = T)
str(CigDay)
CigDayT <- CigDay |> 
  select(randomfamid:age_T6.Tw2, contains("CigDay3L"))
str(CigDayT)
dim(CigDayT)

## Merge
DepSmk_ <- DepSxT |> 
  full_join(CigDayT)

dim(DepSmk_)
str(DepSmk_)


## Sample Size ==========================

DepSmk_ |> 
  filter(!(is.na(DepSxResRN_T1.Tw1) & is.na(DepSxResRN_T1.Tw2) &
             is.na(DepSxResRN_T2.Tw1) & is.na(DepSxResRN_T2.Tw2) &
             is.na(DepSxResRN_T3.Tw1) & is.na(DepSxResRN_T3.Tw2) &
             is.na(DepSxResRN_T4.Tw1) & is.na(DepSxResRN_T4.Tw2) &
             is.na(DepSxResRN_T5.Tw1) & is.na(DepSxResRN_T5.Tw2) &
             is.na(DepSxResRN_T6.Tw1) & is.na(DepSxResRN_T6.Tw2) &
             is.na(CigDay3L_T1.Tw1) & is.na(CigDay3L_T1.Tw2) &
             is.na(CigDay3L_T2.Tw1) & is.na(CigDay3L_T2.Tw2) &
             is.na(CigDay3L_T3.Tw1) & is.na(CigDay3L_T3.Tw2) &
             is.na(CigDay3L_T4.Tw1) & is.na(CigDay3L_T4.Tw2) &
             is.na(CigDay3L_T5.Tw1) & is.na(CigDay3L_T5.Tw2) &
             is.na(CigDay3L_T6.Tw1) & is.na(CigDay3L_T6.Tw2))) |> 
  nrow()

DepSmk_ |> 
  filter(!(is.na(DepSxResRN_T1.Tw1) & is.na(DepSxResRN_T1.Tw2) &
             is.na(DepSxResRN_T2.Tw1) & is.na(DepSxResRN_T2.Tw2) &
             is.na(DepSxResRN_T3.Tw1) & is.na(DepSxResRN_T3.Tw2) &
             is.na(DepSxResRN_T4.Tw1) & is.na(DepSxResRN_T4.Tw2) &
             is.na(DepSxResRN_T5.Tw1) & is.na(DepSxResRN_T5.Tw2) &
             is.na(DepSxResRN_T6.Tw1) & is.na(DepSxResRN_T6.Tw2) &
             is.na(CigDay3L_T1.Tw1) & is.na(CigDay3L_T1.Tw2) &
             is.na(CigDay3L_T2.Tw1) & is.na(CigDay3L_T2.Tw2) &
             is.na(CigDay3L_T3.Tw1) & is.na(CigDay3L_T3.Tw2) &
             is.na(CigDay3L_T4.Tw1) & is.na(CigDay3L_T4.Tw2) &
             is.na(CigDay3L_T5.Tw1) & is.na(CigDay3L_T5.Tw2) &
             is.na(CigDay3L_T6.Tw1) & is.na(CigDay3L_T6.Tw2))) |> 
  count(zygos, sex.Tw1, sex.Tw2)

## N. individuals with at least one DepSx or CigDay
DepSxInd <- read.table(paste0(dataDir,"TEDS_DepSx_byInd_Jun2024.tsv"), header = T)
str(DepSxInd)

DepSxIndT <- DepSxInd |> 
  select(randomfamid:sex, contains("DepSxResRN")) |> 
  tidyr::pivot_longer(cols = -c(randomfamid:sex), 
                      names_to = "wave", 
                      values_to = "DepSxResRN", 
                      names_prefix = "DepSxResRN_") |> 
  filter(!is.na(DepSxResRN))

glimpse(DepSxIndT)
dim(DepSxIndT)

DepSxIDs <- DepSxIndT |> 
  distinct(randomtwinid, .keep_all = T)
nrow(DepSxIDs)

CigDayInd <- read.table(paste0(dataDir,"TEDS_CigDay_byInd_Jun2024.tsv"), header = T)
str(CigDayInd)

CigDayIndT <- CigDayInd |> 
  select(randomfamid:sex, contains("CigDay3L")) |> 
  tidyr::pivot_longer(cols = -c(randomfamid:sex), 
                      names_to = "wave", 
                      values_to = "CigDay3L", 
                      names_prefix = "CigDay3L_") |> 
  filter(!is.na(CigDay3L))

glimpse(CigDayIndT)
dim(CigDayIndT)

CigDayIDs <- CigDayIndT |> 
  distinct(randomtwinid, .keep_all = T)
nrow(CigDayIDs)

## Combined
DepSmkIDs <- DepSxIDs |> 
  select(randomtwinid) |> 
  bind_rows(select(CigDayIDs, randomtwinid)) |> 
  distinct()
nrow(DepSmkIDs)

DepSxIDs |> 
  select(randomtwinid:sex) |> 
  bind_rows(select(CigDayIDs, randomtwinid:sex)) |> 
  distinct() |> 
  count(sex) |> 
  mutate(prop = n/sum(n))

DepSxIDs |> 
  select(randomfamid:sex) |> 
  bind_rows(select(CigDayIDs, randomfamid:sex)) |> 
  distinct() |> 
  nrow()

DepSxIDs |> 
  select(randomfamid:sex) |> 
  bind_rows(select(CigDayIDs, randomfamid:sex)) |> 
  distinct() |> 
  group_by(randomfamid) |> 
  mutate(Twin_nr = row_number()) |> 
  ungroup() |> 
  count(Twin_nr)

## Complete twin pairs
DepSxIDs |> 
  select(randomfamid:sex) |> 
  bind_rows(select(CigDayIDs, randomfamid:sex)) |> 
  distinct() |> 
  group_by(randomfamid) |> 
  mutate(Twin_nr = row_number()) |> 
  ungroup() |> 
  tidyr::pivot_wider(id_cols = c(randomfamid, zygos), names_from = Twin_nr, values_from = sex, names_prefix = "T") |> 
  filter(!is.na(T2)) |> 
  count(zygos, T1, T2)

## Unpaired individuals
DepSxIDs |> 
  select(randomfamid:sex) |> 
  bind_rows(select(CigDayIDs, randomfamid:sex)) |> 
  distinct() |> 
  group_by(randomfamid) |> 
  mutate(Twin_nr = row_number()) |> 
  ungroup() |> 
  tidyr::pivot_wider(id_cols = c(randomfamid, zygos), names_from = Twin_nr, values_from = sex, names_prefix = "T") |> 
  filter(is.na(T2)) |> 
  count(zygos, T1, T2)


#### Two-wave ============

DepSxIndT16 <- DepSxIndT |> 
  filter(wave %in%c("T1","T6"))
dim(DepSxIndT16)

DepSxIDs_T16 <- DepSxIndT16 |> 
  distinct(randomtwinid)
nrow(DepSxIDs_T16)


CigDayIndT16 <- CigDayIndT |> 
  filter(wave %in%c("T1","T6"))
dim(CigDayIndT16)

CigDayIDs_T16 <- CigDayIndT16 |> 
  distinct(randomtwinid)
nrow(CigDayIDs_T16)


DepSmkIDs_T16 <- DepSxIDs_T16 |> 
  bind_rows(CigDayIDs_T16) |> 
  distinct()
nrow(DepSmkIDs_T16)


DepSxIndT16 |> 
  select(randomtwinid:sex) |> 
  bind_rows(select(CigDayIndT16, randomtwinid:sex)) |> 
  distinct() |> 
  count(sex) |> 
  mutate(prop = n/sum(n))

DepSxIndT16 |> 
  select(randomfamid:sex) |> 
  bind_rows(select(CigDayIndT16, randomfamid:sex)) |> 
  distinct() |> 
  nrow()

DepSxIndT16 |> 
  select(randomfamid:sex) |> 
  bind_rows(select(CigDayIndT16, randomfamid:sex)) |> 
  distinct() |> 
  group_by(randomfamid) |> 
  mutate(Twin_nr = row_number()) |> 
  ungroup() |> 
  count(Twin_nr)

## Complete twin pairs
DepSxIndT16 |> 
  select(randomfamid:sex) |> 
  bind_rows(select(CigDayIndT16, randomfamid:sex)) |> 
  distinct() |> 
  group_by(randomfamid) |> 
  mutate(Twin_nr = row_number()) |> 
  ungroup() |> 
  tidyr::pivot_wider(id_cols = c(randomfamid, zygos), names_from = Twin_nr, values_from = sex, names_prefix = "T") |> 
  filter(!is.na(T2)) |> 
  count(zygos, T1, T2)

## Unpaired individuals
DepSxIDs |> 
  select(randomfamid:sex) |> 
  bind_rows(select(CigDayIDs, randomfamid:sex)) |> 
  distinct() |> 
  group_by(randomfamid) |> 
  mutate(Twin_nr = row_number()) |> 
  ungroup() |> 
  tidyr::pivot_wider(id_cols = c(randomfamid, zygos), names_from = Twin_nr, values_from = sex, names_prefix = "T") |> 
  filter(is.na(T2)) |> 
  count(zygos, T1, T2)



## Data ============

## mxData cannot have "." in the variable names
colnames(DepSmk_) <- gsub("\\.","_",colnames(DepSmk_))

head(DepSmk_[,19:54])

DepSmk <- DepSmk_ |> 
  select(-c(contains("random"))) |> 
  select(zygos, sex_Tw1, sex_Tw2, contains("age"), 
         contains("T1_Tw1"), contains("T2_Tw1"), contains("T3_Tw1"), 
         contains("T4_Tw1"), contains("T5_Tw1"), contains("T6_Tw1"), 
         contains("T1_Tw2"), contains("T2_Tw2"), contains("T3_Tw2"), 
         contains("T4_Tw2"), contains("T5_Tw2"), contains("T6_Tw2"))

head(DepSmk)

DepSmk |> 
  count(sex_Tw1, sex_Tw2)


## TEDS21 and TEDS26
DepSmk_T16 <- DepSmk |> 
  select(zygos, # zygSex,
         sex_Tw1, sex_Tw2, age_T16, 
         contains("T1_Tw1"), contains("T6_Tw1"),
         contains("T1_Tw2"), contains("T6_Tw2"))
head(DepSmk_T16)


# TEDS21 and TEDS26 ======================

# Select Continuous Variables
varsc     <- paste0("DepSxResRN_T",c(1,6))            # list of continuous variables names
conVars   <- c(paste0(varsc,"_Tw1"),paste0(varsc,"_Tw2"))
conVars

# Select Ordinal Variables
nth       <- 2                                     # number of thresholds
varso     <- paste0("CigDay3L_T",c(1,6))              # list of ordinal variables names
ordVars   <- c(paste0(varso,"_Tw1"),paste0(varso,"_Tw2"))
ordVars


#### MZ subset ===================

MZDepSmk_T16 <- DepSmk_T16 |> 
  filter(zygos == "MZ") |> 
  select( -c(zygos,sex_Tw1,sex_Tw2,contains("resPRS"),
             contains("age")) )
head(MZDepSmk_T16)

table(rowSums(is.na(MZDepSmk_T16)) < ncol(MZDepSmk_T16))

MZDepSmk_T16 <- MZDepSmk_T16[rowSums(is.na(MZDepSmk_T16)) < ncol(MZDepSmk_T16),]
dim(MZDepSmk_T16)

# Create Data Objects for Multiple Groups
MZDepSmkF_T16   <- cbind( MZDepSmk_T16[,conVars], 
                          mxFactor(x=MZDepSmk_T16[,ordVars], levels=c(0:nth))
) 
str(MZDepSmkF_T16)

## reorder
MZDepSmkF_T16 <- MZDepSmkF_T16 |> 
  select(contains("T1_Tw1"),contains("T6_Tw1"),
         contains("T1_Tw2"),contains("T6_Tw2"))
str(MZDepSmkF_T16)

## cor matrix
MZDepSmkF_T16_cor <- polycor::hetcor(MZDepSmkF_T16, 
                                     ML = T, 
                                     use = "pairwise.complete.obs", 
                                     parallel = T)
MZDepSmkF_T16_cor
save(MZDepSmkF_T16_cor, file = paste0("Out/MZDepSmkF_T16_cor.RData"))


MZdatDepSmk_T16 <- mxData(MZDepSmkF_T16, type = "raw")
str(MZdatDepSmk_T16)
MZdatDepSmk_T16$numObs

## Sum stats
MZsumStatDepSmk_T16 <- omxAugmentDataWithWLSSummary(MZdatDepSmk_T16,
                                                    allContinuousMethod='marginals',
                                                    type='WLS')

MZsumStatDepSmk_T16@observed <- NULL # Delete the raw data from the object

str(MZsumStatDepSmk_T16)
cov2cor(MZsumStatDepSmk_T16$observedStats$cov)

save(MZsumStatDepSmk_T16, file = paste0("Out/MZsumStatDepSmk_T16.RData"))



#### DZ subset ===================

DZDepSmk_T16 <- DepSmk_T16 |> 
  filter(zygos == "DZ") |> 
  select( -c(zygos,sex_Tw1,sex_Tw2,
             contains("age")) )
head(DZDepSmk_T16)

table(rowSums(is.na(DZDepSmk_T16)) < ncol(DZDepSmk_T16))

DZDepSmk_T16 <- DZDepSmk_T16[rowSums(is.na(DZDepSmk_T16)) < ncol(DZDepSmk_T16),]
dim(DZDepSmk_T16)

# Create Data Objects for Multiple Groups
DZDepSmkF_T16   <- cbind( DZDepSmk_T16[,conVars], 
                          mxFactor(x=DZDepSmk_T16[,ordVars], levels=c(0:nth))
) 
str(DZDepSmkF_T16)


## reorder
DZDepSmkF_T16 <- DZDepSmkF_T16 |> 
  select(contains("T1_Tw1"),contains("T6_Tw1"),
         contains("T1_Tw2"),contains("T6_Tw2"))
str(DZDepSmkF_T16)

## cor matrix
DZDepSmkF_T16_cor <- polycor::hetcor(DZDepSmkF_T16, 
                                     ML = T, 
                                     use = "pairwise.complete.obs", 
                                     parallel = T)
DZDepSmkF_T16_cor
save(DZDepSmkF_T16_cor, file = paste0("Out/DZDepSmkF_T16_cor.RData"))

DZdatDepSmk_T16 <- mxData(DZDepSmkF_T16, type = "raw")
str(DZdatDepSmk_T16)
DZdatDepSmk_T16$numObs

## Sum stats
DZsumStatDepSmk_T16 <- omxAugmentDataWithWLSSummary(DZdatDepSmk_T16,
                                                    allContinuousMethod='marginals',
                                                    type='WLS')

DZsumStatDepSmk_T16@observed <- NULL # Delete the raw data from the object

str(DZsumStatDepSmk_T16)
cov2cor(DZsumStatDepSmk_T16$observedStats$cov)

save(DZsumStatDepSmk_T16, file = paste0("Out/DZsumStatDepSmk_T16.RData"))



# All six waves ======================

# Select Continuous Variables
varsc     <- paste0("DepSxResRN_T",1:6)            # list of continuous variables names
conVars   <- c(paste0(varsc,"_Tw1"),paste0(varsc,"_Tw2"))
conVars

# Select Ordinal Variables
nth       <- 2                                     # number of thresholds
varso     <- paste0("CigDay3L_T",1:6)              # list of ordinal variables names
ordVars   <- c(paste0(varso,"_Tw1"),paste0(varso,"_Tw2"))
ordVars


#### MZ subset ===================

MZDepSmk <- DepSmk |> 
  filter(zygos == "MZ") |> 
  select( -c(zygos,sex_Tw1,sex_Tw2,
             contains("age")) )
head(MZDepSmk)

table(rowSums(is.na(MZDepSmk)) < ncol(MZDepSmk))

MZDepSmk <- MZDepSmk[rowSums(is.na(MZDepSmk)) < ncol(MZDepSmk),]
dim(MZDepSmk)

# Create Data Objects for Multiple Groups
MZDepSmkF   <- cbind( MZDepSmk[,conVars], 
                      mxFactor(x=MZDepSmk[,ordVars], levels=c(0:nth)) ) 
str(MZDepSmkF)

## reorder
MZDepSmkF <- MZDepSmkF |> 
  select(contains("T1_Tw1"),contains("T2_Tw1"),contains("T3_Tw1"),
         contains("T4_Tw1"),contains("T5_Tw1"),contains("T6_Tw1"),
         contains("T1_Tw2"),contains("T2_Tw2"),contains("T3_Tw2"),
         contains("T4_Tw2"),contains("T5_Tw2"),contains("T6_Tw2"))
str(MZDepSmkF)

## cor matrix
MZDepSmkF_cor <- polycor::hetcor(MZDepSmkF, 
                                 ML = T, 
                                 use = "pairwise.complete.obs", 
                                 parallel = T)
MZDepSmkF_cor
save(MZDepSmkF_cor, file = paste0("Out/MZDepSmkF_cor.RData"))

MZdatDepSmk <- mxData(MZDepSmkF, type = "raw")
str(MZdatDepSmk)
MZdatDepSmk$numObs

## Sum stats
MZsumStatDepSmk <- omxAugmentDataWithWLSSummary(MZdatDepSmk,
                                                allContinuousMethod='marginals',
                                                type='WLS')

MZsumStatDepSmk@observed <- NULL # Delete the raw data from the object

str(MZsumStatDepSmk)
cov2cor(MZsumStatDepSmk$observedStats$cov)

save(MZsumStatDepSmk, file = paste0("Out/MZsumStatDepSmk.RData"))



#### DZ subset ===================

DZDepSmk <- DepSmk |> 
  filter(zygos == "DZ") |> 
  select( -c(zygos,sex_Tw1,sex_Tw2,
             contains("age")) )
head(DZDepSmk)

table(rowSums(is.na(DZDepSmk)) < ncol(DZDepSmk))

# Create Data Objects for Multiple Groups
DZDepSmkF   <- cbind( DZDepSmk[,conVars], 
                      mxFactor(x=DZDepSmk[,ordVars], levels=c(0:nth)) ) 
str(DZDepSmkF)

## reorder
DZDepSmkF <- DZDepSmkF |> 
  select(contains("T1_Tw1"),contains("T2_Tw1"),contains("T3_Tw1"),
         contains("T4_Tw1"),contains("T5_Tw1"),contains("T6_Tw1"),
         contains("T1_Tw2"),contains("T2_Tw2"),contains("T3_Tw2"),
         contains("T4_Tw2"),contains("T5_Tw2"),contains("T6_Tw2"))
str(DZDepSmkF)

## cor matrix
DZDepSmkF_cor <- polycor::hetcor(DZDepSmkF, 
                                 ML = T, 
                                 use = "pairwise.complete.obs", 
                                 parallel = T)
DZDepSmkF_cor
save(DZDepSmkF_cor, file = paste0("Out/DZDepSmkF_cor.RData"))


DZdatDepSmk <- mxData(DZDepSmkF, type = "raw")
str(DZdatDepSmk)
DZdatDepSmk$numObs

## Sum stats
DZsumStatDepSmk <- omxAugmentDataWithWLSSummary(DZdatDepSmk,
                                                allContinuousMethod='marginals',
                                                type='WLS')

DZsumStatDepSmk@observed <- NULL # Delete the raw data from the object

str(DZsumStatDepSmk)
cov2cor(DZsumStatDepSmk$observedStats$cov)

save(DZsumStatDepSmk, file = paste0("Out/DZsumStatDepSmk.RData"))


# END. =========================================================================