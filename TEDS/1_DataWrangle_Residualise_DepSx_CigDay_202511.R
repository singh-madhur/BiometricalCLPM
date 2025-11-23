## Residualize and transform the variables
## Madhur Singh

library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)

dataDir <- "Data/"
outDir  <- "Out/"
plotDir <- "Plots/"
logDir  <- "Logs/"

## Read in
list.files(dataDir)

teds <- haven::read_spss(paste0(dataDir,"teds_v2_Jun2024.sav"))
dim(teds)
names(teds)

teds |> 
  distinct(randomfamid) |> 
  nrow()

## Background variables =======================
# https://www.teds.ac.uk/datadictionary/studies/variable_lists/background_variables.htm 

# Exclusion
teds |> 
  distinct(randomtwinid, .keep_all = T) |> 
  count(exclude1,exclude2)
## Recommended exclusion from all TEDS analyses

## Parent-reported ethnicity at 1st contact
# 1=white, 0=other.
teds |> 
  distinct(randomfamid, .keep_all = T) |> 
  count(aethnic)
## No further info available on "Other" ehnic groups
## Not appropriate to consider all non-White ethnicities as a single group.

tedsIncl <- teds |> 
  filter(exclude1==0 & exclude2==0 & aethnic==1) |> 
  as.data.frame()

dim(tedsIncl)

# Sex (biological sex)
# 0=female, 1=male (or missing if unknown)
tedsIncl |> 
  distinct(randomtwinid, .keep_all = T) |> 
  count(sex1)


## Dep Sx ========================================================================

### TEDS21 - COVID assessed 8 out of 13 sMFQ items
tedsIncl |> 
  select(contains("u2cmfq")) |> 
  head()

tedsIncl |> 
  select(contains("ucv4mf")) |> 
  head()

### TEDS26 assessed all 13 MFQ items
tedsIncl |> 
  select(contains("zmhmfq")) |> 
  head()

summary(tedsIncl$u2cmfqt1)  ## Range 0-16
summary(tedsIncl$ucv4mfqt1) ## Range 0-16
summary(tedsIncl$zmhmfqt1)  ## Range 0-26

### TEDS 21: https://datadictionary.teds.ac.uk/studies/variable_lists/21yr_variables.htm
# MFQ item 1, miserable (TEDS21 phase 1 twin qnr)
# MFQ item 2, did nothing (TEDS21 phase 1 twin qnr)
# MFQ item 3, restless (TEDS21 phase 1 twin qnr)
# MFQ item 4, cried (TEDS21 phase 1 twin qnr)
# MFQ item 5, hard to concentrate (TEDS21 phase 1 twin qnr)
# MFQ item 6, hated myself (TEDS21 phase 1 twin qnr)
# MFQ item 7, lonely (TEDS21 phase 1 twin qnr)
# MFQ item 8, not as good as others (TEDS21 phase 1 twin qnr)

### TEDS 26: https://datadictionary.teds.ac.uk/studies/variable_lists/26yr_variables.htm
# MFQ item 1: miserable (TEDS26 twin MHQ)
# MFQ item 2: did not enjoy anything (TEDS26 twin MHQ)
# MFQ item 3: tired (TEDS26 twin MHQ), see value labels
# MFQ item 4: restless (TEDS26 twin MHQ)
# MFQ item 5: felt no good (TEDS26 twin MHQ)
# MFQ item 6: cried a lot (TEDS26 twin MHQ)
# MFQ item 7: hard to concentrate (TEDS26 twin MHQ)
# MFQ item 8: hated myself (TEDS26 twin MHQ)
# MFQ item 9: bad person (TEDS26 twin MHQ)
# MFQ item 10: lonely (TEDS26 twin MHQ)
# MFQ item 11: nobody loved me (TEDS26 twin MHQ)
# MFQ item 12: never be as good as others (TEDS26 twin MHQ)
# MFQ item 13: did everything wrong (TEDS26 twin MHQ)

### Overlapping items in TEDS26
# 1, 2, 4, 6, 7, 8, 10, 12


##### Sum of the overlapping 8 items in TEDS26 ==========

### Sum of the overlapping 8 items from the short MFQ in TEDS21 and COVID
# Original Formula = https://datadictionary.teds.ac.uk/studies/derived_variables/26yr_derived_variables.htm#zmhmfqt
# 13 * (Mean of the 13 items), if >50% items are answered (if at least 7 are non-missing)

## Use the new formula
# 8 * (Mean of the 8 items), if at least 4 are non-missing
# At least half the items are required to be non-missing for the scale to be computed.

depT1 <- tedsIncl |> 
  select(randomtwinid, zmhmfq011, zmhmfq021, zmhmfq041, zmhmfq061, 
         zmhmfq071, zmhmfq081, zmhmfq101, zmhmfq121)
dim(depT1)

depT2 <- tedsIncl |> 
  select(randomtwinid, zmhmfq012, zmhmfq022, zmhmfq042, zmhmfq062, 
         zmhmfq072, zmhmfq082, zmhmfq102, zmhmfq122)
dim(depT2)


### Twin 1
depT1 <- depT1 |> 
  rowwise() |> 
  mutate(
    non_missing_count = sum(!is.na(c_across(zmhmfq011:zmhmfq121))),
    zmhmfq_t81 = if_else(
      non_missing_count >= 4,
      8 * mean(c_across(zmhmfq011:zmhmfq121), na.rm = TRUE),
      NA_real_
    )
  ) |> 
  ungroup()

summary(depT1$zmhmfq_t81)

### Twin 2
depT2 <- depT2 |> 
  rowwise() |> 
  mutate(
    non_missing_count = sum(!is.na(c_across(zmhmfq012:zmhmfq122))),
    zmhmfq_t82 = if_else(
      non_missing_count >= 4,
      8 * mean(c_across(zmhmfq012:zmhmfq122), na.rm = TRUE),
      NA_real_
    )
  ) |> 
  ungroup()

summary(depT2$zmhmfq_t82)

## combine the two twin
depTEDS26 <- depT1 |> 
  select(randomtwinid, zmhmfq_t81) |> 
  full_join(depT2 |> select(randomtwinid, zmhmfq_t82))
head(depTEDS26)

table(tedsIncl$randomtwinid == depTEDS26$randomtwinid)


#### Examine using the Total of 8 =======================

depSxT <- tedsIncl |> 
  select(randomfamid,randomtwinid,sex1,sex2,zygos,
         contains("age"),contains("mfqt")) 

depSxT2 <- depSxT |> 
  full_join(depTEDS26)

depSxT2 |> 
  select(u2cmfqt1, zmhmfq_t81) |> 
  rename(TEDS21_DepSx = u2cmfqt1,
         TEDS26_DepSx = zmhmfq_t81) |> 
  pivot_longer(
    cols = TEDS21_DepSx:TEDS26_DepSx,
    # cols_vary = "slowest",
    names_to = c("wave", ".value"),
    names_pattern = "(.*)_(.*)"
  ) |> 
  ggplot(aes(DepSx, fill = wave)) +
  geom_density() + 
  scale_fill_viridis_d(alpha = 0.5, direction = -1) + 
  theme_bw(12) +
  labs(x = "Depressive Symptoms",
       fill = NULL) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85,0.9),
        axis.title.x = element_text(face = "bold"))
ggsave(paste0(plotDir,"DepSx_TEDS21_v_26.jpeg"), 
       width = 6, height = 6, units = "in", dpi = 300)

depSx <- depSxT2 |> 
  select(-c(zmhmfqt1,zmhmfqt2)) |> 
  filter( 
    !( is.na(u2cmfqt1) & is.na(u2cmfqt2) &
         is.na(ucv1mfqt1) & is.na(ucv1mfqt2) &
         is.na(ucv2mfqt1) & is.na(ucv2mfqt2) &
         is.na(ucv3mfqt1) & is.na(ucv3mfqt2) &
         is.na(ucv4mfqt1) & is.na(ucv4mfqt2) &
         is.na(zmhmfq_t81) & is.na(zmhmfq_t82) )
  ) 
dim(depSx)
head(depSx)


## Twin corr
depSxTw <- depSx |> 
  distinct(randomfamid, .keep_all = T)

cor.test(depSxTw$u2cmfqt1, depSxTw$u2cmfqt2)

depSxMZ <- depSxTw |> 
  filter(zygos==1)
depSxDZ <- depSxTw |> 
  filter(zygos==2)

## MZ
cor.test(depSxMZ$u2cmfqt1, depSxMZ$u2cmfqt2)
## DZ
cor.test(depSxDZ$u2cmfqt1, depSxDZ$u2cmfqt2)

##### TEDS26
## MZ
cor.test(depSxMZ$zmhmfq_t81, depSxMZ$zmhmfq_t82)
## DZ
cor.test(depSxDZ$zmhmfq_t81, depSxDZ$zmhmfq_t82)


#### Pivot Long =============================

names(depSx)

depSxWide <- depSx |> 
  select(randomfamid,zygos,randomtwinid,
         ends_with("1")) |> 
  arrange(randomfamid,randomtwinid)
head(depSxWide)  

colnames(depSxWide) <- gsub("1","",colnames(depSxWide))
colnames(depSxWide) <- gsub("u2c","T1_",colnames(depSxWide))
colnames(depSxWide) <- gsub("ucv2","T3_",colnames(depSxWide))
colnames(depSxWide) <- gsub("ucv3","T4_",colnames(depSxWide))
colnames(depSxWide) <- gsub("ucv4","T5_",colnames(depSxWide))
colnames(depSxWide) <- gsub("zmh","T6_",colnames(depSxWide))
colnames(depSxWide) <- gsub("_t8","t",colnames(depSxWide))
colnames(depSxWide) <- gsub("ucv","T2_",colnames(depSxWide))

head(depSxWide)
dim(depSxWide)

depSxLong <- depSxWide %>%
  pivot_longer(
    cols = T1_age:T6_mfqt,
    names_to = c("wave", ".value"),
    names_pattern = "(.*)_(.*)"
  ) |> 
  rename(depSx = mfqt) |> 
  mutate(ID = paste(randomtwinid, wave, sep = "_")) |> 
  as.data.frame()

rownames(depSxLong) <- depSxLong$ID

head(depSxLong, 12)

jpeg(paste0(plotDir,"RawDepSx_Hist.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
hist(depSxLong$depSx,
     xlab = "Raw Depression Sx",
     main = "Depression Sx")
dev.off()


#### Descriptives =============================

## age
depSxLong |> 
  filter(!is.na(depSx)) |> 
  group_by(wave) |> 
  summarise(psych::describe(age))

## sex (N)
with(depSxLong[!is.na(depSxLong$depSx),], table(wave, sex))

## sex (prop)
with(depSxLong[!is.na(depSxLong$depSx),], round(prop.table(table(wave, sex), margin = 1), 3))


#### Residualize and Transform ============================

## Rank-based inverse normal transformation
## Both steps are performed in the long-format data, i.e.,
## across all time points and all individuals.
## So, after rank-normal transformation, the mean of **all** observations is zero
## and the variance of **all** observations is one.
## Pivoting back to wide will allow examining mean and variance differences over time.


lmDepSx <- lm(depSx ~ age + sex, data = depSxLong)
summary(lmDepSx)

regDepSx <- as.data.frame(lmDepSx$residuals)
colnames(regDepSx) <- "DepSxRes"
regDepSx$ID <- rownames(regDepSx)
head(regDepSx)
psych::describe(regDepSx$DepSxRes)
hist(regDepSx$DepSxRes)

jpeg(paste0(plotDir,"ResidDepSx_Hist.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
hist(regDepSx$DepSxRes,  xlab = "Residualized Depression Sx",
     main = "Residualized Depression Sx")
dev.off()

## Rank Normal
DepSxResRN <- RNOmni::RankNorm(regDepSx$DepSxRes)
head(DepSxResRN)
nrow(regDepSx)
length(DepSxResRN)
psych::describe(DepSxResRN)
hist(DepSxResRN)

regDepSx$DepSxResRN <- DepSxResRN

jpeg(paste0(plotDir,"Raw_v_RankNormal_ResidDepSx.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
regDepSx |> 
  ggplot(aes(DepSxRes, DepSxResRN)) +
  geom_point(shape = 1,
             color = "indianred3") +
  theme_bw(14) +
  labs(x = "Residualized Depression Sx",
       y = "Rank-Normalized Residualized Depression Sx")
dev.off()

jpeg(paste0(plotDir,"RankNormalResidDepSx_Hist.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
hist(regDepSx$DepSxResRN, xlab = "Rank-Normalized Residualized Depression Sx",
     main = "Residualized Depression Sx\nRank-Based Inverse-Normal Transformed")
dev.off()


##Save
head(regDepSx)

depSxLong <- depSxLong |> 
  left_join(regDepSx, by = "ID")

str(depSxLong)

## Recode Sex and Zygosity
depSxLong <- depSxLong |> 
  mutate(
    sex = case_when(
      sex == 0 ~ "Female",
      sex == 1 ~ "Male"
    ),
    zygos = case_when(
      zygos == 1 ~ "MZ",
      zygos == 2 ~ "DZ"
    )
  )

str(depSxLong)


depSxLong |> 
  filter(wave %in% c("T1","T6")) |> 
  mutate(wave2 = forcats::fct_recode(wave,
                                     TEDS21 = "T1",
                                     TEDS26 = "T6")) |> 
  ggplot(aes(DepSxRes, fill = wave2)) +
  geom_density() + 
  scale_fill_viridis_d(alpha = 0.5, direction = -1) + 
  theme_bw(12) +
  labs(x = "Residualized Depressive Symptoms",
       fill = NULL) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85,0.9),
        axis.title.x = element_text(face = "bold"))
ggsave(paste0(plotDir,"DepSxRes_TEDS21_v_26.jpeg"), 
       width = 6, height = 6, units = "in", dpi = 300)


depSxLong |> 
  filter(wave %in% c("T1","T6")) |> 
  mutate(wave2 = forcats::fct_recode(wave,
                                     TEDS21 = "T1",
                                     TEDS26 = "T6")) |> 
  ggplot(aes(DepSxResRN, fill = wave2)) +
  geom_density() + 
  scale_fill_viridis_d(alpha = 0.5, direction = -1) + 
  theme_bw(12) +
  labs(x = "Rank-Normalized Residualized Depressive Symptoms",
       fill = NULL) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85,0.9),
        axis.title.x = element_text(face = "bold"))
ggsave(paste0(plotDir,"DepSxResRN_TEDS21_v_26.jpeg"), 
       width = 6, height = 6, units = "in", dpi = 300)



#### Pivot Back to Wide =============================

depSxWideInd <- depSxLong |> 
  pivot_wider(id_cols = c("randomfamid","zygos","randomtwinid","sex"), 
              names_from = wave,
              names_glue = "{.value}_{wave}",
              values_from = c(age, depSx, DepSxRes, DepSxResRN)) |> 
  arrange(randomfamid, sex) |> 
  group_by(randomfamid) |> 
  mutate(twinID = row_number()) |> 
  ungroup()

dim(depSxWideInd)
glimpse(depSxWideInd)

depSxWideFam <- depSxWideInd |> 
  pivot_wider(id_cols = c("randomfamid","zygos"),
              names_from = twinID,
              names_glue = "{.value}.Tw{twinID}",
              values_from = randomtwinid:DepSxResRN_T6)
dim(depSxWideFam)
glimpse(depSxWideFam)



#### Correlations ==================================

#### MZ
depSxMZ <- depSxWideFam |> 
  filter(zygos=="MZ") |> 
  select(contains("DepSxResRN")) |> 
  relocate(ends_with("Tw1"))
colnames(depSxMZ) <- gsub("DepSxResRN_","",colnames(depSxMZ))
head(depSxMZ)

depSxMZ |> 
  ggpairs(color = "lightgrey", shape = 1) +
  labs(title = "MZ Twin and Auto-Correlations: Depression Sx")
ggsave(paste0(plotDir,"DepSx_MZcorr_acrossTime.jpeg"), 
       width = 12, height = 9, units = "in", dpi = 300)


#### DZ
depSxDZ <- depSxWideFam |> 
  filter(zygos=="DZ") |> 
  select(contains("DepSxResRN")) |> 
  relocate(ends_with("Tw1"))
colnames(depSxDZ) <- gsub("DepSxResRN_","",colnames(depSxDZ))
head(depSxDZ)

depSxDZ |> 
  ggpairs(color = "lightgrey", shape = 1) +
  labs(title = "DZ Twin and Auto-Correlations: Depression Sx")
ggsave(paste0(plotDir,"DepSx_DZcorr_acrossTime.jpeg"), 
       width = 12, height = 9, units = "in", dpi = 300)


#### Plot Distribution Over Time ==============

glimpse(depSxWideInd)

depSxLong |> 
  group_by(wave) |> 
  summarise(psych::describe(depSx))

depSxLong |> 
  group_by(wave) |> 
  summarise(psych::describe(DepSxRes))

depSxLong |> 
  group_by(wave) |> 
  summarise(psych::describe(DepSxResRN))

depSxLong |> 
  mutate(study = ifelse(wave %in% c("T1","T6"),
                        "TEDS","COVID")) |> 
  mutate(wave2 = forcats::fct_recode(wave,
                                     "TEDS-21" = "T1",
                                     "COVID-1" = "T2",
                                     "COVID-2" = "T3",
                                     "COVID-3" = "T4",
                                     "COVID-4" = "T5",
                                     "TEDS-26" = "T6")) |> 
  ggplot(aes(wave2, depSx)) +
  geom_violin(aes(fill = study),
              scale = "count",
              # trim = F,
              alpha = 0.5) +
  geom_boxplot(notch = T, width = 0.15, outlier.alpha = 0.5) +
  scale_fill_viridis_d(alpha = 0.5, direction = -1) +
  labs(x = NULL,
       y = "Depressive Symptoms\n(Short Moods & Feelings Questionnaire)") +
  theme_bw(12) + 
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12))
ggsave(paste0(plotDir,"DepSx_acrossTime.jpeg"), 
       width = 9, height = 6, units = "in", dpi = 300)

depSxLong |> 
  mutate(study = ifelse(wave %in% c("T1","T6"),
                        "TEDS","COVID")) |> 
  mutate(wave2 = forcats::fct_recode(wave,
                                     "TEDS-21" = "T1",
                                     "COVID-1" = "T2",
                                     "COVID-2" = "T3",
                                     "COVID-3" = "T4",
                                     "COVID-4" = "T5",
                                     "TEDS-26" = "T6")) |> 
  ggplot(aes(wave2, DepSxRes)) +
  geom_violin(aes(fill = study),
              scale = "count",
              alpha = 0.5) +
  geom_boxplot(notch = T, width = 0.15, outlier.alpha = 0.5) +
  scale_fill_viridis_d(alpha = 0.5, direction = -1) +
  labs(x = NULL,
       y = "Depressive Symptoms\n(Residualized for Age and Sex)") +
  theme_bw(12) + 
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12))
ggsave(paste0(plotDir,"DepSxRes_acrossTime.jpeg"), 
       width = 9, height = 6, units = "in", dpi = 300)

depSxLong |> 
  mutate(study = ifelse(wave %in% c("T1","T6"),
                        "TEDS","COVID")) |> 
  mutate(wave2 = forcats::fct_recode(wave,
                                     "TEDS-21" = "T1",
                                     "COVID-1" = "T2",
                                     "COVID-2" = "T3",
                                     "COVID-3" = "T4",
                                     "COVID-4" = "T5",
                                     "TEDS-26" = "T6")) |> 
  ggplot(aes(wave2, DepSxResRN)) +
  geom_violin(aes(fill = study),
              scale = "count") +
  geom_boxplot(notch = T, width = 0.15, outlier.alpha = 0.5, outlier.color = "grey10") +
  scale_fill_viridis_d(alpha = 0.5, direction = -1) +
  labs(x = NULL,
       y = "Depressive Symptoms\n(Rank-Normalized Residuals)") +
  theme_bw(12) + 
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12))
ggsave(paste0(plotDir,"DepSxResRN_acrossTime.jpeg"), 
       width = 9, height = 6, units = "in", dpi = 300)



#### Zygosity + Sex Variable ==============

depSxWideFam <- depSxWideFam |> 
  mutate(zygSex = case_when(
    (zygos == "MZ" & sex.Tw1 == "Female" & sex.Tw2  == "Female") ~ "MZF",
    (zygos == "MZ" & sex.Tw1 == "Male" & sex.Tw2  == "Male") ~ "MZM",
    (zygos == "DZ" & sex.Tw1 == "Female" & sex.Tw2  == "Female") ~ "DZF",
    (zygos == "DZ" & sex.Tw1 == "Male" & sex.Tw2  == "Male") ~ "DZM",
    (zygos == "DZ" & sex.Tw1 == "Female" & sex.Tw2  == "Male") ~ "DZFM"
  ))

depSxWideFam |> 
  count(zygos,sex.Tw1,sex.Tw2,zygSex)
glimpse(depSxWideFam)


#### Save =======================

write.table(depSxWideInd, paste0(dataDir,"TEDS_DepSx_byInd_Jun2024.tsv"), 
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(depSxWideFam, paste0(dataDir,"TEDS_DepSx_byFam_Jun2024.tsv"), 
            sep = "\t", quote = F, row.names = F, col.names = T)

str(read.table(paste0(dataDir,"TEDS_DepSx_byInd_Jun2024.tsv"), header = T))
str(read.table(paste0(dataDir,"TEDS_DepSx_byFam_Jun2024.tsv"), header = T))




## Cigarettes per Day ========================

#### SmkInt ===========================

## TEDS 21
table(tedsIncl$u2csmok011, useNA = "ifany")

## COVID 1
table(tedsIncl$ucv1smok11, useNA = "ifany")

## COVID 2
table(tedsIncl$ucv2smok11, useNA = "ifany")

## COVID 3
table(tedsIncl$ucv3smok11, useNA = "ifany")

## COVID 4
table(tedsIncl$ucv4smok11, useNA = "ifany")

## TEDS 26
table(tedsIncl$zmhsmo11, useNA = "ifany")


#### Current Smoking =====================

## Not assessed directly in TEDS 21
## Please answer questions 4 to 9 if you currently smoke. (FTND scale)
## Those with NA on all 6 items --> not currently smoking

tedsIncl <- tedsIncl |> 
  mutate(currSmkT211 = ifelse(
    ## initiated smoking
    u2csmok011 == 1 &
      ## but responded to at least one of the FTND items
      !( is.na(u2csmok041) & is.na(u2csmok051) & is.na(u2csmok061) & 
           is.na(u2csmok071) & is.na(u2csmok081) & is.na(u2csmok091) ),
    TRUE,
    FALSE
  ),
  currSmkT212 = ifelse(
    ## initiated smoking
    u2csmok012 == 1 &
      ## but responded to at least one of the FTND items
      !( is.na(u2csmok042) & is.na(u2csmok052) & is.na(u2csmok062) & 
           is.na(u2csmok072) & is.na(u2csmok082) & is.na(u2csmok092) ),
    TRUE,
    FALSE
  )
  )

table(tedsIncl$u2csmok011, tedsIncl$currSmkT211, useNA = "ifany")

tedsIncl |> 
  count(u2csmok011,currSmkT211,u2csmok071)


## In COVID Study, those with SmkInt were asked if they smoked in the past month
## Those with past-month smoking were asked Cigs per day
## So, currSmk == 0 can be recoded as CigDay == 0
table(tedsIncl$ucv1smok21, tedsIncl$ucv1smok31, useNA = "ifany")

## In TEDS 26, those with SmkInt were asked if they were currently smoking
table(tedsIncl$zmhsmocig11, tedsIncl$zmhsmocig51, useNA = "ifany")


#### Cigs per day =========================

## re-code to 0 if SmkInt == 1 and CurrSmk == FALSE
tedsIncl <- tedsIncl |> 
  mutate(
    CigDayT211 = ifelse(u2csmok011 == 1 & currSmkT211 == FALSE,
                        0, u2csmok071),
    CigDayT212 = ifelse(u2csmok012 == 1 & currSmkT212 == FALSE,
                        0, u2csmok072),
    CigDayCV11 = ifelse(ucv1smok21 == 0, 0,
                        ucv1smok31),
    CigDayCV12 = ifelse(ucv1smok22 == 0, 0,
                        ucv1smok32),
    CigDayCV21 = ifelse(ucv2smok21 == 0, 0,
                        ucv2smok31),
    CigDayCV22 = ifelse(ucv2smok22 == 0, 0,
                        ucv2smok32),
    CigDayCV31 = ifelse(ucv3smok21 == 0, 0,
                        ucv3smok31),
    CigDayCV32 = ifelse(ucv3smok22 == 0, 0,
                        ucv3smok32),
    CigDayCV41 = ifelse(ucv4smok21 == 0, 0,
                        ucv4smok31),
    CigDayCV42 = ifelse(ucv4smok22 == 0, 0,
                        ucv4smok32),
    CigDayT261 = ifelse(zmhsmocig11 == 0, 0,
                        zmhsmocig51),
    CigDayT262 = ifelse(zmhsmocig12 == 0, 0,
                        zmhsmocig52)
  )

table(tedsIncl$CigDayT211)

table(tedsIncl$CigDayCV11)

table(tedsIncl$CigDayCV21)

table(tedsIncl$CigDayCV31)

table(tedsIncl$CigDayCV41)

table(tedsIncl$CigDayT261)

## 1 = <=10 cigs
## 2 = 11-20 cigs
## 3 = >20 cigs


#### Transformation ===========================

CigDay <- tedsIncl |> 
  select(randomfamid,randomtwinid,sex1,sex2,zygos,
         contains("age"),contains("CigDay")) |> 
  filter( 
    !( is.na(CigDayT211) & is.na(CigDayT212) &
         is.na(CigDayCV11) & is.na(CigDayCV12) &
         is.na(CigDayCV21) & is.na(CigDayCV22) &
         is.na(CigDayCV31) & is.na(CigDayCV32) &
         is.na(CigDayCV41) & is.na(CigDayCV42) &
         is.na(CigDayT261) & is.na(CigDayT262) )
  ) 
dim(CigDay)
head(CigDay, 10)

## Twin corr (ignoring repeated measures)
CigDayTw <- CigDay |> 
  distinct(randomfamid, .keep_all = T)

polycor::polychor(factor(CigDayTw$CigDayT211, ordered = T),
                  factor(CigDayTw$CigDayT212, ordered = T))

CigDayMZ <- CigDayTw |> 
  filter(zygos==1)
CigDayDZ <- CigDayTw |> 
  filter(zygos==2)

## MZ
polycor::polychor(factor(CigDayMZ$CigDayT211, ordered = T),
                  factor(CigDayMZ$CigDayT212, ordered = T))

## DZ
polycor::polychor(factor(CigDayDZ$CigDayT211, ordered = T),
                  factor(CigDayDZ$CigDayT212, ordered = T))



#### Pivot Long =============================

names(CigDay)

CigDayWide <- CigDay |> 
  select(randomfamid,zygos,randomtwinid,
         ends_with("1")) |> 
  arrange(randomfamid,randomtwinid)
head(CigDayWide)  

colnames(CigDayWide) <- gsub("1","",colnames(CigDayWide))
colnames(CigDayWide) <- gsub("u2c","T1_",colnames(CigDayWide))
colnames(CigDayWide) <- gsub("ucv2","T3_",colnames(CigDayWide))
colnames(CigDayWide) <- gsub("ucv3","T4_",colnames(CigDayWide))
colnames(CigDayWide) <- gsub("ucv4","T5_",colnames(CigDayWide))
colnames(CigDayWide) <- gsub("zmh","T6_",colnames(CigDayWide))
colnames(CigDayWide) <- gsub("ucv","T2_",colnames(CigDayWide))
CigDayWide <- CigDayWide |> 
  rename(
    T1_CigDay = CigDayT2,
    T2_CigDay = CigDayCV,
    T3_CigDay = CigDayCV2,
    T4_CigDay = CigDayCV3,
    T5_CigDay = CigDayCV4,
    T6_CigDay = CigDayT26
  )

head(CigDayWide)
dim(CigDayWide)

CigDayLong <- CigDayWide |> 
  pivot_longer(
    cols = T1_age:T6_CigDay,
    names_to = c("wave", ".value"),
    names_pattern = "(.*)_(.*)"
  ) |> 
  mutate(ID = paste(randomtwinid, wave, sep = "_")) |> 
  as.data.frame()

rownames(CigDayLong) <- CigDayLong$ID

head(CigDayLong, 24) |> 
  tail(12)

jpeg(paste0(plotDir,"RawCigDay_Hist.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
CigDayLong |> 
  ggplot(aes(CigDay)) +
  geom_bar(color = "grey10", fill = "lightgray") +
  theme_bw(14) +
  labs(title = "Cigarettes per Day, TEDS 21")
dev.off()


## Merge CigDay==2 and CigDay==3 (too few CigDay==3 for fitting data to ordinal variables)

CigDayLong <- CigDayLong |> 
  mutate(CigDay3L = ifelse(CigDay==3, 2, CigDay))

CigDayLong |> 
  count(CigDay, CigDay3L)

jpeg(paste0(plotDir,"RawCigDay3L_Hist.jpeg"), 
     width = 6, height = 6, units = "in", res = 300)
CigDayLong |> 
  ggplot(aes(CigDay3L)) +
  geom_bar(color = "grey10", fill = "lightgray") +
  theme_bw(14) +
  labs(title = "Cigarettes per Day, TEDS 21")
dev.off()



#### Descriptives =============================

CigDayLong |> 
  filter(!is.na(CigDay3L)) |> 
  group_by(wave) |> 
  summarise(psych::describe(age))

with(CigDayLong[!is.na(CigDayLong$CigDay3L),], table(wave, sex))

with(CigDayLong[!is.na(CigDayLong$CigDay3L),], round(prop.table(table(wave, sex), margin = 1), 3))

with(CigDayLong[!is.na(CigDayLong$CigDay3L),], table(sex, CigDay3L))

with(CigDayLong[!is.na(CigDayLong$CigDay3L),], round(prop.table(table(sex, CigDay3L), margin = 1), 3))

table(CigDayLong$wave, CigDayLong$CigDay3L)

round(prop.table(table(CigDayLong$wave, CigDayLong$CigDay3L), 1), 2)

CigDayTab <- matrix(table(CigDayLong$wave, CigDayLong$CigDay3L), nrow = 6, 
                    dimnames = list(paste0("T",1:6), 0:2))
rowSums(CigDayTab)

round(prop.table(table(CigDayLong$wave, CigDayLong$CigDay3L), margin = 1), 3)


#### Zygosity + Sex Variable ==============

CigDayWideFam <- CigDayWideFam |> 
  mutate(zygSex = case_when(
    (zygos == "MZ" & sex.Tw1 == "Female" & sex.Tw2  == "Female") ~ "MZF",
    (zygos == "MZ" & sex.Tw1 == "Male" & sex.Tw2  == "Male") ~ "MZM",
    (zygos == "DZ" & sex.Tw1 == "Female" & sex.Tw2  == "Female") ~ "DZF",
    (zygos == "DZ" & sex.Tw1 == "Male" & sex.Tw2  == "Male") ~ "DZM",
    (zygos == "DZ" & sex.Tw1 == "Female" & sex.Tw2  == "Male") ~ "DZFM"
  ))

CigDayWideFam |> 
  count(zygos,sex.Tw1,sex.Tw2,zygSex)
glimpse(CigDayWideFam)


#### Save =======================

write.table(CigDayWideInd, paste0(dataDir,"TEDS_CigDay_byInd_Jun2024.tsv"), 
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(CigDayWideFam, paste0(dataDir,"TEDS_CigDay_byFam_Jun2024.tsv"), 
            sep = "\t", quote = F, row.names = F, col.names = T)

str(read.table(paste0(dataDir,"TEDS_CigDay_byInd_Jun2024.tsv"), header = T))
str(read.table(paste0(dataDir,"TEDS_CigDay_byFam_Jun2024.tsv"), header = T))


# END. =========================================================================