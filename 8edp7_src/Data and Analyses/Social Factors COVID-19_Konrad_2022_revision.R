
---
title: "Konrad et al. Social Factors during COVID-19 pandemic"
author: "Annika Konrad"
---
# This script is an updated verion of the original script, changes were made only during the revision process
  
# install.packages(c("compareGroups", "psych", "lme4", "emmeans", 
#                    "robustlmm", "sjPlot", "performance", "ggeffects", 
#                    "ggplot2", "ggnewscale", "see", "dplyr", "rstatix", "coin", "ggpubr", "effectsize"))

library(compareGroups)
library(psych)
library(lme4)
library(emmeans) # ! Version 1.5.5-1
library(robustlmm) # ! > Version 2.4
library(sjPlot)
library(performance)
library(ggeffects)
library(ggplot2)
library(ggnewscale)
library(see)
library(dplyr)
library(rstatix)
library(coin)
library(ggpubr)
library(effectsize)

#### Load Dataframe and setwd()
# setwd("")
data_long <- read.csv2("data_Konrad_etal.csv")

####################################################################################
####################################################################################
### 1. Data Preparation
####################################################################################
# check string issues
data_long$id <- as.numeric(data_long$id)
data_long$id_measurements <- as.numeric(data_long$id_measurements)
data_long$age <- as.numeric(data_long$age)
data_long[,24:99] <- lapply(data_long[,24:99], as.numeric) # turn integers to numeric
data_long[, 5:16] <- lapply(data_long[,5:16], as.factor) 

# combine treatment setting
data_long$treatment_setting2 <- factor(ifelse(data_long$treatment_setting == "A1: No Psychotherapy", "No Psychotherapy",
                                              ifelse(data_long$treatment_setting == "A2: Subsequent Treatment", "Current Psychotherapy", "Other")))
summary(data_long$treatment_setting2)

# calculate combined score for contact and digital items if necessary
data_long$contact = (data_long$contact_family + data_long$contact_friend)/2
data_long$digital = (data_long$digital_family + data_long$digital_friend)/2

# subset for baseline survey
baseline <- subset(data_long, id_measurements == 1)

# get normality
normality <- baseline %>%
  group_by(diagnosis) %>%
  shapiro_test(age, bsi_gsi, tics_socialisolation, 
               rsa_socialresources, bes_disconnection, household, 
               digital, contact)

tab_df(normality, digits = 3, file = "shapiro.html")

# check densityplots
ggdensity(data = baseline, x = "age", fill = "diagnosis")
ggdensity(data = baseline, x = "bsi_gsi", fill = "diagnosis")
ggdensity(data = baseline, x = "tics_socialisolation", fill = "diagnosis")
ggdensity(data = baseline, x = "rsa_socialresources", fill = "diagnosis")
ggdensity(data = baseline, x = "bes_disconnection", fill = "diagnosis")
ggdensity(data = baseline, x = "household", fill = "diagnosis")
ggdensity(data = baseline, x = "digital", fill = "diagnosis")
ggdensity(data = baseline, x = "contact", fill = "diagnosis")

# get descriptive reports
descr_diag <- compareGroups(diagnosis ~ age + gender + bsi_gsi + tics_socialisolation +
                              rsa_socialresources + bes_disconnection + 
                              household + contact + 
                              digital + work_homeoffice + work_lossofwork + 
                              work_normal + work_morehours + work_shortwork + 
                              diagnosis_category + stat_treatment + 
                              treatment_setting2, method = c(age = NA, household = NA, 
                                                             bsi_gsi = NA, tics_socialisolation = NA,
                                                             rsa_socialresources = NA, bes_disconnection = NA, 
                                                             contact = NA, digital = NA), 
                            data = baseline)

descr_diag
createTable(descr_diag)

# wilcox tests
wilcox_test(age ~ diagnosis, data = baseline, paired = FALSE)
wilcox_test(bsi_gsi ~ diagnosis, data = baseline, paired = FALSE)
wilcox_test(tics_socialisolation ~ diagnosis, data = baseline, paired = FALSE)
wilcox_test(rsa_socialresources ~ diagnosis, data = baseline, paired = FALSE)
wilcox_test(bes_disconnection ~ diagnosis, data = baseline, paired = FALSE)
wilcox_test(household ~ diagnosis, data = baseline, paired = FALSE)
wilcox_test(contact ~ diagnosis, data = baseline, paired = FALSE)
wilcox_test(digital ~ diagnosis, data = baseline, paired = FALSE)

w1 <- wilcox_effsize(age ~ diagnosis, ci = T, data = baseline)
w2 <- wilcox_effsize(bsi_gsi ~ diagnosis, ci = T, data = baseline)
w3 <- wilcox_effsize(tics_socialisolation ~ diagnosis, ci = T, data = baseline)
w4 <- wilcox_effsize(rsa_socialresources ~ diagnosis, ci = T, data = baseline)
w5 <- wilcox_effsize(bes_disconnection ~ diagnosis, ci = T, data = baseline)
w6 <- wilcox_effsize(household ~ diagnosis, ci = T, data = baseline)
w7 <- wilcox_effsize(contact ~ diagnosis, ci = T, data = baseline)
w8 <- wilcox_effsize(digital ~ diagnosis, ci = T, data = baseline)

r <- rbind(w1, w2, w3, w4, w5, w6, w7, w8)
tab_df(r, file = "effectsizes_groups.html")
rm(w1, w2, w3, w4, w5, w6, w7, w8)

# get chisquare-tests for 2x3 tables
chisq.test(baseline$gender, baseline$diagnosis, correct = T) # no yates correction is done

# get Yates Correction since some cells are < 5, and R function does only give yates correction for 2x2 tables
observed <- c(149, 45, 2, 380, 160, 5) # observed cells (female, male, non-binary per group)
expected <- c(chisq.test(baseline$gender, baseline$diagnosis)$expected)  # extract expected values
yates <- data.frame(observed, expected)
yates$o_minus_e <- (yates$observed - yates$expected)
yates$corr <- (abs(yates$o_minus_e)-0.5)
yates$sq <- ((yates$corr)^2)
yates$sq_div_e <- (yates$sq/yates$expected)
sum(yates$sq_div_e) # = Chi-square value: 
pchisq(sum(yates$sq_div_e), df=2, lower.tail=FALSE) # get p value

chisq_to_cramers_v(sum(yates$sq_div_e), 741, nrow = 2, ncol = 3, alternative = "two.sided")

# alternative: fisher.test(baseline$gender, baseline$diagnosis)
chisq.test(baseline$treatment_setting2, baseline$diagnosis)
effectsize::cramers_v(baseline$treatment_setting2, baseline$diagnosis, alternative="two.sided")

# get chisquare-tests for 2x2 tables and effectsizes
chisq.test(baseline$work_homeoffice, baseline$diagnosis, correct = F)
effectsize::phi(baseline$work_homeoffice, baseline$diagnosis, alternative="two.sided")
chisq.test(baseline$work_lossofwork, baseline$diagnosis, correct = F)
effectsize::phi(baseline$work_lossofwork, baseline$diagnosis, alternative="two.sided")
chisq.test(baseline$work_morehours, baseline$diagnosis, correct = F)
effectsize::phi(baseline$work_morehours, baseline$diagnosis, alternative="two.sided")
chisq.test(baseline$stat_treatment, baseline$diagnosis, correct = F)
effectsize::phi(baseline$stat_treatment, baseline$diagnosis, alternative="two.sided")

# Note in the no diagnosis group 10 individuals stated mental health issues (see variable diagnosis_category)
# Since they never received a official diagnosis, we conducted a sensitivity analysis (see in script below)

###################################
### Check Cronbach's alpha at T0
###################################
### Basic Empathy Scale: Subscale Empathic Disconnection
alpha_bes_disconnection <- psych::alpha(baseline[c("bes_01", "bes_07", "bes_08", "bes_13", "bes_18", "bes_19")])
alpha_bes_disconnection$total

### Triers Chronic Stress inventory: Subscale social isolation
alpha_tics_socialisolation <- psych::alpha(baseline[c("tics_11", "tics_29", "tics_34", "tics_42", "tics_51", "tics_56")])
alpha_tics_socialisolation$total

### Resilience scale for adults: subscale social resources
alpha_rsa_socialresources <- psych::alpha(baseline[c("rsa_05", "rsa_11", "rsa_17", "rsa_23", "rsa_28", "rsa_32", "rsa_33")])
alpha_rsa_socialresources$total

### Brief Symptom Inventory 53: Sumscore
alpha_bsi_gsi <- psych::alpha(baseline[c("bsi_01", "bsi_02", "bsi_03", "bsi_04", "bsi_05", "bsi_06", "bsi_07", 
                                 "bsi_08", "bsi_09", "bsi_10", "bsi_11", "bsi_12", "bsi_13", "bsi_14", 
                                 "bsi_15", "bsi_16", "bsi_17", "bsi_18", "bsi_19", "bsi_20", "bsi_21", 
                                 "bsi_22", "bsi_23", "bsi_24", "bsi_25", "bsi_26", "bsi_27", "bsi_28", 
                                 "bsi_29", "bsi_30", "bsi_31", "bsi_32", "bsi_33", "bsi_34", "bsi_35", 
                                 "bsi_36", "bsi_37", "bsi_38", "bsi_39", "bsi_40", "bsi_41", "bsi_42", 
                                 "bsi_43", "bsi_44", "bsi_45", "bsi_46", "bsi_47", "bsi_48", "bsi_49", 
                                 "bsi_50", "bsi_51", "bsi_52", "bsi_53")])
alpha_bsi_gsi$total

####################################################################################
####################################################################################
### MODEL 1
####################################################################################
# due to NAs in age variable (min-/max problems), subset of data_long with no missings
# otherwise comparison of model0 and model0_1 with anova function is not possible
# due to different sample sizes
data_long.out <- subset(data_long, age != "NA")

## 1. Null model
model0 = lmer(bsi_gsi ~ (1|id), 
                data = data_long.out)
summary(model0)

## 2. Null Model with covariates
model0.1 = lmer(bsi_gsi ~ age + gender + (1|id), 
                data = data_long.out)
summary(model0.1)

# Compare models
anova(model0, model0.1)

# refit model0_1 with df "data_long" to compare with other models
model0.1 = lmer(bsi_gsi ~ age + gender + (1|id), 
                data = data_long, 
                na.action = "na.omit")
summary(model0.1)

## 3. Model 1 without interaction effect
model1 = lmer(bsi_gsi ~ daycount + diagnosis + rsa_socialresources + age + gender + (1|id), 
                   data = data_long, 
                   na.action = "na.omit")
summary(model1)

## 4. Add interaction effect
model1_full = lmer(bsi_gsi ~ daycount*diagnosis*rsa_socialresources + age + gender + (1|id), 
                   data = data_long,
                   na.action = "na.omit")
summary(model1_full)

###################################
### Compare models 
## 1. Regarding relevance of fixed effects
# CAVE: models are refitted with ML instead REML to compare likelihoods
anova(model0.1, model1, model1_full) 

## 2. compare model indices (ICC, conditional and margina R?)
tab_model(model0.1, model1, model1_full) 

###################################
### Check model assumptions
check_heteroscedasticity(model1_full)
check_normality(model1_full)
check_outliers(model1_full)
check_autocorrelation(model1_full)

###################################
### Refit robust model
model1_robust <- rlmer(bsi_gsi ~ daycount*diagnosis*rsa_socialresources + age + gender + (1|id), 
                   data = data_long,
                   na.action = "na.omit")
summary(model1_robust)

## Compare models
tab_model(model1_full, model1_robust, 
          show.se = TRUE, digits = 3, title = "Model Comparison")

## get standardized betas
model1_robust_refit <- standardize(model1_robust, two_sd = T)
tab_model(model1_robust_refit)

## check collinearity
check_collinearity((model1_robust_refit))

###################################
### Plot robust model 1
## get predicted Values for different levels of rsa_socialresources
gg_model1 <- ggpredict(model1_robust, c("daycount[0, 25, 50, 75]", "rsa_socialresources[meansd]", "diagnosis"))

plot(gg_model1) +
  labs(title = "Time x Social Resources x Diagnosis", 
       x = "Time [in days]", y= "Psychological Distress",
       colour = "Social\nResources") + 
  guides(colour = guide_legend(reverse = TRUE)) +
  scale_color_hue(labels = c("4.84 (-SD)", "5.77 (Mean)", "6.7 (+SD)")) + 
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.3)) + ylim(0,1.2)

###################################
### Post tests for each group

## Diagnosis group
data_long.affected <- subset(data_long, diagnosis == "Diagnosis")

model1_affected <- rlmer(bsi_gsi ~ daycount*rsa_socialresources + age + gender + (1|id), 
                data = data_long.affected,
                na.action = "na.omit")
summary(model1_affected)

# get standardized betas
model1_affected_refit <- standardize(model1_affected, two_sd = T)
tab_model(model1_affected, model1_affected_refit)

## No Diagnosis group
data_long.unaffected <- subset(data_long, diagnosis == "No Diagnosis")

model1_unaffected <- rlmer(bsi_gsi ~ daycount*rsa_socialresources + age + gender + (1|id), 
                data = data_long.unaffected,
                na.action = "na.omit")
summary(model1_unaffected)

# get standardized betas
model1_unaffected_refit <- standardize(model1_unaffected, two_sd = T)
tab_model(model1_unaffected, model1_unaffected_refit)

## Compare models
tab_model(model1_affected, model1_unaffected, show.se = TRUE, 
          digits = 3, title = "Model Comparison", file = "model1_sep.html")

tab_model(model1_affected_refit, model1_unaffected_refit,
          digits = 2, title = "Model Comparison", file = "model1_refit_sep.html")

## Plot separate models
gg_m1_affected <- ggpredict(model1_affected, c("daycount", "rsa_socialresources[meansd]"))
plot(gg_m1_affected)

gg_m1_unaffected <- ggpredict(model1_unaffected, c("daycount", "rsa_socialresources[meansd]"))
plot(gg_m1_unaffected)

###################################
### Trend Analysis
## 1. Pairwise comparisons of slopes
emtrends_m1 <- emtrends(model1_robust, pairwise ~ rsa_socialresources*diagnosis, 
                         var = c("daycount"), at = list(rsa_socialresources=c(4.84,5.77,6.7)), adjust = "bonf")

tab_df(data.frame(emtrends_m1$emtrends), digits = 3, file = "emtrends1.html")

## 2. Marginal means
emmeans_m1 <- emmeans(model1_robust, pairwise ~ rsa_socialresources*diagnosis,
        at = list(rsa_socialresources=c(4.84,5.77,6.7)), adjust = "bonf")

tab_df(data.frame(emmeans_m1$contrasts), digits = 3, file = "emmeans1.html")

## 3. With standardized models
emtrends_m1_refit <- emtrends(model1_robust_refit, pairwise ~ rsa_socialresources*diagnosis, 
                        var = c("daycount"), at = list(rsa_socialresources=c(-0.5,0,0.5)), adjust = "bonf")

tab_df(data.frame(emtrends_m1_refit$emtrends), digits = 3, file = "emtrends1_refit.html")


emmeans_m1_refit <- emmeans(model1_robust_refit, pairwise ~ rsa_socialresources*diagnosis,
                      at = list(rsa_socialresources=c(-0.5,0,0.5)), adjust = "bonf")

tab_df(data.frame(emmeans_m1_refit$contrasts), digits = 3, file = "emmeans1_refit.html")

####################################################################################
####################################################################################
#### Model 2
####################################################################################

## 1. Model 2 without interaction effect
model2 = lmer(bsi_gsi ~ daycount + diagnosis + bes_disconnection + age + gender + (1|id), 
                data = data_long, 
                na.action = "na.omit")
summary(model2)

# 2. Add interaction effect
model2_full <- lmer(bsi_gsi ~ daycount*diagnosis*bes_disconnection + age + gender + (1|id), 
                    data = data_long,
                    na.action = "na.omit")
summary(model2_full)

###################################
### Compare models 
## 1. Regarding relevance of fixed effects
# CAVE: models are refitted with ML instead REML to compare likelihoods
anova(model0.1, model2, model2_full)

## 2. compare model indices (ICC, conditional and margina R?)
tab_model(model0.1, model2, model2_full)

###################################
### Check model assumptions
check_heteroscedasticity(model2_full)
check_normality(model2_full)
check_outliers(model2_full)
check_autocorrelation(model2_full)

###################################
### Refit robust model
model2_robust <- rlmer(bsi_gsi ~ daycount*diagnosis*bes_disconnection + age + gender + (1 |id), 
              data = data_long,
              na.action = "na.omit")
summary(model2_robust)

## Compare model fit
tab_model(model2_full, model2_robust, 
          show.se = TRUE, digits = 3, title = "Model Comparison")

## get standardized betas
model2_robust_refit <- standardize(model2_robust, two_sd = T)
tab_model(model2_robust_refit)

## check collinearity again
check_collinearity(model2_robust_refit)

###################################
### Plot Model 2
# get estimates at fixed levels 
gg_model2 <- ggpredict(model2_robust, c("daycount[0, 25, 50, 75]", "bes_disconnection[meansd]", "diagnosis"))

plot(gg_model2) +
  labs(title = "Time x Empathic Disconnection x Diagnosis", 
       x = "Time [in days]", y= "Psychological Distress",
       colour = "Empathic\nDisconnection") + 
  guides(colour = guide_legend(reverse = TRUE)) +
  scale_color_hue(labels = c("1.57 (-SD)", "2.15 (Mean)", "2.73 (+SD)")) + 
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.3)) + ylim(0,1.2)

###################################
### Post tests for each group

## Diagnosis group
model2_affected <- rlmer(bsi_gsi ~ daycount*bes_disconnection + age + gender + (1 |id), 
                data = data_long.affected,
                na.action = "na.omit")
summary(model2_affected)

# get standardized betas
model2_affected_refit <- standardize(model2_affected, two_sd = T)
tab_model(model2_affected, model2_affected_refit)

# No Diagnosis group
model2_unaffected <- rlmer(bsi_gsi ~ daycount*bes_disconnection + age + gender + (1 |id), 
                         data = data_long.unaffected,
                         na.action = "na.omit")
summary(model2_unaffected)

# get standardized betas
model2_unaffected_refit <- standardize(model2_unaffected, two_sd = T)
tab_model(model2_unaffected, model2_unaffected_refit)

## Compare models
tab_model(model2_affected, model2_unaffected, 
          show.se = TRUE, digits = 3, title = "Model Comparison Per Group", file = "model2_sep.html")

tab_model(model2_affected_refit, model2_unaffected_refit,
          digits = 2, title = "Model Comparison Per Group", file = "model2_refit_sep.html")

## Plot separate models
gg_m2_affected <- ggpredict(model2_affected, c("daycount", "bes_disconnection[meansd]"))
plot(gg_m2_affected)

gg_m2_unaffected <- ggpredict(model2_unaffected, c("daycount", "bes_disconnection[meansd]"))
plot(gg_m2_unaffected)

###################################
### Trend Analysis
## 1. Pairwise comparison of slopes
emtrends_m2 <- emtrends(model2_robust, pairwise ~ bes_disconnection*diagnosis, 
                        var = c("daycount"), 
                        at = list(bes_disconnection=c(1.57,2.15,2.73)), 
                        adjust = "bonf")

tab_df(data.frame(emtrends_m2$emtrends), digits = 3, file = "emtrends2.html")

### 2. Marginal means
emmeans_m2 <- emmeans(model2_robust, pairwise ~ bes_disconnection*diagnosis,
                        at = list(bes_disconnection=c(1.57,2.15,2.73)), 
                      adjust = "bonf")

tab_df(data.frame(emmeans_m2$contrasts), digits = 3, file = "emmeans2.html")

## 3. With standardized models
emtrends_m2_refit <- emtrends(model2_robust_refit, pairwise ~ bes_disconnection*diagnosis, 
                        var = c("daycount"), 
                        at = list(bes_disconnection=c(-0.5,0,0.5)), 
                        adjust = "bonf")

tab_df(data.frame(emtrends_m2_refit$emtrends), digits = 3, file = "emtrends2_refit.html")

emmeans_m2_refit <- emmeans(model2_robust_refit, pairwise ~ bes_disconnection*diagnosis,
                      at = list(bes_disconnection=c(-0.5,0,0.5)), 
                      adjust = "bonf")

tab_df(data.frame(emtrends_m2_refit$contrasts), digits = 3, file = "emmeans2_refit.html")

####################################################################################
####################################################################################
### MODEL 3
####################################################################################

## 1. Model 3 with Only covariates
model0.2 <- lmer(bsi_gsi ~ age + gender + 
                      household + digital + contact + (1|id), 
                    data = data_long,
                    na.action = "na.omit")
summary(model0.2)

## 2. Model 3 without interaction effect
model3 <- lmer(bsi_gsi ~ daycount + diagnosis + tics_socialisolation + age + gender + 
                      household + digital + contact + (1|id), 
                    data = data_long,
                    na.action = "na.omit")
summary(model3)

## 3. Add interaction effect
model3_full <- lmer(bsi_gsi ~ daycount*diagnosis*tics_socialisolation + age + gender + 
                 household + digital + contact + (1|id), 
               data = data_long,
               na.action = "na.omit")
summary(model3_full)

###################################
### Compare models 
## 1. Regarding relevance of fixed effects
# CAVE: models are refitted with ML instead REML to compare likelihoods
anova(model0.2, model3, model3_full)

## 2. compare model indices (ICC, conditional and marginal R)
tab_model(model0.2, model3, model3_full)

###################################
### Check model assumptions
check_heteroscedasticity(model3_full)
check_normality(model3_full)
check_outliers(model3_full)
check_autocorrelation(model3_full)

###################################
### Refit robust model
model3_robust <- rlmer(bsi_gsi ~ daycount*diagnosis*tics_socialisolation + 
                         age + gender + household + digital + contact + 
                         (1|id), 
                       data = data_long,
                       na.action = "na.omit")
summary(model3_robust)

## Compare models
tab_model(model3_full, model3_robust, 
          show.se = TRUE, digits = 3, title = "Model Comparison")

## get standardized betas
model3_robust_refit <- standardize(model3_robust, two_sd = T)
tab_model(model3_robust_refit)

## check collinearity again
check_collinearity(model3_robust_refit)

###################################
# Plot robust model 3
gg_model3 <- ggpredict(model3_robust, c("daycount[0, 25, 50, 75]", "tics_socialisolation[meansd]", "diagnosis"))

plot(gg_model3) +
  labs(title = "Time x Social Isolation x Diagnosis", x = "Time [in days]", y= "Psychological Distress",
                     colour = "Social\nIsolation") + 
  guides(colour = guide_legend(reverse = TRUE)) +
  scale_color_hue(labels = c("0.53 (-SD)", "1.52 (Mean)", "2.51 (+SD)")) + 
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.3)) + ylim(0,1.2)

###################################
### Post tests for each group

## Diagnosis group
model3_affected <- rlmer(bsi_gsi ~ daycount*tics_socialisolation + age + gender + household + 
                  digital + contact + (1 |id), 
                data = data_long.affected,
                na.action = "na.omit")
summary(model3_affected)

# get standardized betas
model3_affected_refit <- standardize(model3_affected, two_sd = T)
tab_model(model3_affected, model3_affected_refit)

## No Diagnosis group
model3_unaffected <- rlmer(bsi_gsi ~ daycount*tics_socialisolation + age + gender + household + 
                           digital + contact + (1 |id), 
                         data = data_long.unaffected,
                         na.action = "na.omit")
summary(model3_unaffected)

# get standardized betas
model3_unaffected_refit <- standardize(model3_unaffected, two_sd = T)
tab_model(model3_unaffected, model3_unaffected_refit)

## Compare models
tab_model(model3_affected, model3_unaffected, show.se = TRUE, digits = 3,
          itle = "Model Comparison Per Group", file = "model3_sep.html")

# standardized betas         
tab_model(model3_affected_refit, model3_unaffected_refit, digits = 2, 
          title = "Model Comparison Per Group", file = "model3_refit_sep.html")

## Plot separate models
gg_m3_affected <- ggpredict(model3_affected, c("daycount", "tics_socialisolation[meansd]"))
plot(gg_m3_affected)

gg_m3_unaffected <- ggpredict(model3_unaffected, c("daycount", "tics_socialisolation[meansd]"))
plot(gg_m3_unaffected)

###################################
### Trend Analysis Model 3
## 1. Pairwise comparisons of slopes
emtrends_m3 <- emtrends(model3_robust, pairwise ~ tics_socialisolation*diagnosis, 
                        var = c("daycount"), 
                        at = list(tics_socialisolation=c(0.53,1.52,2.51)), adjust = "bonf")

tab_df(data.frame(emtrends_m3$emtrends), digits = 3, file = "emtrends3.html")

## 2. Marginal means
emmeans_m3 <- emmeans(model3_robust, pairwise ~ tics_socialisolation*diagnosis,
                        at = list(tics_socialisolation=c(0.53,1.52,2.51)), adjust = "bonf")

tab_df(data.frame(emmeans_m3$contrasts), digits = 3, file = "emmeans3.html")

## 3. With standardized models
emtrends_m3_refit <- emtrends(model3_robust_refit, pairwise ~ tics_socialisolation*diagnosis, 
                        var = c("daycount"), 
                        at = list(tics_socialisolation=c(-0.5, 0, 0.5)), adjust = "bonf")

tab_df(data.frame(emtrends_m3_refit$emtrends), digits = 3, file = "emtrends3_refit.html")

emmeans_m3_refit <- emmeans(model3_robust_refit, pairwise ~ tics_socialisolation*diagnosis,
                      at = list(tics_socialisolation=c(-0.5, 0, 0.5)), adjust = "bonf")

tab_df(data.frame(emmeans_m3_refit$contrasts), digits = 3, file = "emmeans3_refit.html")


####################################################################################
# summary all robust models
tab_model(model1_robust, model2_robust, model3_robust,  
          show.stat = TRUE, digits = 3, file = "Allrobustmodels.html")

# summary all standardized models
tab_model(model1_robust_refit, model2_robust_refit, model3_robust_refit,  
          show.stat = TRUE, digits = 2, file = "Allrobustmodels_refit.html")

###################################################################################
# Figure 1
###################################################################################
data_long$Measurements <- ifelse(data_long$id_measurements == 1, "T0", "Follow-Ups")
data_long$Measurements <- factor(data_long$Measurements)
summary(data_long$Measurements)

data_long$`Pre-existing Diagnosis` <- factor(data_long$diagnosis, labels = c("Yes", "No"))

### Plot 1
plot1 <- ggplot(data_long, mapping = aes(daycount, bsi_gsi)) +
  labs(x = "Date", y = "Psychological Distress") +
  ylim(0,4) +
  xlim(-25,100) +
  geom_point2(aes(colour = Measurements), alpha = 0.5) +
  scale_color_manual(values=c("grey40", "orange"), 
                     guide = guide_legend(reverse=TRUE, 
                                          override.aes = list(size=2, stroke=1.5)), 
                     labels = c("Follow-Ups:\n - State: Psychological Distress,\n               Social Isolation,\n               Contact Items", 
                                "First Survey:\n - Trait: Demographics,\n              Social Resources,\n              Empathic Disconnection\n - State: Psychological Distress,\n               Social Isolation,\n               Contact Items\n")) +
  new_scale_color() +
  geom_smooth(aes(colour = `Pre-existing Diagnosis`, linetype = `Pre-existing Diagnosis`), method = "lm") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_color_see() +
  theme_classic()
plot1

### Plot 1.1
plot1.1 <- plot1 +
  scale_x_continuous(limits = c(-25,100), breaks = c(-25, 0 , 25, 50, 75, 100), 
                     labels= c("10th March","3rd April", "28th April", "23th May", "17th June", "14th July")) +
  theme_minimal(base_size = 8)
plot1.1

ggsave(filename = "Fig_1.tiff", 
       plot = last_plot(), width = 15, height = 8, dpi = 300, units = "cm")

###################################################################################
# Figure 2
###################################################################################
gg_model1$facet <- factor(gg_model1$facet, 
                          labels = c("With Pre-existing Diagnosis", "Without Pre-existing Diagnosis"))

### Panel A
P1 <- plot(gg_model1) + 
  labs(title = "Time x Social Resources x Diagnosis", x = "Time [in days]", 
       y= "Psychological Distress",
       colour = "Social\nResources", tag = "(A)") +
  scale_color_hue(labels = c("4.84 (-SD)", "5.77 (Mean)", "6.7 (+SD)")) + 
  scale_y_continuous(breaks = c(0, 0.4, 0.8, 1.2)) + ylim(0,1.2) + 
  theme_minimal(base_size = 8)

### Panel B
gg_model2$facet <- factor(gg_model2$facet, 
                          labels = c("With Pre-existing Diagnosis", "Without Pre-existing Diagnosis"))

P2 <- plot(gg_model2) + 
  labs(title = "Time x Empathic Disconnection x Diagnosis", x = "Time [in days]", 
       y= "Psychological Distress",
       colour = "Empathic\nDisconnection", tag = "(B)") +
  scale_color_hue(labels = c("1.57 (-SD)", "2.15 (Mean)", "2.73 (+SD)")) + 
  scale_y_continuous(breaks = c(0, 0.4, 0.8, 1.2)) + ylim(0,1.2) + 
  theme_minimal(base_size = 8)

### Panel C
gg_model3$facet <- factor(gg_model3$facet, 
                          labels = c("With Pre-existing Diagnosis", "Without Pre-existing Diagnosis"))

P3 <- plot(gg_model3) + 
  labs(title = "Time x Social Isolation x Diagnosis", x = "Time [in days]", 
       y= "Psychological Distress",
       colour = "Social\nIsolation", tag = "(C)") +
  scale_color_hue(labels = c("0.53 (-SD)", "1.52 (Mean)", "2.51 (+SD)")) + 
  scale_y_continuous(breaks = c(0, 0.4, 0.8, 1.2)) + ylim(0,1.2) + 
  theme_minimal(base_size = 8)

# Combine Plots
P_all <- ggarrange(P1, P2, P3, ncol = 3, legend = "bottom") 

# Save
ggsave(filename = "Fig_2.tiff" , 
       plot = last_plot(), width = 25, height = 10, dpi = 300, units = "cm")

####################################################################################
####################################################################################
# SENSITIVITY ANALYSIS
####################################################################################

data_long_noout <- subset(data_long, diagnosis != "No Diagnosis" | stat_treatment != "Yes") # exclude 7 individuals with 
data_long_outl <- subset(data_long, diagnosis == "No Diagnosis" & stat_treatment == "Yes")
View(data_long_outl)


###################################
### Robust model 1
model1_noout <- rlmer(bsi_gsi ~ daycount*diagnosis*rsa_socialresources + age + gender + (1|id), 
                      data = data_long_noout,
                      na.action = "na.omit")
summary(model1_noout)

tab_model(model1_robust, model1_noout, digits = 3, file = "model1_noout.html")

###################################
### Robust model 2
model2_noout <- rlmer(bsi_gsi ~ daycount*diagnosis*bes_disconnection + age + gender + (1 |id), 
                      data = data_long_noout,
                      na.action = "na.omit")
summary(model2_noout)

tab_model(model2_robust, model2_noout, digits = 3, file = "model2_noout.html")

###################################
### Robust model 3
model3_noout <- rlmer(bsi_gsi ~ daycount*diagnosis*tics_socialisolation + age + gender + household + 
                        digital + contact + (1 |id), 
                      data = data_long_noout,
                      na.action = "na.omit")
summary(model3_noout)

tab_model(model3_robust, model3_noout, digits = 3, file = "model3_noout.html")

####################################################################################
# 2ND SENSITIVITY ANALYSIS
####################################################################################
data_long_noout_curtreat <- subset(data_long, 
                                   diagnosis != "No Diagnosis" | treatment_setting2 != "Current Psychotherapy") # exclude 28 individuals with 
data_long_outl_curtreat <- subset(data_long, 
                                  diagnosis == "No Diagnosis" & treatment_setting2 == "Current Psychotherapy")

###################################
### Robust model 1
model1_noout_curtreat <- rlmer(bsi_gsi ~ daycount*diagnosis*rsa_socialresources + age + gender + (1|id), 
                               data = data_long_noout_curtreat,
                               na.action = "na.omit")
summary(model1_noout_curtreat)

tab_model(model1_robust, model1_noout_curtreat, digits = 3, file = "model1_noout_curtreat.html")

###################################
### Robust model 2
model2_noout_curtreat <- rlmer(bsi_gsi ~ daycount*diagnosis*bes_disconnection + age + gender + (1 |id), 
                               data = data_long_noout_curtreat,
                               na.action = "na.omit")
summary(model2_noout_curtreat)

tab_model(model2_robust, model2_noout_curtreat, digits = 3, file = "model2_noout_curtreat.html")

###################################
### Robust model 3
model3_noout_curtreat <- rlmer(bsi_gsi ~ daycount*diagnosis*tics_socialisolation + age + gender + household + 
                                 digital + contact + (1 |id), 
                               data = data_long_noout_curtreat,
                               na.action = "na.omit")
summary(model3_noout_curtreat)
tab_model(model3_robust, model3_noout_curtreat, digits = 3, file = "model3_noout_curtreat.html")

####################################################################################
# 3RD SENSITIVITY ANALYSIS
####################################################################################
data_long_noout_diagcat <- subset(data_long, diagnosis != "No Diagnosis" | is.na(data_long$diagnosis_category)) # exclude 44 obs. 
data_long_outl_diagcat <- subset(data_long, diagnosis == "No Diagnosis" & diagnosis_category != "NA")


###################################
### Robust model 1
model1_noout_diagcat <- rlmer(bsi_gsi ~ daycount*diagnosis*rsa_socialresources + age + gender + (1|id), 
                              data = data_long_noout_diagcat,
                              na.action = "na.omit")
summary(model1_noout_diagcat)

tab_model(model1_robust, model1_noout_diagcat, digits = 3, file = "model1_noout_diagcat.html")

###################################
### Robust model 2
model2_noout_diagcat <- rlmer(bsi_gsi ~ daycount*diagnosis*bes_disconnection + age + gender + (1 |id), 
                              data = data_long_noout_diagcat,
                              na.action = "na.omit")
summary(model2_noout_diagcat)

tab_model(model2_robust, model2_noout_diagcat, digits = 3, file = "model2_noout_diagcat.html")

###################################
### Robust model 3
model3_noout_diagcat <- rlmer(bsi_gsi ~ daycount*diagnosis*tics_socialisolation + age + gender + household + 
                                digital + contact + (1 |id), 
                              data = data_long_noout_diagcat,
                              na.action = "na.omit")
summary(model3_noout_diagcat)
tab_model(model3_robust, model3_noout_diagcat, digits = 3, file = "model3_noout_diagcat.html")
