## Load Packages ---------------------------------------------------------------

library(survival)
library(survminer)
library(tidyverse)
library(patchwork)

## LAKI Data -------------------------------------------------------------------

## Survival Analysis of LAKI mice treated with
## DMSO or UMK57. First, I tried to analyse first
## without considering the sex of the animals.

df <- read.csv("input_data\\survival\\survival_laki.csv", 
               header = TRUE, sep = ";")
head(df, 20)


# Create the survival object:

surv_object <- Surv(time = df$Death, event = rep(1, nrow(df)))  # 1 = All died

# Fit the Cox model:

cox_model <- coxph(surv_object ~ Condition, data = df)

## coxph() fits a Cox regression model with Condition 
## (DMSO vs UMK57) as the predictor. This model estimates 
## the hazard ratio (HR) — how much more (or less) likely a mouse 
## in one group is to die at any given moment compared to the other group. 
## A HR < 1 for UMK57 would suggest a protective effect.

summary(cox_model) ## The p in LogRank test here is not the correct one.

## Create a survival plot:

fit <- survfit(Surv(Death, rep(1, nrow(df))) ~ Condition, data = df)

summary(fit)$table

p1 <- ggsurvplot(fit, data = df, pval = TRUE,
                 surv.median.line = "hv") ## P-value is Log-Rank Test
p1

## survfit() generates Kaplan-Meier survival curves per condition
## ggsurvplot() plots them. The pval = TRUE argument automatically 
## overlays the log-rank test p-value, which tests whether the two survival 
## curves are statistically different.

s <- summary(cox_model)
s

pval <- s$coefficients[,"Pr(>|z|)"] ## Get the Wald Test pval

## Get the correct LogRank test p-value:

res <- survdiff(Surv(Death, rep(1, nrow(df))) ~ Condition, data = df)
res

pval_logrank <- 1 - pchisq(res$chisq, length(res$n) - 1)
pval_logrank

## Get the Wald Test pval again:

pval_wald_test <- s$coefficients[,"Pr(>|z|)"]

HR <- s$coefficients[,"exp(coef)"]

pval_logrank
pval_wald_test

summary_cox_hazards <- data.frame(
  "logrank_pval" = pval_logrank,
  "wald_test_pval" = pval_wald_test,
  "hazard_ratio" = HR
)

summary_cox_hazards

summary_cox_hazards$logrank_pval

signif(summary_cox_hazards$logrank_pval, digits = 2)

p1$plot + coord_cartesian(xlim = c(12, 20)) + scale_x_continuous(breaks=seq(12, 20, 1)) +
  annotate("text", x = 12, y = 0.25, label = paste("Log-Rank Test: P = ", signif(summary_cox_hazards$logrank_pval, digits = 2)), hjust = 0) +
  annotate("text", x = 12, y = 0.2, label = paste("Wald Test (Cox Prop. Hazards): P = ", signif(summary_cox_hazards$wald_test_pval, digits = 2)), hjust = 0) +
  annotate("text", x = 12, y = 0.15, label = paste("Hazard Ratio UMK57 vs DMSO (HR) = ", signif(summary_cox_hazards$hazard_ratio, digits = 2)), hjust = 0) +
  labs(y = "Survival Probability", x = "Time (Weeks)") + 
  theme_bw(base_size = 14)


## LAKI Survival Data ----------------------------------------------------------

## Survival Analysis of LAKI mice treated with
## DMSO or UMK57, adjusted for sex.

## Load survival data:

df <- read.csv("input_data\\survival\\survival_laki.csv", 
               header = TRUE, sep = ";")


## Load information about sex of the animals:

sex_info <- read.csv("input_data\\survival\\animal_sex_info.csv", sep = ";")
length(unique(sex_info$Code))


##Ensures columns in df and sex_info2 are both named "Animal":

sex_info2 <- sex_info %>% distinct() %>% rename(Animal = Code)


##Merge survival and sex data:

df_sex <- merge(df, sex_info2, by.x = "Animal", by.y = "Animal")


## Survival analysis inlcuding sex as covariate:

## coxph() fits a Cox regression model with Condition adjusted for Sez
## (DMSO vs UMK57) as the predictor. This model estimates 
## the hazard ratio (HR) — how much more (or less) likely a mouse 
## in one group is to die at any given moment compared to the other group. 
## A HR < 1 for UMK57 would suggest a protective effect.

cox_model <- coxph(
  Surv(time = Death, event = rep(1, nrow(df_sex))) ## 1 = Dead
  ~ Sex + Condition, data = df_sex) ## Adjusted for Sex

summary_results <- summary(cox_model) 

## Reveals p-value and effect size for each variable (sex or condition)

summary_results

## Confirm if proportional hazards change over time. 
## They shouldn't. Should be p > 0.05

cox.zph(cox_model) 

## No significant change over time, we can proceed.


## LogRank, Cox and Plot -------------------------------------------------------

## Get LogRank test results as comparison:

res <- survdiff(Surv(Death, rep(1, nrow(df_sex))) ~ Condition, data = df_sex)
res
pval_logrank <- 1 - pchisq(res$chisq, length(res$n) - 1)
pval_logrank


pval <- summary_results$coefficients
pval


pval_sex <- pval["SexM",5]
pval_condition <- pval["ConditionUMK57",5]

umk57_HR <- summary_results$coefficients["ConditionUMK57","exp(coef)"]


## Create a plot, first survfit then ggsurvplot:

fit_long_rank <- survfit(Surv(Death, rep(1, nrow(df_sex))) ~ Condition, 
                         data = df_sex)


p1 <- ggsurvplot(fit_long_rank, data = df_sex, pval = T, conf.int = F)
p1

p1$plot + coord_cartesian(xlim = c(12, 20)) + 
  scale_x_continuous(breaks=seq(12, 20, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  annotate(
    "text", 
    x = 12, 
    y = 0.25, 
    label = paste("Log-Rank Test: P = ", 
                  signif(pval_logrank, digits = 1)), hjust = 0) +
  annotate(
    "text", 
    x = 12, 
    y = 0.2, 
    label = paste("Cox Proportional Hazards (Condition): P = ", 
                  signif(pval_condition, digits = 2), "**"), hjust = 0) +
  annotate(
    "text", 
    x = 12, 
    y = 0.15, 
    label = paste("Cox Proportional Hazards (Sex): P = ", 
                  signif(pval_sex, digits = 2)), hjust = 0) +
  annotate(
    "text", 
    x = 12, 
    y = 0.10, 
    label = paste("Hazard Ratio UMK57 vs DMSO (HR) = ", 
                  signif(umk57_HR, digits = 2)), hjust = 0) +
  labs(y = "Survival Probability", x = "Time (Weeks)") + 
  theme_bw(base_size = 14)

ggsave(
  "umk57_survival_cph_plot.pdf",
  path = "output_data\\plots\\main_figure",
  width = 12,
  height = 7,
  unit = "in"
)

writeLines(capture.output(sessionInfo()), "sessionInfo_survival.txt")





