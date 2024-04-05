library(dplyr)
library(ggplot2)
library(ggsurvfit)
library(tidycmprsk)
library(ggpubr)
library(survival)
library(survminer)
library(survivalAnalysis)
library(phenoTest)
rppa_data <- read.table("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Code/Data/RPPA_Data.csv", sep = "\t")
rppa_data$poi_median <- factor(rppa_data$poi_median, levels = c("LowASNS", "HighASNS"))
rppa_data$OS.Wks.SMK.yr <- rppa_data$OS.Wks.SMK/52
## 2-year survival
rppa_data$OS_2year <- ifelse(rppa_data$OS.Wks.SMK.yr <= 2, rppa_data$OS.Wks.SMK.yr, 2.01)
rppa_data$OS_Status <- ifelse(rppa_data$OS_2year <= 2, 1, 0)
rppa_data$OS_Status <- as.numeric(rppa_data$OS_Status)
median_fit_2 <- survfit(Surv(OS_2year, OS_Status) ~ poi_median, data = rppa_data)
cox_median_2 <- coxph(Surv(OS_2year, OS_Status) ~ poi_median, data = rppa_data)
cp <- coxph(Surv(OS_2year, OS_Status) ~ ASNS, data = rppa_data)
summary(cp)

pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure 3a.pdf")
ggsurvplot(median_fit_2, data=rppa_data, pval=T, palette = c("steelblue4", "red3"), size = 1,
           legend.labs = c("Low", "High"),ylab = "Overall Survival Probability", xlab= "Survival Time (years)", 
          legend=c(0.9,0.9), pval.coord = c(0.1,0.05), legend.title="ASNS", title="ASNS Median Survival Curve", conf.int = FALSE,
           risk.table = TRUE) 
dev.off()

pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure 3b.pdf")
smoothCoxph(rppa_data$OS_2year, rppa_data$OS_Status, rppa_data$ASNS, xlab="Log2 Normalized ASNS Expression Levels", ylab="Log Hazard Ratio") 
abline(v=median(rppa_data$ASNS),lty=2, col = "grey")
title(main = "Change in HR with ASNS expression")
text(5, 0.65, "HR = 1.11 ± 0.04 \n Logrank p = 0.02")
text(1, 0.4, "H=410", col="red3")
text(-0.5, 0.4, "L=410", col="steelblue4")
dev.off()
