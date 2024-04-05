library(dplyr)
library(ggplot2)
library(ggsurvfit)
library(tidycmprsk)
library(ggpubr)
library(survival)
library(survminer)
library(survivalAnalysis)
library(phenoTest)
rppa_data <- read.table("RPPA_Data.csv", sep = "\t")
rppa_data$poi_median <- factor(rppa_data$poi_median, levels = c("LowASNS", "HighASNS"))
rppa_data$poi_tertile <- factor(rppa_data$poi_tertile, levels = c("LowASNS", "MediumASNS","HighASNS"))
rppa_data$poi_quartile <- factor(rppa_data$poi_quartile, levels = c("Q1", "Q2", "Q3", "Q4"))
rppa_data$OS.Wks.SMK.yr <- rppa_data$OS.Wks.SMK/52
rppa_data$low <- ifelse(rppa_data$TX1.Venetoclax == "Y" & rppa_data$poi_median == "LowASNS", "Ven:LowASNS", ifelse(rppa_data$TX1.Venetoclax == "N" & rppa_data$poi_median == "LowASNS", "NoVen:LowASNS", NA))
rppa_data$low <- factor(rppa_data$low, levels = c("NoVen:LowASNS", "Ven:LowASNS"))
rppa_data$high <- ifelse(rppa_data$TX1.Venetoclax == "N" & rppa_data$poi_median == "HighASNS", "NoVen:HighASNS", ifelse(rppa_data$TX1.Venetoclax == "Y" & rppa_data$poi_median == "HighASNS", "Ven:HighASNS", NA))
rppa_data$high <- factor(rppa_data$high, levels = c("NoVen:HighASNS", "Ven:HighASNS"))

## 2-year survival
rppa_data$OS_2year <- ifelse(rppa_data$OS.Wks.SMK.yr <= 2, rppa_data$OS.Wks.SMK.yr, 2.01)
rppa_data$OS_Status <- ifelse(rppa_data$OS_2year <= 2, 1, 0)
rppa_data$OS_Status <- as.numeric(rppa_data$OS_Status)
ven <- survfit(Surv(OS_2year, OS_Status) ~ TX1.Venetoclax, data = rppa_data)
cox <- coxph(Surv(OS_2year, OS_Status) ~ TX1.Venetoclax, data = rppa_data)
fit_ven <- survfit(Surv(OS_2year, OS_Status) ~ low, data = rppa_data)
cox_ven <- coxph(Surv(OS_2year, OS_Status) ~ low, data = rppa_data)
fit_ven_high <- survfit(Surv(OS_2year, OS_Status) ~ high, data = rppa_data)
cox_ven_high <- coxph(Surv(OS_2year, OS_Status) ~ high, data = rppa_data)
pdf("Figure 3c.pdf", width = 8)
ggsurvplot(fit_ven_high, data=rppa_data, pval=T,palette = c("red3", "steelblue4"),
           legend.labs = c("No", "Yes"),xlab= "Survival Time (years)",ylab = "Overall Survival Probability",
           legend=c(0.9,0.9), pval.coord = c(0.1,0.05),legend.title="Venetoclax", 
           title="Venetoclax Survival Curve of High ASNS Patients", conf.int = FALSE,
           risk.table = TRUE) 
dev.off()
pdf("Figure S4a.pdf", width = 8)
ggsurvplot(ven, data=rppa_data, pval=T, palette = c("red3", "steelblue4"),
           legend.labs = c("No", "Yes"), xlab= "Survival Time (years)",ylab = "Overall Survival Probability",
           legend=c(0.9,0.9), pval.coord = c(0.1,0.05),legend.title="Venetoclax", title="Venetoclax Survival Curve", conf.int = FALSE,
           risk.table = TRUE) 
dev.off()
pdf("Figure S5a.pdf", width = 8)
ggsurvplot(fit_ven, data=rppa_data, pval=T, palette = c("red3","steelblue4"),           
           legend.labs = c("No", "Yes"), xlab= "Survival Time (years)",ylab = "Overall Survival Probability",
          legend=c(0.9,0.9), pval.coord = c(0.1,0.05),legend.title="Venetoclax", title="Venetoclax Survival Curve of Low ASNS Patients", conf.int = FALSE,
           risk.table = TRUE) 
dev.off()

## Ven multivariate
rppa_data$Cytogenetics.Risk.Group <- factor(rppa_data$Cytogenetics.Risk.Group, levels = c("Favorable", "Intermediate", "Adverse"))
rppa_data$Age.Dx <- ifelse(rppa_data$Age.Dx  > 60, "Greater", "Less")
rppa_data$Age.Dx  <- factor(rppa_data$Age.Dx , levels=c("Less", "Greater"))
rppa_data$ven_asns <- ifelse(rppa_data$TX1.Venetoclax == "N" & rppa_data$poi_median == "LowASNS", "NoVen:LowASNS", ifelse(rppa_data$TX1.Venetoclax == "N" & rppa_data$poi_median == "HighASNS", "NoVen:HighASNS", ifelse(rppa_data$TX1.Venetoclax == "Y" & rppa_data$poi_median == "LowASNS", "Ven:LowASNS", ifelse(rppa_data$TX1.Venetoclax == "Y" & rppa_data$poi_median == "HighASNS", "Ven:HighASNS", NA))))
rppa_data$ven_asns <- factor(rppa_data$ven_asns, levels = c("NoVen:LowASNS", "NoVen:HighASNS", "Ven:LowASNS", "Ven:HighASNS"))
rppa_data$TX1.Venetoclax <- factor(rppa_data$TX1.Venetoclax, levels = c("N", "Y"))
rppa_data$poi_median <- factor(rppa_data$poi_median, levels = c("LowASNS", "HighASNS"))
covariate_names_0 <- c(`Age.Dx:Greater` ="Age > 60", 
                     `poi_median:HighASNS`="High ASNS", 
                     `poi_median:LowASNS`="Low ASNS", 
                     `TX1.Venetoclax:Y`="Venetoclax Treated", 
                     `TX1.Venetoclax:N`="No Venetoclax Treatment", 
                     `Cytogenetics.Risk.Group:Favorable` = "ELN2017 Favorable Cytogenetics",
                     `Cytogenetics.Risk.Group:Intermediate` = "ELN 2017 Intermediate Cytogenetics",
                     `Cytogenetics.Risk.Group:Adverse` = "ELN2017 Adverse Cytogenetics")

covariate_names_2 <- c(`Age.Dx:Greater` ="Age > 60", 
                     `low:NoVen:LowASNS`="No Venetoclax and Low ASNS", 
                     `low:Ven:LowASNS`="Venetoclax and Low ASNS", 
                     `high:NoVen:HighASNS`="No Venetoclax and High ASNS", 
                     `high:Ven:HighASNS`="Venetoclax and High ASNS", 
                     `Cytogenetics.Risk.Group:Favorable` = "ELN2017 Favorable Cytogenetics",
                     `Cytogenetics.Risk.Group:Intermediate` = "ELN 2017 Intermediate Cytogenetics",
                     `Cytogenetics.Risk.Group:Adverse` = "ELN2017 Adverse Cytogenetics")

result3 <- rppa_data %>%
  analyse_multivariate(vars(OS.Wks.SMK, status),
                       covariates = vars(Age.Dx, Cytogenetics.Risk.Group, TX1.Venetoclax, poi_median), covariate_name_dict = covariate_names_0)
result4 <- rppa_data %>%
  analyse_multivariate(vars(OS.Wks.SMK, status),
                       covariates = vars(Age.Dx, Cytogenetics.Risk.Group, high), covariate_name_dict = covariate_names_2)
result5 <- rppa_data %>%
  analyse_multivariate(vars(OS.Wks.SMK, status),
                       covariates = vars(Age.Dx, Cytogenetics.Risk.Group, low), covariate_name_dict = covariate_names_2)

pdf("Figure 3d.pdf")
forest_plot(result4,
            factor_labeller = covariate_names_2,
            endpoint_labeller = c(time="OS_Weeks"),
            orderer = ~order(HR),
            labels_displayed = c("factor"),
            ggtheme = ggplot2::theme_bw(base_size = 8),
            relative_widths = c(1, 1.5, 1),
            HR_x_breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2,3,4,5,7,9))
dev.off()
pdf("Figure S4b.pdf")
forest_plot(result3,
            factor_labeller = covariate_names_0,
            endpoint_labeller = c(time="OS_Weeks"),
            orderer = ~order(HR),
            labels_displayed = c("factor"),
            ggtheme = ggplot2::theme_bw(base_size = 8),
            relative_widths = c(1, 1.5, 1),
            HR_x_breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2,3,4))
dev.off()
pdf("Figure S5b.pdf")
forest_plot(result5,
            factor_labeller = covariate_names_2,
            endpoint_labeller = c(time="OS_Weeks"),
            orderer = ~order(HR),
            labels_displayed = c("factor"),
            ggtheme = ggplot2::theme_bw(base_size = 8),
            relative_widths = c(1, 1.5, 1),
            HR_x_breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2,3,4,5,7,8))
dev.off()

