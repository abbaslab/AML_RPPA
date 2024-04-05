library(ggplot2)
rppa_data <- read.table("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/ASNS/RPPA_Data.csv", sep = "\t")

pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure 1b.pdf")
ggplot(rppa_data, aes(x = ASNS)) + 
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.5, color = "black",fill = "steelblue3") + 
  geom_density(color="red", lwd=1)  +
  xlab("Log2 Normalized ASNS Expression Values") +
  ylab("Density") +
  theme_classic()
dev.off()
mean(rppa_data$ASNS)
median(rppa_data$ASNS)
shapiro.test(rppa_data$ASNS)
pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure S1.pdf")
qqnorm(rppa_data$ASNS)
dev.off()

rppa_data$poi_median <- factor(rppa_data$poi_median, levels = c("LowASNS", "HighASNS"))
rppa_data$OS.Wks.SMK.yr <- rppa_data$OS.Wks.SMK/52
rppa_data$status <- ifelse(rppa_data$vital.status=="D", 0, 1)
library(survival)
library(survminer)
fit <- survfit(Surv(OS.Wks.SMK, status==0) ~ 1, data = rppa_data)
surv_median(fit)