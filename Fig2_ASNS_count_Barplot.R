library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggpubr)
rppa_data <- read.table("RPPA_Data.csv", sep = "\t")
rppa_data$Cytogenetics.Risk.Group <- factor(rppa_data$Cytogenetics.Risk.Group, c("Favorable", "Intermediate", "Adverse"))
rppa_data$FAB <- rppa_data$PseudoFAB
rppa_data$FAB <- gsub("M4e", "M4", rppa_data$FAB)
rppa_data$FAB <- gsub("M5a", "M5", rppa_data$FAB)
rppa_data$FAB <- gsub("M5b", "M5", rppa_data$FAB)
rppa_data$FAB <- gsub("M6a", "M6", rppa_data$FAB)
rppa_data$FAB <- gsub("M6b", "M6", rppa_data$FAB)
rppa_data$FAB <- gsub("Unclassifiable", NA, rppa_data$FAB)
rppa_data$FAB <- factor(rppa_data$FAB, levels = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "RAEB-T"))
rppa_data$Cytogenetics.Cat.Summary <- gsub("COMPLEX", "Complex", rppa_data$Cytogenetics.Cat.Summary)
rppa_data$Cytogenetics.Cat.Summary <- gsub("DIPLOID", "Diploid", rppa_data$Cytogenetics.Cat.Summary)
rppa_data$Cytogenetics.Cat.Summary <- gsub("Misc", "Others", rppa_data$Cytogenetics.Cat.Summary)
rppa_data$Cytogenetics.Cat.Summary <- gsub("trisomy8", "Trisomy8", rppa_data$Cytogenetics.Cat.Summary)
rppa_data$Cytogenetics.Cat.Summary <- factor(rppa_data$Cytogenetics.Cat.Summary, c("t(8;21)", "inv16", "Diploid", "Trisomy8", "t(9;11)", "Others", "Complex",  "-5/5q-", "-7/7q-", "11q23", "t(6;9)"))
rppa_data$X...7.7q <- factor(rppa_data$X...7.7q)
tapply(rppa_data$ASNS , rppa_data$FAB, median)
tapply(rppa_data$ASNS , rppa_data$X...7.7q, median) 
tapply(rppa_data$ASNS, rppa_data$Cytogenetics.Cat.Summary, median)
bar_data <- rppa_data %>% dplyr::select(Cytogenetics.Risk.Group, poi_median)
bar_data <- bar_data %>% group_by(Cytogenetics.Risk.Group,poi_median) %>%   mutate(count_name_occurr = n())
bar_data_fab <- rppa_data %>% dplyr::select(FAB, poi_median)
bar_data_fab <- bar_data_fab %>% group_by(FAB, poi_median) %>%   mutate(count_name_occurr = n())
bar_data_cyto <- rppa_data %>% dplyr::select(Cytogenetics.Cat.Summary, poi_median)
bar_data_cyto <- bar_data_cyto %>% group_by(Cytogenetics.Cat.Summary, poi_median) %>%   mutate(count_name_occurr = n())
bar_data_7 <- rppa_data %>% dplyr::select(X...7.7q, poi_median)
bar_data_7 <- bar_data_7 %>% group_by(X...7.7q, poi_median) %>%   mutate(count_name_occurr = n())
bar_data_7 <- na.omit(bar_data_7)
chi_cyto_poi<- chisq.test(bar_data$poi_median, bar_data$Cytogenetics.Risk.Group)
chi_fab_poi <- chisq.test(bar_data_fab$poi_median, bar_data_fab$FAB)
chi_poi <- chisq.test(bar_data_cyto$poi_median, bar_data_cyto$Cytogenetics.Cat.Summary)
chi_7_poi<- chisq.test(bar_data_7$poi_median, bar_data_7$X...7.7q)
bar_data_fab <- unique(bar_data_fab)
bar_data <- unique(bar_data)
bar_data_cyto <- unique(bar_data_cyto)
bar_data_7 <- unique(bar_data_7)
bar_data_7$X...7.7q <- ifelse(bar_data_7$X...7.7q == "N", "Absent", "Present")
pdf("Figure 2a.pdf")
ggplot(bar_data_fab, aes(x = FAB, y = count_name_occurr, fill = poi_median )) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("High ASNS", "Low ASNS"), values = c("#D5695D", "#5D8CA8")) +
  xlab("FAB Classification") +
  #scale_x_discrete(guide = guide_axis(angle = 90)) +
  ylab("Percentage of Patients") +
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  geom_text(aes(label = count_name_occurr), position = position_fill(vjust = 0.6)) +
  labs(fill = "ASNS Expression Groups", caption = paste0("p = ",round(chi_fab_poi$p.value,3))) +
  theme_classic()
dev.off()
pdf("Figure 2c.pdf")
ggplot(bar_data, aes(x = Cytogenetics.Risk.Group, y = count_name_occurr,fill=poi_median)) +
  geom_bar(width=0.5,position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("High ASNS", "Low ASNS"), values = c("#D5695D", "#5D8CA8")) +
  xlab("ELN 2017 Risk Stratification") +
  labs(fill = "ASNS Expression Groups", caption = paste0("p = ",round(chi_cyto_poi$p.value,3))) +
  ylab("Percentage of Patients") +
  geom_text(aes(label = count_name_occurr), position = position_fill(vjust = 0.6)) +
  theme_classic()
dev.off()
pdf("Figure S1a.pdf")
ggplot(bar_data_cyto, aes(x = Cytogenetics.Cat.Summary, y = count_name_occurr, fill=poi_median)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("High ASNS", "Low ASNS"), values = c("#D5695D", "#5D8CA8")) +
  xlab("AML Cytogenetics") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ylab("Percentage of Patients") +
  labs(fill = "ASNS Expression Groups", caption = paste0("p = ",round(chi_poi$p.value,3))) +
  geom_text(aes(label = count_name_occurr), position = position_fill(vjust = 0.6)) +
  theme_classic()
dev.off()
pdf("Figure S1c.pdf")
ggplot(bar_data_7, aes(x = X...7.7q, y = count_name_occurr,fill=poi_median)) +
  geom_bar(width=0.3,position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("High ASNS", "Low ASNS"), values = c("#D5695D", "#5D8CA8")) +
  xlab("Chromsome 7/7q Deletion") +
  labs(fill = "ASNS Expression Groups", caption = paste0("p = ",round(chi_cyto_poi$p.value,3))) +
  ylab("Percentage of Patients") +
  geom_text(aes(label = count_name_occurr), position = position_fill(vjust = 0.6)) +
  theme_classic()
dev.off()
