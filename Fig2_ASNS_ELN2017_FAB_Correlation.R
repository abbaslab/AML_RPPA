library(ggplot2)
library(rstatix)
library(RColorBrewer)
library(colorspace)
rppa_data <- read.table("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Code/Data/RPPA_Data.csv", sep = "\t")
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

tapply(rppa_data$ASNS , rppa_data$FAB, median) 
results_list <- list()
anovalist <- list()
n <- list()
for(var_name in c("Cytogenetics.Risk.Group", "Cytogenetics.Cat.Summary", "FAB")) {
  df1 <- rppa_data %>% dplyr::select(ASNS, var_name)
  colnames(df1) <- c("ASNS", "Variable")
  df <- na.omit(df1) #omiting NA values as it is not required for statistical analysis
  n[[var_name]] = nrow(df)
  fit <- df %>% rstatix::anova_test(ASNS ~ Variable) %>% adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test <- df %>% 
    rstatix::tukey_hsd(ASNS ~ Variable) %>% 
    rstatix::add_significance() %>% 
    rstatix::add_xy_position()
  results_list[[var_name]] <- stat.test
  anovalist[[var_name]] <- fit

} 
results_list[["Cytogenetics.Risk.Group"]]$y.position[1] = 5
results_list[["Cytogenetics.Risk.Group"]]$y.position[2] = 6
results_list[["Cytogenetics.Risk.Group"]]$y.position[3] = 4
results_list[["Cytogenetics.Cat.Summary"]]$y.position[15] = 5.5
results_list[["FAB"]]$y.position[6] = 3
results_list[["FAB"]]$y.position[8] = 3.5
results_list[["FAB"]]$y.position[9] = 4
results_list[["FAB"]]$y.position[13] = 4.5
results_list[["FAB"]]$y.position[15] = 5
results_list[["FAB"]]$y.position[19] = 5.5
results_list[["FAB"]]$y.position[28] = 6
results_list[["FAB"]]$y.position[31] = 6.3


color_pal <- c("#5D8CA8", "#D5695D", "lightpink", "cornflowerblue", "#F15C80", "khaki","#6C71C4", "#FC8D62", "orchid","tan", "#C24841", "lightgray")
pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure 2b.pdf", width=8)
ggpubr::ggboxplot(subset(rppa_data, !is.na(FAB)), x = "FAB", y = "ASNS", fill = "FAB",  ylim = c(-2,8)) +
  scale_fill_manual(values = color_pal) +
  ggpubr::stat_pvalue_manual(results_list[["FAB"]], label = "p.adj.signif", hide.ns=T) +
  geom_text(data=rppa_data %>%
              group_by(FAB) %>% 
              dplyr::summarise(top = max(ASNS), n=dplyr::n()), 
            aes(x=FAB, y=7, label= paste0("n = ", n)), 
            nudge_y=1) +
  ggplot2::labs(fill="FAB", subtitle=rstatix::get_test_label(anovalist[["FAB"]], detailed = FALSE), caption = paste0("n = ", n[["FAB"]])) + 
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  
  xlab("FAB Classification") +
  scale_x_discrete(limits = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "RAEB-T"))+
  ylab("Log2 Normalized ASNS Expression") + theme(legend.position = "none")
dev.off()

pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure 2d.pdf", width=8)
ggpubr::ggboxplot(subset(rppa_data, !is.na(Cytogenetics.Risk.Group)), x = "Cytogenetics.Risk.Group", y = "ASNS",  fill = "Cytogenetics.Risk.Group", ylim = c(-2,8)) +
  scale_fill_manual(values = color_pal) +
  ggpubr::stat_pvalue_manual(results_list[["Cytogenetics.Risk.Group"]], label = "p.adj.signif", hide.ns=F) +
  geom_text(data=rppa_data %>%
              group_by(Cytogenetics.Risk.Group) %>% 
              dplyr::summarise(top = max(ASNS), n=dplyr::n()), 
            aes(x=Cytogenetics.Risk.Group, y=6, label= paste0("n = ", n)), 
            nudge_y=1) +
  ggplot2::labs(fill="ELN 2017", subtitle = rstatix::get_test_label(anovalist[["Cytogenetics.Risk.Group"]], detailed = FALSE)) +
  xlab("ELN 2017 Risk Stratification") +
  scale_x_discrete(limits = c("Favorable", "Intermediate", "Adverse")) +
  ylab("Log2 Normalized ASNS Expression") + theme(legend.position = "none")
dev.off()
pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure S1b.pdf", width=8)
ggpubr::ggboxplot(subset(rppa_data, !is.na(Cytogenetics.Cat.Summary)), x = "Cytogenetics.Cat.Summary", y = "ASNS",  fill = "Cytogenetics.Cat.Summary", ylim = c(-2,8)) +
  scale_fill_manual(values = color_pal) +
  ggpubr::stat_pvalue_manual(results_list[["Cytogenetics.Cat.Summary"]], label = "p.adj.signif", hide.ns=T) +
  geom_text(data=rppa_data %>%
              group_by(Cytogenetics.Cat.Summary) %>% 
              dplyr::summarise(top = max(ASNS), n=dplyr::n()), 
            aes(x=Cytogenetics.Cat.Summary, y=6, label= paste0("n = ", n)), 
            nudge_y=1) +
  ggplot2::labs(fill="Cytogenetics", subtitle = rstatix::get_test_label(anovalist[["Cytogenetics.Cat.Summary"]], detailed = FALSE),
                caption = paste0("n = ", n[["Cytogenetics.Cat.Summary"]])) + 
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  xlab("AML Cytogenetics") +
  scale_x_discrete(limits = c("t(8;21)", "inv16", "Diploid", "Trisomy8", "t(9;11)", "Others", "Complex",  "-5/5q-", "-7/7q-", "11q23", "t(6;9)"))+
  ylab("Log2 Normalized ASNS Expression") + theme(legend.position = "none")
dev.off()
