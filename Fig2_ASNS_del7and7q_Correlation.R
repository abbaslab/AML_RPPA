library(ggplot2)
library(rstatix)
rppa_data <- read.table("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Code/Data/RPPA_Data.csv", sep = "\t")
rppa_data$X...7.7q <- ifelse(rppa_data$X...7.7q=="Y", "-7/7q", "No -7/7q")
rppa_data$X...7.7q <- factor(rppa_data$X...7.7q, levels = c("-7/7q", "No -7/7q"))
rppa_data$ch7del <- ifelse(rppa_data$Cytogenetics.Cat.Summary %in% c("-5/5q-", "11q23", "COMPLEX", "DIPLOID", "inv16", "Misc", "t(6;9)", "t(8;21)", "t(9;11)", "trisomy8"), "Other", ifelse(rppa_data$Cytogenetics.Cat.Summary == "-7/7q-", "-7/7q-", NA))
rppa_data$ch7del <- factor(rppa_data$ch7del, levels = c("-7/7q-", "Other"))
rppa_data$del7complex <- ifelse(rppa_data$Cytogenetics.Cat.Summary == "-7/7q-", "Isolated -7/7q", ifelse((rppa_data$Cytogenetics.Cat.Summary == "COMPLEX" & (rppa_data$AML_Cat_Cyto == "-7/7q-"| rppa_data$AML_Cat_Cyto == "-5/5q- and -7/7q-")), "CK with -7/7q", ifelse((rppa_data$Cytogenetics.Cat.Summary == "COMPLEX" & (rppa_data$AML_Cat_Cyto != "-7/7q-"| rppa_data$AML_Cat_Cyto != "-5/5q- and -7/7q-")), "CK without -7/7q", NA)))
rppa_data$del7complex <- factor(rppa_data$del7complex, levels = c("Isolated -7/7q", "CK with -7/7q", "CK without -7/7q", NA))
tapply(rppa_data$ASNS , rppa_data$del7complex, median)
tapply(rppa_data$ASNS , rppa_data$Cytogenetics.Cat.Summary, median)
tapply(rppa_data$ASNS , rppa_data$X...7.7q, median)
results_list <- list()
anovalist <- list()
n <- list()
for(var_name in c("X...7.7q", "ch7del", "del7complex")){
  df1 <- rppa_data %>% dplyr::select(ASNS, var_name)
  colnames(df1) <- c("ASNS", "Variable")
  df <- na.omit(df1) #omiting NA values as it is not required for statistical analysis
  n[[var_name]] = nrow(df)
  if(nlevels(rppa_data[[var_name]]) <=2) {
    #print(var_name)
    stat.test <- df %>% 
      t_test(ASNS ~ Variable, paired = FALSE) 
  }
  else{
    fit <- df %>% rstatix::anova_test(ASNS ~ Variable)
    stat.test <- df %>% 
      rstatix::tukey_hsd(ASNS ~ Variable) %>% 
      rstatix::add_significance() %>% 
      rstatix::add_xy_position(step.increase = 0.05)
    anovalist[[var_name]] <- fit
  }
  results_list[[var_name]] <- stat.test
}
#results_list[["del7complex"]]$y.position[2] = 5
pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure 2e.pdf", width=8)
ggpubr::ggboxplot(subset(rppa_data, !is.na(X...7.7q)), x = "X...7.7q", y = "ASNS", 
                        fill = "X...7.7q",  ylim = c(-3,8)) +
        scale_fill_manual(labels = c("-7/7q", "No -7/7q"),values = c("#5D8CA8", "#D5695D", "lightpink")) +
        ggpubr::stat_pvalue_manual(results_list[["X...7.7q"]], label = "p", y.position = 5) +  
        geom_text(data=rppa_data %>%
                    group_by(X...7.7q) %>% 
                    dplyr::summarise(top = max(ASNS), n=dplyr::n()), 
                  aes(x=X...7.7q, y=5.5, label= paste0("n = ", n)), 
                  nudge_y=1) +
        ggplot2::labs(fill = "Deletion 7/7q",title = paste0("ASNS ~ Del 7/7q")) + 
  xlab("") +
  scale_x_discrete(limits = c("-7/7q", "No -7/7q")) +
  ylab("Log2 Normalized ASNS Expression") +theme(legend.position = "none")
dev.off()

pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure 2f.pdf", width=8)
ggpubr::ggboxplot(subset(rppa_data, !is.na(del7complex)), x = "del7complex", y = "ASNS",  fill = "del7complex",  ylim = c(-2,8)) +
  scale_fill_manual(labels = c("Isolated -7/7q", "CK with -7/7q", "CK without -7/7q", "Other"),values = c("#5D8CA8", "#D5695D", "lightpink")) +
  ggpubr::stat_pvalue_manual(results_list[["del7complex"]], label = "p.adj.signif", hide.ns=T) +
  geom_text(data=rppa_data %>%
              group_by(del7complex) %>% 
              dplyr::summarise(top = max(ASNS), n=dplyr::n()), 
            aes(x=del7complex, y=5.5, label= paste0("n = ", n)), 
            nudge_y=1) +
  ggplot2::labs(fill = "Cytogenetics",title = paste0("ASNS ~ del7complex"), subtitle = rstatix::get_test_label(fit, detailed = FALSE),) +
  xlab("") +
  scale_x_discrete(limits = c("Isolated -7/7q","CK with -7/7q","CK without -7/7q"))+
  ylab("Log2 Normalized ASNS Expression")+ theme(legend.position = "none")
dev.off()

  

