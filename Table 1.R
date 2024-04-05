library(dplyr)
library(gtsummary)
rppa_data <- read.table("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/ASNS/RPPA_Data.csv", sep = "\t")
rppa_data$FAB <- rppa_data$PseudoFAB
rppa_data$FAB <- gsub("M4e", "M4", rppa_data$FAB)
rppa_data$FAB <- gsub("M5a", "M5", rppa_data$FAB)
rppa_data$FAB <- gsub("M5b", "M5", rppa_data$FAB)
rppa_data$FAB <- gsub("M6a", "M6", rppa_data$FAB)
rppa_data$FAB <- gsub("M6b", "M6", rppa_data$FAB)
rppa_data$FAB <- gsub("Unclassifiable", NA, rppa_data$FAB)
rppa_data$FAB <- factor(rppa_data$FAB, levels = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "RAEB-T"))
rppa_data$TX1.BestResponse1 <- rppa_data$TX1.BestResponse
rppa_data$TX1.BestResponse1 <- gsub("CR|CRi", "CR/CRi", rppa_data$TX1.BestResponse1)
rppa_data$TX1.BestResponse1 <- gsub("PR|HI|SD", "PR/HI/SD", rppa_data$TX1.BestResponse1)
rppa_data$TX1.BestResponse1 <- gsub("ED|NED", "ED/NED", rppa_data$TX1.BestResponse1)
rppa_data$TX1.BestResponse1 <- gsub("NR", "NR", rppa_data$TX1.BestResponse1)
unique(rppa_data$TX1.BestResponse1)
clin_data <- rppa_data %>% dplyr::select(Diagnosis, FAB, Age.Dx, Race.Dx, Gender, WBC, Blasts.absolute, Blasts, BM.BlastPercent, Mono...209, Mono...550, Cytogenetics.Risk.Group, TX1.BestResponse1, TX1.Venetoclax, TX1.Remission, vital.status, OS.Wks.SMK, ASNS, poi_quartile, poi_tertile, poi_median, status, Cytogenetics.Cat.Summary)
colnames(clin_data)[1:17] <- c("Diagnosis", "FAB", "Age","Race", "Gender",  "WBC", "Absolute_Blast", "PB_Blast", "BM_Blast", "PB_Mono", "BM_Mono", "Cytogenetics", "Response", "Venetoclax_Treated", "Remission_Status", "Vital_Status", "OS_Weeks")

table0 <- tbl_summary(clin_data, 
            include=c(Age, Gender, Race, FAB, Cytogenetics, WBC, Absolute_Blast, PB_Blast, BM_Blast, PB_Mono, BM_Mono, Venetoclax_Treated, Response, Vital_Status, OS_Weeks,  ASNS),
            by = poi_median,
            missing = "no") %>%
  add_overall() %>%
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.val=TRUE)) %>%
  modify_header(label = "**Variable**") %>%
  bold_labels() 
modify_0 <- table0 %>% modify_caption("**Demographic and clinical characteristics of AML patients in the study set** (N = {N})")
t0 <- add_q(modify_0, method = "fdr", pvalue_fun = NULL, quiet = NULL) %>% as_gt()
#gt::gtsave(t0, filename = "Clinical_Table_by_Medians.png")
gt::gtsave(t0, filename = "/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Code/Table 1/Demographic_Clinical_Table.docx") 

