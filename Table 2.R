library(dplyr)
library(gtsummary)
rppa_data <- read.table("RPPA_Data.csv", sep = "\t")
rppa_data$trisomy8 <- ifelse(is.na(rppa_data$trisomy8), NA, ifelse(rppa_data$trisomy8=="N", "No", "Trisomy8"))
rppa_data$MD.SF3B1.mutated <- ifelse(rppa_data$MD.SF3B1 == "NA", NA, ifelse(rppa_data$MD.SF3B1=="N" | rppa_data$MD.SF3B1=="NEG", "N", "Y"))
rppa_data$MD.DNMT3A.mutated <- ifelse(rppa_data$MD.DNMT3A == "NA", NA, ifelse(rppa_data$MD.DNMT3A=="N" | rppa_data$MD.DNMT3A=="NEG", "N", "Y"))
rppa_data$Karyotype <- ifelse(rppa_data$Cytogenetics.Cat.Summary == "DIPLOID", "Diploid", ifelse(rppa_data$Cytogenetics.Cat.Summary=="COMPLEX", "Complex", ifelse(is.na(rppa_data$Cytogenetics.Cat.Summary), NA, "Not Complex")))
rppa_data$del7_7q <- ifelse(rppa_data$X...7.7q == "Y", "Del7_7q", ifelse(rppa_data$X...7.7q == "N", "No", NA))
rppa_data$del5_5q <- ifelse(rppa_data$X...5.5q == "Y", "Del5_5q", ifelse(rppa_data$X...5.5q == "N", "No", NA))
rppa_data$del17_17p <- ifelse(rppa_data$X...17...abn.17p. == "Y", "Del17_17p", ifelse(rppa_data$X...17...abn.17p. == "N", "No", NA))

oncoprint <- rppa_data[,colnames(rppa_data) %in% c("Karyotype","del5_5q","del7_7q","trisomy8","del17_17p", "MD.EZH2.mutated","MD.TP53.mutated","MD.FLT3.ITD.mutated","MD.NPM1.mutated","MD.RAS.mutated","MD.SF3B1.mutated","MD.IDH1.mutated","MD.IDH2.mutated","MD.DNMT3A.mutated","MD.TET2.mutated","MD.ASXL1.mutated","Cytogenetics.Risk.Group")]
colnames(oncoprint)
colnames(oncoprint)=c("ELN2017", "Trisomy8", "ASXL1", "EZH2", "FLT3", "IDH1", "IDH2", "NPM1", "RAS", "TET2", "TP53", "SF3B1", "DNMT3A", "Karyotype", "Del7/7q", "Del5/5q", "Del17/17p") 
col_order = c("Karyotype","Del5/5q","Del7/7q","Trisomy8","Del17/17p","EZH2","TP53","FLT3","NPM1","RAS","SF3B1","IDH1","IDH2","DNMT3A","TET2","ASXL1","ELN2017")
oncoprint <- oncoprint[,col_order]
oncoprint$poi_median <- rppa_data$poi_median

table1 <- tbl_summary(oncoprint, by = poi_median, type = all_dichotomous() ~ "categorical") %>%
  add_overall() %>%
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.val=TRUE)) %>%
  modify_header(label = "**Variable**") %>%
  bold_labels() 
modify_1 <- table1 %>% modify_caption("**Mutation characteristics of AML patients in the study set** (N = {N})")
t1 <- add_q(modify_1, method = "fdr", pvalue_fun = NULL, quiet = NULL) %>% as_gt()

gt::gtsave(t1, filename = "Mutation_Table.docx") 

