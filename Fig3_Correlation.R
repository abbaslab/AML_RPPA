library(corrplot)
library(dplyr)
library(ggpubr)
rppa_data <- read.table("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Code/Data/RPPA_Data.csv", sep = "\t")
corr_matrix <- as.matrix(rppa_data[,609:1037]) #extract the protein data into a matrix
source("http://www.sthda.com/upload/rquery_cormat.r")
cormat <- rquery.cormat(corr_matrix, type="full", graph=FALSE)
#saveRDS(cormat, "Protein_Correlation.RData")
corr_mat <- as.data.frame(cormat[["r"]])
corr_mat <- corr_mat %>% dplyr::select(ASNS)
colnames(corr_mat) <- "Correlation"
p_mat <- as.data.frame(cormat[["p"]])
p_mat <- p_mat %>% dplyr::select(ASNS)
colnames(p_mat) <- "p"
corrmatrix <- cbind(corr_mat, p_mat)
class(corrmatrix$p) <- "numeric"
corrmatrix$Protein <- rownames(corrmatrix)
corrmatrix$padj <- p.adjust(corrmatrix$p,method = "BH", n=length(corrmatrix$p))
corrmatrix_padj <- corrmatrix %>% dplyr::filter(padj <= 0.0001) #padj cutoff = 0.05/428 proteins

write.table(corrmatrix_padj, "Correlations_Significant.txt", sep = "\t")

corrmatrix_plot <- corrmatrix_padj %>% filter(Correlation>=0.2 | Correlation<=-0.2) #extract positive correlation values
corrmatrix_plot <- corrmatrix_plot %>% filter(Correlation != 1) %>% arrange(desc(Correlation))
corrmatrix_plot$corr <- ifelse(corrmatrix_plot$Correlation >=0, "Positive", "Negative")
pos <- corrmatrix_plot %>% filter(corr=="Positive")
neg <- corrmatrix_plot %>% filter(corr=="Negative")
corrmatrix_plot$Protein <- rownames(corrmatrix_plot)
corrmatrix_plot$corr <- factor(corrmatrix_plot$corr, levels=c("Positive", "Negative"))
write.table(corrmatrix_plot, "/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Code/Figure 4/Correlations_Significant_ASNS.txt", sep = "\t")

pdf("/Users/nnarayanan/Library/CloudStorage/OneDrive-InsideMDAnderson/NishaNarayanan/ASNS/Figures/Figure 3a.pdf", width = 14)
ggbarplot(corrmatrix_plot, x = "Protein", y = "Correlation", fill = "corr", color = NA, palette = c("red3", "steelblue4")) +  
  xlab("") +
  ylab("Pearson Correlation with ASNS") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none") + font("y.text", size = 12, face = "bold")
dev.off()