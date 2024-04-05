library(ggplot2)
rppa_data <- read.table("/RPPA_Data.csv", sep = "\t")

pdf("Figure 1b.pdf")
ggplot(rppa_data, aes(x = ASNS)) + 
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.5, color = "black",fill = "steelblue3") + 
  geom_density(color="red", lwd=1)  +
  xlab("Log2 Normalized ASNS Expression Values") +
  ylab("Density") +
  theme_classic()
dev.off()

