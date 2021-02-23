library(ggrepel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv('SKO_ALL_Studies.csv',row.names = 1)

ggplot(data, aes(Phenotype_Essentiality, Phenotype_Safety, label = rownames(data))) +
  geom_text_repel() +
  geom_point(color = 'red') +
 # geom_label_repel()+
  theme_classic(base_size = 15)
ggsave("Figs/SKO_Scatter_Plot.png",dpi = 300)

#plot(Phenotype_Essentiality ~Phenotype_Safety, col="lightblue", pch=19, cex=2,data=data)
#text(Phenotype_Essentiality ~Phenotype_Safety, labels=rownames(data),data=data, cex=0.9, font=1)
