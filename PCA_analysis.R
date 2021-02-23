
setwd(".")

library(readxl)
library(tibble)
library("FactoMineR")
library("factoextra")
library(DESeq2)
library(knitr)
library(arrayQualityMetrics)
library(devtools)
#install_github("vqv/ggbiplot")
data = read.csv('./data/severity_study/RPKM_3_all.csv',row.names = 1)
#data = t(data)
rownames(data)
# read metadata
target = read.csv('./data/severity_study/Metadata_3.csv')
# remove samples not in SARS-CoV-2
series_selected = unique(target[target$Type=='SARS-CoV-2','Series'])
target =target[target$Series %in% series_selected,]
rownames(target) = target$Sample
# remove phenotypes with less than 6 samples ( 3 replicates per case)
tab <- table(target$Series)
target =target[target$Series %in% names(tab)[tab>5],]
# remove sample with drug pertubagen
target = target[- grep("Rux", target$Sample),]

#change sample names in the count table such as in the metadata
colnames(data) = gsub('SARS.CoV.2','Infected',colnames(data))
colnames(data) = gsub('A549.ACE2','A549_ACE2',colnames(data))
target$Sample = gsub('SARS.CoV.2','Infected',target$Sample)
target$Sample = gsub('A549.ACE2','A549_ACE2',target$Sample)
target$Type = gsub('SARS-CoV-2','Infected',target$Type)

# select samples only in the metadata after filteration
data = data[,colnames(data) %in% as.character(target$Sample)]
data = data[,order(colnames(data))]
target = target[order(target$Sample),]
#define the condition
condition = target$Series
genotype = as.factor(target$Type)

meta = DataFrame(condition,genotype)


# On Normalized Counts

norm <- data

#norm <- column_to_rownames(norm,'Gene')
norm <- norm+0.001 #to remove -inf

norm_log <- log2(norm)
rownames <- rownames(norm_log)
colnames <- colnames(norm_log)
new_labels = c("NBHE_2","A549_2_ACE2","A549_0.02","A549_2","A549_0.2_ACE2","Calu3_2")
new_labels = rep(new_labels,each=6)
conditions = target$Series
type <- target$Type


pal_m=c("greenyellow", "pink", "deepskyblue","red", "brown", "violet")#,"turquoise", "brown")

PC <- PCA(t(norm_log), graph = FALSE)
svg(file="./Figs/PC1_vs_PC2_condition.svg")
fviz_pca_ind(PC, 
             col.ind = factor(new_labels),palette = pal_m,
             #palette = "jco",
             col.ind.sup =factor(target$Type),
             pointshape = 20, pointsize = 2,
             geom =  c("point"), 
             #addEllipses = TRUE,
             #ellipse.type = "confidence",  ellipse.level=0.95,
             title='log2(normalized counts)')
dev.off()

svg(file="Figs/PC1_vs_PC3_condition.svg")

fviz_pca_ind(PC, axes=c(1,3),palette = pal_m,
             col.ind = factor(new_labels),
             #palette = "jco",
             col.ind.sup =factor(target$Type),
             pointshape = 20, pointsize = 2,
             geom =  c("point"), 
             #addEllipses = TRUE,
             #ellipse.type = "confidence",  ellipse.level=0.95,
             title='log2(normalized counts)')
dev.off()


# Cluster by infected vs mock
svg(file="Figs/PC1_vs_PC2_type.svg")

fviz_pca_ind(PC, 
             col.ind = factor(type),
             palette = "jco",
             col.ind.sup =factor(target$Type),
             pointshape = 20, pointsize = 2,
             geom =  c("point"), 
             #addEllipses = TRUE,
             #ellipse.type = "confidence",  ellipse.level=0.95,
             title='log2(normalized counts)')
dev.off()
svg(file="Figs/PC1_vs_PC3_type.svg")

fviz_pca_ind(PC, axes=c(1,3),
             col.ind = factor(type),
             palette = "jco",
             col.ind.sup =factor(target$Type),
             pointshape = 20, pointsize = 2,
             geom =  c("point"), 
             #addEllipses = TRUE,
             #ellipse.type = "confidence",  ellipse.level=0.95,
             title='log2(normalized counts)')
dev.off()
