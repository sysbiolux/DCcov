setwd(".")
library(DESeq2)
library(knitr)
library(arrayQualityMetrics)
library(devtools)
#install_github("vqv/ggbiplot")
data = read.csv('./data/severity_study/GSE147507_RawReadCounts_Human.tsv',sep = '\t',row.names = 1)
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

# Create quality report
#minimalSet <- ExpressionSet(assayData=as.matrix(data))
#arrayQualityMetrics(expressionset = minimalSet, outdir = "Quality_Report_GSE147507", force = TRUE, do.logtransform = FALSE)

# detect outliers via hierarchical clustering
#distmat = dist(t(data))
#hcldat = hclust(distmat, method="average")
#plot(hcldat)

#BiocManager::install('arrayQualityMetrics')
# Create quality report for Moran et al. dataset
library(arrayQualityMetrics)
minimalSet <- ExpressionSet(assayData=as.matrix(data))
arrayQualityMetrics(expressionset = minimalSet, outdir = "Quality_Report_Moran", force = TRUE, do.logtransform = FALSE)


#We filter out lowly expressed genes.

library(edgeR)
keep <- filterByExpr(data)
table(keep)
keep

data <- data[keep,]

#genotype + condition + genotype:condition
dds <- DESeqDataSetFromMatrix(countData = data,colData = meta,~condition +genotype:condition)

rld <- rlog(dds)
plotPCA(rld)

# also possible to perform custom transformation:
dds <- estimateSizeFactors(dds)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
plotPCA( DESeqTransform( se ) )

boxplot(log10(counts(dds)+1))

dds <- estimateSizeFactors(dds)
boxplot(log10(counts(dds,normalized=TRUE)+1))

vsd <- vst(dds)
class(vsd)
assay(vsd)[1:3,1:3]

plotPCA(vsd, "condition")

## PCA
library("ggplot2")
pcaData <- plotPCA(vsd, intgroup = c( "genotype", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = genotype, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()


#dds$genotype <- target$Tyoe
dds$condition 
dds$genotype
dds$genotype = relevel( dds$genotype, "Mock")

#design(dds) <- ~ genotype + condition + genotype:conditionx 

assay(dds)[1:10,]

dds <- DESeq(dds)
resultsNames(dds)

#res <- results(dds, contrast=c("condition","Series1","Series2"))
#res <- results(dds, contrast=c("genotype","Infected","Mock"))
res <- results(dds, name=c("Intercept"))
res <- subset(res, padj < 0.05)
res <- subset(res, abs(log2FoldChange) >1)
write.table(res,paste0('./data/severity_study/DEGS_Intercept.txt'))
#res <- results(dds, list(c("condition_Infected_vs_Mock","conditionSeries16.genotypeInfected")))
#c("condition_B_vs_A","genotypeIII.conditionB")
#res <- subset(res, padj < 0.05)

conditions_unq = unique(condition)
for (i in 1:length(conditions_unq)){
  condition_ = conditions_unq[i]
  res <- results(dds, name=c(paste0("condition",condition_,".genotypeInfected")))
  res <- subset(res, padj < 0.05)
  res <- subset(res, abs(log2FoldChange) >1)
  write.table(res,paste0('./data/severity_study/DEGs/',condition_,'.txt'))
}

design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype","Infected","Mock"))
res <- subset(res, padj < 0.05)
res <- subset(res, abs(log2FoldChange) >1)
res
write.table(res,paste0('./data/severity_study/DEGS_Inf_vs_Mock.txt'))
