## Doing DEG analysis on the Mock samples between Series5 and Series16 (then with sars-cov-2 infected samples)
# To find the effect of ACE2 transfection
setwd(".")
library(DESeq2)
library(knitr)
library(devtools)
require(data.table) # v1.9.0+
library(edgeR)

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

## Select Mock samples only in Series5 and Series16
target = target[target$Type %in% c('Mock'),]

# select samples only in the metadata after filteration
data = data[,colnames(data) %in% as.character(target$Sample)]
data = data[,order(colnames(data))]
target = target[order(target$Sample),]
#define the condition
condition = as.factor(target$Series)
genotype = as.factor(target$Type)
meta = DataFrame(condition,genotype)

#We filter out lowly expressed genes.
keep <- filterByExpr(data)
data <- data[keep,]
dds <- DESeqDataSetFromMatrix(countData = data,colData = meta,~condition)
dds <- estimateSizeFactors(dds)
design(dds) <- ~  condition

## The effect of ACVE2 transfection in mock samples only
dds$'condition' <- relevel( dds$condition, 'Series5')
dds <- DESeq(dds)
res <- results(dds, name=c('condition_Series16_vs_Series5'))
res <- subset(res, padj < 0.05)
res1 <- subset(res, abs(log2FoldChange) >1)
res1
write.table(res1,paste0('./data/severity_study/DEGs/ACE2_Mock.txt'))



## The effect of ACE2 transfection in infected samples only
#install_github("vqv/ggbiplot")
data = read.csv('./data/severity_study/GSE147507_RawReadCounts_Human.tsv',sep = '\t',row.names = 1)
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

## Select infected samples only in Series5 and Series16
target = target[target$Type %in% c('Infected'),]

# select samples only in the metadata after filteration
data = data[,colnames(data) %in% as.character(target$Sample)]
data = data[,order(colnames(data))]
target = target[order(target$Sample),]
condition = as.factor(target$Series)
genotype = as.factor(target$Type)
meta = DataFrame(condition,genotype)

#We filter out lowly expressed genes.
keep <- filterByExpr(data)
table(keep)
data <- data[keep,]
dds <- DESeqDataSetFromMatrix(countData = data,colData = meta,~condition)
dds <- estimateSizeFactors(dds)

design(dds) <- ~ condition 
dds$'condition' <- relevel( dds$condition, 'Series5')
dds <- DESeq(dds)#modelMatrixType = 
resultsNames(dds)
#res <- results(dds, contrast=c('condition:genotype','Series16:Infected','Series5:Infected'))
res <- results(dds, name=c('condition_Series16_vs_Series5'))
res <- subset(res, padj < 0.05)
res2 <- subset(res, abs(log2FoldChange) >1)
res2
write.table(res2,paste0('./data/severity_study/DEGs/ACE2_Infected.txt'))