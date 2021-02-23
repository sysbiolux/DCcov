setwd("./") #~/Desktop/ISB/COVID19/

file.readcounts <- "data/severity_study/GSE147507_RawReadCounts_Human.tsv"
#file.annotations <- "data/severity_study/Biomart.annotations.hg38.txt"

# Import the read count matrix data into R.
counts <- read.csv(file.readcounts, sep="\t",row.names = 1)

# Import feature annotations. 
# Assign feature lenght into a numeric vector.
#annotations <- read.table(file.annotations, sep="\t", header=TRUE)
library(DESeq2)

# Convert from gene.symbol to ensembl.gene
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v79")

library(EnsDb.Hsapiens.v79)

geneSymbols <-  rownames(counts)
geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

# Get the length for all genes
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EDASeq")
library (EDASeq)
length = getGeneLengthAndGCContent(geneIDs2$GENEID, "hsa")
geneIDs2['len'] = length[,1]

geneIDs2<- geneIDs2[complete.cases(geneIDs2),]
geneIDs2 = remove_empty(geneIDs2,"rows") 


write.csv(geneIDs2,'data/severity_study/genes_length.csv')

