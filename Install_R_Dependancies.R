if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR",version ="3.30.3" )
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2",version ="1.28.1")

if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!requireNamespace("FactoMineR", quietly = TRUE))
  install.packages("FactoMineR",)
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
library('devtools')
#networkD3==0.4
devtools::install_github('christophergandrud/networkD3')




  