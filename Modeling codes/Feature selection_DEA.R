# Note: This code has been taken form an existing research on effective methods of selecting relevant genes for a particular cancer type

assignInNamespace("version_info", c(devtools:::version_info, list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), "devtools")
find_rtools() # is TRUE now

devtools::install_github("r-lib/devtools")
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
library(TCGAbiolinks)

install.packages("devtools.zip", repos = NULL, type = "source")
library(devtools) 
build_github_devtools()

install_github("BioinformaticsFMRP/TCGAbiolinks")
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks", dependencies = T,force = TRUE)

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "miRNA expression array",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  #barcode = listSamples, 
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCAmiRNAseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCAmiRNAseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")

# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
BRCAmiRnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCAmiRnaseqSE)



# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo =  geneInfo)
# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])