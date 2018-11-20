install.packages("biclust")
library("biclust")
library(Binarize)

selected_genes <- c(hsa-mir-22, hsa-mir-30a, hsa-mir-31, hsa-let-7b, hsa-mir-3622a, hsa-mir-618, hsa-mir-874, hsa-mir-1307, hsa-mir-1224)

merged_dat <- read.csv('C:/Users/garim/Documents/Capstone/Full data/Normalized/merged_rna_mirna_cnv_log.csv')
rownames(merged_dat) <- merged_dat[,1]
selected_dat <- t(as.matrix(merged_dat[,-c(1,which(nacount>0))]))
get_mean <- apply(selected_dat,2, function(x) mean(x))

plaid_clus <- biclust(selected_dat, method=BCPlaid(), cluster="b")
plotclust(plaid_clus,selected_dat)

bccc_clus <- biclust(selected_dat, method=BCCC(), number=100)
bccc_clus

bccc_genes <- as.data.frame(bccc_clus@RowxNumber)
rownames(bccc_genes) <- colnames(merged_dat[,-c(1,which(nacount>0))])
colnames(bccc_genes) <- paste0("cluster-", 1:100)
write.csv(bccc_genes,'C:/Users/garim/Documents/Capstone/Codes/Flexmix/bccc_genes.csv')

bccc_pat <-as.data.frame(bccc_clus@NumberxCol)
colnames(bccc_pat) <- rownames(merged_dat[,-c(1,which(nacount>0))])
rownames(bccc_pat) <- paste0("cluster-", 1:100)
write.csv(bccc_pat,'C:/Users/garim/Documents/Capstone/Codes/Flexmix/bccc_pat.csv')


binarized_dat <- binarizeMatrix(selected_dat,  method = c("BASCA"))

for (i in 1:100){
  
  
  
}

motifs_clus <- biclust(selected_dat, method=BCXmotifs(), number=100)
motifs_clus

spectral_clus <- biclust(selected_dat,  method=BCSpectral(), number=100)
spectral_clus

bcquest_clust <- biclust(selected_dat,   method=BCQuest(), number=100)
bcquest_clust

nancount <- apply(merged_dat[,-c(1,which(nacount>0))],2, function(x) mean(x))
nancnancountsum(nacount)                 
length(nacount[which(nacount>0)])
datatype <- sapply(merged_dat[,-c(1,which(nacount>0))], class)
table(datatype)

sumdat <- apply(merged_dat[,-c(1,which(nacount>0))],2, function(x) sum(x))



merged_dat[, which(colnames(merged_dat) %in% c('chr1.120000001.121000000'))]
nacount[which(colnames(merged_dat) %in% c('chr1.120000001.121000000'))]
