
# Data load
newdat_mirna <- read.csv('C:/Users/garim/Documents/Capstone/Codes/Flexmix/predictive fxn/mirna_dat.csv')
newdat_rna <- read.csv('C:/Users/garim/Documents/Capstone/Codes/Flexmix/predictive fxn/rna_dat.csv')
newdat_meth <- read.csv('C:/Users/garim/Documents/Capstone/Codes/Flexmix/predictive fxn/meth_dat.csv')
mirna_gm <- read.csv('C:/Users/garim/Documents/Capstone/Codes/Flexmix/predictive fxn/mirna_gm.csv')
rna_gm <- read.csv('C:/Users/garim/Documents/Capstone/Codes/Flexmix/predictive fxn/rna_gm.csv')
drug_clinical <- read.csv('C:/Users/garim/Documents/Capstone/Codes/Flexmix/predictive fxn/drug_clin.csv')

# Load model
model_obj <- readRDS('C:/Users/garim/Documents/Capstone/Codes/Flexmix/predictive fxn/mymodel.rds')


#Process raw data files to have patient IDs as column names and gene names as row names

rownames(mirna_gm) <- mirna_gm[,1]
mirna_gm <- mirna_gm[,-1]
rownames(rna_gm) <- rna_gm[,1]
rna_gm <- rna_gm[,-1]

newdat_mirna_num <-as.data.frame(-t(newdat_mirna[2:nrow(newdat_mirna), 2:ncol(newdat_mirna)]))
rownames(newdat_mirna_num) <-newdat_mirna[1,]
colnames(newdat_mirna_num) <-newdat_mirna[,1]

newdat_rna_num <-as.data.frame(-t(newdat_rna[2:nrow(newdat_rna), 2:ncol(newdat_rna)]))
rownames(newdat_rna_num) <-newdat_rna[1,]
colnames(newdat_rna_num) <-newdat_rna[,1]

#Normalize the mirna, rna, and methylation datasets

#Input transposed mirna/ rna dataframe with genes as row names and patient ID as column names
#The Geometric mean map dataframe should have one column with geometric means and gene names as row names

normalize <- function(data, gm){
  dat_num <- merge(data,gm, by = "row.names")
  colnames(dat_num) <- dat_num[,1]
  dat_num <- dat_num[,-1]
  dat_ratio <- dat_num[,1:(ncol(dat_num)-1)]/dat_num$geom_mean
  dat_med <- apply(dat_ratio, 2, function(x) median(x))
  
  for (i in 1:ncol(dat_ratio)){
    dat_ratio[,i] <- dat_ratio[,i]/dat_med[i]
  }
  
  dat_ratio <- as.data.frame(dat_ratio)
  rownames(dat_ratio) <- rownames(data)
  colnames(dat_ratio) <- colnames(data)
  return(dat_ratio)
  
}

mirna_norm <- normalize(newdat_mirna_num, mirna_gm)
rna_norm <- normalize(newdat_rna_num, rna_gm)

meth_probes <- c('cg05337753','cg02282892')
meth_norm <- apply(meth_probes[,c(colnames(newdat_meth) %in% meth_probes)], 2, function(x) log2(x/(1-x)))
rownames(meth_norm) <- newdat_meth[,1]

#process the drug data
rownames(drug_clinical) <- drug_clinical[,1]
drug_clinical <- drug_clinical[,-1]

final_data <- merge(mirna_norm, rna_norm, by = 'row.names')
colnames(final_data) <- final_data[,1]
final_data <- final_data[,-1]

final_data <- merge(final_data, meth_norm, by = 'row.names')
colnames(final_data) <- final_data[,1]
final_data <- final_data[,-1]

final_data <- merge(final_data, drug_clinical, by = 'row.names')
colnames(final_data) <- final_data[,1]
final_data <- final_data[,-1]

# Prediction for all possible drugs
drug_names <- c("Cyclophosphamide","Doxorubicin", "Tamoxifen", "Anastrozole", "Paclitaxel", "Docetaxel", "Fluorouracil", "Trastuzumab", "Letrozole", "Exemestane", "Epirubicin", "Carboplatin","Methotrexate")
PFI_save <- as.data.frame(matrix(NA, length(drug_names), nrow(final_data)))
rownames(PFI_save)<- drug_names
colnames(PFI_save)<- rownames(final_data)

for (i in 1:length(drug_names)){
    mod_final_data <- final_data
    mod_final_data[,grep(drug_names[i], colnames(mod_final_data))]<- mod_final_data[,grep(drug_names[i], colnames(mod_final_data))]+1
    PFI_save[,grep(drug_names[i], colnames(PFI_save))] <- predict(model_obj, data=mod_final_data, type = 'class')
  }

#return(PFI_save)

