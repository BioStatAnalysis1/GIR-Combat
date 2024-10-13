
rm(list=ls())
library(sva)
library(Seurat)
packageVersion('Seurat')
library(matrixStats)

data_dir <- 'D:/Users/lenovo/Desktop/小论文/小论文真实数据2/Batch-effect-removal-benchmarking-master/Data/dataset1/'

utils_dir <- 'D:/Users/lenovo/Desktop/小论文/小论文真实数据2/Batch-effect-removal-benchmarking-master/Script/Combat/'

source(paste0(utils_dir,'combat_functions.R'))


TPM_file <- 'dataset1_sm_uc3.txt'  # replace by link to dataset
sample_file <- 'sample_sm_uc3.txt' # replace by link to dataset


myData <- read.table(paste0(data_dir,TPM_file),sep="\t",header=T,row.names=1,check.names=F)


mySample <- read.table(paste0(data_dir,sample_file),sep="\t",header=T,row.names=1,check.names=F)


opt1 <- 'D:/Users/lenovo/Desktop/小论文/小论文真实数据2/Batch-effect-removal-benchmarking-master/Script/Combat/parametric/'



myFilteredData <- filter_data_mtx(myData, base_name, is_filter_cells=TRUE,min_genes=300, 
                                  is_filter_genes=TRUE, min_cells=10)
myFilteredData <- as.matrix(myFilteredData)
row_max = matrixStats::rowMaxs(myFilteredData)
row_min = matrixStats::rowMins(myFilteredData)

if(sum(row_max==row_min)>0){
  myFilteredData = myFilteredData[row_max != row_min,]
}

cells_use <- colnames(myFilteredData) # 565(细胞) 细胞名
mySample <- mySample[cells_use,] # 565 (细胞)  *  3( "cell" "celltype" "batch")

########################
X<-t(myFilteredData) # 565(细胞)*15074(基因)
# batch1
r.datah1_final <- t(X[which(mySample$batch =="Batch1"),])

unc.merge_batch1 <- mySample[which(mySample$batch =="Batch1"),] 
table(unc.merge_batch1$celltype)


# batch1中的细胞类型1
result1_1=NULL
for(i in ((unc.merge_batch1[which(unc.merge_batch1$celltype=="CD141"),])$cell) ) { 
  # for(i in test2) {            # Head of for-loop
  # i="D101_51"
  result11111111 <- r.datah1_final[,which(colnames(r.datah1_final)==i)]
  print(which(colnames(r.datah1_final)==i))
  result1_1 = cbind(result1_1,result11111111)
}
sd(result1_1) 

# batch1中的细胞类型2
result1_2=NULL
for(i in ((unc.merge_batch1[which(unc.merge_batch1$celltype=="DoubleNeg"),])$cell) ) { 
  # for(i in test2) {            # Head of for-loop
  # i="D101_51"
  result11111111 <- r.datah1_final[,which(colnames(r.datah1_final)==i)]
  print(which(colnames(r.datah1_final)==i))
  result1_2 = cbind(result1_2,result11111111)
}
sd(result1_2)  

# batch1中的细胞类型3
result1_3=NULL
for(i in ((unc.merge_batch1[which(unc.merge_batch1$celltype=="pDC"),])$cell) ) { 
  # for(i in test2) {            # Head of for-loop
  # i="D101_51"
  result11111111 <- r.datah1_final[,which(colnames(r.datah1_final)==i)]
  print(which(colnames(r.datah1_final)==i))
  result1_3 = cbind(result1_3,result11111111)
}
sd(result1_3)  

r.datah2_final <- t(X[which(mySample$batch =="Batch2"),]) 

unc.merge_batch2 <- mySample[which(mySample$batch =="Batch2"),] 
table(unc.merge_batch2$celltype)

result2_1=NULL
for(i in ((unc.merge_batch2[which(unc.merge_batch2$celltype=="CD1C"),])$cell) ) { 
  
  result11111111 <- r.datah2_final[,which(colnames(r.datah2_final)==i)]
  print(which(colnames(r.datah2_final)==i))
  result2_1 = cbind(result2_1,result11111111)
}
sd(result2_1) 

# batch2中的细胞类型2
result2_2=NULL
for(i in ((unc.merge_batch2[which(unc.merge_batch2$celltype=="DoubleNeg"),])$cell) ) { 
 
  result11111111 <- r.datah2_final[,which(colnames(r.datah2_final)==i)]
  print(which(colnames(r.datah2_final)==i))
  result2_2 = cbind(result2_2,result11111111)
}
sd(result2_2) 


result2_3=NULL
for(i in ((unc.merge_batch2[which(unc.merge_batch2$celltype=="pDC"),])$cell) ) { 
  # for(i in test2) {            # Head of for-loop
  # i="D101_51"
  result11111111 <- r.datah2_final[,which(colnames(r.datah2_final)==i)]
  print(which(colnames(r.datah2_final)==i))
  result2_3 = cbind(result2_3,result11111111)
}
sd(result2_3)

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization


data <- myFilteredData# 15074   (基因)*565(细胞)
######################################
repeat_measure <- as.numeric(as.factor(mySample$batch)) # 相当于原来的batch
people_kind <- as.numeric(as.factor(mySample$celltype))
people_kind_factor <- as.factor(people_kind)
old_pheno <-data.frame(num = 1:ncol(data), measure=repeat_measure, people_kind=people_kind_factor)

combat_data= ComBat(dat=data, batch=old_pheno$measure, mod=NULL, par.prior=TRUE)

addmod_real_8 = model.matrix(~as.factor(people_kind),data=old_pheno)#  704(细胞)*3("(Intercept)"  "as.factor(people_kind)2" "as.factor(people_kind)3")

combat_edata_addmod_real_noref = ComBat(dat=data, batch=old_pheno$measure, mod=addmod_real_8, par.prior=TRUE)

# 打印 batch1 中每种细胞类型的数量
# cat("Batch1 cell type counts:\n")
# batch1_counts <- table(mySample[mySample$batch == "Batch1", "celltype"])
# print(batch1_counts[c("CD141", "CD1C", "DoubleNeg", "pDC")])

# 打印 batch2 中每种细胞类型的数量
# cat("Batch2 cell type counts:\n")
# batch2_counts <- table(mySample[mySample$batch == "Batch2", "celltype"])
# print(batch2_counts[c("CD141", "CD1C", "DoubleNeg", "pDC")])
