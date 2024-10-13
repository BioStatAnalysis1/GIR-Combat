rm(list=ls())
require(devtools)

require(exploBATCH)
require(exploBATCHcolon)   
data(Colon) # Colon 110 细胞*20534 基因                       
data(batchColon)            
data(bioCL)   


library("RColorBrewer")
library("Matrix")
library("ggplot2")
library(RISC)
library(data.table)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
X<-Colon #  110 细胞*20534 基因  
rownames(X)<-c(1:nrow(X))
pheno<-data.frame(sample=c(1:length(batchColon)),
                  batch=batchColon,
                  cancer=bioCL
)


X_middle_batch2_2 = X[53:86,]
cell_name<-rownames(X_middle_batch2_2)
set.seed(1)
cell_name_new <- sample(cell_name)
ord <- as.vector(cell_name_new) 
indx <- factor(rownames(X_middle_batch2_2), levels = ord) 
X_middle_batch2_2<-X_middle_batch2_2[order(indx),]

pheno_middle_batch2_2 = pheno[53:86,]
cell_namepheno_middle_batch2_2<-rownames(pheno_middle_batch2_2)
set.seed(1)
cell_name_newpheno_middle_batch2_2 <- sample(cell_namepheno_middle_batch2_2)
ord1 <- as.vector(cell_name_newpheno_middle_batch2_2) 
phenoindx <- factor(rownames(pheno_middle_batch2_2), levels = ord1) 
pheno_middle_batch2_2<-pheno_middle_batch2_2[order(phenoindx),]
#####


pheno_part1 =pheno[(1:52),] 
pheno_part2 =pheno[(87:110),]
pheno = rbind(pheno_part1,pheno_middle_batch2_2)
pheno = rbind(pheno,pheno_part2) 

pheno_part1 =pheno[(1:52),] 
pheno_part2 =pheno[(87:110),]
pheno = rbind(pheno_part1,pheno_middle_batch2_2)
pheno = rbind(pheno,pheno_part2) 

X<-X[-(58:61),]
X<-X[-(75:80),]
pheno<-pheno[-(58:61),]
pheno<-pheno[-(75:80),]
#####3
# batch1
r.datah1_final <- t(X[which(pheno$batch =="1"),]) 
# colnames(r.datah1_final) <- colnames(r.datah1_final) 
unc.merge_batch1 <- pheno[which(pheno$batch =="1"),] 
table(unc.merge_batch1$cancer)



result1_1=NULL
for(i in ((unc.merge_batch1[which(unc.merge_batch1$cancer=="2"),])$sample) ) {

  result11111111 <- r.datah1_final[,which(colnames(r.datah1_final)==i)]
  print(which(colnames(r.datah1_final)==i))
  result1_1 = cbind(result1_1,result11111111)
}
sd(result1_1)  # 批次1 细胞类型1方差


########################
# batch2

r.datah2_final <- t(X[which(pheno$batch =="2"),]) 

unc.merge_batch2 <- pheno[which(pheno$batch =="2"),] 
table(unc.merge_batch2$cancer)



middle_r.datah2_final<-r.datah2_final
# batch2中的细胞类型1
result2_1=NULL
for(i in ((unc.merge_batch2[which(unc.merge_batch2$cancer=="1"),])$sample) ) {
  # for(i in test2) {            # Head of for-loop
  # i="87"
  result11111111 <- r.datah2_final[,which(colnames(r.datah2_final)==i)]
  print(which(colnames(r.datah2_final)==i))
  result2_1 = cbind(result2_1,result11111111)
}
sd(result2_1) 

result2_2=NULL
for(i in ((unc.merge_batch2[which(unc.merge_batch2$cancer=="2"),])$sample) ) {
 
  result11111111 <- r.datah2_final[,which(colnames(r.datah2_final)==i)]
  print(which(colnames(r.datah2_final)==i))
  result2_2 = cbind(result2_2,result11111111)
}
sd(result2_2)  

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(sva)


data <- t(X)# 20534   (基因)*100(细胞)

mySample <-pheno


repeat_measure <- as.numeric(as.factor(mySample$batch)) 
people_kind <- as.numeric(as.factor(mySample$cancer))
people_kind_factor <- as.factor(people_kind)
old_pheno <-data.frame(num = 1:ncol(data), measure=repeat_measure, people_kind=people_kind_factor)

combat_data= ComBat(dat=data, batch=old_pheno$measure, mod=NULL, par.prior=TRUE)

addmod_real_8 = model.matrix(~as.factor(people_kind),data=old_pheno)#  704(细胞)*3("(Intercept)"  "as.factor(people_kind)2" "as.factor(people_kind)3")


combat_edata_addmod_real_noref = ComBat(dat=data, batch=old_pheno$measure, mod=addmod_real_8, par.prior=TRUE)

# 计算每个批次中 normal 和 cancer 的数量
cat("Batch 1 counts:\n")
print(table(pheno[pheno$batch == "1", "cancer"]))

cat("Batch 2 counts:\n")
print(table(pheno[pheno$batch == "2", "cancer"]))
