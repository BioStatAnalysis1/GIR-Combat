rm(list=ls())
x <- c("ggplot2", "reshape2", "plyr", "moments", "sva", "gridExtra", "Biobase", "RColorBrewer", "BatchQC", "gridExtra")
lapply(x, require, character.only = TRUE)

batch_sum = 3

npeople = 150
npeople_type = npeople

npeople_batch3 = 150
nindex = 100

xmus=c(0,   5,    5)
xsds=c(1,   1,    1)

ymus=c(5,   5,    0)
ysds=c(1,   1,    1)

prop1=c(0.3,  0.5,  0.2)# 批次1细胞类型比例
prop2=c(0.65, 0.3,  0.05)# 批次2细胞类型比例
prop3=c(0.2,  0.2,  0.6)# 批次3细胞类型比例



set.seed(1)
comp1 <- sample(1:3, prob=prop1, size=npeople_type, replace=TRUE)

set.seed(1)
samples1 <- cbind(rnorm(n=npeople_type, mean=xmus[comp1],sd=xsds[comp1]),
                  rnorm(n=npeople_type, mean=ymus[comp1],sd=ysds[comp1]))

set.seed(1)
proj <- matrix(rnorm(nindex*npeople_batch3), nrow=nindex, ncol=2)#100*2

A1 <- samples1 %*% t(proj)

A1_raw <- samples1 %*% t(proj)
rownames(A1_raw) <- paste0("People", seq_len(npeople_type), "-1")
colnames(A1_raw) <- paste0("Index", seq_len(nindex))

set.seed(1)

A1 <- A1 + rnorm(nindex*npeople_type,mean = 0.2, sd = 1)

rownames(A1) <- paste0("People", seq_len(npeople_type), "-1")
colnames(A1) <- paste0("Index", seq_len(nindex))

set.seed(2000)
comp2 <- sample(1:3, prob=prop2, size=npeople_type, replace=TRUE)

set.seed(2000)
samples2 <- cbind(rnorm(n=npeople_type, mean=xmus[comp2], sd=xsds[comp2]),
                  rnorm(n=npeople_type, mean=ymus[comp2], sd=ysds[comp2]))

A2 <- samples2 %*% t(proj)

A2_raw <- samples2 %*% t(proj)
rownames(A2_raw) <- paste0("People", seq_len(npeople_type), "-2")
colnames(A2_raw) <- paste0("Index", seq_len(nindex))

set.seed(2000)

A2 <- A2 + matrix(rep(rnorm(nindex), each=npeople_type), ncol=nindex) 
A2 <- A2 + rnorm(nindex*npeople_type,mean = 0.3, sd = 1) # noise
rownames(A2) <- paste0("People", seq_len(npeople_type), "-2")
colnames(A2) <- paste0("Index", seq_len(nindex))

set.seed(34343)
comp3 <- sample(1:3, prob=prop3, size=npeople_batch3, replace=TRUE)

set.seed(33434)
samples3 <- cbind(rnorm(n=npeople_batch3, mean=xmus[comp3],sd=xsds[comp3]),
                  rnorm(n=npeople_batch3, mean=ymus[comp3],sd=ysds[comp3]))#产生2000个细胞 1000*2


A3 <- samples3 %*% t(proj)

A3_raw <- samples3 %*% t(proj)
rownames(A3_raw) <- paste0("People", seq_len(npeople_batch3), "-3")
colnames(A3_raw) <- paste0("Index", seq_len(nindex))

set.seed(33434)

A3 <- A3 + matrix(rep(rnorm(nindex), each=npeople_batch3), ncol=nindex) 
A3 <- A3 + rnorm(nindex*npeople_batch3,mean = 0.4, sd = 1) # noise
rownames(A3) <- paste0("People", seq_len(npeople_batch3), "-3")
colnames(A3) <- paste0("Index", seq_len(nindex))

####################
B1 <- t(A1)#100(指标)*150(人)
B2 <- t(A2)#100(指标)*150(人)
B3 <- t(A3)#100(指标)*70(人)
B1_raw <- t(A1_raw)#100(指标)*150(人)
B2_raw <- t(A2_raw)#100(指标)*150(人)
B3_raw <- t(A3_raw)#100(指标)*70(人)

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization


comp1 # B1的细胞类型
comp2 # B2的细胞类型
comp3 # B3的细胞类型

data_raw <- cbind(B1_raw,B2_raw,B3_raw)# 100(基因)*450(细胞)

data <- cbind(B1,B2,B3)# 100(基因)*450(细胞)

repeat_measure <- c(rep(1,150),rep(2,150),rep(3,150))# 相当于原来的batch
people_kind <- as.numeric(c(comp1,comp2,comp3))
people_kind_factor <- as.factor(people_kind)

old_pheno <-data.frame(num = 1:ncol(data), measure=repeat_measure, people_kind=people_kind_factor)

combat_data= ComBat(dat=data, batch=old_pheno$measure, mod=NULL, par.prior=TRUE)

addmod_real_8 = model.matrix(~as.factor(people_kind),data=old_pheno)


combat_edata_addmod_real_noref = ComBat(dat=data, batch=old_pheno$measure, mod=addmod_real_8, par.prior=TRUE)
