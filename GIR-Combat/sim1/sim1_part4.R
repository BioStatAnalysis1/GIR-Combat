Xmnn = SummarizedExperiment::assays(Xmnn)
Xmnn = Xmnn$'corrected'
#*********************************************4 first_correct_matrix 
first_correct_matrix <-cbind(data, Xmnn)
dim(first_correct_matrix) #100(指标)*600(人)
#*********************************************6 last_combat_edata 
new_pheno <- rbind(old_pheno, old_pheno)# 4000 *3
new_pheno[(nrow(old_pheno)+1):nrow(new_pheno),]$measure <- (batch_sum + 1)
combat_edata33 = ComBat(dat=first_correct_matrix, batch=new_pheno$measure, mod=NULL, par.prior=TRUE, ref.batch=(batch_sum + 1))


last_combat_edata<- combat_edata33[,1:nrow(old_pheno)]

data_for_kmeans <- t(data)# 300(人)*100(指标)
df <- scale(data_for_kmeans)# 标准化数据

set.seed(1000002324)
# 球形kmeans聚类
library("skmeans")
library("slam")
data_for_SphericalKmeans <- t(data)# 300(人)*100(指标)
data_for_SphericalKmeans_simplematrix <- as.simple_triplet_matrix(data_for_SphericalKmeans)
hparty <- skmeans(data_for_SphericalKmeans_simplematrix, 3, control = list(verbose = TRUE)) 
# save(hparty,file = "hparty.RData")
# load(file = "hparty.RData")
skmeans_comp_cluster <- hparty$cluster
comp_cluster <- as.numeric(skmeans_comp_cluster)

##########
pheno_addmod <- cbind(new_pheno,comp_cluster)
addmod <- model.matrix(~as.factor(comp_cluster),data=pheno_addmod)

combat_edata_addmod = ComBat(dat=first_correct_matrix,
                             batch=pheno_addmod$measure,
                             mod=addmod, par.prior=TRUE, ref.batch=(batch_sum + 1))


last_combat_edata_addmod<- combat_edata_addmod[,1:nrow(old_pheno)]
#*********************************************
addmod_real <- model.matrix(~as.factor(people_kind),data=new_pheno)

combat_edata_addmod_real = ComBat(dat=first_correct_matrix, batch=new_pheno$measure, mod=addmod_real, par.prior=TRUE, ref.batch=(batch_sum + 1))

last_combat_edata_addmod_real<- combat_edata_addmod_real[,1:nrow(old_pheno)]
#*********************************************
# combat_edata_addmod_noref
old_pheno_9 = cbind(old_pheno,comp_cluster)# old_pheno_9 704(细胞)*4("num" "measure" "people_kind" "comp_cluster")
addmod_9 <- model.matrix(~as.factor(comp_cluster),data=old_pheno_9)# addmod_9  704(细胞)* 3("(Intercept)" "as.factor(comp_cluster)2" "as.factor(comp_cluster)3")

combat_edata_addmod_noref = ComBat(dat=data, batch=old_pheno$measure, mod=addmod_9, par.prior=TRUE)
# 
# ######################################data_raw
###################################################################

Scatter_Density <- function(data = data, batch = batch, trt = NULL,
                            batch.legend.title = 'Batch',
                            trt.legend.title = 'Treatment', density.lwd = 0.2,
                            title = people_kind, title.cex = 10, legend.cex = 2.5, legend.title.cex =2.5){
  
  color_value = brewer.pal(8,'Set1')
  data = as.data.frame(data)
  xlim = c(min(data[ ,1])*1.1,max(data[ ,1])*1.1)
  ylim = c(min(data[ ,2])*1.1,max(data[ ,2])*1.1)
  # xlim = c(-10,13)
  # ylim = c(-10,13)
  batch = as.factor(batch)
  trt = as.factor(trt)
  if(nlevels(trt) >= 2){
    pMain <- ggplot(data = data, aes(x = data[ ,1], y = data[ ,2], colour = batch, shape = trt)) +
      geom_point(size= 4) + 
      xlab(paste0('UMAP1')) +
      ylab(paste0('UMAP2')) +
      # axisLabSize = 16 +
      scale_color_manual(values = color_value) + theme_bw() + xlim(xlim[1], xlim[2]) +
      ylim(ylim[1], ylim[2]) + labs(colour = batch.legend.title, shape = trt.legend.title)
    # +coord_fixed(ratio = 0.5)# 横纵比
    
  }else{
    pMain <- ggplot(data = data, aes(x = data[ ,1], y=data[ ,2], colour = batch)) +
      geom_point(shape = 16) + xlab(paste0('UMAP1: ')) +
      ylab(paste0('UMAP2: ')) +
      scale_color_manual(values = color_value) + theme_bw() + xlim(xlim[1], xlim[2]) +
      ylim(ylim[1], ylim[2]) + labs(colour = batch.legend.title)
    
  }
  
  pTop <- labs(title = title)
  g <- ggplotGrob(pMain + theme(legend.position = 'right', legend.box = 'vertical',
                                legend.direction = 'vertical',
                                legend.key.height = unit(0.2, 'cm'),
                                legend.key.width = unit(0.1, 'cm'),
                                legend.title = element_text(size = rel(legend.title.cex)),
                                legend.spacing.x = unit(0.1, 'cm'),
                                legend.spacing.y = unit(0.1, 'cm'),
                                legend.margin = margin(t = 0, r = 0, b = 0, l = 20, unit = "pt"),
                                legend.background = element_rect(fill = "transparent"),
                                legend.text = element_text(size = rel(legend.cex))))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  
  grid.arrange(pMain + theme(legend.position = 'none',axis.title.x = element_text(size = 23),axis.title.y = element_text(size = 23),plot.margin = margin(2, 0, 0, 2, "cm")) +# margin 上右下左
                 pTop+ theme(plot.title = element_text(hjust = 0.5,size = 30)),legend,
               ncol = 2, nrow = 1, widths = c(8, 1))
  
}


#----------------------------------------------
library("FactoMineR")
library(ggplot2)

pca_mat = PCA(as.data.frame(t(data_raw)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比


plot1 <- Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure,
                         trt = people_kind,
                         batch.legend.title = 'batch',
                         trt.legend.title = 'label',
                         title = '(a) No Batch Effect')



######################################
###################################################################改
#----------------------------------------------
set.seed(0)
pca_mat = PCA(as.data.frame(t(data)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot2 <-Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure,
                        trt = people_kind,
                        batch.legend.title = 'batch',
                        trt.legend.title = 'label',
                        title = '(b) Raw')

#------------------------------------------------------
#----------------------------------------------
pca_mat = PCA(as.data.frame(t(combat_data)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot3 <- Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure,
                         trt = people_kind,
                         batch.legend.title = 'batch',
                         trt.legend.title = 'label',
                         title = '(c) Combat')

pca_mat = PCA(as.data.frame(t(combat_edata_addmod_real_noref)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot9 <-Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure,
                        trt = people_kind,
                        batch.legend.title = 'batch',
                        trt.legend.title = 'label',
                        title = '9无参考②有mod①(真实)仅combat校正过的矩阵UMAP图')


pca_mat = PCA(as.data.frame(t(Xmnn)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot4 <-Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure,
                        trt = people_kind,
                        batch.legend.title = 'batch',
                        trt.legend.title = 'label',
                        title = '(d) Mnn')


pca_mat = PCA(as.data.frame(t(first_correct_matrix)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot5 <-Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = new_pheno$measure,
                        trt = new_pheno$people_kind,
                        batch.legend.title = 'batch',
                        trt.legend.title = 'label',
                        title = '5 UMAP图')


pca_mat = PCA(as.data.frame(t(last_combat_edata)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot6 <-Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure,
                        trt = people_kind,
                        batch.legend.title = 'batch',
                        trt.legend.title = 'label',
                        title = '6 做参考的combat校正过的矩阵 UMAP图')

pca_mat = PCA(as.data.frame(t(last_combat_edata_addmod)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot7 <-Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure,
                        trt = people_kind,
                        batch.legend.title = 'batch',
                        trt.legend.title = 'label',
                        title = '(e) GIR-Combat')


#------------------------------------------------------
#----------------------------------------------
pca_mat = PCA(as.data.frame(t(last_combat_edata_addmod_real)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot8 <-Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure,
                        trt = people_kind,
                        batch.legend.title = 'batch',
                        trt.legend.title = 'label',
                        title = '8 有参考有mod的combat校正过的矩阵(真实当mod) UMAP图')


pca_mat = PCA(as.data.frame(t(combat_edata_addmod_noref)),ncp = 5,graph = FALSE)
eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot10 <-Scatter_Density(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure,
                         trt = people_kind,
                         batch.legend.title = 'batch',
                         trt.legend.title = 'label',
                         title = '10  无参考 ②有mod(聚的)①仅combat校正过的矩阵UMAP图')



# ######################
# png(
#   filename = "name1.png", # 文件名称 # 原来的测试，每张图片的名字还是中文
#   width = 1500,           # 宽
#   height = 1000,          # 高
#   units = "px",          # 单位
#   bg = "white",          # 背景颜色
#   res = 100)              # 分辨率
# # 2. 绘图
# grid.arrange(plot1,
#              plot2,
#              plot3,
#              plot4,
#              plot7,
#              ncol=3,nrow=2)
# # 3. 关闭画布
# dev.off()

######################
png(
  filename = "experment1.png", # 文件名称 # 原来的测试，每张图片的名字还是中文
  width = 2700,           # 宽
  height = 2250,          # 高
  units = "px",          # 单位
  bg = "white",          # 背景颜色
  res = 100)              # 分辨率
# 2. 绘图
grid.arrange(plot1,
             plot2,
             plot3,
             plot4,
             plot7,
             ncol=2,
             nrow=3)
# 3. 关闭画布
dev.off()
#------------------------------------------------
#用指标评估校正的前奏
data            
combat_data       
Xmnn              
last_combat_edata  
###################################################################改
last_combat_edata_addmod
last_combat_edata_addmod_real  
###################################################################改

uncorrected <- data
Xcom <- combat_data
MNN_correct_t <- Xmnn
RefMnn_t <- last_combat_edata

people_kind 
old_pheno$measure
cluster.id <- people_kind
batch.id <- old_pheno$measure
MNN_correct <- t(MNN_correct_t)
comp_cluster 
###################################################################

uncorrected                     
Xcom                            
combat_edata_addmod_real_noref  
MNN_correct_t = t(MNN_correct)
RefMnn_t = last_combat_edata    
###################################################################改
last_combat_edata_addmod        
last_combat_edata_addmod_real  
combat_edata_addmod_noref       
###################################################################改


data_raw_t <- t(data_raw) # 1 # 300(人)*100(指标)
sponge.lm.data_raw <- lm(data_raw_t[ ,9] ~  comp_cluster + batch.id)# 第9个指标
summary(sponge.lm.data_raw)

uncorrected_t <- t(uncorrected) # 2 # 300(人)*100(指标)
sponge.lm.before <- lm(uncorrected_t[ ,9] ~  comp_cluster + batch.id)# 第9个指标
summary(sponge.lm.before)

Xcom_t<- t(Xcom) # 3 # 300(人)*100(指标)
sponge.lm.combat <- lm(Xcom_t[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.combat)

combat_edata_addmod_real_noref_t<- t(combat_edata_addmod_real_noref) # 9 # 704(人)*1045(指标)
sponge.lm.combat_edata_addmod_real_noref <- lm(combat_edata_addmod_real_noref_t[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.combat_edata_addmod_real_noref)

MNN_correct # 4 # 300(人)*50(主成分)
sponge.lm.MNN <- lm(MNN_correct[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.MNN)

RefMnn<- t(RefMnn_t) # 6 # 300(人)*100(指标)
sponge.lm.RefMnn <- lm(RefMnn[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.RefMnn)

last_combat_edata_addmod_t<- t(last_combat_edata_addmod) # 7 # 300(人)*100(指标)
sponge.lm.last_combat_edata_addmod <- lm(last_combat_edata_addmod_t[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.last_combat_edata_addmod)

last_combat_edata_addmod_real_t<- t(last_combat_edata_addmod_real) # 8 # 300(人)*100(指标)
sponge.lm.last_combat_edata_addmod_real <- lm(last_combat_edata_addmod_real_t[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.last_combat_edata_addmod_real)

combat_edata_addmod_noref_t<- t(combat_edata_addmod_noref) # 10 # 704(人)*1045(指标)
sponge.lm.combat_edata_addmod_noref <- lm(combat_edata_addmod_noref_t[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.combat_edata_addmod_noref)

# Sponge data
sponge.form <- ~  people_kind +  batch.id
sponge.info <- as.data.frame(cbind(rownames(uncorrected_t), people_kind,batch.id))
rownames(sponge.info) <- rownames(uncorrected_t)

# before
library(variancePartition)
sponge.varPart.data_raw <- fitExtractVarPartModel(exprObj = t(data_raw_t), # 1
                                                  formula = sponge.form,
                                                  data = sponge.info)

sponge.varPart.before <- fitExtractVarPartModel(exprObj = t(uncorrected_t), # 2
                                                formula = sponge.form,
                                                data = sponge.info)
# combat
sponge.varPart.combat <- fitExtractVarPartModel(exprObj = t(Xcom_t), # 3
                                                formula = sponge.form,
                                                data = sponge.info)
# combat
sponge.varPart.combat_edata_addmod_real_noref <- fitExtractVarPartModel(exprObj = t(combat_edata_addmod_real_noref_t), # 9
                                                                        formula = sponge.form,
                                                                        data = sponge.info)

#mnn
sponge.varPart.mnn <- fitExtractVarPartModel(exprObj = t(MNN_correct), # 4
                                             formula = sponge.form,
                                             data = sponge.info)


sponge.varPart.RefMnn <- fitExtractVarPartModel(exprObj =  t(RefMnn), # 6
                                                formula = sponge.form,
                                                data = sponge.info)

sponge.varPart.last_combat_edata_addmod <- fitExtractVarPartModel(exprObj = t(last_combat_edata_addmod_t), # 7
                                                                  formula = sponge.form,
                                                                  data = sponge.info)

sponge.varPart.last_combat_edata_addmod_real <- fitExtractVarPartModel(exprObj = t(last_combat_edata_addmod_real_t), # 8
                                                                       formula = sponge.form,
                                                                       data = sponge.info)

sponge.varPart.combat_edata_addmod_noref <- fitExtractVarPartModel(exprObj = t(combat_edata_addmod_noref_t), # 10
                                                                   formula = sponge.form,
                                                                   data = sponge.info)

# extract the variance of trt and batch
# before
sponge.varmat.data_raw <- as.matrix(sponge.varPart.data_raw[ ,1:2]) # 1

# before
sponge.varmat.before <- as.matrix(sponge.varPart.before[ ,1:2]) # 2
# ComBat
sponge.varmat.combat <- as.matrix(sponge.varPart.combat[ ,1:2]) # 3

sponge.varmat.combat_edata_addmod_real_noref <- as.matrix(sponge.varPart.combat_edata_addmod_real_noref[ ,1:2]) # 9
# mnn
sponge.varmat.mnn <- as.matrix(sponge.varPart.mnn[ ,1:2]) # 4

sponge.varmat.RefMnn <- as.matrix(sponge.varPart.RefMnn[ ,1:2]) # 6

sponge.varmat.last_combat_edata_addmod <- as.matrix(sponge.varPart.last_combat_edata_addmod[ ,1:2]) # 7

sponge.varmat.last_combat_edata_addmod_real <- as.matrix(sponge.varPart.last_combat_edata_addmod_real[ ,1:2]) # 8

sponge.varPart.combat_edata_addmod_noref <- as.matrix(sponge.varPart.combat_edata_addmod_noref[ ,1:2]) # 10

sponge.variance <- c(as.vector(sponge.varmat.data_raw), # 1
                     as.vector(sponge.varmat.before), # 2
                     as.vector(sponge.varmat.combat), # 3
                     as.vector(sponge.varmat.mnn), # 4
                     as.vector(sponge.varmat.last_combat_edata_addmod) # 7
)# 10
###############################

sponge.variance <- cbind(variance = sponge.variance,
                         Type = rep(c('Label','Batch'), each = ncol(uncorrected_t)),
                         method = rep(c('No Batch Effect',
                                        'Raw',
                                        'ComBat',
                                        # 'ComBat havemod_real',
                                        'Mnn',
                                        # 'RefMnn',
                                        'GIR-Combat'
                                        # 'last_combat_edata_addmod_real',
                                        # 'Orign ComBat havemod_cluster'
                         ), each = 2*ncol(uncorrected_t)))

###############################

# reorder levels
sponge.variance <- as.data.frame(sponge.variance)
sponge.variance$method <- factor(sponge.variance$method,
                                 levels = unique(sponge.variance$method))
sponge.variance$variance <- as.numeric(as.character(sponge.variance$variance))

png(
  filename = "experment11.png", # 文件名称 # 原来的测试，每张图片的名字还是中文
  width = 1875,           # 宽
  height = 1000,          # 高
  units = "px",          # 单位
  bg = "white",          # 背景颜色
  res = 100)              # 分辨率

ggplot(sponge.variance, aes(x = Type, y = variance, fill= Type)) +
  geom_boxplot() + facet_grid(cols = vars(method)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 25),
        strip.text = element_text(size = 30), panel.grid = element_blank(),# strip.text 标题的字
        axis.text = element_text(size = 25), axis.title = element_text(size = 30),
        legend.title = element_text(size = 25), legend.text = element_text(size = 25)) +
  labs(x = 'Type', y = 'Proportion Variance', name = 'Type') + ylim(0,1)+scale_fill_brewer(palette="Blues")+ geom_jitter(aes(fill=Type),width =0.2,shape = 21,size=1)
# 3. 关闭画布
dev.off()



data_raw                      
Xcom                            
combat_edata_addmod_real_noref  
MNN_correct_t = t(MNN_correct)  
RefMnn_t = last_combat_edata    
last_combat_edata_addmod       
last_combat_edata_addmod_real   
combat_edata_addmod_noref       

library(FactoMineR)
library(factoextra)

################################################################################
################################################################################

################################################################################
################################################################################
library(lisi)#原始数据与整合后的数据的UMAP值输入到LISI中

lisi_score <- function(Matrix0,Batch_type,comp_cluster){
 
  BC <- cbind(Batch_type,comp_cluster)
  BC = as.data.frame(BC)
  colnames(BC) = c("Batch_type","Cell_type")
  BC$Cell_type = as.factor(BC$Cell_type)
  BC$Batch_type = as.factor(BC$Batch_type)
  F1 = 0
  for(i in 1:20){
   
    extract_samples = sample( c(1:nrow(Matrix0)),size = 0.8 * nrow(Matrix0) ,replace = FALSE )
    
    lisi0 <- compute_lisi(X = Matrix0[extract_samples,],meta_data = BC[extract_samples,],label_colnames = c("Batch_type","Cell_type"))
    
    score_Cell = mean(lisi0$Cell_type)/nlevels(BC$Cell_type)
    score_batch = mean(lisi0$Batch_type)/nlevels(BC$Batch_type)
    F1 = F1 + 2*(1-score_Cell)*score_batch/(1-score_Cell+score_batch)
  }
  return(F1/20)
}
library(FactoMineR)
library(factoextra)

dat_raw = as.data.frame(t(data_raw))
dat.pca_raw <- PCA(dat_raw,ncp = 20,graph = F)
umap_raw = umap::umap(dat.pca_raw$ind$coord)$layout
E = lisi_score(umap_raw ,batch.id,comp_cluster)
E #  t(data_raw)

dat_uncorrected = as.data.frame(t(uncorrected))
dat.pca_uncorrected <- PCA(dat_uncorrected,ncp = 20,graph = F)
umap_uncorrected = umap::umap(dat.pca_uncorrected$ind$coord)$layout

B = lisi_score(umap_uncorrected ,batch.id,people_kind)


dat_Xcom = as.data.frame(t(Xcom))
dat.pca_Xcom <- PCA(dat_Xcom,ncp = 20,graph = F)
umap_Xcom = umap::umap(dat.pca_Xcom$ind$coord)$layout
C = lisi_score(umap_Xcom ,batch.id,people_kind)


dat_combat_edata_addmod_real_noref = as.data.frame(t(combat_edata_addmod_real_noref))
dat.pca_combat_edata_addmod_real_noref <- PCA(dat_combat_edata_addmod_real_noref,ncp = 20,graph = F)
umap_combat_edata_addmod_real_noref = umap::umap(dat.pca_combat_edata_addmod_real_noref$ind$coord)$layout

Q = lisi_score(umap_combat_edata_addmod_real_noref ,batch.id,people_kind)
Q # 9 #  t(combat_edata_addmod_real_noref)

dat_MNN_correct_t = as.data.frame(t(MNN_correct_t)) 
dat.pca_MNN_correct_t <- PCA(dat_MNN_correct_t,ncp = 20,graph = F)
umap_MNN_correct_t = umap::umap(dat.pca_MNN_correct_t$ind$coord)$layout

D = lisi_score(umap_MNN_correct_t ,batch.id,people_kind)
D # 4 # t(MNN_correct_t)

dat_RefMnn = as.data.frame(t(RefMnn_t))
dat.pca_RefMnn <- PCA(dat_RefMnn,ncp = 20,graph = F)
umap_RefMnn = umap::umap(dat.pca_RefMnn$ind$coord)$layout

A = lisi_score(umap_RefMnn ,batch.id,people_kind)
A # 6 # t(RefMnn_t)

dat_last_combat_edata_addmod = as.data.frame(t(last_combat_edata_addmod))
dat.pca_last_combat_edata_addmod <- PCA(dat_last_combat_edata_addmod,ncp = 20,graph = F)
umap_last_combat_edata_addmod = umap::umap(dat.pca_last_combat_edata_addmod$ind$coord)$layout

G = lisi_score(umap_last_combat_edata_addmod ,batch.id,people_kind)
G # 7 #  t(last_combat_edata_addmod)

dat_last_combat_edata_addmod_real = as.data.frame(t(last_combat_edata_addmod_real))
dat.pca_last_combat_edata_addmod_real <- PCA(dat_last_combat_edata_addmod_real,ncp = 20,graph = F)
umap_last_combat_edata_addmod_real = umap::umap(dat.pca_last_combat_edata_addmod_real$ind$coord)$layout

H = lisi_score(umap_last_combat_edata_addmod_real ,batch.id,people_kind)
H # 8 #  t(last_combat_edata_addmod_real)

dat_combat_edata_addmod_noref = as.data.frame(t(combat_edata_addmod_noref))
dat.pca_dat_combat_edata_addmod_noref <- PCA(dat_combat_edata_addmod_noref,ncp = 20,graph = F)
umap_combat_edata_addmod_noref = umap::umap(dat.pca_dat_combat_edata_addmod_noref$ind$coord)$layout

W = lisi_score(umap_combat_edata_addmod_noref ,batch.id,people_kind)
W # 10  t(combat_edata_addmod_noref)
################################################################################
################################################################################

################################################################################
################################################################################
library(cluster)

asw_score <- function(Matrix0,Batch_type,comp_cluster){
  # Matrix0 <- dat.pca_uncorrected$ind$coord # 300(人)*20(主成分)
  # Batch_type <- batch.id
  
  
  BC <- cbind(Batch_type,comp_cluster)
  BC = as.data.frame(BC)
  colnames(BC) = c("Batch_type","Cell_type")
  BC$Cell_type = as.factor(BC$Cell_type)
  BC$Batch_type = as.factor(BC$Batch_type)
  F1 = 0
  F1_sum =0
  for( i in 1:20){
    # i=1
    extract_samples = sample(c(1:nrow(Matrix0)),size = 0.8 * nrow(Matrix0) ,replace = FALSE)
    dist.obs = dist(Matrix0[extract_samples,])
    
    sih.Cell <- silhouette(as.integer(BC$Cell_type[extract_samples]), dist = dist.obs)
    sih.Cell = as.numeric(summary(sih.Cell)[[2]], 1)
    sih.Cell1 = round(sih.Cell, 2)
    
    sih.batch <- silhouette(as.integer(BC$Batch_type[extract_samples]), dist = dist.obs)
    sih.batch = as.numeric(summary(sih.batch)[[2]], 1)
    sih.batch1 = round(sih.batch, 2)
    
    sih.Cell1 =(sih.Cell1+1)/2 
    sih.batch1 = (sih.batch1+1)/2 
    score_Cell = length(sih.Cell1)/sum(1/sih.Cell1)
    score_batch = length(sih.batch1)/sum(1/sih.batch1)
    
    F1 = 2*(1-score_batch)*(score_Cell)/(1-score_batch+score_Cell)
    cat(F1,"第",i,"次","\n")
    F1_sum = F1_sum + F1
  }
  
  score_F1 = F1_sum/20
  return(score_F1)
}
a11 =asw_score(dat.pca_raw$ind$coord,batch.id,comp_cluster)
a11

a2 =asw_score(dat.pca_uncorrected$ind$coord,batch.id,people_kind) 
a2 

a3 =asw_score(dat.pca_Xcom$ind$coord,batch.id,people_kind)
a3

a9 =asw_score(dat.pca_combat_edata_addmod_real_noref$ind$coord,batch.id,people_kind) 
a9

a4 =asw_score(dat.pca_MNN_correct_t$ind$coord,batch.id,people_kind) 
a4 

a6 =asw_score(dat.pca_RefMnn$ind$coord,batch.id,people_kind) 
a6

a7 =asw_score(dat.pca_last_combat_edata_addmod$ind$coord,batch.id,people_kind) 
a7

a8 =asw_score(dat.pca_last_combat_edata_addmod_real$ind$coord,batch.id,people_kind) 
a8 

a10 =asw_score(dat.pca_dat_combat_edata_addmod_noref$ind$coord,batch.id,people_kind)
a10 

a11 
a2  
a3 
a9
a4 
a6
a7
a8 
a10 

################################################################################
################################################################################

################################################################################

ari_score <- function(Matrix0,Batch_type,comp_cluster){
  # Matrix0 <- dat.pca_uncorrected$ind$coord
  # Batch_type <- batch.id
  # comp_cluster
  
  require(mclust)
  
  ARI_Cell = 0
  M <- droplevels(as.factor(comp_cluster))
  n.M  <- nlevels(M) 
  Ms <- lapply(levels(M), function(x) which(M == x))
  n.Ms <- sapply(Ms, length)
  

  select_samples = c()
  for(j in 1: n.M){
    # j=3
    sample.n = ceiling(0.8 * n.Ms[j])
    extract_samples = sample(Ms[[j]],size = sample.n,replace = FALSE)
    select_samples = c(select_samples,extract_samples)
  }
  select_samples = sort(select_samples)
  kms <- kmeans(Matrix0[select_samples,],n.M)
  ARI_Cell = adjustedRandIndex(kms$cluster,M[select_samples])
  
  
  
  ARI_batch = 0
  batch <- droplevels(as.factor(Batch_type))
  n.batch  <- nlevels(batch)
  batches <- lapply(levels(batch), function(x) which(batch == x))
  n.batches <- sapply(batches, length)
  

  select_samples = c()
  for(j in 1: n.batch){
    # j=2
    sample.n = ceiling(0.8 * n.batches[j])
    extract_samples = sample(batches[[j]],size = sample.n,replace = FALSE)
    select_samples = c(select_samples,extract_samples)
  }
  select_samples = sort(select_samples)
  kms <- kmeans(Matrix0[select_samples,],n.batch)
  ARI_batch = adjustedRandIndex(kms$cluster,batch[select_samples])
  
  
  min_val <- -1
  max_val <- 1
  ARI_batch_norm <- (ARI_batch-min_val)/(max_val-min_val)
  ARI_Cell_norm <- (ARI_Cell-min_val)/(max_val-min_val)
  
  cat(ARI_batch_norm,"score_批次","\n")
  cat(ARI_Cell_norm,"score_细胞","\n")
  # produce final fscore ARI, similar to scMerge paper
  fscoreARI <- (2 * (1 - ARI_batch_norm)*(ARI_Cell_norm))/(1 - ARI_batch_norm + ARI_Cell_norm)
  
  
  

  return(fscoreARI)
}
b11 =ari_score(dat.pca_raw$ind$coord,batch.id,comp_cluster) 
b11

b2 =ari_score(dat.pca_uncorrected$ind$coord,batch.id,people_kind) 
b2 

b3 =ari_score(dat.pca_Xcom$ind$coord,batch.id,people_kind) 
b3 

b9 =ari_score(dat.pca_combat_edata_addmod_real_noref$ind$coord,batch.id,people_kind) 
b9

b4 =ari_score(dat.pca_MNN_correct_t$ind$coord,batch.id,people_kind) 
b4 

b6 =ari_score(dat.pca_RefMnn$ind$coord,batch.id,people_kind)
b6

b7 =ari_score(dat.pca_last_combat_edata_addmod$ind$coord,batch.id,people_kind)
b7 

b8 =ari_score(dat.pca_last_combat_edata_addmod_real$ind$coord,batch.id,people_kind) 
b8

b10 =ari_score(dat.pca_dat_combat_edata_addmod_noref$ind$coord,batch.id,people_kind)
b10 


