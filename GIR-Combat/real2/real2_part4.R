Xmnn = SummarizedExperiment::assays(Xmnn)
Xmnn = Xmnn$'corrected'


first_correct_matrix <-cbind(data, Xmnn)

new_pheno <- rbind(old_pheno, old_pheno)# (34917(细胞)*2) *3("num"  "measure" "people_kind")
new_pheno[(nrow(old_pheno)+1):nrow(new_pheno),]$measure <- 3

combat_edata33 = ComBat(dat=first_correct_matrix, batch=new_pheno$measure, mod=NULL, par.prior=TRUE, ref.batch=3)


last_combat_edata<- combat_edata33[,1:nrow(old_pheno)]

data_for_kmeans <- t(data)# 704(细胞)*1045(基因)
df <- scale(data_for_kmeans)# 标准化数据



library("skmeans")
library("slam")
data_for_SphericalKmeans <- t(data)# 5498(细胞)*7066(基因)
data_for_SphericalKmeans_simplematrix <- as.simple_triplet_matrix(data_for_SphericalKmeans)
set.seed(10000000)

#load(file = "F:/Rcode/benmark/benmark_docker_data/dataset1/20220601_2_20230306.RData")
hparty <- skmeans(data_for_SphericalKmeans_simplematrix, 2, control = list(verbose = TRUE))
skmeans_comp_cluster <- hparty$cluster
comp_cluster <- as.numeric(skmeans_comp_cluster)


pheno_addmod <- cbind(new_pheno,comp_cluster)
addmod <- model.matrix(~as.factor(comp_cluster),data=pheno_addmod)

combat_edata_addmod = ComBat(dat=first_correct_matrix,
                             batch=pheno_addmod$measure,
                             mod=addmod, par.prior=TRUE, ref.batch=3)

last_combat_edata_addmod<- combat_edata_addmod[,1:nrow(old_pheno)]

gc()

addmod_real <- model.matrix(~as.factor(people_kind),data=new_pheno)# 1408(704(细胞)*2) * 3("(Intercept)" "as.factor(people_kind)2" "as.factor(people_kind)3")

combat_edata_addmod_real = ComBat(dat=first_correct_matrix, batch=new_pheno$measure, mod=addmod_real, par.prior=TRUE, ref.batch=3)


last_combat_edata_addmod_real<- combat_edata_addmod_real[,1:nrow(old_pheno)]

gc()

old_pheno_9 = cbind(old_pheno,comp_cluster)
addmod_9 <- model.matrix(~as.factor(comp_cluster),data=old_pheno_9)

combat_edata_addmod_noref = ComBat(dat=data, batch=old_pheno$measure, mod=addmod_9, par.prior=TRUE)


Scatter_Density_batch <- function(data = data, batch = batch, 
                                  batch.legend.title = 'Batch', 
                                  density.lwd = 0.2,
                                  title = people_kind, title.cex = 10, legend.cex = 2.5, legend.title.cex =2.5){
  
  
  color_value = brewer.pal(8,'Set1')
  data = as.data.frame(data)
  xlim = c(min(data[ ,1])*1.2,max(data[ ,1])*1.2)
  ylim = c(min(data[ ,2])*1.2,max(data[ ,2])*1.2)
  batch = as.factor(batch)
  
  pMain <- ggplot(data = data, aes(x = data[ ,1], y = data[ ,2], colour = batch)) + 
    geom_point(size= 4) + xlab(paste0('UMAP1')) + 
    ylab(paste0('UMAP2')) + 
    scale_color_manual(values = color_value) + theme_bw() + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) + labs(colour = batch.legend.title)
  
  
  
  pTop <- labs(title = title)
  g <- ggplotGrob(pMain + theme(legend.position = 'right', legend.box = 'horizontal',
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

Scatter_Density_cell <- function(data = data, cell_type = cell_type, 
                                 cell_type.legend.title = 'cell_type', 
                                 density.lwd = 0.2,
                                 title = people_kind, title.cex = 10, legend.cex = 2.5, legend.title.cex =2.5){
  
  
  color_value = brewer.pal(5, "Set2")
  data = as.data.frame(data)
  xlim = c(min(data[ ,1])*1.2,max(data[ ,1])*1.2)
  ylim = c(min(data[ ,2])*1.2,max(data[ ,2])*1.2)
  cell_type = as.factor(cell_type)
  
  pMain <- ggplot(data = data, aes(x = data[ ,1], y = data[ ,2], colour = cell_type)) + 
    geom_point(size= 4) + xlab(paste0('UMAP1')) + 
    ylab(paste0('UMAP2')) + 
    scale_color_manual(values = color_value) + theme_bw() + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) + labs(colour = cell_type.legend.title)
  
  
  
  pTop <- labs(title = title)
  g <- ggplotGrob(pMain + theme(legend.position = 'right', legend.box = 'horizontal',
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


library(FactoMineR)
library(gridExtra)
library(RColorBrewer)

set.seed(123)
pca_mat2 = PCA(as.data.frame(t(data)),ncp = 5,graph = FALSE)

plot2_batch <- Scatter_Density_batch(data = umap::umap(pca_mat2$ind$coord)$layout[,1:2], batch = old_pheno$measure, 
                                     
                                     batch.legend.title = 'batch', 
                                     
                                     title = '(a1) Raw')
plot2_cell <- Scatter_Density_cell(data = umap::umap(pca_mat2$ind$coord)$layout[,1:2], cell_type = people_kind, 
                                   
                                   cell_type.legend.title = 'cell_type', 
                                   
                                   title = '(a2) Raw')

set.seed(123)
pca_mat3 = PCA(as.data.frame(t(combat_data)),ncp = 5,graph = FALSE)
# eig_df = as.data.frame(pca_mat$eig)#eig_df pca分解的特征值的属性 100*3("eigenvalue" "percentage of variance" "cumulative percentage of variance")
# eig_df$`percentage of variance`# 特征值的方差百分比

plot3_batch <- Scatter_Density_batch(data = umap::umap(pca_mat3$ind$coord)$layout[,1:2], batch = old_pheno$measure, 
                                     
                                     batch.legend.title = 'batch', 
                                     
                                     title = '(b1) Combat')
plot3_cell <- Scatter_Density_cell(data = umap::umap(pca_mat3$ind$coord)$layout[,1:2], cell_type = people_kind, 
                                   
                                   cell_type.legend.title = 'cell_type', 
                                   
                                   title = '(b2) Combat')



set.seed(123)
pca_mat4 = PCA(as.data.frame(t(Xmnn)),ncp = 5,graph = FALSE)

plot4_batch <- Scatter_Density_batch(data = umap::umap(pca_mat4$ind$coord)$layout[,1:2], batch = old_pheno$measure, 
                                     
                                     batch.legend.title = 'batch', 
                                     
                                     title = '(c1) Mnn')
plot4_cell <- Scatter_Density_cell(data = umap::umap(pca_mat4$ind$coord)$layout[,1:2], cell_type = people_kind, 
                                   
                                   cell_type.legend.title = 'cell_type', 
                                   
                                   title = '(c2) Mnn')



set.seed(123)
pca_mat = PCA(as.data.frame(t(last_combat_edata_addmod)),ncp = 5,graph = FALSE)

plot7_batch <- Scatter_Density_batch(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], batch = old_pheno$measure, 
                                     
                                     batch.legend.title = 'batch', 
                                     
                                     title = '(d1) GIR-Combat')

plot7_cell <- Scatter_Density_cell(data = umap::umap(pca_mat$ind$coord)$layout[,1:2], cell_type = people_kind, 
                                   
                                   cell_type.legend.title = 'cell_type', 
                                   
                                   title = '(d2) GIR-Combat')





set.seed(123)
pca_mat8 = PCA(as.data.frame(t(last_combat_edata_addmod_real)),ncp = 5,graph = FALSE)


plot8_batch <- Scatter_Density_batch(data = umap::umap(pca_mat8$ind$coord)$layout[,1:2], batch = old_pheno$measure, 
                                     
                                     batch.legend.title = 'batch', 
                                     
                                     title = '(d) GIR-Combat')
plot8_cell <- Scatter_Density_cell(data = umap::umap(pca_mat8$ind$coord)$layout[,1:2], cell_type = people_kind, 
                                   
                                   cell_type.legend.title = 'cell_type', 
                                   
                                   title = '(d) GIR-Combat')

png( 
  filename = "realdata2.png", # 文件名称 # 原来的测试，每张图片的名字还是中文
  width = 2700,           # 宽
  height = 2667,          # 高
  units = "px",          # 单位
  bg = "white",          # 背景颜色
  res = 100)              # 分辨率

grid.arrange(
  # plot1,
  plot2_batch, plot2_cell,
  plot3_batch, plot3_cell,
  plot4_batch ,plot4_cell,
  plot7_batch ,plot7_cell,
  # plot8_batch ,plot8_cell,
  ncol=2,nrow=4)
# 3. 关闭画布
dev.off()

data                        
combat_data                  
combat_edata_addmod_real_noref 
Xmnn                         
last_combat_edata            
###################################################################改
last_combat_edata_addmod      
last_combat_edata_addmod_real 
combat_edata_addmod_noref     
###################################################################改
uncorrected <- data 
Xcom <- combat_data 
MNN_correct_t <- Xmnn 
RefMnn_t <- last_combat_edata 


gc()

people_kind 
old_pheno$measure
cluster.id <- people_kind
batch.id <- old_pheno$measure
MNN_correct <- t(MNN_correct_t)
comp_cluster 




uncorrected                    
Xcom                           

MNN_correct_t = t(MNN_correct) 
RefMnn_t = last_combat_edata    

last_combat_edata_addmod      
last_combat_edata_addmod_real  

uncorrected_t <- t(uncorrected) # 2 # 704(人)*1045(指标)
sponge.lm.before <- lm(uncorrected_t[ ,9] ~  comp_cluster + batch.id)# 第9个指标
summary(sponge.lm.before)

Xcom_t<- t(Xcom) # 3 # 704(人)*1045(指标)
sponge.lm.combat <- lm(Xcom_t[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.combat)


MNN_correct # 4 # 300(人)*50(主成分)
sponge.lm.MNN <- lm(MNN_correct[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.MNN)



last_combat_edata_addmod_t<- t(last_combat_edata_addmod) # 7 # 704(人)*1045(指标)
sponge.lm.last_combat_edata_addmod <- lm(last_combat_edata_addmod_t[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.last_combat_edata_addmod)

last_combat_edata_addmod_real_t<- t(last_combat_edata_addmod_real) # 8 # 704(人)*1045(指标)
sponge.lm.last_combat_edata_addmod_real <- lm(last_combat_edata_addmod_real_t[ ,9] ~ comp_cluster +  batch.id)
summary(sponge.lm.last_combat_edata_addmod_real)




sponge.form <- ~  people_kind +  batch.id

sponge.info <- as.data.frame(cbind(rownames(uncorrected_t), people_kind, batch.id))
rownames(sponge.info) <- rownames(uncorrected_t)


library(variancePartition)

sponge.varPart.before <- fitExtractVarPartModel(exprObj = t(uncorrected_t), # 2
                                                formula = sponge.form,
                                                data = sponge.info)

sponge.varPart.combat <- fitExtractVarPartModel(exprObj = t(Xcom_t), # 3
                                                formula = sponge.form,
                                                data = sponge.info)


sponge.varPart.mnn <- fitExtractVarPartModel(exprObj = t(MNN_correct), # 4
                                             formula = sponge.form,
                                             data = sponge.info)

sponge.varPart.last_combat_edata_addmod <- fitExtractVarPartModel(exprObj = t(last_combat_edata_addmod_t), # 7
                                                                  formula = sponge.form,
                                                                  data = sponge.info)

sponge.varPart.last_combat_edata_addmod_real <- fitExtractVarPartModel(exprObj = t(last_combat_edata_addmod_real_t), # 8
                                                                       formula = sponge.form,
                                                                       data = sponge.info)

sponge.varmat.before <- as.matrix(sponge.varPart.before[ ,1:2]) # 2

sponge.varmat.combat <- as.matrix(sponge.varPart.combat[ ,1:2]) # 3


sponge.varmat.mnn <- as.matrix(sponge.varPart.mnn[ ,1:2]) # 4

sponge.varmat.last_combat_edata_addmod <- as.matrix(sponge.varPart.last_combat_edata_addmod[ ,1:2]) # 7

sponge.varmat.last_combat_edata_addmod_real <- as.matrix(sponge.varPart.last_combat_edata_addmod_real[ ,1:2]) # 8



sponge.variance <- c(

  as.vector(sponge.varmat.before), # 2
  as.vector(sponge.varmat.combat), # 3
 
  as.vector(sponge.varmat.mnn), # 4

  as.vector(sponge.varmat.last_combat_edata_addmod) # 7
  
)# 10


sponge.variance <- cbind(variance = sponge.variance,
                         Type = rep(c('CellType','Batch'), each = ncol(uncorrected_t)),
                         method = rep(c(
                         
                           'Raw',
                           'ComBat',
                          
                           'Mnn',
                         
                           'GIR-Combat'
       
                         ), each = 2*ncol(uncorrected_t)))


sponge.variance <- as.data.frame(sponge.variance)
sponge.variance$method <- factor(sponge.variance$method,
                                 levels = unique(sponge.variance$method))
sponge.variance$variance <- as.numeric(as.character(sponge.variance$variance))


png( 
  filename = "realdata21.png", # 文件名称 # 原来的测试，每张图片的名字还是中文
  width = 2000,           # 宽
  height = 1000,          # 高
  units = "px",          # 单位
  bg = "white",          # 背景颜色
  res = 100)              # 分辨率
# 2. 绘图

ggplot(sponge.variance, aes(x = Type, y = variance, fill = Type)) +
  geom_boxplot() + facet_grid(cols = vars(method)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 25),
        strip.text = element_text(size = 30), panel.grid = element_blank(),
        axis.text = element_text(size = 25), axis.title = element_text(size = 30),
        legend.title = element_text(size = 25), legend.text = element_text(size = 25)) +
  labs(x = 'Type', y = 'Proportion Variance', name = 'Type') + ylim(0,1)+scale_fill_brewer(palette="Blues")+ geom_jitter(aes(fill=Type),width =0.2,shape = 21,size=1)

dev.off()


uncorrected                   
Xcom                        

MNN_correct_t = t(MNN_correct) 

last_combat_edata_addmod       
last_combat_edata_addmod_real   
     

library(FactoMineR)
library(factoextra)



library(kBET)

kBET_score <- function(Matrix0, Batch_type ){

  Matrix0 = as.matrix(Matrix0)
  n_repeat = 10
  rejection = kBET(df = Matrix0, batch = batch.id,do.pca = F, addTest = F,  n_repeat = n_repeat, verbose = F)
  
  plot.data <- data.frame(class=rep(c('observed', 'expected'),
                                    each=length(rejection$stats$kBET.observed)),
                          data =  c(rejection$stats$kBET.observed,
                                    rejection$stats$kBET.expected))
  

  avg_score <- sum(rejection$stats$kBET.observed)/n_repeat
  avg_score
}

B = kBET_score(t(uncorrected),batch.id)# 2  t(uncorrected) 57(样本)* 100(指标)
B

C = kBET_score(t(Xcom),batch.id)# 3 #  t(Xcom) 57(样本)* 100(指标)
C



D = kBET_score(t(MNN_correct_t),batch.id)# 4 # t(MNN_correct_t) 57(样本)* 50(主成分)
D


G = kBET_score(t(last_combat_edata_addmod),batch.id)# 7  t(last_combat_edata_addmod) 57(样本)* 100(指标)
G

H = kBET_score(t(last_combat_edata_addmod_real),batch.id)# 8  t(last_combat_edata_addmod_real) 57(样本)* 100(指标)
H


library(lisi)



lisi_score <- function(Matrix0,Batch_type,comp_cluster){
  # Matrix0 <- umap_combat
  # Batch_type <-  batch.id
  BC <- cbind(Batch_type,comp_cluster)
  BC = as.data.frame(BC)
  colnames(BC) = c("Batch_type","Cell_type")
  BC$Cell_type = as.factor(BC$Cell_type)
  BC$Batch_type = as.factor(BC$Batch_type)
  F1 = 0
  for(i in 1:20){
    # i=1
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


dat_uncorrected = as.data.frame(t(uncorrected))# 2 #  t(uncorrected) 57(样本)* 100(指标)
dat.pca_uncorrected <- PCA(dat_uncorrected,ncp = 20,graph = F)
umap_uncorrected = umap::umap(dat.pca_uncorrected$ind$coord)$layout

B = lisi_score(umap_uncorrected ,batch.id,people_kind)
B

dat_Xcom = as.data.frame(t(Xcom))# 3 #  t(Xcom) 57(样本)* 100(指标)
dat.pca_Xcom <- PCA(dat_Xcom,ncp = 20,graph = F)
umap_Xcom = umap::umap(dat.pca_Xcom$ind$coord)$layout

C = lisi_score(umap_Xcom ,batch.id,people_kind)
C



dat_MNN_correct_t = as.data.frame(t(MNN_correct_t)) # 4 # t(MNN_correct_t) 57(样本)* 50(主成分)
dat.pca_MNN_correct_t <- PCA(dat_MNN_correct_t,ncp = 20,graph = F)
umap_MNN_correct_t = umap::umap(dat.pca_MNN_correct_t$ind$coord)$layout

D = lisi_score(umap_MNN_correct_t ,batch.id,people_kind)
D


dat_last_combat_edata_addmod = as.data.frame(t(last_combat_edata_addmod))# 7 #  t(last_combat_edata_addmod) 57(样本)* 100(指标)
dat.pca_last_combat_edata_addmod <- PCA(dat_last_combat_edata_addmod,ncp = 20,graph = F)
umap_last_combat_edata_addmod = umap::umap(dat.pca_last_combat_edata_addmod$ind$coord)$layout

G = lisi_score(umap_last_combat_edata_addmod ,batch.id,people_kind)
G

dat_last_combat_edata_addmod_real = as.data.frame(t(last_combat_edata_addmod_real))# 8 #  t(last_combat_edata_addmod_real) 57(样本)* 100(指标)
dat.pca_last_combat_edata_addmod_real <- PCA(dat_last_combat_edata_addmod_real,ncp = 20,graph = F)
umap_last_combat_edata_addmod_real = umap::umap(dat.pca_last_combat_edata_addmod_real$ind$coord)$layout

H = lisi_score(umap_last_combat_edata_addmod_real ,batch.id,people_kind)
H



library(cluster)



asw_score <- function(Matrix0,Batch_type,comp_cluster){

  
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
    dist.obs = dist(Matrix0[extract_samples,])# Matrix0[extract_samples,]的维数：240*20
    
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


a2 =asw_score(dat.pca_uncorrected$ind$coord,batch.id,people_kind) 
a2  

a3 =asw_score(dat.pca_Xcom$ind$coord,batch.id,people_kind) 
a3 


a4 =asw_score(dat.pca_MNN_correct_t$ind$coord,batch.id,people_kind)
a4



a7 =asw_score(dat.pca_last_combat_edata_addmod$ind$coord,batch.id,people_kind)
a7

a8 =asw_score(dat.pca_last_combat_edata_addmod_real$ind$coord,batch.id,people_kind) 
a8


a2
a3 
a4
a7
a8 

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

  fscoreARI <- (2 * (1 - ARI_batch_norm)*(ARI_Cell_norm))/(1 - ARI_batch_norm + ARI_Cell_norm)
  
  
  

  return(fscoreARI)
}


b2 =ari_score(dat.pca_uncorrected$ind$coord,batch.id,people_kind)
b2 

b3 =ari_score(dat.pca_Xcom$ind$coord,batch.id,people_kind) 
b3


b4 =ari_score(dat.pca_MNN_correct_t$ind$coord,batch.id,people_kind)
b4 


b7 =ari_score(dat.pca_last_combat_edata_addmod$ind$coord,batch.id,people_kind) 
b7

b8 =ari_score(dat.pca_last_combat_edata_addmod_real$ind$coord,batch.id,people_kind)
b8


