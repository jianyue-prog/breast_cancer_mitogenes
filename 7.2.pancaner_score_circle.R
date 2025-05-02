
rm(list=ls())
setwd("~//R//my_projects//metabolism_sig//10.泛癌分数圈图")

#setwd('D:\\科研\\肺癌乳酸代谢\\10.泛癌分数圈图')
if(T){
  library(circlize)
  library(ComplexHeatmap)
  library(dplyr)
  library(tibble)
  library(data.table)
  ##load #META
  rawAnno <- read.delim("merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) 
  rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
  samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
  samAnno <- samAnno[which(samAnno$`cancer type` != ""),]
  # 输出样本分组信息
  #write.table(samAnno,"easy_input_sample_annotation.txt",sep = "\t",row.names = F,col.names = T,quote = F)
  # 直接读取样本信息，这里是每个sample对应的肿瘤类型
  samAnno <- read.table("easy_input_sample_annotation.txt", sep = "\t",header = T, check.names = F)
  # 快速读取表达谱数据并做数据预处理
  expr <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
  #expr[1:3,1:3] #20531
  expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
  gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1)
  expr$gene <- gene
  expr <- expr[!duplicated(expr$gene),]#20502
  rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]
  #expr[expr < 0] <- 0 # 对于这份泛癌数据，将略小于0的数值拉到0，否则不能取log（其他途径下载的泛癌数据可能不需要此操作）
  colnames(expr) <- substr(colnames(expr),1,15) # 非TCGA数据可忽略这行
  gc()
  expr<-expr[-1,]
  #dat <- t(as.matrix(log2(expr + 1)))#11069-20501
  dat<-as.data.frame( log2(expr + 1) )
  dat<-as.data.frame(t(dat))#as.data.frame(t(dat))
  #dat<-na.omit(dat)#20501 -->9318 lose gene
  rownames(dat)<-gsub("\\.","-",rownames(dat))
}
#colnames(dat)[grep("TIMM",colnames(dat) )] #无TIMM23  TIMM13替换！ 
#step_rsf[["xvar.names"]][grep("TIMM",step_rsf[["xvar.names"]] )]


#生存数据
if(T){
  cli <- read.csv("~//R//my_projects//metabolism_sig//12.泛癌生存曲线//Survival_SupplementalTable_S1_20171025_xena_sp.csv",header = T,row.names = 1)
  jj <- intersect(rownames(dat),rownames(cli))
  cli <- cli[jj,] #9565
  #合并
  cli<-cli[25:26]
  cli$ID <-rownames(cli)
  
  dat$ID<-rownames(dat)
  dat<-merge(dat,cli,by="ID")
  dat<-dat[!duplicated(dat$ID),] #删除重复
  rownames(dat)<-dat$ID;dat<-dat[,-1]
  library(tidyverse)
  dat<-dat %>% select(c("OS","OS.time") , everything())
  colnames(dat)[1]<-"fustat"
  colnames(dat)[2]<-"futime"
  dat$futime<-dat$futime/365
}
plot(dat$futime)#0-30

#输出文件夹
folder_path <- "./out" # 检查文件夹是否存在
if (!dir.exists(folder_path)) {   # 如果文件夹不存在，则创建它
  dir.create(folder_path)
  print(paste("Folder created at", folder_path))
} else {
  print(paste("Folder already exists at", folder_path))
}

######加载我的得分
# best.iter <- readRDS("best.iter.rds")#500L
# gbm_fit <- readRDS("gbm_fit.rds")#GBM_FIT
# key_gene <- readRDS("key_gene.rds")#29
# #head(expr)[,1:6]
# library(gbm)

# my_model_name <-"StepCox[forward] + RSF"
# my_model <-readRDS("~//R//my_projects//metabolism_sig//5.1-101建模/out/101_results.rds")
# #取出my_model
# my_model <-my_model$ml.res[my_model_name] #自己单独计算？F2查看Mime1::ML.Dev.Prog.Sig
# my_model$`StepCox[forward] + RSF` #比较是否一致？

#重新运行模型
# step_fit <- step(coxph(Surv(futime,fustat)~.,data = trainset ),dirtion="forward") #逐步回归_both
# for (y in 1:length(cohort)) {#id <- paste0("step_",i)
#   cohort[[y]] <- as.data.frame(cohort[[y]])
#   risk_score <-  as.numeric(predict(step_fit, newdata = cohort[[y]][,key_gene]))
#   id=names(predict(step_fit, newdata = cohort[[y]][,key_gene]))
#   risk_score_all <- rbind( risk_score_all,
#                            cbind(algorithm=rep("Step (forward)",length(risk_score)),
#                                  data=rep(y,length(risk_score)),
#                                  risk_score=risk_score,
#                                  id=id,
#                                  futime= cohort[[y]][,"futime"],
#                                  fustat= cohort[[y]][,"fustat"])
#   )
#   
#   
# }
#LOAD cohort
if(T){
  id <- list.files(path = "~//R//my_projects//metabolism_sig//5.确定基因//data",pattern = "csv")
  id
  cohort <- list()
  library("data.table");library(dplyr);library(tibble)
  for (i in 1:length(id)) {
    cohort[[i]] <- fread(paste0("~//R//my_projects//metabolism_sig//5.确定基因//data//",id[i]),#id[i],
                         header = T,data.table = F)
    cohort[[i]]  <-na.omit( column_to_rownames(cohort[[i]] ,var = colnames(cohort[[i]] )[1]))
    #print(max(cohort[[i]][,3:ncol(cohort[[i]])]))
    
    cohort[[i]] <- cohort[[i]] %>%
      rownames_to_column("ID")
    colnames(cohort[[i]] )[2:3] <- c("futime", "fustat")
    #删除生存=0
    cohort[[i]] <- cohort[[i]][cohort[[i]]$futime != 0,]
  }
  names(cohort) <- as.character(lapply(strsplit(id,"_"),function(i){i <- i[1]}))
  cohort <- cohort[7:1]
  names(cohort) <- toupper(names(cohort) )

}
#提取队列+基因表达

#key_gene <-read.csv("~/R/my_projects/metabolism_sig/5.确定基因/out/2024-08-18-确定基因.csv")
#key_gene <-as.character(key_gene$id)
#提取出模型基因表达-重新建模应该输入COX基因
#my_model <-readRDS("~//R//my_projects//metabolism_sig//5.1-101建模/out/101_results_178input.rds")
#key_gene <-my_model$Sig.genes #202 #colnames(my_model[["StepCox[forward] + RSF"]][["xvar"]])
#rm(my_model)
#输入线粒体基因重新建模
genelist <- readxl::read_xls("~/R/my_projects/metabolism_sig/1.scRNA/Human.MitoCarta3.0.xls",sheet = 2)
genelist <-unique(genelist$Symbol)#1136
df<-read.csv("~/R/my_projects/metabolism_sig/2.enrichment/out/差异基因.csv",row.names = 1)
gene <-intersect(genelist,df[df$Group !="N.S.",]$id )#244
gene <-intersect(gene,colnames(cohort$TCGA) ) #与表达矩阵取交集 176

for (y in names(cohort) ) {  #1:length(cohort)
  
  cohort[[y]] <- cohort[[y]][,c("futime","fustat",gene)]
  cohort[[y]] <- cohort[[y]][cohort[[y]]$futime>0,]
} #提出keygene表达矩阵


#run mymodel----
if(T){
  names(cohort)[1]<-"train_set"
  #
  trainset <- cohort$train_set#建模.R
  library(survival)
  library(BiocParallel)
  register(MulticoreParam(workers = 20, progressbar = TRUE))#20 ~1H 30min
  step_fit <- step(coxph(Surv(futime,fustat)~.,data = trainset ),dirtion="forward")#slow
  key_gene1 <-  names(step_fit$coefficients )#96   63
  trainset <- cohort$train_set [,c("futime","fustat",key_gene1 )]
  #提取基因再次建模
  library(obliqueRSF);library(randomForestSRC)
  seed=123
  step_rsf <- rfsrc(Surv(futime,fustat)~.,data = trainset,           #随机森林
                    ntree = 1000,  
                    splitrule = 'logrank',
                    importance = T,
                    proximity = T,
                    forest = T,
                    seed = seed)#
  step_rsf_tc <- tune.rfsrc(Surv(futime, fustat) ~ ., data = trainset, ntreeTry = 1000, stepFactor = 1.5) #23:26
  step_rsf <- rfsrc(Surv(futime,fustat)~.,data = trainset,
                    ntree = which.min(step_rsf$err.rate),nodesize = step_rsf_tc$optimal[1], 
                    mtry=step_rsf_tc$optimal[2],
                    splitrule = 'logrank',
                    importance = T,
                    proximity = T,
                    forest = T,
                    seed = seed)#随机森林
  risk_score_all<-data.frame()
  for (y in 1:length(cohort)) {
    cohort[[y]] <- as.data.frame(cohort[[y]])
    risk_score <-  as.numeric( predict(step_rsf,cohort[[y]][,key_gene1])$predicted)
    
    risk_score_all <- rbind( risk_score_all,
                             cbind(algorithm=rep("RSF",length(risk_score)),
                                   data=rep(y,length(risk_score)),
                                   risk_score =risk_score ,
                                   id=rownames(cohort[[y]]),
                                   futime= cohort[[y]][,"futime"],
                                   fustat= cohort[[y]][,"fustat"])
    )
  }
  #saveRDS(step_fit,"./out/step_fit.rds")#20240914覆盖  63
  #saveRDS(step_rsf,"./out/step_rsf.rds")
  #all_step_both <- mlr(cohort)
  #all_step_both$algorithm <- paste("Step (both) +",all_step_both$algorithm)
  #save.image("./out/Step (both).rdata")
}


#apply to pancancer----
#dat<-as.data.frame(dat)
#best.iter=500
#计算样本风险得分
#colnames(step_rsf$xvar) #96
step_rsf<-readRDS("./out/step_rsf.rds")
library(obliqueRSF);library(randomForestSRC)
#dat<-dat[intersect(colnames(cohort$TCGA),colnames(dat))[colnames(cohort$TCGA)] ] #203 少一个基因
dat <-dat[is.na(dat$futime) !=T,]#11007 #7
dat <-dat[is.na(dat$fustat) !=T,]#10952

# predict.rfsrc(step_rsf,as.data.frame(cohort$TCGA)# dat #cohort$TCGA
#               #dat[colnames(cohort$TCGA)]
# )

#predict.rfsrc(step_rsf,dat[,colnames(cohort$train_set)])

# table(colnames(dat) %in%  colnames(cohort$TCGA) ) 
# all(colnames(dat) %in%  colnames(cohort$TCGA) )
# which(colnames(dat) %in%  colnames(cohort$TCGA))
# intersect(colnames(dat),step_rsf[["xvar.names"]] ) #96-->95 差一个
# 
# which(colnames(dat) %in%  step_rsf[["xvar.names"]] )
# colnames(dat)[which(colnames(dat) %in%  step_rsf[["xvar.names"]] )]
# y=step_rsf[["xvar.names"]]

#dat1<-dat[,colnames(cohort$TCGA)]
# predict.rfsrc(step_rsf,dat[,c("futime","fustat",step_rsf[["xvar.names"]]) ]
# )
# #"TIMM23" %in%  colnames(dat) 
# 
# predict.rfsrc(step_rsf,dat[,
#                            intersect(step_rsf[["xvar.names"]],colnames(dat)) ] 
# )
# predict(step_rsf,#gbm_fit,
#         dat[,c("futime","fustat",step_rsf[["xvar.names"]]) ]
#         #as.data.frame(cohort[[2]])
#         #dat[,intersect(colnames(dat),colnames(step_rsf$xvar) ) ] #[,key_gene] #intersect(colnames(dat),colnames(step_rsf$xvar) )#95 #key_gene
#         #Error in eval(predvars, data, env) : object 'TIMM23' not found
#         #,best.iter
# )
#undefined columns selected
#key_gene1 <-step_rsf[["xvar.names"]]#96

#step_rsf$xvar.names
#x-variables in test data do not match original training data 生存数据-输入格式！
# riskscore <- predict(step_rsf,#gbm_fit,
#                      as.data.frame(cohort[[2]])
# )
# riskscore <-predict.rfsrc(step_rsf,cohort[[2]])
# riskscore<-as.data.frame(riskscore$predicted)
#输入格式！！！
risk_score<-predict(step_rsf,#gbm_fit,
                    dat[,c("futime","fustat",step_rsf[["xvar.names"]]) ]  )
riskscore <- as.numeric(risk_score[["predicted"]] )

es=data.frame(id=rownames(risk_score[["xvar"]]),riskscore=riskscore)
me=merge(samAnno,es,by.x = 2,by.y = 1)
# 假设你的数据框叫做data
average_risk_scores <- me %>%
  group_by(`cancer type`) %>%
  summarise(average_risk_score = mean(riskscore, na.rm = TRUE))
average_risk_scores=column_to_rownames(average_risk_scores,var = "cancer type")
average_risk_scores=na.omit(average_risk_scores)
##进行0-1标准化
#定义矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对score进行矫正
average_risk_scores$average_risk_score=normalize(average_risk_scores$average_risk_score)
## 输入队列信息（队列名称、各队列样本数、需要绘制的得分）
average_risk_scores=average_risk_scores[c("ACC", "PCPG", "PAAD", "THCA", "LGG", "GBM", "UVM", "MESO",
                                          "SKCM", "SARC", "UCS", "UCEC", "OV", "CESC", "BRCA", "TGCT",
                                          "PRAD", "BLCA", "DLBC", "THYM", "KICH", "KIRP", "KIRC",
                                          "ESCA", "STAD", "HNSC", "LIHC", "CHOL", "LUAD", "LUSC", "READ", "COAD"),,drop=F]
CohortInfo <- data.frame(
  "Cohort" = c("ACC", "PCPG", "PAAD", "THCA", "LGG", "GBM", "UVM", "MESO",
               "SKCM", "SARC", "UCS", "UCEC", "OV", "CESC", "BRCA", "TGCT",
               "PRAD", "BLCA", "DLBC", "THYM", "KICH", "KIRP", "KIRC",
               "ESCA", "STAD", "HNSC", "LIHC", "CHOL", "LUAD", "LUSC", "READ", "COAD"),
  "Number" = c(77, 182, 179, 512, 522, 165, 79, 87,
               469, 262, 57, 181, 426, 306, 1098, 137,
               496, 407, 47, 119, 66, 289, 531,
               182, 414, 520, 371, 36, 515, 498, 92, 288),
  "CS.score" = average_risk_scores$average_risk_score # 随机生成的CS得分
)
Class1 <- c("ACC" = "Secretory gland", "PCPG" = "Secretory gland", "PAAD" = "Secretory gland", "THCA" = "Secretory gland",
            "MESO" = "Soft Tissue", "SKCM" = "Soft Tissue", "SARC" = "Soft Tissue", 
            "UCS" = "Female\nreproductive\norgan", "UCEC" = "Female\nreproductive\norgan",
            "OV" = "Female\nreproductive\norgan", "CESC" = "Female reproductive\norgan",
            "BRCA" = "Female\nreproductive organ", "TGCT" = "Male genitourinary system",
            "PRAD" = "Male genitourinary system", "BLCA" = "Male genitourinary system",
            #"LAML" = "Hematopoietic\nand immune", 
            "DLBC" = "Hematopoietic\nand immune", 
            "THYM" = "Hematopoietic\nand immune", "ESCA" = "Digestive\nsystem", 
            "STAD" = "Digestive\nsystem", "HNSC" = "Head and neck", "LIHC" = "Hepatobiliary system")
Class2 <- c("ACC" = "Adrenal\ngland", "PCPG" = "Adrenal\ngland", "THCA" = "Thyroid",
            "LGG" = "Brain", "GBM" = "Brain", "UVM" = "Eye", "SKCM" = "Skin", 
            "UCS" = "Uterus", "UCEC" = "Uterus", "OV" = "Ovary", 
            "BRCA" = "Breast", "PRAD" = "Prostate", "KIRP" = "Kidney", "KIRC" = "Kidney",
            "ESCA" = "Stomach", "STAD" = "Stomach", "HNSC" = "Head\nand neck", "LIHC" = "Liver",
            "LUAD" = "Lung", "LUSC" = "Lung", "READ" = "Colon", "COAD" = "Colon")

## 转换为样本信息（每个样本为一行）
## 注意安排好队列顺序，保证同类队列相邻
## Cohort、Class1和Class2用于区分不同分区
## CohortLabel、Class1Label和Class2Label用于绘制标识文字
SampleInfo <- data.frame("Cohort" = rep(CohortInfo$Cohort, CohortInfo$Number))
rownames(SampleInfo) = paste0("Sample", 1:nrow(SampleInfo))
SampleInfo$Class1 <- plyr::mapvalues(x = SampleInfo$Cohort, from = names(Class1), to = Class1)
SampleInfo$Class1 <- ifelse(test = SampleInfo$Class1%in%Class1, yes = SampleInfo$Class1, no = NA)
SampleInfo$Class2 <- plyr::mapvalues(x = SampleInfo$Cohort, from = names(Class2), to = Class2)
SampleInfo$CohortLabel <- paste0("n=", CohortInfo$Number[match(SampleInfo$Cohort, CohortInfo$Cohort)], "\n", SampleInfo$Cohort)
SampleInfo$Class1Label <- ifelse(test = SampleInfo$Class1%in%Class1, yes = SampleInfo$Class1, no = NA)
SampleInfo$Class2Label <- ifelse(test = SampleInfo$Class2%in%Class2, yes = SampleInfo$Class2, no = NA)
SampleInfo$CS.score <- CohortInfo$CS.score[match(SampleInfo$Cohort, CohortInfo$Cohort)]
openxlsx::write.xlsx(SampleInfo, "easy_input.xlsx")
SampleInfo <- openxlsx::read.xlsx("easy_input.xlsx")
head(SampleInfo)
## 指示各环对应颜色
Cohort.Col <- c("ACC" = "#F5CCDC", "PCPG" = "#FED9EA", "PAAD" = "#F6B5B6", "THCA" = "#EEDEED",
                "LGG" = "#AE8780", "GBM" = "#C59B95", "UVM" = "#7B7C75", "MESO" = "#F7F7BC",
                "SKCM" = "#D0D36F", "SARC" = "#E5E6AB", "UCS" = "#88C4A8", "UCEC" = "#79D9E4",
                "OV" = "#D8F3FA", "CESC" = "#83C3E4", "BRCA" = "#91BBE5", "TGCT" = "#71A5D0",
                "PRAD" = "#76A3CC", "BLCA" = "#A0C2DD", 
                #"LAML" = "#F7908F", 
                "DLBC" = "#ED5C64", 
                "THYM" = "#A55A5E", "KICH" = "#CE8BC0", "KIRP" = "#D2ACD1", "KIRC" = "#D8B4D5",
                "ESCA" = "#DECEE8", "STAD" = "#CAB4D8", "HNSC" = "#FBF7B8", "LIHC" = "#F0A363",
                "CHOL" = "#AA5F39", "LUAD" = "#D7E7E7", "LUSC" = "#CCD4D6", "READ" = "#98D38F",
                "COAD" = "#8BCF9D")

Class1.Col <- c("Secretory gland" = "#D079AB",
                "Soft Tissue" = "#B4B350",
                "Female\nreproductive\norgan" = "#A3C4E3",
                "Male genitourinary system" = "#38679A",
                "Hematopoietic\nand immune" = "#EA7C67",
                "Digestive\nsystem" = "#8F65A9",
                "Hepatobiliary system" = "#EC7F29")

Class2.Col <- c("Adrenal\ngland" = "#F6B7D0", "Thyroid" = "#E6CEE0",
                "Brain" = "#C28C83", "Eye" = "#585D5B",
                "Skin" = "#BCBD30", "Uterus" = "#21C2CE",
                "Ovary" = "#B6EDF9", "Breast" = "#7FA7D4",
                "Prostate" = "#5193BF", "Kidney" = "#EDD5EB",
                "Stomach" = "#C6B0D6", "Head\nand neck" = "#FDEA7E",
                "Liver" = '#FCC88F', "Lung" = "#BBC0C3",
                "Colon" = "#76C48A")

Col = c(Cohort.Col, Class1.Col, Class2.Col)
color_mine=c("#074481","#2D7DB8","#6DAFD2","#B7D8E9","#EAF1F4","#FBE4D6",
             "#F6B394","#DD6D56","#B81C2C","#6D0019")
CS.Col = colorRamp2(breaks = seq(0, 1, length.out = 10), color_mine)#连续变量所使用的颜色
# 用于批量标注颜色和标签的函数
SectorColor <- function(sector, col, track.index, sampleInfo, currentCol, standardCol = "Cohort", ...){
  ## sector：需要绘制的区域（可为多个区域）
  ## col：区域所绘制的颜色
  ## track.index：在第几环（从外向内数）绘制
  ## sampleInfo：样本信息
  ## currentCol：绘制区域（sector）使用的哪一列的分类方法
  ## standardCol：样本最基础的分类列（通常即为队列来源）
  cohorts.to.plot = unique(SampleInfo[[standardCol]][SampleInfo[[currentCol]] %in% sector])
  st.degree = max(unlist(lapply(cohorts.to.plot, get.cell.meta.data, name = "cell.start.degree")))
  ed.degree = min(unlist(lapply(cohorts.to.plot, get.cell.meta.data, name = "cell.end.degree")))
  
  draw.sector(start.degree = st.degree, end.degree = ed.degree,
              rou1 = get.cell.meta.data("cell.top.radius", track.index = track.index),
              rou2 = get.cell.meta.data("cell.bottom.radius", track.index = track.index),
              col = col[sector], ...)
}

SectorLabel <- function(sector, col, track.index, sampleInfo, currentCol, standardCol = "Cohort", niceFacing = T, ...){
  ## sector：需要绘制的区域（可为多个区域）
  ## col：区域所绘制的颜色
  ## track.index：在第几环（从外向内数）绘制
  ## sampleInfo：样本信息
  ## currentCol：绘制区域（sector）使用的哪一列的分类方法
  ## standardCol：样本最基础的分类列（通常即为队列来源）
  ## niceFacing：是否自动调整样本角度
  plot.data = sampleInfo[sampleInfo[[currentCol]] %in% sector, ]
  plot.data = plot.data[round(0.5*nrow(plot.data)), ]
  circos.text(x = plot.data$ID, y = 0.5, labels = sector, sector.index = plot.data$Cohort, track.index = track.index, col = col, niceFacing = niceFacing, ...)
}
SampleInfo$Cohort <- factor(SampleInfo$Cohort, levels = unique(SampleInfo$Cohort)) # 将Cohort转为因子变量，保证顺序不乱
SampleInfo$ID <- 1:nrow(SampleInfo) # 手动编序号，以便后续打标签时确定位置

### 初始化
circos.par("cell.padding" = c(0.02, 0.00, 0.02, 0.00), "start.degree" = 90) # 设置起始角度
circos.initialize(sectors = SampleInfo$Cohort, x = SampleInfo$ID) # 设置弦图分区

### 第一环：连续变量，对每一个样本分别绘制热图，通过col参数调整所绘制的列和需要的配色
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, cell.padding = c(0, 1))
for (i in 1:nrow(SampleInfo)) circos.rect(xleft = SampleInfo$ID[i], xright = SampleInfo$ID[i]+1, ybottom = 0, ytop = 1, 
                                          sector.index = SampleInfo$Cohort[i], track.index = 1, border = NA, col = CS.Col(SampleInfo$CS.score[i]))
## 第二环：分类变量，为最基础的队列信息，所使用的颜色存储在Col中
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2)
for (i in unique(SampleInfo$Cohort)) SectorColor(sector = i, col = Col, track.index = 2, sampleInfo = SampleInfo, currentCol = "Cohort") #绘制队列颜色
for (i in unique(SampleInfo$CohortLabel)) SectorLabel(sector = i, col = "black", track.index = 2, sampleInfo = SampleInfo, currentCol = "CohortLabel") # 绘制队列标签

### 第三环：分类变量，所使用的颜色存储在Col中
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2)
for (i in unique(SampleInfo$Class2)) SectorColor(sector = i, col = Col, track.index = 3, sampleInfo = SampleInfo, currentCol = "Class2") #绘制队列颜色
for (i in unique(SampleInfo$Class2Label)) SectorLabel(sector = i, col = "white", track.index = 3, sampleInfo = SampleInfo, currentCol = "Class2Label") # 绘制队列标签

### 第四环：分类变量，所使用的颜色存储在Col中
#### 此处仿照原图对文字进行了调整
#### niceFacing = F：不自动调整文字角度（自动调整会使得部分文字向环外延申）
#### facing = "clockwise"：使文字垂直于环形
#### pos = 2：文字右对齐
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.lty = "blank")
for (i in unique(SampleInfo$Class1)) SectorColor(sector = i, col = Col, track.index = 4, sampleInfo = SampleInfo, currentCol = "Class1", lty = "blank") #绘制队列颜色
for (i in unique(SampleInfo$Class1Label)) SectorLabel(sector = i, col = "black", track.index = 4, sampleInfo = SampleInfo, currentCol = "Class1Label", niceFacing = F, facing = "clockwise", pos = 2) # 绘制队列标签
circos.clear()

#### 添加连续变量的图例
lgd = Legend(col_fun = CS.Col, title = "Score")
draw(lgd, x = unit(0.9, "npc"), y = unit(0.1, "npc"))



#2.pheatmap----
library(pheatmap)
# save_pheatmap_pdf <- function(x, filename, width, height) {
#   library(grid)
#   x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.02)#修改聚类树线条宽度
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height,family = "serif")
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }#保存热图的函数
#注释
SampleInfo$Class2Label <-gsub("\n"," ",SampleInfo$Class2Label)

my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', 
                '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', 
                '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8',
                '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755',
                '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3',
                '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD',
                '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', 
                '#968175') #颜色设置 
ann_colors<-list(
  Cohort=c("ACC"="#E5D2DD","PCPG"='#53A85F',"PAAD"='#F1BB72',"THCA"='#F3B1A0',
           "LGG"='#D6E7A3', "GBM"="#57C3F3","UVM"='#476D87',"MESO"='#E95C59',
           "SKCM"='#E59CC4', "SARC"='#AB3282', "UCS"='#23452F', "UCEC"='#BD956A',
           "OV"='#8C549C',"CESC"='#585658', "BRCA"='#9FA3A8',"TGCT"='#E4C755',
           "PRAD"='#5F3D69', "BLCA"='#C5DEBA',"DLBC"='#58A4C3', "THYM"="#C1E6F3",
           "KICH"="#F7F398","KIRP"="#AA9A59","KIRC"="#E63863","ESCA"="#712820",
           "STAD"="#6778AE","HNSC"="#91D0BE","LIHC"="#B53E2B","CHOL"="#DCC1DD",
           "LUAD"="#CCE0F5","LUSC"="#CCC9E6","READ"="#625D9E","COAD"="#3A6963"),
  
  Class2Label=c("Adrenal gland"="#A6CEE3","NA"="skyblue","Thyroid"="#B2DF8A","Brain"="#33A02C","Eye"="#4DAF4A","Skin"="#FB9A99",
                "Uterus"="#E31A1C","Ovary"="#FDBF6F","Breast"="#FF7F00","Prostate"="#CAB2D6","Kidney"="#984EA3","Stomach"="#6A3D9A",
                "Head and neck"="#FFFF99","Liver"="#B10928","Lung"="#A65628","Colon"="#F781BF"
                )
) #32
colnames(SampleInfo[,c(1,6)])
unique(SampleInfo[,c(6)]);unique(SampleInfo[,c(1)])
#RColorBrewer::brewer.pal(12,"Paired");RColorBrewer::brewer.pal(9,"Set1")[1:9]
#[1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#4DAF4A" "#FB9A99"http://biotrainee.vip:21043/graphics/plot_zoom_png?width=272&height=577 "#E31A1C" "#FDBF6F" "#FF7F00"
#[9] "#CAB2D6" "#984EA3" "#6A3D9A" "#FFFF99""#FFFF33"  "#B10928" "#A65628""#F781BF" "#999999"
library(RColorBrewer)
pdf("./out/pheatmap_pancancer_score-2-240914.pdf",family = "serif",width = 2.8)#,width = 6.66,15 #20/1.5
p1<-pheatmap(as.matrix(SampleInfo[,7]),cluster_rows = F,cluster_cols = F,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(10),#RdYlBu
         annotation_row = SampleInfo[c(1,6)]#,annotation_colors=ann_colors
        # width = 0.1#,height = 5,
         ) #太长-合并-到平均值?
p1
dev.off()
#样本数？3
#匹配其他TCGA数值数据？https://guolab.wchscu.cn/GSCA/#/download-page
length(unique(SampleInfo$Cohort))#32 配色

#3.ggplot2 mean----
library(ggplot2)
SampleInfo %>% ggplot(aes(x=reorder(Cohort,-CS.score),y=CS.score)) +
  #geom_boxplot()+
  geom_point(aes(size = CS.score,color=CS.score))+
  labs(x="",y="Score",color="Score",size="Score")+theme_bw()+theme(panel.grid.minor  = element_blank(),axis.text.x =element_text(angle = 90) )+
  scale_color_viridis_c()#+facet_grid(~Class1,scales = "free")
ggsave("./out/dot_pancancer_score.pdf",width = 6,height = 3.5,family="serif")
#length(unique(SampleInfo$Class1Label)) 
#注释癌症类型？