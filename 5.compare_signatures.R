
rm(list = ls())
setwd("~//R//my_projects//metabolism_sig//7.sig比较")
if(T){
  library(data.table)
  library(survival)
  library(caret)
  library(glmnet)
  library(survminer)
  library(timeROC)
  library(tibble)
  library(gbm)
}



if(T){
  #加载队列数据
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
    colnames(cohort[[i]] )[2:3] <- c("OS.time", "OS")
    #删除生存=0
    cohort[[i]] <- cohort[[i]][cohort[[i]]$OS.time != 0,]
  }
  names(cohort) <- as.character(lapply(strsplit(id,"_"),function(i){i <- i[1]}))
  cohort <- cohort[7:1]
  names(cohort) <- toupper(names(cohort) )
  
  ###
  sig <- read.csv("sig1.csv",header = T)
  sig <- sig[sig$Gene_id%in%colnames(cohort[[1]]),]
}

###########
sig$Coef <- as.numeric(sig$Coef)
length(unique(sig$PMID))##共101个sig

tab<- data.frame()
for (i in names(table(sig$PMID))) {
  #i=30152019
  dat_sig <- sig[sig$PMID==i,]
  for (y in names(cohort)){
    #y="TCGA"
    dat <- cohort[[y]]
    dat <- dat[,dat_sig$Gene_id,drop=F]
    myFun=function(x){crossprod(as.numeric(x),dat$Coef)}
    riskscore <- apply(dat ,1,myFun)
    tab <- rbind(tab,cbind(PMID=rep(i,length(riskscore)),
                           data=rep(y,length(riskscore)),
                           RS=as.numeric(riskscore ),
                           ID=names(riskscore),
                           futime=cohort[[y]][,"futime"],
                           fustat=cohort[[y]][,"fustat"])) #all score sample
  }
  
}

######我的得分
library(glmnet)
my_model=readRDS("gbm_fit.rds")
best.iter <- readRDS("best.iter.rds")
#summary(my_model)

#替换？加载之前模型
# my_model_name <-"StepCox[forward] + RSF"
# my_model <-readRDS("~//R//my_projects//metabolism_sig//5.1-101建模/out/101_results.rds")
# #取出my_model
# my_model <-my_model$ml.res[my_model_name] #自己单独计算？F2查看Mime1::ML.Dev.Prog.Sig
#放弃此部分代码-使用mime-101结果：~/R/my_projects/metabolism_sig/5.1-101建模/out

#自己单独计算？使用Mine1代码----
if (T) {
  library(Matrix)
  library(survival)
  library(randomForestSRC)
  library(glmnet)
  library(plsRcox)
  library(superpc)
  library(gbm)
  library(CoxBoost)
  library(survivalsvm)
  library(dplyr)
  library(tibble)
  library(BART)
  library(miscTools)
  library(compareC)
  library(ggplot2)
  library(ggsci)
  library(tidyr)
  library(ggbreak)
  library(mixOmics)
  library(data.table)
  Sys.setenv(LANGUAGE = "en")
  options(stringsAsFactors = FALSE)
}
seed=123
set.seed(seed)

candidate_genes <- readxl::read_xls("~/R/my_projects/metabolism_sig/1.scRNA/Human.MitoCarta3.0.xls",sheet = 2)
candidate_genes <-unique(candidate_genes$Symbol)
train_data<-cohort$TCGA
#colnames(train_data)[2:3]<-c("OS.time", "OS")
common_feature = c("ID","OS.time", "OS", candidate_genes)# "OS.time", "OS"

for ( i in names(cohort) ){
  print(i)
  common_feature = intersect(common_feature, colnames(cohort[[i]]))
}

#candidate_genes<-intersect(candidate_genes,colnames(train_data))#769



est_dd <- as.data.frame(train_data)[, common_feature[-1]]#common_feature[-1]

fit <- rfsrc(Surv(OS.time, OS) ~ ., #OS.time, OS
             data = est_dd, #cohort$TCGA,#est_dd, 
             ntree = 1000, nodesize = 5, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed) #slow
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- train_data[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(cohort,#list_train_vali_Data, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                  alpha = 0, type.measure = "class")
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "response", newx = as.matrix(x[, 
                                                                                  -c(1, 2)]), s = fit$lambda.min)))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("RSF + ", "Ridge")
  #result <- rbind(result, cc)
  ml.res = list();riskscore = list()
  ml.res[[paste0("RSF + ", "Ridge")]] = fit
  
  returnIDtoRS = function(rs.table.list, rawtableID) {
    for (i in names(rs.table.list)) {
      rs.table.list[[i]]$ID = rawtableID[[i]]$ID
      rs.table.list[[i]] = rs.table.list[[i]] %>% dplyr::select("ID", 
                                                                everything())
    }
    return(rs.table.list)
  }
  
  rs = returnIDtoRS(rs.table.list = rs, rawtableID = cohort)#list_train_vali_Data
  riskscore[[paste0("RSF + ", "Ridge")]] = rs
}
 # else {
 #  warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
 # } 




# #加载自己的key_gene!
# #key_gene <- readRDS("key_gene.rds") #
# key_gene <-read.csv("~/R/my_projects/metabolism_sig/5.确定基因/out/2024-08-18-确定基因.csv")
# key_gene <-as.character(key_gene$id) #421
# best.iter <- readRDS("best.iter.rds")
# for (y in names(cohort) ) {  #1:length(cohort)
#   cohort[[y]] <- cohort[[y]][,c("futime","fustat",key_gene)]
#   cohort[[y]] <- cohort[[y]][cohort[[y]]$futime>0,]
# } #提出keygene表达矩阵-基因数目减少，其他模型无法应用


myscore <- data.frame()
for (i in names(cohort) ) {
  #i="TCGA"
  dat <- cohort[[i]]
  riskscore <- as.numeric( predict(fit,#my_model,
                                   dat,#  [,key_gene]
                                   best.iter
                                   ))      #SBS
  myscore<- rbind(myscore,cbind(PMID=rep("My_model",length(riskscore)),   
                                ####换成自己sig的名字
                                data=rep(i,length(riskscore)),
                                RS=as.numeric(riskscore ),
                                ID=rownames(dat),
                                futime=dat[,"futime"],
                                fustat=dat[,"fustat"]))
}

tab <- rbind(myscore,tab)
########################################
if(T){
  library(survival)
  library(paletteer)
  library(ggplot2)
  library(xgboost)
  library(ggsci)
  col <-RColorBrewer::brewer.pal(9,"Paired") #c( '#D6E7A3', '#57C3F3', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  
}
# col <- c('#313c63','#ebc03e','#377b4c',
#          '#7bc7cd','#5d84a4',"#EC7232",
#          pal_lancet(palette = c("lanonc"))(9),
#          pal_npg("nrc")(9),
#          pal_nejm("default", alpha = 0.6)(8),
#          pal_d3("category20")(20))

pic <- list()
for (i in names(table(tab$data))) {
  #i="GSE13213"
  dat <- tab[tab$data==i,]
  colo <- which(names(table(tab$data))==i)
  colo <- col[colo]
  dat$RS <- as.numeric(dat$RS)
  rt_tcga <- data.frame()
  for (y in names(table(dat$PMID))) {
    #y="30152019"
    dat_rs <- dat[dat$PMID==y,]
    dat_rs$futime <- as.numeric(dat_rs$futime)
    #dat_rs[dat_rs$futime==0,]$futime=0.001
    dat_rs$fustat <- as.numeric(dat_rs$fustat)
    mylog=summary(coxph(Surv(futime, fustat)~dat_rs$RS, data=dat_rs))
    
    c_index <-  mylog$concordance[1]
    se <- mylog$concordance[2]
    rt_tcga <- rbind(rt_tcga,cbind(
      pmc_id=paste0("PMID: ",y),
      c_index=c_index,
      c_index_high=c_index+1.96*se,
      c_index_low =c_index-1.96*se ))
  }

  rt_tcga$pmc_id <-gsub('PMID: Prol',"Prol",rt_tcga$pmc_id)
  rt_tcga$c_index <- as.numeric(rt_tcga$c_index)
  rt_tcga$c_index_high <- as.numeric(rt_tcga$c_index_high )
  rt_tcga$c_index_low <- as.numeric(rt_tcga$c_index_low )
  rt_tcga <- rt_tcga[order(rt_tcga$c_index,decreasing = F),]
  rt_tcga$pmc_id <- factor(rt_tcga$pmc_id,levels = rt_tcga$pmc_id)
  which(rt_tcga$pmc_id=="Prol")
  sig_camp <-  ggplot(rt_tcga ) +
    geom_segment(aes(x=pmc_id, xend=pmc_id, y=c_index_high, yend=c_index_low), color="grey",size=1)+
    geom_point( aes(x=pmc_id, y=c_index), color=colo, size=3 )+
    geom_hline(yintercept=0.6, linetype = "dashed" ,color="grey")+
    #theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #axis.text.y= element_text(face = "bold.italic",colour = "#441718",size = 10,hjust = 0.5),
          panel.grid = element_blank(),
          panel.border = element_rect(fill=NA,
                                      color="black", 
                                      size=1.5, 
                                      linetype="solid"),
          panel.background = element_rect(fill = "white"))+
    labs(title=i)+ylab("C-Index")+xlab("")+ #xlab="C-Index",ylab = "",
    coord_flip() 
pic[[i]] <-  sig_camp
  
}
library(patchwork)
pic$TCGA+pic$GSE42127+pic$GSE31210+pic$GSE30219+pic$GSE29016+pic$GSE26939+pic$GSE13213+plot_layout(nrow = 1)
ggsave("campare-1.pdf",width = 210*2,height = 180*2,units = "mm",family="serif")
