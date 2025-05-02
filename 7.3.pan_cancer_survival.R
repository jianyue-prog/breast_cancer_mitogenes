
if(T){
  rm(list=ls())
  setwd("~/R/my_projects/metabolism_sig/12.泛癌生存曲线")
  #setwd("D:\\科研\\肺癌乳酸代谢\\12.泛癌生存曲线")
  #输出文件夹
  folder_path <- "./out" # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {   # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
}

#load exp----
if(T){
  library(data.table)
  ann <- fread("probeMap_gencode.v23.annotation.gene.probemap",header = T,data.table = F)
  rt <- fread("tcga_RSEM_gene_tpm",header = T,data.table = F)
  ann <- ann[,c(1,2)]
  colnames(rt)[1] <- "id"
  rt <- dplyr::inner_join(ann,rt,by="id")
  rt <- rt[,-1]
  rt <- limma::avereps(rt[,-1],rt$gene)
  rt <-as.data.frame(rt)
}

#load model----
# best.iter <- readRDS("best.iter.rds")
# gbm_fit <- readRDS("gbm_fit.rds")
# key_gene <- readRDS("key_gene.rds")
if(T){
  #加载模型
  step_rsf<-readRDS("~//R//my_projects//metabolism_sig//10.泛癌分数圈图//out/step_rsf.rds")
  library(obliqueRSF);library(randomForestSRC)
  library(survival)
  key_gene <-step_rsf[["xvar.names"]]
}

# risk_score<-predict(step_rsf,#gbm_fit,
#                     dat[,c("futime","fustat",step_rsf[["xvar.names"]]) ]  )
# riskscore <- as.numeric(risk_score[["predicted"]] ) #10125
# 
# #riskscore <- as.numeric( predict(gbm_fit,as.data.frame(scale(dat[,key_gene])),best.iter))
# es=data.frame(id=rownames(risk_score[["xvar"]]), #rownames(dat)
#               riskscore=riskscore)



  rt <- rt[key_gene,]
  rt <- t(rt)
  cli <- read.csv("Survival_SupplementalTable_S1_20171025_xena_sp.csv",header = T,row.names = 1)
  jj <- intersect(rownames(rt),rownames(cli))
  cli <- cli[jj,]
  rt <- rt[jj,]
  rt <- (2^rt)-0.001
  rt <- log2(rt+1)
  
  cli <- cli[,c(2,25,26)]
  cli <- cli[jj,]
  rt <- rt[jj,]
  rt <- cbind(cli,rt)
  colnames(rt)[1:3] <- c("type","fustat","futime")
  #
  library(survival);library(survminer);library(randomForestSRC);library(gbm)
  rt$futime <- rt$futime/365
  #
  
  type<- c()
  for (i in names(table(rt$type))) {dat <- rt[rt$type==i,]
  dat[,4:ncol(dat)] <- scale(dat[,4:ncol(dat)] )
  dat <- na.omit(dat)
  #dat$riskscore <- as.numeric( predict(gbm_fit,as.data.frame(dat[,key_gene]),best.iter))
  risk_score<-predict(step_rsf,#gbm_fit,
                      dat[,c("futime","fustat",step_rsf[["xvar.names"]]) ]  )
  dat$riskscore <- as.numeric(risk_score[["predicted"]] )
  
  cut_value <- surv_cutpoint(dat, #数据集
                             time = "futime", #生存状态
                             event = "fustat", #生存时间
                             variables = c("riskscore") )$cutpoint[,1] 
  dat$Risk <- ifelse(dat$riskscore>cut_value,"High","Low")
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data =dat)
  pval <- survdiff(Surv(futime, fustat) ~ Risk, rho = 0, data =dat)
  pval <- pval$pvalue
  hr <- summary(coxph(Surv(futime, fustat) ~ Risk, data =dat))
  hr <-hr$conf.int[,"exp(coef)"]
  if (pval<0.05&hr<1) {pdf(paste0("./out/",i,"_km.pdf"),3,3)
    print(ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                                            conf.int = F,### 是否添加置信区间
                                            legend.title = "Risk", # 设置图例标题
                                            legend.labs = c("High", "Low"), # 指定图例分组标签
                                            risk.table = F, # 是否添加风险表
                                            risk.table.col = "strata",
                                            title = i, # 添加标题

                                            font.main = c(16), # 设置标题字体大小、格式和颜色, "bold.italic", "black"
                                            font.x = c(14), # 设置x轴字体大小、格式和颜色
                                            font.y = c(14), # 设置y轴字体大小、格式和颜色
                                            font.tickslab = c(12),# 设置坐标轴刻度字体大小、格式和颜色
                                            #linetype = "solid",#改变线条的渐变
                                            legend=c(0.8,0.7),
                                            ###linetype = "strata",
                                            #surv.median.line = "hv", # 是否添加中位生存线
                                            risk.table.y.text.col = F,risk.table.y.text = FALSE,
                                            ggtheme = theme_bw()+#theme(legend.text = element_text(colour = c("red", "blue")))
                                            theme(legend.position = "right", panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                                            +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                                            +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                                            +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                                            +theme(legend.title = element_text(family = "Times",colour = "black",size = 12))
                                            +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                                            palette = c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                                            #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                                            xlab = "Years")##随访的时间时天，就写Days，月份就写Months

          
    )
    
    type<- c(type,i)
    dev.off()
    
    
    
  }
  
  
  }
saveRDS(type,"具有预后价值的肿瘤ID.rds")
#   for (i in names(table(rt$type)) ) {
#     dat <- rt[rt$type==i,]
#   dat[,4:ncol(dat)] <- scale(dat[,4:ncol(dat)] )
#   #dat <- na.omit(dat)
# 
# #
# #dat$riskscore <- as.numeric( predict(gbm_fit,as.data.frame(dat[,key_gene]),best.iter))
# risk_score<-predict(step_rsf,#gbm_fit,
#                     dat[,c("futime","fustat",step_rsf[["xvar.names"]]) ]  )
# dat$riskscore <- as.numeric(risk_score[["predicted"]] )
# 
# cut_value <- surv_cutpoint(dat, #数据集
#                            time = "futime", #生存状态
#                            event = "fustat", #生存时间
#                            variables = c("riskscore") )$cutpoint[,1] 
# dat$Risk <- ifelse(dat$riskscore>cut_value,"High","Low")
# 
# fit <- survfit(Surv(futime, fustat) ~ Risk, data =dat)
# pval <- survdiff(Surv(futime, fustat) ~ Risk, rho = 0, data =dat)
# pval <- pval$pvalue
# hr <- summary(coxph(Surv(futime, fustat) ~ Risk, data =dat))
# hr <-hr$conf.int[,"exp(coef)"]
# if (pval<0.05&hr<1) {pdf(paste0("./out/",i,"_km.pdf"),3,3) #paste0(i,"_km.pdf")
#   ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
#                   conf.int = F,### 是否添加置信区间
#                   legend.title = "Risk", # 设置图例标题
#                   legend.labs = c("High", "Low"), # 指定图例分组标签
#                   risk.table = F, # 是否添加风险表
#                   risk.table.col = "strata", 
#                   title = i, # 添加标题
#                   
#                   font.main = c(16), # 设置标题字体大小、格式和颜色, "bold.italic", "black"
#                   font.x = c(14), # 设置x轴字体大小、格式和颜色
#                   font.y = c(14), # 设置y轴字体大小、格式和颜色
#                   font.tickslab = c(12),# 设置坐标轴刻度字体大小、格式和颜色
#                   #linetype = "solid",#改变线条的渐变
#                   legend=c(0.8,0.8),
#                   ###linetype = "strata",
#                   #surv.median.line = "hv", # 是否添加中位生存线
#                   risk.table.y.text.col = F,risk.table.y.text = FALSE,
#                   ggtheme = theme_bw()+#theme(legend.text = element_text(colour = c("red", "blue")))
#                   theme(legend.position = "right", panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#                   +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
#                   +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
#                   +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
#                   +theme(legend.title = element_text(family = "Times",colour = "black",size = 12))
#                   +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
#                   palette = c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
#                   #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
#                   xlab = "Years")##随访的时间时天，就写Days，月份就写Months
#   
#  
#   # print(ggsurvplot(fit, 
#   #                  #surv.median.line = "hv", # 添加中位数生存时间线
#   #                  
#   #                  # Change legends: title & labels
#   #                  legend.title = "Risk", # 设置图例标题
#   #                  legend.labs = c("High", "Low"), # 指定图例分组标签
#   #                  
#   #                  
#   #                  # Add p-value and tervals
#   #                  pval = TRUE, # 设置添加P值
#   #                  pval.method = TRUE, #设置添加P值计算方法
#   #                  conf.int = F, # 设置添加置信区间
#   #                  size=1.8,
#   #                  # Add risk table
#   #                  #risk.table = TRUE, # 设置添加风险因子表
#   #                  #tables.height = 0.2, # 设置风险表的高度
#   #                  #tables.theme = theme_dark(), # 设置风险表的主题
#   #                  
#   #                  # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
#   #                  # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
#   #                  palette = c("#C31820","#106C61"), # 设置颜色画板
#   #                  ggtheme = theme(panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
#   #                                  panel.background = element_rect(fill = "white"),
#   #                                  panel.grid = element_blank(),
#   #                                  legend.key = element_rect(fill = NA),
#   #                                  legend.background = element_rect(fill = NA),
#   #                                  legend.text = element_text(face ="bold.italic" ),
#   #                                  legend.title = element_text(face ="bold.italic" )),
#   #                  title = i, # 添加标题
#   #                  font.main = c(16, "bold.italic", "black"), # 设置标题字体大小、格式和颜色
#   #                  font.x = c(14, "bold.italic", "black"), # 设置x轴字体大小、格式和颜色
#   #                  font.y = c(14, "bold.italic", "black"), # 设置y轴字体大小、格式和颜色
#   #                  font.tickslab = c(12, "plain", "black"),# 设置坐标轴刻度字体大小、格式和颜色
#   #                  #linetype = "solid",#改变线条的渐变
#   #                  legend=c(0.8,0.8))+ylab("Overall survival")
#   #       
#   # )
#   
#   type<- c(type,i)
#   dev.off()
# 
#   
#   
# }
#   }

saveRDS(type,"./out/具有预后价值的肿瘤ID.rds")

