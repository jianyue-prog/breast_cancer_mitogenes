#见5.1 和6.0

#ROC+KM
#手动显示101最佳模型的ROC和KM-riskscore-median
#保存关键基因表达+预后+风险分组！

rm(list = ls())
setwd("~/R/my_projects/metabolism_sig/6.生存曲线与ROC")
#加载之前的模型结果
#setwd("~//R//my_projects//metabolism_sig//5.1-101建模");getwd()
#saveRDS(res,"./out/101_results.rds")
my_model <-"StepCox[forward] + RSF"
res <-readRDS("~//R//my_projects//metabolism_sig//5.1-101建模/out/101_results_178input.rds")
#取出my_model
#res1 <-res$ml.res[my_model] #自己单独计算？F2查看Mime1::ML.Dev.Prog.Sig
riskscore<-res[["riskscore"]][[my_model]] # 提出riskscore
saveRDS(riskscore,"各个队列评分.rds")
#write.csv(riskscore,"tcga_risk.csv")
#riskscore<- readRDS("各个队列评分.rds")
#加载队列数据
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

#准备数据
if(T){
  #key_gene <-read.csv("~/R/my_projects/metabolism_sig/5.确定基因/out/2024-08-18-确定基因.csv")
  #key_gene <-as.character(key_gene$id)
  key_gene <-res$Sig.genes
  for (y in names(cohort) ) {  #1:length(cohort)
    cohort[[y]] <- cohort[[y]][,c("futime","fustat",key_gene)]
    cohort[[y]] <- cohort[[y]][cohort[[y]]$futime>0,]
  } #提出keygene表达矩阵
}



if(T){
  library(survival)
  library(survminer)
  library(timeROC)
  library(ggsci)#配色
  colors <-c( '#D6E7A3', '#57C3F3', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  
}
#' colors=c(#'#313c63','#b42e20','#ebc03e','#377b4c',
#'          #'#7bc7cd','#5d84a4','#4B4B5D',"#EC7232",
#'          pal_lancet(palette = c("lanonc"))(9),
#'          pal_npg("nrc")(9),
#'          pal_nejm("default", alpha = 0.6)(8),
#'          pal_d3("category20")(20),
#'          pal_igv("default")(51),
#'          pal_uchicago("default", alpha = 0.6)(9),
#'          pal_uchicago("dark", alpha = 0.6)(9))

##

for (i in names(cohort)) {
  dat <- riskscore[[i]]#cohort[[i]]
  dat$riskscore <- dat$RS#risk_score
  colnames(dat)[c(2,3)]<-c("futime","fustat")
# dat$Risk <- ifelse(dat$riskscore>surv_cutpoint(dat, #数据集
#                                                time = "futime", #生存状态
#                                                event = "fustat", #生存时间
#                                               variables = c("riskscore") )$cutpoint[,1] ,"High","Low")

dat$Risk <- ifelse(dat$riskscore>median(dat$riskscore) ,"High","Low")


fit <- survfit(Surv(futime, fustat) ~ Risk, data =dat)
# print(ggsurvplot(fit, 
#                  #surv.median.line = "hv", # 添加中位数生存时间线
#                  # Change legends: title & labels
#                  legend.title = "Risk", # 设置图例标题
#                  legend.labs = c("High", "Low"), # 指定图例分组标签
#                  # Add p-value and tervals
#                  pval = TRUE, # 设置添加P值
#                  pval.method = TRUE, #设置添加P值计算方法
#                  conf.int = F, # 设置添加置信区间
#                  size=1.8,
#                  # Add risk table
#                  #risk.table = TRUE, # 设置添加风险因子表
#                  #tables.height = 0.2, # 设置风险表的高度
#                  #tables.theme = theme_dark(), # 设置风险表的主题
#                  # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
#                  # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
#                  palette = c("#C31820","#106C61"), # 设置颜色画板
#                  ggtheme = theme(panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
#                                  panel.background = element_rect(fill = "white"),
#                                  panel.grid = element_blank(),
#                                  legend.key = element_rect(fill = NA),
#                                  legend.background = element_rect(fill = NA),
#                                  legend.text = element_text(face ="bold.italic" ),
#                                  legend.title = element_text(face ="bold.italic" )),
#                  title = i, # 添加标题
#                  font.main = c(16, "bold.italic", "black"), # 设置标题字体大小、格式和颜色
#                  font.x = c(14, "bold.italic", "black"), # 设置x轴字体大小、格式和颜色
#                  font.y = c(14, "bold.italic", "black"), # 设置y轴字体大小、格式和颜色
#                  font.tickslab = c(12, "plain", "black"),# 设置坐标轴刻度字体大小、格式和颜色
#                  #linetype = "solid",#改变线条的渐变
#                  legend=c(0.8,0.8))+ylab("Overall survival")
#       
# )
print(
  
  ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                  conf.int = F,### 是否添加置信区间
                  legend.title = "Risk", # 设置图例标题
                  legend.labs = c("High", "Low"), # 指定图例分组标签
                  risk.table = F, # 是否添加风险表
                  risk.table.col = "strata", 
                  ###linetype = "strata",
                  #surv.median.line = "hv", # 是否添加中位生存线
                  risk.table.y.text.col = F,risk.table.y.text = FALSE,
                  ggtheme = theme_bw()+theme(legend.text = element_text(colour = c("red", "blue")))
                  +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                  +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                  +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                  +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                  +theme(legend.title = element_text(family = "Times",colour = "black",size = 12))
                  +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                  palette = c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                  #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                  xlab = "Years")##随访的时间时天，就写Days，月份就写Months
  
)

ggsave(filename = paste0(i,"_生存曲线.pdf"), height = 3.5, width = 3.5)
mk=c(1,3,5)
ROC_rt=timeROC(T=dat$futime,delta=dat$fustat,
               marker=dat$riskscore,cause=1,
               weighting='aalen',
               times=mk,ROC=TRUE)

p.dat=data.frame()
for(x in which(ROC_rt$times>0)){
  lbs=paste0('AUC at ',mk[x],' years: ',round(ROC_rt$AUC[x],2))
  p.dat=rbind.data.frame(p.dat,data.frame(V1=ROC_rt$FP[,x],V2=ROC_rt$TP[,x],Type=lbs))
}

ggplot(p.dat, aes(x=V1,y=V2, fill=Type))+
  stat_smooth(aes(colour=Type),se = FALSE, size = 1)+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  ylab('True positive fraction')+
  theme(axis.text=element_text(face = "bold.italic",colour = "#441718",size = 9), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        # axis.text.x=element_blank(),
        axis.title=element_text(face = "bold.italic",colour = "#441718",size = 12), #设置y轴标题的字体属性
        #plot.title=element_blank(),
        legend.text=element_text(face = "bold.italic",colour = "#441718",size = 12),  #设置图例的子标题的字体属性
        legend.title=element_text(face = "bold.italic",colour = "#441718",size = 12),#设置图例的总标题的字体属性
        legend.justification=c(1,1),
        legend.position=c(1,0.5),
        legend.background = element_rect(fill = NA, colour = NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#E6F0F7"),
        plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"),
        plot.title=element_blank())+xlab('False positive fraction')

ggplot(p.dat, aes(x=V1,y=V2, color=Type))+
  #stat_smooth(aes(colour=Type),se = FALSE, size = 1)+
  geom_line(aes(x=V1,y=V2), size = 0.7)+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  #geom_abline(slope = 1, intercept = 0, color = "grey", size = 0.5, linetype = 2)+
  theme_bw() +
  labs(x = "1-Specificity (FPR)", y = "Sensitivity (TPR)",size=16,family="serif")+
  theme(
    legend.position=c(0.68,0.20),legend.text = element_text(size = 11),
    axis.text = element_text(size = 15, color = "black",family="serif"),#坐标轴字体
    axis.title.x = element_text( size = 20, color = "black",family="serif"),
    axis.title.y = element_text( size = 20, color = "black",family="serif"),panel.grid = element_blank()#删除网格
  )+
  scale_x_continuous(expand = c(0,0.05)) +scale_y_continuous(expand = c(0,0.01))


#丑修改
ggsave(filename =paste0(i,"_ROC.pdf"),width=3.7,height=3.5,family="serif")

}

