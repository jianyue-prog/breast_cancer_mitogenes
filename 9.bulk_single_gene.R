
rm(list = ls())
#setwd("~//R//my_projects/metabolism_sig//8.PCA和risk热图//risk热图")
setwd("~/R/my_projects/metabolism_sig/19.bulk_single_gene")

#设置输出文件夹
folder_path <- "./out"
# 检查文件夹是否存在
if (!dir.exists(folder_path)) {
  # 如果文件夹不存在，则创建它
  dir.create(folder_path)
  print(paste("Folder created at", folder_path))
} else {
  print(paste("Folder already exists at", folder_path))
}


#0.load data----
#risk分组+表达+临床
#library(limma)
#riskFile="tcga_risk.csv"      #风险+表达文件 热图展示建模关键基因


if(T){
  my_model <-"StepCox[forward] + RSF"
  res <-readRDS("~//R//my_projects//metabolism_sig//5.1-101建模/out/101_results_178input.rds")
  #risk<-res[["riskscore"]][[my_model]]$TCGA
  #基因重要性
  impor<-as.data.frame(res[["ml.res"]][[my_model]]$importance)  #[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
  colnames(impor)<-"importance"
  plot(density(impor$importance))
  #impor[impor$importance >mean(impor$importance),] #~100
  library(tidyverse)
  impor<-arrange(impor,desc(importance) )#impor[order(impor$importance),] #降序
  #length(impor[impor$importance>0,]) #45
  
  #A.加载风险得分
  riskscore<-res[["riskscore"]][[my_model]]$TCGA
  colnames(riskscore)[4]<-"riskscore"
  riskscore$risk <- ifelse(riskscore$riskscore>median(riskscore$riskscore) ,"High","Low")#分组
  #B.合并关键基因表达
  risk<-cbind(riskscore,
              res[["ml.res"]][[my_model]]$xvar#503 samples -202
  ) 
  rownames(risk)<-risk$ID;risk<-risk[,-1]
  
  
  
  cliFile="~//R//my_projects/metabolism_sig//8.PCA和risk热图//risk热图//clinical.txt"        #临床数据文件
  #读取临床数据文件
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1) #读取临床数据
  # "fustat" "age"    "gender" "stage"  "T"      "M"      "N" 
  
  #读取风险-表达文件
  # # risk=read.csv(riskFile, row.names=1) #得分和riskscore分组
  # # risk$risk=factor(risk$risk, levels=c("high", "low"))#low
  # risk<-read.csv("~//R//my_projects//metabolism_sig//8.PCA和risk热图//pca//meta-cohort.csv",row.names = 1)
  # risk <- risk[risk$cohort=="TCGA",] #训练集
  # #risk <- risk[!duplicated(risk$ID),]
  # # risk <- risk %>%
  # #   rownames_to_column("ID")
  # rownames(risk)<-risk$ID;risk <-risk[,-1]
  
  #合并数据
  samSample=intersect(row.names(risk), row.names(cli))#316  合并临床
  cli=cli[samSample,,drop=F]
  risk=risk[samSample,,drop=F]
  data=cbind( cli,risk)
  data=data[order(data$riskscore),,drop=F]      #根据病人的风险得分对样品排序
  exp=data[,(12:(ncol(data) ))]     #         #提取模型基因的表达量
  #显示前top个基因
  #exp=exp[,colnames(exp) %in%  rownames(impor)[1:50] ] #
  
  # match(colnames(exp),rownames(impor) ) #407 202
  # length(intersect( colnames(exp),rownames(impor) ) )#143
  
  Type=data[,1:11]   #(ncol(risk):ncol(data))        #提取临床信息，作为热图注释文件
  Type=Type[-c(8,10)]
  Type <- Type[c(9,1:8)]#Type[c(ncol(Type),(ncol(Type)-1):1)]
  
  #对临床性状进行循环，得到显著性标记
  sigVec=c("risk")
  for(clinical in colnames(Type[,1:ncol(Type)])){
    #clinical="fustat"
    data=Type[c("risk", clinical)]
    colnames(data)=c("risk", "clinical")
    data=data[(data[,"clinical"]!="unknow"),]
    tableStat=table(data)
    stat=chisq.test(tableStat)
    pvalue=stat$p.value
    Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    sigVec=c(sigVec, paste0(clinical, Sig))
  }
  
 # table(Type$risk)#160 156
}


#1.0 批量单基因KM.ROC----
colnames(risk)
dat<-risk
colnames(dat)[c(1,2)]<-c("futime","fustat")

if(T){
  library(survival)
  library(survminer)
  library(timeROC)
  library(ggsci)#配色
  colors <-c( '#D6E7A3', '#57C3F3', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  
}

for (i in colnames(dat)[5:68] ) {
  dat<-risk
  colnames(dat)[c(1,2)]<-c("futime","fustat")
  #dat <- risk[[i]]#cohort[[i]]
  dat$riskscore <- dat[,i]#risk_score --> gene exp
  #colnames(dat)[c(2,3)]<-c("futime","fustat")
  # dat$Risk <- ifelse(dat$riskscore>surv_cutpoint(dat, #数据集
  #                                                 time = "futime", #生存状态
  #                                                 event = "fustat", #生存时间
  #                                                variables = c("riskscore") )$cutpoint[,1] ,"High","Low")

  dat$Risk <- ifelse(dat$riskscore>median(dat$riskscore) ,"High","Low")

  fit <- survfit(Surv(futime, fustat) ~ Risk, data = dat) #dat

  print(
    
    ggsurvplot_list(fit,data = dat,
                    pval = T,pval.method = TRUE,##是否添加P值
                    conf.int = F,### 是否添加置信区间
                    legend.title = i,#"Risk", # 设置图例标题
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
                    palette = c("pink", "skyblue"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                    #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                    xlab = "Years")##随访的时间时天，就写Days，月份就写Months
    
  )
  
  #ggsave(filename = paste0("./out/",i,"_km.pdf"), height = 3.5, width = 3.5,family="serif")
  ggsave(filename = paste0("./out/",i,"_km_median.pdf"), height = 3.5, width = 3.5,family="serif")
  
}

#ROC
for (i in colnames(dat)[5:68] ) {
  mk=c(1,3,5)
  ROC_rt=timeROC(T=dat$futime,delta=dat$fustat,
                 marker=dat[,i],cause=1,
                 weighting='aalen',
                 times=mk,ROC=TRUE)
  
  p.dat=data.frame()
  for(x in which(ROC_rt$times>0)){
    lbs=paste0('AUC at ',mk[x],' years: ',round(ROC_rt$AUC[x],2))
    p.dat=rbind.data.frame(p.dat,data.frame(V1=ROC_rt$FP[,x],V2=ROC_rt$TP[,x],Type=lbs))
  }
  
  
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
  
  
  #ggsave(filename =paste0("./out/",i,"_ROC.pdf"),width=3.7,height=3.5,family="serif")
  
}

#2.0 stage single gene ----
Type$id <-rownames(Type)
risk$id <- rownames(risk)
df <- merge(Type,risk[-c(2,4)],by="id")
# for (j in colnames(df)[13:ncol(df)] ){
#   print(j)
#   gene<- j
#   ggplot(df,aes(x=N,y=gene ))+geom_boxplot()+geom_point() #`T`
# }  
#转为长数据
#df1 <-reshape2::melt(df, id.vars = "group", measure.vars = c("var1", "var2", "var3"))
##2.1 RISKSCORE N T M----
ggplot(df[df$N != "N3",],aes(x=N,y=riskscore,color=N))+
  geom_point(alpha=0.5,size=2,
             position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
               size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
               outlier.stroke = 0.5)+
  theme_bw(base_size = 16)+theme(legend.position = "none")+
  ggpubr::stat_compare_means()+labs(y="MitoScore")+
  scale_color_manual(values =  c("#EDA065","#66CCFF","#7EC7A7"))
#ggsave(paste0("./out/boxplot_Ttest_",i,".pdf"),width = 2,height = 3,family="serif")
#ggsave("./out/boxplot_riskscore_N-1.pdf",width = 3,height = 3,family="serif")#4 3

ggplot(df,aes(x=df$T,y=riskscore,color=df$T))+
  geom_point(alpha=0.5,size=2,
             position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
               size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
               outlier.stroke = 0.5)+
  theme_bw(base_size = 16)+theme(legend.position = "none")+
  ggpubr::stat_compare_means()+labs(y="MitoScore",x="T")+
  scale_color_manual(values =  colors[4:16] )
#ggsave("./out/boxplot_riskscore_T-1.pdf",width = 3.3,height = 3,family="serif")#4 3

ggplot(df,aes(x=stage,y=riskscore,color=stage))+
  geom_point(alpha=0.5,size=2,
             position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
               size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
               outlier.stroke = 0.5)+
  theme_bw(base_size = 16)+theme(legend.position = "none")+
  ggpubr::stat_compare_means()+labs(y="MitoScore",x="Pathologic Stage")+
  scale_color_manual(values =  colors[4:16] )
ggsave("./out/boxplot_riskscore_stage.pdf",width = 4,height = 3,family="serif")#4 3


for (k in c("age","gender","M") ){
  print(k)
  df1 <-df
  df1$type <-df1[,k]
  ggplot(df1,aes(x=type,y=riskscore,color=type ))+
    geom_point(alpha=0.5,size=2,
               position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
    geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
                 size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    theme_bw(base_size = 16)+theme(legend.position = "none")+
    ggpubr::stat_compare_means()+labs(y="MitoScore",x=k)+
    scale_color_manual(values =  colors[4:16] )
  ggsave(paste0("./out/boxplot_riskscore_",k,".pdf"),width = 3,height = 3,family="serif")
  
}

#2.2 single gene stage----
#colnames(df)[13:ncol(df)] #64

for (gene in  colnames(df)[13:ncol(df)] ){
  print(gene)
  for (i in colnames(df)[4:9] ){
    print(i)
    df1 <-df[df$N !="N3",]
    df1$type <-df1[,i]
    df1$gene <-df1[,gene]
    ggplot(df1,aes(x=type,y=gene,color=type ))+
      geom_point(alpha=0.5,size=2,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.8))+
      geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
                   size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      theme_bw(base_size = 16)+theme(legend.position = "none")+
      ggpubr::stat_compare_means()+labs(y=gene,x="")+ #,x=i
      scale_color_manual(values =  colors[4:16] )
    ggsave(paste0("./out/boxplot_",i,"_",gene,".pdf"),width = 3.5,height = 3,family="serif")
    
  }
}
#gene<-colnames(df)[13:ncol(df)][1]

#3.0 tumor vs. normal exp genes----
expr <- read.table("~/R/my_projects/metabolism_sig/2.enrichment/luad_mRNA.txt",header=T, sep="\t", check.names=F, row.names=1)
#expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
#gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1)
#expr$gene <- gene
# expr <- expr[!duplicated(expr$gene),]
# rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]
# expr[expr < 0] <- 0 # 对于这份泛癌数据，将略小于0的数值拉到0，否则不能取log（其他途径下载的泛癌数据可能不需要此操作）
# colnames(expr) <- substr(colnames(expr),1,15) # 非TCGA数据可忽略这行
# expr<-expr[-1,]

expr<-as.data.frame( t(as.matrix(log2(expr + 1))) )
genes<- res[["Sig.genes"]]#colnames(dat)[5:68]  #兴趣基因
expr<-expr[genes]
table(ifelse(as.numeric(substr(rownames(expr),14,15))<10,"Tumor","Normal" )) #N 59 T 515
expr$Sample <-ifelse(as.numeric(substr(rownames(expr),14,15))<10,"Tumor","Normal" )
expr$id <-substr(rownames(expr),1,12)
library(tidyverse)
expr<-expr %>% select(c("id","Sample") , everything())

for (i in genes){
  print(i)
  library(ggplot2);library(ggpubr)
  expr1<-expr[,c("id","Sample",i)]
  #expr1$Gene <-expr[,i]
  colnames(expr1)[3]<-"Gene"
  ggplot(expr1,#
         aes(Sample,Gene,color=Sample))+  #age_group
    geom_point(alpha=0.2,size=1.5,
               position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
    geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
                 size=0.4,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    scale_color_manual(values =  c("#EDA065","#66CCFF","#7EC7A7"))+ #blue c("#0066CC","#66CCFF","#9DC3E6")
    stat_compare_means(method = "wilcox.test",#t.test
                       size=5,show.legend= F,label = "p.signif",#"p.format",
                       label.x =0.75*2,
                       label.y =max(expr1$Gene)-0.06)+ #,label.y =max(df_gene$Gene)-0.55,label.x =1.5,label.y = 1.3  0.8 ,label.y = 4 #非参数检验kruskal.test not anova
  ylim(NA,max(expr1$Gene)+0.5) +
    labs(color="",title = i,y="mRNA expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
    theme_bw(base_family = "serif") +
    theme(text = element_text(size=16,family = "serif"),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
          panel.grid = element_blank(),legend.key = element_blank(),legend.position="none" ) 
    
  ggsave(paste0("./out/boxplot_TvsN_WILCOX_",i,".pdf"),width = 2,height = 3,family="serif")
}

# library(showtext) # Using Fonts More Easily in R Graphs
# #showtext_auto() 
# font_families()
