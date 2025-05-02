rm(list = ls())
setwd("~//R//my_projects/metabolism_sig//8.PCA和risk热图//risk热图")
#setwd("D:\\科研\\肺癌乳酸代谢\\8.PCA和risk热图\\risk热图")

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


#data----
#library(limma)
#riskFile="tcga_risk.csv"      #风险+表达文件 热图展示建模关键基因
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



cliFile="clinical.txt"        #临床数据文件
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
samSample=intersect(row.names(risk), row.names(cli))#316
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
#colnames(Type)=sigVec

#plot----
if(T){
  library(ComplexHeatmap)
  library(circlize)
  library(tidyverse)
  library(data.table)
  library(RColorBrewer)
}
#brewer.pal(8,"Dark2")
#brewer.pal(9,"Paired")[3:9]
#display.brewer.all()

#"#8DD3C7FF" "#BFE5BFFF" "#ECF7B7FF" "#F2EFBDFF" "#D8D3CDFF" "#C2B8D6FF"
# 颜色和选项配置
if(T){
  nejm <- c("pink","skyblue","#E18727FF","#20854EFF","#D8D3CDFF", "#C2B8D6FF","#FFDC91FF","#EE4C97FF")
  npg <- brewer.pal(11,"Spectral")#brewer.pal(9,"Set3")#(8,"Dark2")#c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF")
  nj=brewer.pal(9,"Paired")[3:9] #c('#313c63','#ebc03e','#5d84a4','#4B4B5D',"#EC7232","#F8F4EB")
  #定义热图注释的颜色
  colorList=list()
  #Type=Type[apply(Type,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
  bioCol=c(nejm,npg,nj)
  j=0
  for(cli in colnames(Type[,1:ncol(Type)])){
    cliLength=length(levels(factor(Type[,cli])))
    cliCol=bioCol[(j+1):(j+cliLength)]
    j=j+cliLength
    names(cliCol)=levels(factor(Type[,cli]))
    if("unknow" %in% levels(factor(Type[,cli]))){
      cliCol["unknow"]="grey75"}
    colorList[[cli]]=cliCol
  }
}

#ht_opt(legend_border='black', heatmap_border = TRUE, annotation_border = TRUE)
#colnames(Type)
# 创建热图注释
ha <- HeatmapAnnotation(df = Type,
                        annotation_name_gp = gpar(fontsize = 10, # 列标题大小
                                                  #fontface = "bold.italic", # 列标题字体
                                                  #fill = "#0E434F", # 列标题背景色
                                                  col = "black"),gap = unit(0.5, "mm"),
                        border = F, 
                        col = colorList)
#draw(ha)
# 由于示例缺少实际的热图数据，这里只提供了如何创建注释。你需要添加热图数据部分。
#range(t(exp))
#hist(as.matrix(exp), main="直方图：矩阵m中的数值分布", xlab="数值", ylab="频率", col="lightblue", border="black")
da=Heatmap(t(exp), 
           name = "Expression",  #热图1的图例名称
           column_split = Type$risk,
           cluster_rows = T, #关闭行聚类
           cluster_columns = F, #打开列聚类
           
           row_names_gp =gpar(fontsize = 10, # 列标题大小 6
                              #fontface = "bold.italic", # 列标题字体
                              #fill = "#0E434F", # 列标题背景色
                              col = "black" # 列标题颜色
                            
                              #border = "black"# 列标题边框色
           ),  #聚类树的位置
           row_labels = rownames(t(exp)),   #设置行名，默认为输入矩阵mat的行名，也可以自行修改
           show_row_names = T, 
           show_column_names = F,  #是否显示行名与列名
           col = colorRamp2(c(-2, 0, 2), c("#61AACF", "white", "#DA9599")),
             #circlize::colorRamp2(c(-2, 0, 2), c("#404B77", "#5F96B9", "#DAEACD")), #颜色范围及类型，超过此范围的值的颜色都将使用阈值所对应的颜色,也支持颜色向量：c('#7b3294','#c2a5cf','#f7f7f7','#a6dba0','#008837'
           border_gp = gpar(col = "black", lwd = 0.05), #设置热图的边框与粗度
           #row_gap = unit(1, "mm"), 
           column_gap = unit(1, "mm"), #设置行列分割宽度
           top_annotation = ha,
           ht_opt(legend_border='black',
                  heatmap_border = TRUE,
                  annotation_border = TRUE) 
) 

#output_file <- "heatmap-top10.pdf" #基因数目太多了！407
# 打开PDF图形设备
pdf(file = "./out/heatmap-240914.pdf", width = 10, height = 10,family="serif")#12 6
draw(da)
# 关闭图形设备
dev.off()

#重要性排序----


df <- impor %>% #read_delim(file = "DEG_result.txt", col_names = T, delim = "\t") %>%
  dplyr::mutate(rank = rank(importance, ties.method= "max")) #
df$Symbol<-rownames(df)
library(ggrepel)
df  %>% ggplot() + 
  geom_point(aes(x = rank, y = importance, color = importance, size =rank))+
  #scale_color_gradient2(low = "#fb8072",high = "#80b1d3",mid = "grey70", midpoint = mean(df$importance), name = "importance") 
 scale_color_gradient(low = "#80b1d3",high = "#fb8072", name = "importance") +
  geom_text_repel(data = df, #[1:10,]
                  aes(x = rank , y = importance, label = Symbol),
                  box.padding = 0.5,
                  #nudge_x = 10,
                  #nudge_y = -0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  direction = "y", 
                  hjust = "right") +
  theme_bw(base_size = 18) + 
  theme(
    panel.border = element_rect(linewidth = 1.25),
    #legend.background = element_roundrect(color = "#808080",linetype = 1),#边框
    axis.text = element_text(color = "#000000", size = 12.5),
    axis.title = element_text(color = "#000000", size = 15)
  )

ggsave(filename = "./out/rank-240914.pdf",height = 4.5,width = 6,family="serif")

#HR-重要性排序热图？----
#兴趣基因KM,ROC #单细胞 得分 高低 免疫
#HR重要性排序加入热图注释中#所有队列
#1.0 计算HR----
#准备输入数据
colnames(risk)[1:2]<-c("futime", "fustat")
risk<-risk[,-c(3:4)]
#cox
library(survival)
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(risk[,3:ncol(risk)])){ #trainset[,3:ncol(trainset)]
  cox <- coxph(Surv(futime, fustat) ~ risk[,i], data = risk)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(T){ #coxP<0.2
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR= coxSummary$conf.int[,"exp(coef)"] ,
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
colnames(outTab);str(outTab)
outTab[2:5]<-apply(outTab[2:5],2,as.numeric)
outTab[2:4]<-round(outTab[2:4],2)
outTab$Sign. <-ifelse(outTab$pvalue<0.001,"***",ifelse(outTab$pvalue<0.01,"**",ifelse(outTab$pvalue<0.05,"*","ns")))
outTab$pvalue<-round(outTab$pvalue,3)

#merge importance score
colnames(df)[3]<-"id"
#impor<-as.data.frame(res[["ml.res"]][[my_model]]$importance)  #[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
#colnames(impor)<-"importance"
impor<-arrange(impor,desc(importance) )#impor[order(impor$importance),] #降序
impor$rank<-c(1:nrow(impor))
impor$id <-rownames(impor)
#merge
outTab_impor <-merge(outTab,impor,by="id")
##annotation for pheatmap----
rownames(outTab_impor)<-outTab_impor$id #2:8
# va <- HeatmapAnnotation(df = outTab_impor[c(2,5:8)],
#                         annotation_name_gp = gpar(fontsize = 10, # 列标题大小
#                                                   #fontface = "bold.italic", # 列标题字体
#                                                   #fill = "#0E434F", # 列标题背景色
#                                                   col = "black"),gap = unit(0.5, "mm"),
#                         border = F #,col = colorList
#                         )
outTab_impor$HR.<-ifelse(outTab_impor$HR >1,"> 1","<1")


colorList2<-list("HR."= c("> 1"="#E84D94","<1"="#F4A7C1"),
    "pvalue"=c("*"="#BEE8E8","**"="#9FD1EE","***"="#4C9BE6","ns"="#98A9D0")
                                 )
outTab_impor<-outTab_impor[row.names(t(exp)),] #保持和表达矩阵一样顺序
va <- rowAnnotation(HR. = outTab_impor$HR., 
                    pvalue = outTab_impor$Sign.,
                    Importance =outTab_impor$importance,
                    col= colorList2
                    #Rank=outTab_impor$rank
)

db=Heatmap(t(exp), 
           name = "Expression",  #热图1的图例名称
           column_split = Type$risk,
           cluster_rows = F, #关闭行聚类!!!!!!!!!!!!!!!!!!不聚类否则对不上！
           cluster_columns = F, #打开列聚类
           row_names_gp =gpar(fontsize = 10, # 列标题大小 6
                              #fontface = "bold.italic", # 列标题字体
                              #fill = "#0E434F", # 列标题背景色
                              col = "black" # 列标题颜色
                              
                              #border = "black"# 列标题边框色
           ),  #聚类树的位置
           row_labels = rownames(t(exp)),   #设置行名，默认为输入矩阵mat的行名，也可以自行修改
           show_row_names = T, 
           show_column_names = F,  #是否显示行名与列名
           col = colorRamp2(c(-2, 0, 2), c("#61AACF", "white", "#DA9599")),
           #circlize::colorRamp2(c(-2, 0, 2), c("#404B77", "#5F96B9", "#DAEACD")), #颜色范围及类型，超过此范围的值的颜色都将使用阈值所对应的颜色,也支持颜色向量：c('#7b3294','#c2a5cf','#f7f7f7','#a6dba0','#008837'
           border_gp = gpar(col = "black", lwd = 0.05), #设置热图的边框与粗度
           #row_gap = unit(1, "mm"), 
           column_gap = unit(1, "mm"), #设置行列分割宽度
           top_annotation = ha,left_annotation =va,
           ht_opt(legend_border='black',
                  heatmap_border = TRUE,
                  annotation_border = TRUE) 
) 
pdf(file = "./out/heatmap-240914-HR-p.pdf", width = 12, height = 10,family="serif")#12 6
draw(db)
# 关闭图形设备
dev.off()

##annotation ggplot2----
outTab_impor  %>% ggplot() + 
  geom_vline(xintercept = 1,color="grey")+geom_hline(yintercept = 0,color="grey")+
  geom_point(aes(x = HR, y = importance, color = importance, size =HR ))+
  #scale_color_gradient2(low = "#fb8072",high = "#80b1d3",mid = "grey70", midpoint = mean(df$importance), name = "importance") 
  #scale_color_gradient(low = "#80b1d3",high = "#fb8072", name = "importance") +
  scale_color_viridis_c(alpha = 0.6)+
  geom_text_repel(data = outTab_impor , #[1:10,]
                  aes(x = HR , y = importance, label = id),
                  box.padding = 0.5,
                  #nudge_x = 10,
                  #nudge_y = -0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  direction = "y", 
                  hjust = "right") +
  
  theme_bw(base_size = 18) + 
  theme(
    panel.border = element_rect(linewidth = 1.25),
    #legend.background = element_roundrect(color = "#808080",linetype = 1),#边框
    panel.grid  = element_blank(),#panel.grid.minor
    axis.text = element_text(color = "#000000", size = 12.5),
    axis.title = element_text(color = "#000000", size = 15)
  )+ scale_x_continuous(limits = c(0.6, 1.5), breaks = seq(0.8, 1.5, 0.2))
#plot(outTab_impor$HR)
ggsave(filename = "./out/rank-240914-HR-impor.pdf",height = 4,width = 5.5,family="serif")
