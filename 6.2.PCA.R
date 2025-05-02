rm(list = ls())
setwd("~//R//my_projects//metabolism_sig//8.PCA和risk热图//pca")

#setwd("D:\\科研\\肺癌乳酸代谢\\8.PCA和risk热图\\pca")
if(T){
  library(Rtsne)
  library(ggplot2)
  library(patchwork)
  col=c("pink","skyblue")
}
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


##A.加载自己的数据----

#####合并 之前的关键基因表达信息和risk分组信息
#1.0 加载队列数据
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


#2.0加载riskscore
my_model <-"StepCox[forward] + RSF"
res <-readRDS("~//R//my_projects//metabolism_sig//5.1-101建模/out/101_results_178input.rds")
#取出my_model
#res1 <-res$ml.res[my_model] #自己单独计算？F2查看Mime1::ML.Dev.Prog.Sig
riskscore<-res[["riskscore"]][[my_model]] 

#3.0提出keygene表达矩阵,应该选择建模基因！
#key_gene <-read.csv("~/R/my_projects/metabolism_sig/5.确定基因/out/2024-08-18-确定基因.csv")
#key_gene <-as.character(key_gene$id)

key_gene <-res$Sig.genes #64
for (y in names(cohort) ) {  #1:length(cohort)
  cohort[[y]] <- cohort[[y]][,c("ID","futime","fustat",key_gene)] #"ID"
  cohort[[y]] <- cohort[[y]][cohort[[y]]$futime>0,]
} 

#4.0 riskscore 分组！
alldat <- data.frame()
for (i in names(cohort)) {
  dat <- riskscore[[i]]#cohort[[i]]
  dat$riskscore <- dat$RS#risk_score
  colnames(dat)[c(2,3)]<-c("futime","fustat")
  # dat$Risk <- ifelse(dat$riskscore>surv_cutpoint(dat, #数据集
  #                                                time = "futime", #生存状态
  #                                                event = "fustat", #生存时间
  #                                               variables = c("riskscore") )$cutpoint[,1] ,"High","Low")
  
  dat$risk <- ifelse(dat$riskscore>median(dat$riskscore) ,"High","Low")
  #加入队列信息
  dat$cohort <-i
  dat<-merge(dat,cohort[[i]][,-c(2:3)],by="ID")  #合并表达矩阵
  alldat <- rbind(alldat,dat)
  bioPCA(data=dat, pcaFile=paste0("./out/",i,"_PCA.pdf") ) #PCA函数
}
write.csv(alldat,"meta-cohort.csv")#保存
table(alldat$cohort)

##B.plot function 重新写----
bioPCA=function(data=null, pcaFile=null){
  #读取输入文件,提取数据
  #inputFile="gse13213_sva_risk.csv"
  #rt=read.csv(inputFile, row.names=1)
  #data=dat;rownames(data)<-dat$ID
  data=dat[,c(8:(ncol(dat)))];rownames(data)<-dat$ID
  risk=dat[,"risk"]
  
  #PCA分析
  data.pca=prcomp(data, center = TRUE,scale. = TRUE)
  summ<-summary(data.pca)
  #提取PC score；
  df1<-data.pca$x
  df1<-data.frame(df1,risk=risk)
  #提取主成分的方差贡献率，生成坐标轴标题；
  xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
  ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
  #计算质心，也就是“中心点”，centroid (average) position for each group；
  centroid <- aggregate(df1[,1:4],by=list(df1$risk), FUN = mean)
  #centroid
  #合并两个数据框；
  coords <- merge(df1, centroid,
                  by.x = "risk",
                  by.y ="Group.1",
                  all.x = T,
                  suffixes = c("",".centroid"),)
  #head(coords)
  #使用ggplot2绘制散点图；
  pdf(file=pcaFile, height=6/1.5, width=7/1.5)       #保存输入出文件
  
  p <- ggplot(data = coords, aes(x = PC1, y = PC2))+
    labs(x = xlab, y = ylab)+
    geom_hline(yintercept = 0,lty="dashed",colour="grey60",size=0.4)+
    geom_vline(xintercept = 0,lty="dashed",colour="grey60",size=0.4)+
    stat_ellipse(aes(fill=risk),
                 type = "norm", geom ="polygon",
                 alpha=0.2,color=NA,show.legend = F)+
    # geom_curve(aes(xend= PC1.centroid, yend= PC2.centroid, colour = risk),
    #            angle=150,curvature = -0.5,show.legend = FALSE)+
    geom_point(aes(fill = risk,colour = risk),
               shape = 21, size = 1.5, #colour="black",
               show.legend = FALSE)+
    scale_fill_manual(values = col,name="Group")+
    scale_colour_manual(values = col,guide=F)+
    theme_bw(base_size = 14)+
    theme(#axis.text= element_text(face = "bold.italic",colour = "#441718",size = 12,hjust = 0.5),
          #axis.line = element_blank(),
          #plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16,hjust = 0.5),
          #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
          #legend.text = element_text(face ="bold"),
          #panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
          #panel.background = element_rect(fill = "#F8F4EB"),
          panel.grid= element_blank(),
          #legend.title = element_text(face ="bold.italic",size = 13),
          legend.position = "right"
          #axis.title = element_blank() ,
    )+
    geom_label(data = centroid, aes(label = Group.1, fill = Group.1),
               size = 4,hjust=0.5, color="#F8F4EB",show.legend = FALSE)
  
  p1 <- ggplot(data = coords) +
    geom_density(aes(x = PC1, fill=risk,colour = risk), #同样，颜色区分实验/对照，线型区分队列地区；
                  alpha = 0.5, position = 'identity') +
    scale_fill_manual(values =col[1:2])+theme_bw()+
    theme(axis.text= element_blank(),
          axis.line = element_blank(),
          axis.ticks =element_blank(), 
          #plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16,hjust = 0.5),
          title = element_blank(),
          #legend.text = element_text(face ="bold.italic"),
          #panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
          #panel.background = element_rect(fill = "#F8F4EB"),
          panel.grid= element_blank(),
          #legend.title = element_text(face ="bold.italic",size = 13),
          legend.position = "none" #right
          #axis.title = element_blank() ,
    )
  p2 <- ggplot(data = coords) +
    geom_density(aes(x = PC2, fill=risk,colour = risk), #同样，颜色区分实验/对照，线型区分队列地区；
                 alpha = 0.5, position = 'identity') +
    scale_fill_manual(values =col[1:2])+theme_bw()+
    theme(axis.text= element_blank(),
          axis.line = element_blank(),
          axis.ticks =element_blank(), 
          #plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16,hjust = 0.5),
          title = element_blank(),
          #legend.text = element_text(face ="bold.italic"),
          #panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
          #panel.background = element_rect(fill = "#F8F4EB"),
          panel.grid= element_blank(),
          #legend.title = element_text(face ="bold.italic",size = 13),
          legend.position = "none" #right
          #axis.title = element_blank() ,
    )+coord_flip()
  # p2 <- ggplot(data = coords) +
  #   geom_density(aes(x = PC2, fill=risk,colour = risk), #同样，颜色区分实验/对照，线型区分队列地区；
  #                color = 'black', alpha = 1, position = 'identity') +
  #   scale_fill_manual(values =col[1:2])+
  #   theme(axis.text= element_blank(),
  #         axis.line = element_blank(),
  #         axis.ticks =element_blank(), 
  #         #plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16,hjust = 0.5),
  #         title = element_blank(),
  #         #legend.text = element_text(face ="bold.italic"),
  #         #panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
  #         #panel.background = element_rect(fill = "#F8F4EB"),
  #         panel.grid= element_blank(),
  #         #legend.title = element_text(face ="bold.italic",size = 13),
  #         legend.position = "right"
  #         #axis.title = element_blank() ,
  #   )+coord_flip() #这张图要拼在pca图的右侧，因此需要翻转坐标轴
  pp <- p1 + plot_spacer() + p + p2 +
    plot_layout(guides = 'collect') &
    theme(legend.position='right') #合并图例，并将其放在右侧
  f=pp + plot_layout(widths = c(4, 1), heights = c(1,4))
  print(f)
  dev.off()
}



# file=list.files(pattern = ".csv")
# alldat <- data.frame()
# for (i in file) {
#   dat <- read.csv(i,header = T,row.names = 1)
#   alldat <- rbind(alldat,dat)
# }
#write.csv(alldat,"meta-cohort.csv")



# bioPCA=function(inputFile=null, pcaFile=null){
#   #读取输入文件,提取数据
#   #inputFile="gse13213_sva_risk.csv"
#   rt=read.csv(inputFile, row.names=1)
#   data=rt[,c(3:(ncol(rt)-2))]
#   risk=rt[,"risk"]
#   
#   #PCA分析
#   data.pca=prcomp(data, center = TRUE,scale. = TRUE)
#   summ<-summary(data.pca)
#   #提取PC score；
#   df1<-data.pca$x
#   df1<-data.frame(df1,risk=risk)
#   #提取主成分的方差贡献率，生成坐标轴标题；
#   xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
#   ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
#   #计算质心，也就是“中心点”，centroid (average) position for each group；
#   centroid <- aggregate(df1[,1:4],by=list(df1$risk), FUN = mean)
#   #centroid
#   #合并两个数据框；
#   coords <- merge(df1, centroid,
#                   by.x = "risk",
#                   by.y ="Group.1",
#                   all.x = T,
#                   suffixes = c("",".centroid"),)
#   #head(coords)
#   #使用ggplot2绘制散点图；
#   pdf(file=pcaFile, height=6/1.5, width=7/1.5)       #保存输入出文件
#   p <- ggplot(data = coords, aes(x = PC1, y = PC2))+
#     labs(x = xlab, y = ylab)+
#     geom_hline(yintercept = 0,lty="dashed",colour="grey60",size=0.4)+
#     geom_vline(xintercept = 0,lty="dashed",colour="grey60",size=0.4)+
#     stat_ellipse(aes(fill=risk),
#                  type = "norm", geom ="polygon",
#                  alpha=0.0,color=NA,show.legend = F)+
#     geom_curve(aes(xend= PC1.centroid, yend= PC2.centroid, colour = risk),
#                angle=150,curvature = -0.5,show.legend = FALSE)+
#     geom_point(aes(fill = risk),
#                colour="black",shape = 21, size = 2.5,
#                show.legend = FALSE)+
#     scale_fill_manual(values = col,name="Group")+
#     scale_colour_manual(values = col,guide=F)+
#     theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 12,hjust = 0.5),
#           axis.line = element_blank(),
#           plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16,hjust = 0.5),
#           title = element_text(face = "bold.italic",colour = "#441718",size = 13),
#           legend.text = element_text(face ="bold.italic"),
#           panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
#           panel.background = element_rect(fill = "#F8F4EB"),
#           panel.grid= element_blank(),
#           legend.title = element_text(face ="bold.italic",size = 13),
#           legend.position = "right"
#           #axis.title = element_blank() ,
#     )+
#     geom_label(data = centroid, aes(label = Group.1, fill = Group.1),
#                size = 4,hjust=0.5, color="#F8F4EB",show.legend = FALSE)
#   
#   p1 <- ggplot(data = coords) +
#     geom_density(aes(x = PC1, fill=risk), #同样，颜色区分实验/对照，线型区分队列地区；
#                  color = 'black', alpha = 1, position = 'identity') +
#     scale_fill_manual(values =col[1:2])+
#     theme(axis.text= element_blank(),
#           axis.line = element_blank(),
#           axis.ticks =element_blank(), 
#           plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16,hjust = 0.5),
#           title = element_blank(),
#           legend.text = element_text(face ="bold.italic"),
#           panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
#           panel.background = element_rect(fill = "#F8F4EB"),
#           panel.grid= element_blank(),
#           legend.title = element_text(face ="bold.italic",size = 13),
#           legend.position = "right"
#           #axis.title = element_blank() ,
#     )
#   p2 <- ggplot(data = coords) +
#     geom_density(aes(x = PC2, fill=risk), #同样，颜色区分实验/对照，线型区分队列地区；
#                  color = 'black', alpha = 1, position = 'identity') +
#     scale_fill_manual(values =col[1:2])+
#     theme(axis.text= element_blank(),
#           axis.line = element_blank(),
#           axis.ticks =element_blank(), 
#           #plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16,hjust = 0.5),
#           title = element_blank(),
#           legend.text = element_text(face ="bold.italic"),
#           panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
#           panel.background = element_rect(fill = "#F8F4EB"),
#           panel.grid= element_blank(),
#           legend.title = element_text(face ="bold.italic",size = 13),
#           legend.position = "right"
#           #axis.title = element_blank() ,
#     )+
#     coord_flip() #这张图要拼在pca图的右侧，因此需要翻转坐标轴
#   pp <- p1 + plot_spacer() + p + p2 +
#     plot_layout(guides = 'collect') &
#     theme(legend.position='right') #合并图例，并将其放在右侧
#   f=pp + plot_layout(widths = c(4, 1), heights = c(1,4))
#   print(f)
#   dev.off()
# }
# file=list.files(pattern = ".csv")
# for(i in file){
#   #i="gse13213_sva_risk.csv"
#   name=sapply(strsplit(i,"\\."), "[",1)
#   bioPCA(inputFile=i, pcaFile=paste(name,"PCA.pdf"))
# }

