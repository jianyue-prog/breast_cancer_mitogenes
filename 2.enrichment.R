
if(T){
  rm(list = ls())
  gc()
  setwd("~//R//my_projects//metabolism_sig//2.enrichment")
  library(data.table)
  library(limma)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
}

#输出文件夹
folder_path <- "./out" # 检查文件夹是否存在
if (!dir.exists(folder_path)) {   # 如果文件夹不存在，则创建它
  dir.create(folder_path)
  print(paste("Folder created at", folder_path))
} else {
  print(paste("Folder already exists at", folder_path))
}

genelist <- readxl::read_xls("~/R/my_projects/metabolism_sig/1.scRNA/Human.MitoCarta3.0.xls",sheet = 2) #read.csv("all染色质基因.csv",header =T,fill = T)
#genelist <- genelist$x
#genelist <- genelist [!genelist=="" ]
genelist <- unique(genelist$Symbol )#1136
rt <- fread("luad_mRNA.txt",header = T,data.table = F)
rt <- column_to_rownames(rt,var = "Tag")
type <- ifelse(as.numeric(substr(colnames(rt),14,15))<10,"Tumor","Normal")
max(rt)
#########
rt <- log2(rt+1)
list <- factor(type)
list <- model.matrix(~factor(list)+0)
colnames(list) <- c("Normal","Tumor")
df.fit <- lmFit(rt, list)
df.matrix <- makeContrasts(Tumor - Normal , levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
tempOutput$id <- rownames(tempOutput)
tempOutput <- tempOutput[tempOutput$id%in%genelist,]
###########
df <- tempOutput 
df$Group <- factor(ifelse(df$adj.P.Val < 0.05 & abs(df$logFC) >=1,
                          ifelse(df$logFC >=1, 'Up','Down'),'N.S.'))
write.csv(df,"./out/差异基因.csv")
############
gene <- rownames(df[abs(df$logFC) >2,])#c("NPM3","DCAF13","RPP40")
p1 <- ggplot(df,aes(logFC, -log10(adj.P.Val)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "skyblue")+
  #geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "skyblue")+
  geom_point(aes(size = AveExpr , #-(adj.P.Val)
                 color = AveExpr , #-log10(adj.P.Val)
                 shape=Group
                 ))+
  #scale_color_gradientn(values = seq(0,1,0.2), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_colour_viridis_c(direction = -1)+ #option=
  scale_size_continuous(range = c(0,1))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.justification = c(0,1))+
  # 设置图例
  guides(col = 
           guide_colorbar(title = "AveExpr" ,#"-Log10 adj.Pval",
                          ticks.colour = NA,
                          reverse = T,
                          title.vjust = 0.8,
                          barheight = 8,
                          barwidth = 1),
         size = "none") +
  # 添加标签：
  geom_text_repel(data = df[df$id%in%gene,],
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  # 这里的filter很关键，筛选你想要标记的基因
                  aes(label =id),
                  size = 3, 
                  #box.padding = 0.5, #字到点的距离 
                  #point.padding = 0.8, #字到点的距离，点周围的空白宽度 
                  min.segment.length = 0.1, 
                  segment.size=0.3, 
                  color = 'black') +
  labs(x="Log2(FC)",y="-Log10(adj.Pval)")

p1
ggsave("./out/volcano.pdf",family="serif",width = 5,height = 4)
#标记基因类型数目？

library(ggplot2)
df<-read.csv("./out/差异基因.csv",row.names = 1)
aflogFC <- 1 ;adjPval <- 0.05    
#取基因
genelist <- readxl::read_xls("~/R/my_projects/metabolism_sig/1.scRNA/Human.MitoCarta3.0.xls",sheet = 2)
genelist <-unique(genelist$Symbol) #1132
#标记线粒体基因
#df$Mitogene <-ifelse((df$id %in% genelist)==T,"Yes","No")
#差异表达线粒体基因
Symbol<-intersect(genelist,df[df$Group !="N.S.",]$id)#244
for_label_2<-df[df$id %in% Symbol,]

ggplot(df,aes(logFC, -log10(adj.P.Val),shape = Group ) )+ #
  geom_point(aes(col=Group, ),alpha=0.8,size=3)+ #, size= abs(log2FoldChange) 
  #scale_color_manual(values=c(ggsci::pal_nejm()(2)[2], "#838B8B", ggsci::pal_nejm()(2)[1]))+
  scale_color_manual(values=c("skyblue","grey","pink") )+
  labs(title = " ",x="log2(FoldChange)",y="-log10(adj.pvalue)" )+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),legend.position = "top")+
  geom_hline(aes(yintercept=-log10(adjPval)), colour="black", linetype=2,size=0.5)+
  geom_vline(aes(xintercept=aflogFC), colour="black", linetype=2,size=0.3)+
  geom_vline(aes(xintercept=-aflogFC), colour="black", linetype=2,size=0.3)+
  #geom_vline(aes(xintercept=-2), colour="black", linetype=2,size=0.3)+
  #geom_vline(aes(xintercept=2), colour="black", linetype=2,size=0.3)+
  theme_bw(base_size = 14)+#ylim(0,4.5)+
  
  ggrepel::geom_text_repel(
    aes(label = Symbol), #Symbol
    data = for_label_2, #[abs(for_label_2$logFC) >2,]
    color="black",max.overlaps=80,
    #label.size =0.01, 
    #max.overlaps = Inf,
    segment.size=0.3, box.padding = 0.6#,force = 10
  )
#ggsave("./out/volcano_diff_mito.pdf",family="serif",width = 5,height = 4)

#KEGG----
id <- df$id[!df$Group=="N.S."]
library(clusterProfiler)
trans = bitr(id, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

ego <- enrichKEGG(
  gene          = trans$ENTREZID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)

ego <-setReadable(ego,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
ego <- ego@result
ego <- ego [order(ego$Count,decreasing = T),]
#ego$Count

#ego <- head(ego,15)
#ego <- ego[,c(2,5,7,9)]
ego <-ego[ego$pvalue < 0.01,]
library(ggplot2);library(ggsci);library("scales")
ego <-ego[-c(1:11),]
# show_col(pal_npg("nrc")(3))
# colorRampPalette(c("#778997","red"),space = "rgb")(20)
# ego$Description <- factor(ego$Description,levels = c(ego$Description))
# windowsFonts()
# Fon <- 'serif'
# kegg_plot <- ggplot(data=ego, aes(x=Count, y=Description))+
#   geom_col(aes(color=qvalue))+ #fill = "#9c393a",color="black"
#   theme_classic()+ylab("")+ 
#   theme(axis.title.y = element_text(size = 12,colour = "black", face='bold'))
# kegg_plot

#
mytheme<- theme(axis.title = element_text(size = 13),
                axis.text = element_text(size = 11),
                plot.title = element_text(size = 14,
                                          hjust= 0.5,
                                          face= "bold"),
                legend.title = element_text(size = 13),
                legend.text = element_text(size = 11))

p<- ggplot(data = ego,
           aes(x = Count, y = reorder(Description,-Count), fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "RdPu",direction = 1) + #RdPu YlOrRd
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "Number of Genes", y = "", title = "KEGG pathway enrichment") + 
  theme_bw()+ mytheme
p
ggsave("./out/KEGG_DEG_my.pdf",width = 6,height = 5,family="serif")

#线粒体基因/基因集富集KEGG

  library(clusterProfiler)
  trans1 = bitr(genelist, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  ego1 <- enrichKEGG(
    gene          = trans1$ENTREZID,
    keyType     = "kegg",
    organism   = 'hsa',
    pvalueCutoff      = 0.05,
    pAdjustMethod     = "BH",
    qvalueCutoff  = 0.05
  )
  ego1 <- ego1@result
  ego1 <- ego1 [order(ego$Count,decreasing = T),]
  #ego1 <- ego1[,c(2,5,7,9)]
  ego1 <-ego1[ego1$pvalue < 0.001,]
  library(ggplot2)#;library(ggsci);library("scales")
  mytheme<- theme(axis.title = element_text(size = 13),
                  axis.text = element_text(size = 11),
                  plot.title = element_text(size = 14,
                                            hjust= 0.5,
                                            face= "bold"),
                  legend.title = element_text(size = 13),
                  legend.text = element_text(size = 11))
  
  p<- ggplot(data = ego1,
             aes(x = Count, y = reorder(Description,-Count), fill = -log10(pvalue)))+
    scale_fill_distiller(palette = "BuGn",direction = 1) + #RdPu YlOrRd
    #RColorBrewer::display.brewer.all()
    geom_bar(stat = "identity", width = 0.8) +
    labs(x = "Number of Genes", y = "", title = "KEGG pathway enrichment") + 
    theme_bw()+ mytheme
  p
  ggsave("./out/KEGG_genelist_1.pdf",width = 6.5,height = 5,family="serif")


#取交集
intersect(ego$Description,ego1$Description)
#50

#GO----
#GO<-enrichGO(gene=trans$ENTREZID,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",ont="ALL",qvalueCutoff = 0.05,readable = T) 
# GO1 <- GO@result

#丑！ 
# GO <- GO1[,c(1,3,8,10)]
# library(dplyr)
# GO <- GO[order(GO$Count,decreasing = T),]
# eGoBP <- GO  %>% 
#   filter(ONTOLOGY=="BP") %>%
#   filter(row_number() >= 1,row_number() <= 5)
# eGoCC <- GO  %>% 
#   filter(ONTOLOGY=="CC") %>%
#   filter(row_number() >= 1,row_number() <= 5)
# eGoMF <- GO  %>% 
#   filter(ONTOLOGY=="MF") %>%
#   filter(row_number() >= 1,row_number() <= 5)
# eGo10 <- rbind(eGoBP,eGoMF,eGoCC)
# eGo10$Description <- factor(eGo10$Description,levels = c(eGo10$Description))
# go_plot <-  ggplot(data=eGo10, aes(x=Count, y=Description))+geom_col(aes(fill = ONTOLOGY),color="black")+scale_fill_manual(values=c("#ea8c7e","#9ddae9","#5eab9e"))+theme_classic()+theme(legend.position = "top")+ylab("")
# kegg_plot <- kegg_plot+theme(legend.position = "none")
# library(patchwork)
# 
# 
# pdf("富集分析.pdf",width =8/1.5,height =12/1.5)
# kegg_plot +go_plot+plot_layout(nrow = 2)
# dev.off()

##流星图
#https://mp.weixin.qq.com/s/XVSOR0Gr43TrCylYrdo5PQ

library(clusterProfiler)
library(dplyr)
#GO<-enrichGO(gene=trans$ENTREZID,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",ont="ALL",qvalueCutoff = 0.05,readable = T) 
eego <- enrichGO(gene          = trans$ENTREZID,
                 #universe      = names(geneList),
                 OrgDb         = "org.Hs.eg.db",
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 #pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

df_go <- data.frame(eego) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 6) %>%
  arrange(desc(pvalue))
ratio <- lapply(df_go$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df_go$ratio <- ratio
df_go$Description <- factor(df_go$Description,levels = df_go$Description)
#plot
ggplot(df_go) +
  ggforce::geom_link(aes(x = 0,y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = stat(index),
                         color = ONTOLOGY,
                         size = 10),#after_stat(index)),#流星拖尾效果
                     n = 500,
                     #color = "#FF0033",
                     show.legend = F) +
  geom_point(aes(x = -log10(pvalue),y = Description),
             color = "black",
             fill = "white",size = 4,shape = 21) +
  geom_line(aes(x = ratio*100,y = Description,group = 1),
            orientation = "y",linewidth = 1,color = "#FFCC00") + #线条-比例 #FFCC00 #E59CC4
  scale_x_continuous(sec.axis = sec_axis(~./100,
                                         labels = scales::label_percent(),
                                         name = "Percent of geneRatio")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        axis.text = element_text(color = "black")) +
  ylab("") + xlab("-log10 Pvalue") +
  facet_wrap(~ONTOLOGY,scales = "free",ncol = 1) + #3
  scale_color_brewer(palette = "Set1")#OrRd Set1 Paired
  #scale_color_manual(values = c('#CCE0F5', '#CCC9E6', '#625D9E') )
my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 

ggsave("./out/GO.pdf",family="serif",width = 4,height = 6)
#Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
  
#热图-1----
rt <- rt[rownames(rt)%in%id,]
Type=type
names(Type)=colnames(rt)
Type=as.data.frame(Type)

data=rt
data=na.omit(data)
data1=scale(t(data), center = TRUE, scale = TRUE)##对矩阵的列进行标准化
data1=t(data1)
##
# normalize=function(x){
#   return((x-min(x))/(max(x)-min(x)))}
# #对score进行矫正
# data2=normalize(data1)
# 
# 
# library(dendextend)
# library(circlize)
# column_dend = as.dendrogram(hclust(dist(data2)))
# column_dend = color_branches(column_dend, k = 3) # `color_branches()` returns a dendrogram object
# column_dend = dendrapply(column_dend, function(d) {
#   if(runif(1) > 0.5) attr(d, "nodePar") = list(cex = 0.5, pch = sample(20, 1), col = rand_color(1))
#   return(d)
# })
# hist(data1)
# range(data1)

###丑
#library(ComplexHeatmap)
#pdf("热图.pdf",width = 10/1.5,20/1.5)
# Heatmap(data1, 
#         name = "Expression",  #热图1的图例名称
#         cluster_rows = T, #关闭行聚类
#         cluster_columns = F, #打开列聚类
#         show_row_dend = F,
#         #row_km = 2,
#         row_split = 2,
#         #row_dend_side = "left",
#         # row_names_gp =gpar(fontsize = 6, # 列标题大小
#         #                    fontface = "bold.italic" # 列标题字体
#         #                    #fill = "#0E434F", # 列标题背景色
#         #                    #col = "black", border = "black"# 列标题边框色
#         # ),  #聚类树的位置
#         show_row_names = T, 
#         show_column_names = F,  #是否显示行名与列名
#         col = circlize::colorRamp2(c(-1, 0, 1), c("#404B77", "#5F96B9", "#DAEACD")), #颜色范围及类型，超过此范围的值的颜色都将使用阈值所对应的颜色,也支持颜色向量：c('#7b3294','#c2a5cf','#f7f7f7','#a6dba0','#008837')
#         na_col = "black", #若矩阵中存在NA值，特定设置为黑色
#         border_gp = gpar(col = "black", lwd = 2), #设置热图的边框与粗度
#         #rect_gp = gpar(col = "white", lwd = 1), #设置热图中格子的边框与粗度
#         #row_split = afsplit$Type,        #按照Type对热图行进行分割
#         column_split = Type$Type,#生成24个模拟分类数据对热图列进行分割
#         row_gap = unit(1, "mm"), 
#         column_gap = unit(1, "mm"), #设置行列分割宽度
# 
#         top_annotation = HeatmapAnnotation(typet1 = anno_block(gp = gpar(fill = c(2,7)), 
#                                                                labels = c("Normal", "Tumor"),#
#                                                                labels_gp = gpar(col = "white", fontface = "bold.italic",
#                                                                                 fontsize = 10))) 
# ) 

#dev.off()


#热图-2----
#自己画的

library(pheatmap);library(RColorBrewer);library("gplots")
#pheatmap(data1[,1:6])
ann_colors=list(#FDR=c('[0,0.05]'='#FDDCA9','(0.05,0.5]'="#DDC9E3"),
                Type=c("Normal"="skyblue","Tumor"="pink" ))

# hallmark_p <-pheatmap(data1,cluster_rows = T,cluster_cols  = F, show_colnames = F,show_rownames = T,
#                       #scale="row",#标准化样本
#                       cutree_rows = 2, cutree_cols = NA, 
#                       border_color = NA,#"white",#cellwidth=10,cellheight = 10,
#                       gaps_col= as.numeric(table(Type$Type)[2]), #table(exp_marker1$response)
#                       color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(10) ), #colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
#                       ##PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 Spectral PRGn RdYlBu
#                       #color=bluered(50),#library("gplots")
#                       annotation_col=Type,annotation_names_col = F,annotation_colors =ann_colors #,annotation_row = anno_row
# ) #slow

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

#pdf("./out/hallmark_pheatmap.pdf",family = "serif",width = 6.66,15) #20/1.5
pheatmap(data1,cluster_rows = T,cluster_cols  = F, show_colnames = F,show_rownames = T,
         #scale="row",#标准化样本
         cutree_rows = 2, cutree_cols = NA, fontsize =6 ,
         border_color = NA,#"white",#cellwidth=10,cellheight = 10,
         gaps_col= as.numeric(table(Type$Type)[2]), #table(exp_marker1$response)
         color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(10) ), #colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
         ##PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 Spectral PRGn RdYlBu
         #color=bluered(50),#library("gplots")
         annotation_col=Type,annotation_names_col = F,annotation_colors =ann_colors #,annotation_row = anno_row
)
#dev.off()

