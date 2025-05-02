#renv::deactivate()



if(T){
  rm(list = ls())
  gc()
  
  setwd("~//R//my_projects//metabolism_sig//1.scRNA")
  
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(clustree)
  library(cowplot)
  library(dplyr)
  library(data.table)
  library(tibble)
  library(SCP)
}

#LOAD RDATA！
sceList <-qread("./sceList.qs")#


##0.单细胞读取----
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
rt <- fread("GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",header = T,data.table = F)
ann <- fread("GSE131907_Lung_Cancer_cell_annotation.txt.gz",data.table = F,header = T)
length(table(ann$Cell_subtype))
rt <- column_to_rownames(rt,var="Index")
ann  <- column_to_rownames(ann ,var="Index")
sceList=CreateSeuratObject(counts = rt  ,
                           project ="LUAD",
                           min.cells = 0,
                           min.features =0,
                           meta.data = ann)

# library(BiocParallel)
# register(MulticoreParam(workers = 30, progressbar = TRUE))

##1.0 单细胞数据处理----
sceList <- Standard_SCP(srt = sceList) # 降维聚类 slow 15:07 多核？  1h; 30 cores  1h无法限制？


table(sceList$Cell_type)
sceList$Cell_type <- ifelse(sceList$Cell_type=="Epithelial cells"&sceList$Cell_subtype=="Malignant cells","Malignant cells",sceList$Cell_type)
#?CellDimPlot

p1 <- CellDimPlot(
  srt = sceList, group.by = c("Standardclusters"),
  reduction = "StandardUMAP2D", theme_use = "theme_blank"
  ,show_stat =F)
p2 <- CellDimPlot(
  srt = sceList, group.by = c( "Cell_type"),
  reduction = "StandardUMAP2D", theme_use = "theme_blank"
  ,show_stat =F)  #重新画，点大小？

##打分
#genelist <- read.csv("all染色质基因.csv",header =T,fill = T)
#genelist <-  genelist$x
#替换基因集1136 线粒体基因 https://personal.broadinstitute.org/scalvo/MitoCarta3.0/human.mitocarta3.0.html
genelist <- readxl::read_xls("Human.MitoCarta3.0.xls",sheet = 2)
genelist <-  genelist$Symbol
genelist <- intersect(rownames(sceList),genelist) 
#sceList <-qread("./sceList.qs")#read_sav("./sceList.qs") #加载之前处理好的
genelist <- list(genelist)#1053
# 
# genelist <- genelist [!genelist=="" ]
# genelist <- unique(genelist )
# genelist <- intersect(rownames(sceList),genelist)  #与单细胞基因取交集
# gene <- as.list(genelist) #1053

#LOAD
#library(haven);library(readr)
sceList <- AddModuleScore(object = sceList,features = genelist,#gene,
                          name = "score" #ctrl = 100,
                          ) #slow 慢的要死



colnames(sceList@meta.data)[14] <- 'Score'
#sceList$score1
p3 <- FeatureDimPlot(sceList, features = "Score",show_stat = F,theme_use = "theme_blank")+theme(legend.position = "right")
p1+p2+p3+patchwork::plot_layout(ncol = 3,widths = c(0.4,0.4,0.4))
#ggsave("乳酸评分-20240805.pdf",height = 60*2,width = 210*2,units = "mm")
ggsave(paste0("score-",substr(Sys.time(),1,10),".pdf"),height = 60*2,width = 210*2,units = "mm")


#qsave(sceList,"sceList.qs")

folder_path <- "./out"
# 检查文件夹是否存在
if (!dir.exists(folder_path)) {
  # 如果文件夹不存在，则创建它
  dir.create(folder_path)
  print(paste("Folder created at", folder_path))
} else {
  print(paste("Folder already exists at", folder_path))
}

##？
my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
color_used = c("#B84D64","#864A68","#E32D32","#5E549A","#8952A0","#384B97","#911310",
               "#7CA878","#35A132","#6B70B0","#20ACBD","#959897",'skyblue',
               "#F4A2A3","#B6CCD7","#AF98B5","#E01516","#FFFFFF")

# colnames(sceList@meta.data)
# pdf("./out/umap_clutser.pdf",width = 4,height = 4)
# p1<- DimPlot(sceList, reduction = "StandardUMAP2D",raster=T,label = T,group.by = "Cell_type",
#              cols = my36colors ) +
#   NoLegend()+NoAxes()
# p1
# dev.off()
# 
# ##dotplot
# table(sceList@meta.data$Cell_type) #11
# #sceList@meta.data$ident <-sceList@meta.data$Cell_type
# Idents(sceList)<-sceList$Cell_type#@meta.data$Cell_type Cell_type有NA?
# #替换NA
# 
# DotPlot(sceList , features =c("Score"),group.by  ="Cell_type",# "Sample_Origin
#         cols =dittoColors()# c("skyblue", "pink")
# )
# 
# #所有基因太多，线粒体marker基因
# TOMM <-unlist(genelist)[grep("TOMM", unlist(genelist) )] #tomm开头基因
# DotPlot(sceList , features =TOMM) #unlist(genelist)
# 
# library(dittoSeq)
# DotPlot(sceList, 
#         features = c("Score"), 
#         split.by = 'Cell_type',
#         cols = dittoColors())
# 
# library(RColorBrewer)
# DotPlot(sceList , features ="Score",group.by ="Cell_type")+
#   scale_color_viridis(8000, option = "G",alpha = 0.7)
#   #scale_colour_viridis_c(direction = -1)
# 
# DotPlot(sceList , features =c("Score"),split.by ="Cell_type",# "Cluster",#"Major.celltype"#,
#         cols = my36colors#c("skyblue", "pink")
# )+
#   #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
#   #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
#   labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
#   theme_bw()+
#   theme(#panel.grid = element_blank(), 
#     axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
# #ggsave("dotplot_major.pdf",family="serif",width = 3,height = 3)#保存3,3  4-12
# 
# rev(colorRampPalette(brewer.pal(11, "Spectral"))(100) )


##3.0 20240815重新整理细胞类型，有NA！----
#Cell_type 中NAT替换为undetermined
meta<-sceList@meta.data
table(meta$Cell_type)
table(is.na(meta$Cell_type)) #1628
meta$Cell_type <-ifelse(is.na(meta$Cell_type) == T,"Undetermined",meta$Cell_type)
#Cell_type.refined
table(is.na(meta$Cell_type.refined)) #28437
table(meta$Cell_type.refined)
meta$Cell_type.refined<-ifelse(is.na(meta$Cell_type.refined) == T,"Undetermined",meta$Cell_type.refined)

sceList@meta.data <-meta
table(sceList@meta.data$Cell_type)
table(sceList@meta.data$Cell_type.refined)
#删除undetermined细胞？上皮细胞！
#DimPlot(sceList, reduction = "StandardUMAP2D",raster=T,label = T,group.by = "Cell_type",cols = my36colors )

##是否恶性？合并后识别恶性细胞！----

#样本信息热图？
table(meta$Sample)
table(meta$Sample_Origin)

##marker DOTPLOT----
colnames(meta)
markers <-c("EPCAM","CDH1","KRT7","KRT19",
            "PECAM1","VWF",
            "PDGFRA","LUM","COL1A1",
            "MYLK","ACTA2","PDGFRB",
            "MKI67","STMN1",
            "LYZ","CD1C","CD163","CD14",
            "S100A9","S100A8","CSF3R",
            "KIT","GATA2","CPA3","MS4A2",
            "MS4A1","CD79A","CD19",
            "IGHG1","MZB1","JCHAIN","SDC1",
            "IGHD","FCER2","TCL1A",
            "NKG7","KLRD1","KLRB1","NCR1",
            "CD3D","CD3E","CD4","CD8A",
            "FOXP3","TNFRSF4","IKZF2","PDCD1","TIGIT"
            ) #原文https://onlinelibrary.wiley.com/doi/10.1111/jcmm.18516

DotPlot(sceList , features =markers,group.by ="Standardclusters")+RotatedAxis() #,group.by ="Cell_type"
#Standardclusters 分群较好？
DotPlot(sceList , features =markers,group.by ="Cell_type")+RotatedAxis() #,group.by ="Cell_type"
DotPlot(sceList , features =markers,group.by ="Cell_type.refined")+RotatedAxis()# +theme_bw()#,group.by ="Cell_type"
##找一篇文献?
markers_0 <-c("PTPRC",#IMMUNE CELL
              "EPCAM",#epithelial
              "PECAM1"  ) #stromal #"MME",无差异
DotPlot(sceList , features =markers_0,group.by ="Cell_type")+RotatedAxis() #,group.by ="Cell_type"
#
DotPlot(sceList ,  
        features =markers_0,group.by ="Cell_type" )+ 
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  scale_color_viridis(7, option = "D",alpha = 0.7)+ #library(viridis) #C D G 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw()+RotatedAxis()
  #theme(axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5))
ggsave("./out/dotplot_marker_3.pdf",family="serif",width = 3.5,height = 3.5)#保存4,3.5 

DotPlot(sceList ,  
        features =markers,group.by ="Cell_type" )+ 
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  scale_color_viridis(7, option = "D",alpha = 0.7)+ #library(viridis) #C D G 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw()+RotatedAxis()
ggsave("./out/dotplot_marker.pdf",family="serif",width = 9,height = 3.5)#保存4,3.5 

#markers_1 #根据上图调整
table(meta$Cell_type)
markers_1 <-c("EPCAM","CDH1","KRT7","KRT19",#Epithelial cells
"PECAM1","VWF",#Endothelial cells
"PDGFRA","LUM","COL1A1",#Fibroblasts
#"MYLK","ACTA2","PDGFRB", #
#"MKI67","STMN1",
"LYZ","CD163","CD14",#Myeloid cells #"CD1C",
#"S100A9","S100A8","CSF3R",#neutrophil
"KIT","GATA2","CPA3","MS4A2",#mast
"MS4A1","CD79A","CD19",#B cells B lymphocytes
#"IGHG1","MZB1","JCHAIN","SDC1",#Plasma
#"IGHD","FCER2","TCL1A", #navie B
"NKG7","KLRD1","KLRB1",#"NCR1",
"CD3D","CD3E","CD4","CD8A"#,"FOXP3","TNFRSF4","IKZF2","PDCD1","TIGIT"
) 
DotPlot(sceList ,  
        features =markers_1,group.by ="Cell_type" )+ 
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  scale_color_viridis(7, option = "D",alpha = 0.7)+ #library(viridis) #C D G 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw()+RotatedAxis()
ggsave("./out/dotplot_marker_1.pdf",family="serif",width = 7,height = 3.5)#保存4,3.5 

##Cell_type umap
#pdf("./out/UMAP_CellDimPlot-1.pdf",width = 8,height = 6,family="serif")
CellDimPlot(
  srt = sceList, group.by = c("Cell_type"), reduction = "umap",
  #title = "beforeStandard", 
  label.size = 20,pt.size = 1,show_stat=F,#add_density=T,
  theme_use = "theme_blank")
ggsave("./out/UMAP_CellDimPlot-1.pdf",width = 5,height = 4,family="serif")
#dev.off()
##score tumor vs. normal ----

df_epi_score<-meta[meta$Cell_type %in% c("Epithelial cells","Malignant cells"),]

library(ggplot2);library(EnvStats);library(ggpubr)
table(all_merge$Efficacy);colnames(all_merge)
ggplot(df_epi_score,
       aes(x = Cell_type, y =Score,fill =Cell_type) )+ #,fill =Efficacy y=Srps
  geom_point(aes(color =Cell_type),alpha=0.5,size=0.2,
             position=position_jitterdodge(jitter.width = 0.45,
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  geom_boxplot(alpha=0.5,width=0.55,
               position=position_dodge(width=0.8),
               size=0.05,outlier.colour = NA)+
  labs(y="Score",x= NULL)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.25)+
  stat_n_text()+ #library(EnvStats) # 显示样本数NR 38  N 31
  scale_color_manual(values = c("#CC79A7", "#56B4E9") )+scale_fill_manual(values = c("#CC79A7", "#56B4E9") )+
  #scale_color_manual(values = cols)+
  #facet_grid(.~Treatment,space = "free",scales = "free")+ #根据treatment分页
  stat_compare_means(#aes(group = Group) ,   #library(ggpubr)
    #comparisons=list(c("Non-responder","Responder")),#my_comparisons,
    #label = "p.format",#"p.format p.signif 1.1e-06
    method = "wilcox.test", #wilcox t
    show.legend= F,#删除图例中的"a"
    label.x=1.2,bracket.size=0.1,vjust=0.1,#label.y=0.35,
    #label.y = max(log2(GDSC2$Camptothecin_1003)),
    hide.ns = T,size=4)+
  #ylim(c(0.1,4))+ #y轴范围0.36
  theme_bw()+#facet_grid(~Treatment)+
  theme(axis.title  = element_text(size=12), axis.text = element_text(size=12),legend.position ="none",  
        panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )#+ 
#scale_x_discrete(labels= c("PD-L1+Chemo","Chemo"))
#ggsave("./out/boxplot_score_epi.pdf",width = 4,height = 4,family="serif")

#UMAP-score----
FeatureDimPlot(sceList, features = "Score",show_stat = F,xlab = "UMAP_1",ylab = "UMAP_2",
               theme_use = "theme_blank")+theme(legend.position = "right")
ggsave(paste0("./out/UMAP_score-",substr(Sys.time(),1,10),".pdf"),height = 5,width = 4)
#scp
FeatureDimPlot(
  srt = sceList, features = c("Score"), #,
  palette ="Reds",raster = T,xlab = "UMAP_1",ylab = "UMAP_2",show_stat=F,
  #compare_features = TRUE, label = TRUE, label_insitu = TRUE,
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave(paste0("./out/UMAP_score-",substr(Sys.time(),1,10),"-1.pdf"),height = 5,width = 4)

#SLOW
# pdf("./out/UMAP_score_TOMM.pdf",width = 5,height = 5,family="serif")
# FeatureDimPlot(
#   srt = sceList, features = c("Score", "TOMM20"), #TOMM
#   #palette ="Reds",
#   raster = T,
#   compare_features = TRUE, label_insitu = TRUE,#label = TRUE, 
#   reduction = "UMAP", theme_use = "theme_blank"
# ) 
# dev.off()


#dotplot score----
DotPlot(sceList , features =c("Score"),group.by ="Cell_type"#, "Cluster",#"Major.celltype"#,
        #cols = c("skyblue", "pink")
)+
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
  scale_color_viridis(8, option = "C",alpha = 0.7)+ #library(viridis) G D C
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
ggsave("./out/dotplot_score.pdf",family="serif",width = 3,height = 3)#保存3,3  4-12

