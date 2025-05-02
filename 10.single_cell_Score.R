#score single cell

rm(list=ls())
setwd("~/R/my_projects/metabolism_sig/20.single_cell_Score")
#load Rdata qs  qsave(sceList,"sceList.qs")
sceList <-qread("./sceList.qs")# #加载打分完成的2.5G

if(T){

  library(qs)
  library(Seurat)
  library(ggplot2)
  library(clustree)
  library(cowplot)
  library(dplyr)
  library(data.table)
  library(tibble)
  library(SCP)
  folder_path <- "./out"
  # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {
    # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
  my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  color_used = c("#B84D64","#864A68","#E32D32","#5E549A","#8952A0","#384B97","#911310",
                 "#7CA878","#35A132","#6B70B0","#20ACBD","#959897",'skyblue',
                 "#F4A2A3","#B6CCD7","#AF98B5","#E01516","#FFFFFF")
}




#0. LOAD sce----
#LOAD RDATA！QS
sceList <-qread("~//R//my_projects//metabolism_sig//1.scRNA/sceList.qs")#

my_model <-"StepCox[forward] + RSF"
res <-readRDS("~//R//my_projects//metabolism_sig//5.1-101建模/out/101_results_178input.rds")

#1.0 AddModuleScore
genelist <- list(res$Sig.genes)
sceList <- AddModuleScore(object = sceList,features = genelist,#gene,
                          name = "score" #ctrl = 100,
) #slow 慢的要死
colnames(sceList@meta.data)[15] <- 'MitoScore'
#Cell_type NA
meta<-sceList@meta.data

table(is.na(sceList@meta.data$Cell_type))# T 1628

sceList@meta.data$Cell_type <-ifelse(is.na(sceList@meta.data$Cell_type)==T,"Undetermined",sceList@meta.data$Cell_type)

#1.0 plot----
##score tumor vs. normal ----

df_epi_score<-meta[meta$Cell_type %in% c("Epithelial cells","Malignant cells"),]

library(ggplot2);library(EnvStats);library(ggpubr)
#table(all_merge$Efficacy);colnames(all_merge)
ggplot(df_epi_score,
       aes(x = Cell_type, y =MitoScore,fill =Cell_type) )+ #,fill =Efficacy y=Srps
  # geom_point(aes(color =Cell_type),alpha=0.5,size=0.2,
  #            position=position_jitterdodge(jitter.width = 0.45,
  #                                          jitter.height = 0,
  #                                          dodge.width = 0.8))+
  geom_boxplot(alpha=0.7,width=0.55,
               position=position_dodge(width=0.8),
               size=0.05,outlier.colour = NA)+
  labs(y="Score",x= NULL)+
  geom_violin(alpha=0.4,width=0.9,
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
  theme_bw(base_size = 14,base_family = "serif")+#facet_grid(~Treatment)+
  theme(axis.title  = element_text(size=14,family = "serif"), axis.text = element_text(size=12),legend.position ="none",  
        panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )#+ 
#scale_x_discrete(labels= c("PD-L1+Chemo","Chemo"))
#ggsave("./out/boxplot_mitoscore_epi.pdf",width = 4,height = 4,family="serif")

#UMAP-score----
FeatureDimPlot(sceList, features = "MitoScore",show_stat = F,xlab = "UMAP_1",ylab = "UMAP_2",
               theme_use = "theme_blank")+theme(legend.position = "right")
#ggsave(paste0("./out/UMAP_MitoScore-",substr(Sys.time(),1,10),".pdf"),height = 5,width = 4)
#scp
FeatureDimPlot(
  srt = sceList, features = c("MitoScore"), #,
  palette ="Reds",raster = T,xlab = "UMAP_1",ylab = "UMAP_2",show_stat=F,
  #compare_features = TRUE, label = TRUE, label_insitu = TRUE,
  reduction = "UMAP", theme_use = "theme_blank"
)
#ggsave(paste0("./out/UMAP_score-",substr(Sys.time(),1,10),"-1.pdf"),height = 5,width = 4)


#dotplot score----

DotPlot(sceList , features =c("MitoScore"),group.by ="Cell_type"#, "Cluster",#"Major.celltype"#,
        #cols = c("skyblue", "pink")
)+
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置
  scale_color_viridis(8, option = "D",alpha = 0.7)+ #library(viridis) G D C
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw(base_family = "serif")+
  theme(#panel.grid = element_blank(),
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5))
#ggsave("./out/dotplot_mitoscore_D.pdf",family="serif",width = 3.2,height = 3.5)#保存3,3  4-12

# Error in `$<-.data.frame`(`*tmp*`, "id", value = NA_character_) : 
#   replacement has 1 row, data has 0  替换NA

impor<-as.data.frame(res[["ml.res"]][[my_model]]$importance)  #[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
colnames(impor)<-"importance"
plot(density(impor$importance))
#impor[impor$importance >mean(impor$importance),] #~100
library(tidyverse)
impor<-arrange(impor,desc(importance) )#impor[order(impor$importance),] #降序

#res$Sig.genes model gene
DotPlot(sceList , features =c("MitoScore",rownames(impor)),
        group.by ="Cell_type" 
)+scale_color_viridis(8, option = "C",alpha = 0.7)+ #library(viridis) G D C
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw(base_family = "serif")+
  theme(#panel.grid = element_blank(),
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5))
#ggsave("./out/dotplot_mitoscore_64_C.pdf",family="serif",width = 12,height =3.8)#保存3,3  4-12

#single gene exp in celltype
for (i in rownames(impor) ){
  print(i)
  DotPlot(sceList , features =c(i),group.by ="Cell_type"#, "Cluster",#"Major.celltype"#,
          #cols = c("skyblue", "pink")
  )+scale_color_viridis(8, option = "D",alpha = 0.7)+ #library(viridis) G D C
    labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
    theme_bw(base_family = "serif")#+
    #theme(axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5))
  ggsave(paste0("./out/dotplot_",i,"_D.pdf"),family="serif",width = 3,height =3)#保存3,3  4-12
  
  ##tumor vs epi boxplot
  df_epi_score<-meta[meta$Cell_type %in% c("Epithelial cells","Malignant cells"),]
  df_epi_score$id <-rownames(df_epi_score)
  #提取基因表达
  data_df <- Seurat::FetchData(sceList,vars = i);colnames(data_df)<-"gene"
  data_df$id <-rownames(data_df)
  df_epi_score<-merge(df_epi_score,data_df,by="id")   #合并
  library(ggplot2);library(EnvStats);library(ggpubr)
  #table(all_merge$Efficacy);colnames(all_merge)
  ggplot(df_epi_score,
         aes(x = Cell_type, y =gene,fill =Cell_type) )+ #,fill =Efficacy y=Srps
    # geom_point(aes(color =Cell_type),alpha=0.5,size=0.2,
    #            position=position_jitterdodge(jitter.width = 0.45,
    #                                          jitter.height = 0,
    #                                          dodge.width = 0.8))+
    geom_boxplot(alpha=0.8,width=0.55,
                 position=position_dodge(width=0.8),
                 size=0.05,outlier.colour = NA)+
    labs(y=i,x= NULL)+
    geom_violin(alpha=0.2,width=0.9,
                position=position_dodge(width=0.8),
                size=0.25)+
    stat_n_text()+ #library(EnvStats) # 显示样本数NR 38  N 31
    scale_color_manual(values = c("#CC79A7", "#56B4E9"))+ # c("#EDA065","#66CCFF","#7EC7A7")
    scale_fill_manual(values = c("#CC79A7", "#56B4E9"))+
    #scale_color_manual(values = cols)+
    #facet_grid(.~Treatment,space = "free",scales = "free")+ #根据treatment分页
    stat_compare_means(#aes(group = Group) ,   #library(ggpubr)
      #comparisons=list(c("Non-responder","Responder")),#my_comparisons,
      label = "p.signif",#"p.format p.signif 1.1e-06
      method = "wilcox.test", #wilcox t
      show.legend= F,#删除图例中的"a"
      label.x=1.5,bracket.size=0.1,vjust=0.1,#label.y=0.35,
      #label.y = max(log2(GDSC2$Camptothecin_1003)),
      hide.ns = T,size=5)+
    #ylim(c(0.1,4))+ #y轴范围0.36
    theme_bw(base_size = 16,base_family = "serif")+#facet_grid(~Treatment)+
    theme(axis.title  = element_text(size=14,family = "serif"), axis.text = element_text(size=14),legend.position ="none",  
          panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )#+ 
  #scale_x_discrete(labels= c("PD-L1+Chemo","Chemo"))
  ggsave(paste0("./out/dotplot_boxplit_TvsN_",i,"_D.pdf"),family="serif",width = 3.5,height =3.5)#保存3,3  4-12
   
}

#2.0 monocle2+cytotrace----

##2.1 DEG----
#SCE Differential expression analysis
#https://github.com/zhanghao-njmu/SCP?tab=readme-ov-file#differential-expression-analysis
#limit worker!
##ALL DEG
library(BiocParallel)
register(MulticoreParam(workers = 24, progressbar = TRUE))
sceList <- RunDEtest(srt = sceList, group_by = "Cell_type", fc.threshold = 1, only.pos = FALSE)
#slow 4h   +/- tumor DEG enrichment
VolcanoPlot(srt =sceList, group_by = "Cell_type",pt.size = 0.1)
#save  qsave(sceList,"sceList.qs") 12G

DEGs <- sceList@tools$DEtest_Cell_type$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
## Annotate features with transcription factors and surface proteins-enrichment
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

#sceList <- AnnotateFeatures(sceList, species = "Homo_sapiens", db = c("TF", "CSPA"))
ht <- FeatureHeatmap(
  srt = sceList, group.by = "Cell_type", features = DEGs$gene, feature_split = DEGs$group1,
  topTerm = 5,species = "Homo_sapiens", db = c("GO_BP"), anno_terms = TRUE, #, "KEGG"
  heatmap_palette = "PiYG", #RColorBrewer::display.brewer.all()
  #feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  cell_split_palette = "npg",feature_split_palette = "npg", #show_palettes()
  height = 5, width = 4
)

pdf("./out/DEG_all.pdf",width = 14,height = 8,family="serif") #save
print(ht$plot)
dev.off()

##Enrichment analysis(over-representation)
sceList <- RunEnrichment(
  srt = sceList, group_by = "Cell_type", db = "GO_BP", #KEGG
  species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
)

table(sceList$Cell_type)
pdf("./out/DEG_GO_BP.pdf",width = 8,height = 4,family="serif") 
EnrichmentPlot(
  srt = sceList, group_by = "Cell_type", group_use = c("Epithelial cells", "Malignant cells"),
  plot_type = "comparison"  
)
dev.off()

pdf("./out/DEG_GO_BP_lollipop.pdf",width = 12,height = 8,family="serif") 
EnrichmentPlot(
  srt = sceList, group_by = "Cell_type", group_use = c("Epithelial cells", "Malignant cells"),
  topTerm = 5,
  plot_type = "lollipop"  #, ncol = 1
)&theme(panel.background = element_rect(fill=NA) )
dev.off()

##GSEA-KEGG-Malignant cells
sceList <- RunGSEA(
  srt = sceList, group_by = "Cell_type", db = "KEGG",#"KEGG",  GO_BP
  species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05"
)
#SHOW
GSEA_KEGG<-sceList@tools[["GSEA_Cell_type_wilcox"]][["results"]][["Malignant cells-KEGG"]]@result
GSEA_KEGG<-GSEA_KEGG[GSEA_KEGG$p.adjust <0.01,]
write.csv(GSEA_KEGG,"Malignant cells-KEGG.csv")  #SAVE

pdf("./out/GSEA_KEGG-Malignant cells_Oxidative phosphorylation.pdf",width = 4,height = 4,family="serif") 
GSEAPlot(srt = sceList, db = "KEGG", palette = "Paired",only_sig = TRUE,
         group_by = "Cell_type", group_use = "Malignant cells", subplots = 1:2,
         #topTerm = 3,
         n_coregene =0,base_size = 14,#plot_type = c("line"), #, "comparison"
         geneSetID = c("hsa00190")
         )
dev.off()

##2.2 Trajectory inference----
library(BiocParallel)
register(MulticoreParam(workers = 24, progressbar = TRUE))
sceList <- RunSlingshot(srt = sceList, group.by = "Cell_type", reduction = "UMAP")#slow
#FeatureDimPlot(sceList, features = paste0("Lineage", 1:3), reduction = "UMAP", theme_use = "theme_blank")
#SLOW-all select epi !
CellDimPlot(sceList, group.by = "Cell_type", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1)
FeatureDimPlot(sceList, features = paste0("Lineage", 1:3), reduction = "UMAP", theme_use = "theme_blank")#qsave(sceList,"sceList.qs")

#Dynamic features
sceList <- RunDynamicFeatures(srt = sceList, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
ht <- DynamicHeatmap(
  srt = sceList, lineages = c("Lineage1", "Lineage2"),#c("Lineage1", "Lineage2"),
  use_fitted = TRUE, n_split = 6, reverse_ht = "Lineage1",
  species = "Homo_sapiens", db = "GO_BP", anno_terms = TRUE, anno_keys = F, anno_features = F,
  heatmap_palette = "viridis", cell_annotation = "Cell_type",
  #separate_annotation = list("Cell_type", c("Sample_Origin", "Sample")), separate_annotation_palette = c("Paired", "Set1"),
  #feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
) #slow
print(ht$plot)#失败

#genes 拟时序点图 slow
#for (i in genelist[[1]] ){
DynamicPlot(
  srt = sceList, lineages = c("Lineage1", "Lineage2"), group.by = "Cell_type",
  features =i,#markers_1[1:5],# c("Plk1", "Hes1", "Neurod2", "Ghrl", "Gcg", "Ins2"),
  compare_lineages = TRUE, compare_features = FALSE)
#}

#基因亚群表达小提琴图-exp in subcelltype violin
for (i in genelist[[1]] ){
  print(i)
  FeatureStatPlot(
    srt = sceList, group.by = "Cell_type", bg.by = "Cell_type",
    stat.by = i,#c("Sox9", "Neurod2", "Isl1", "Rbp4"), 
    add_box = TRUE,legend.position = "",plot_type = c("violin", "box", "bar", "dot", "col"),
    comparisons = list(
      c("Epithelial cells", "Malignant cells") )
  )
  ggsave(paste0("./out/violin_Cell_type_",i,".pdf"),width = 5,height = 3,family="serif")
  

}


##上皮 恶性细胞 单独拟时序
sce<-sceList[,sceList$Cell_type %in% c("Epithelial cells", "Malignant cells")]
sce<- RunSlingshot(srt = sce, group.by = "Cell_type", reduction = "UMAP")#slow  不符合趋势
#FeatureDimPlot(sceList, features = paste0("Lineage", 1:3), reduction = "UMAP", theme_use = "theme_blank")
#SLOW-all select epi !
CellDimPlot(sce, group.by = "Cell_type", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1)
FeatureDimPlot(sce, features = paste0("Lineage", 1:3), reduction = "UMAP", theme_use = "theme_blank")#qsave(sceList,"sceList.qs")



sce <- RunDynamicFeatures(srt = sce, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)#20min
ht <- DynamicHeatmap(
  srt = sce, lineages = c("Lineage1", "Lineage2"),#c("Lineage1", "Lineage2"),
  use_fitted = TRUE, n_split = 6, reverse_ht = "Lineage1",
  species = "Homo_sapiens", db = "GO_BP", anno_terms = TRUE, anno_keys = F, anno_features = F,
  heatmap_palette = "viridis", cell_annotation = "Cell_type",
  #separate_annotation = list(c("Sample_Origin", "Sample")), separate_annotation_palette = c("Paired", "Set1"),
  #feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
) #slow
print(ht$plot)
#slow 1min
DynamicPlot(
  srt = sce, lineages = c("Lineage1", "Lineage2"), group.by = "Cell_type",
  features =i,#markers_1[1:5],# c("Plk1", "Hes1", "Neurod2", "Ghrl", "Gcg", "Ins2"),
  compare_lineages = TRUE, compare_features = FALSE) #丑

#plot(sceList@meta.data[["MitoScore"]])  +- score -->DEG-->hallmark,cell chat
