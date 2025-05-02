
if(T){
  rm(list = ls())
  setwd("~/R/my_projects/metabolism_sig/14.免疫检查点")
  #setwd("D:\\科研\\肺癌乳酸代谢\\14.免疫检查点")
  folder_path <- "./out" # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {   # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
  options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
  library(ComplexHeatmap) # 用于绘制热图
  library(circlize) # 用于热图颜色设置
  library(ChAMPdata) # 用于提供甲基化注释文件
  library(data.table) # 用于读取大文件
  library(tibble)
  library(limma)
  data("probe.features")
  Sys.setenv(LANGUAGE = "en") #显示英文报错信息
  options(stringsAsFactors = FALSE) #禁止chr转成factor
  # 设置热图颜色
  heatmap.BlWtRd <- c("#1F66AC", "#75AFD3", "grey90", "#FAB99B", "#B2192B")
  
  # 设置感兴趣基因集
  immunomodulator <- read.table("immunomodulator.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  
  # 数据处理 #
  # 该部分是为产生最终用于绘图的文件
  ## 表达谱
  expr <- fread("convert_exp.txt",data.table = F,header = T,check.names = F)
  expr <- avereps(expr[,-1],expr$Tag)
  table(substr(colnames(expr),14,15))
  expr <- expr [,as.numeric(substr(colnames(expr ),14,15))<10]
  colnames(expr) <- substr(colnames(expr) ,1,12)
  expr <- t(scale(t(expr)))
  #is.element(rownames(immunomodulator),rownames(expr)) # 所有基因都在表达谱内
  
}

if(T){
  #分型的结果文件-训练集riskScore和risk分组+生存-表达矩阵
  #sinfo <- read.table("TCGA.csv",sep = ",",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  my_model <-"StepCox[forward] + RSF"
  res <-readRDS("~//R//my_projects//metabolism_sig//5.1-101建模/out/101_results_178input.rds")
  riskscore<-res[["riskscore"]][[my_model]]$TCGA
  colnames(riskscore)[4]<-"riskscore"
  riskscore$risk <- ifelse(riskscore$riskscore>median(riskscore$riskscore) ,"High","Low")#分组
  rownames(riskscore)<-riskscore$ID #exp is in res
  rm(res)
  
  sinfo=riskscore[,"risk",drop=F]
  colnames(sinfo)="subtype"
  jj <- intersect(rownames(immunomodulator),rownames(expr))
  expr <- expr[jj,]
}

if(T){
  ## 甲基化谱
  meth <- readRDS("甲基化肺腺癌.rds")
  meth <- na.omit(meth) #library(ChAMPdata) #500MB 手动安装
  probeOfInterest <- probe.features[which(probe.features$gene %in% rownames(immunomodulator)),]
  probeOfInterest <- probeOfInterest[intersect(rownames(probeOfInterest), rownames(meth)),]
  #is.element(rownames(immunomodulator), probeOfInterest$gene)
  meth <- meth[rownames(probeOfInterest),]
  meth$gene <- probeOfInterest$gene
  #获取每个基因的平均甲基化值
  meth <- as.data.frame(apply(meth[,setdiff(colnames(meth), "gene")], 2, function(x) tapply(x, INDEX=factor(meth$gene), FUN=median, na.rm=TRUE)))
  
  ## 拷贝数变异 (-2,-1,0,1,2: 2 copy del, 1 copy del, no change, amplification, high-amplification)
  cna <- fread("TCGA.PANCAN.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz",sep = "\t",data.table = F,header = T)
  cna <- column_to_rownames(cna,var = "Sample")
  cna$gene <- sapply(strsplit(rownames(cna),"|",fixed = T),"[",1)
  cna <- cna[!duplicated(cna$gene),]; cna <- cna[,setdiff(colnames(cna),"gene")]
  is.element(rownames(immunomodulator),rownames(cna)) # 有些基因没有对应的拷贝数结果
  cna <- cna[intersect(rownames(cna),rownames(immunomodulator)),]
  cna[cna > 1] <- 1 # 统一扩增
  cna[cna < -1] <- -1 # 统一缺失
  cna <- cna [,as.numeric(substr(colnames(cna ),14,15))<10]
  colnames(cna) <- substr(colnames(cna) ,1,12)
  
  
  ## 提取共同样本
  comsam <- intersect(colnames(expr), colnames(meth))
  comsam <- intersect(comsam, colnames(cna))
  sinfo <- sinfo[comsam,,drop = F]
  expr <- expr[,comsam]
  meth <- meth[,comsam]
  cna <- cna[,comsam]
}

write.table(sinfo[,"subtype",drop = F], file = "easy_input_subtype.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(expr, file = "easy_input_expr.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(meth, file = "easy_input_meth.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(cna, file = "easy_input_cna.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

if(T){
  # 读取数据 #
  sinfo <- read.table(file = "easy_input_subtype.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  expr <- read.table(file = "easy_input_expr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  meth <- read.table(file = "easy_input_meth.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  cna <- read.table(file = "easy_input_cna.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  
  # 亚型数目（注意，亚型必须以C1，C2，C3...等命名）
  (n.subt <- length(unique(sinfo$subtype))) # 获取亚型数目
  subt <- unique(sinfo$subtype) # 获取亚型名
  
  # 初始化绘图矩阵
  expMat <- as.data.frame(t(expr[rownames(immunomodulator),]))
  expMat$subtype <- sinfo[rownames(expMat), "subtype"]
  expMat <- as.data.frame(t(apply(expMat[,setdiff(colnames(expMat), "subtype")], 2, 
                                  function(x) 
                                    tapply(x, 
                                           INDEX = factor(expMat$subtype), 
                                           FUN = median, 
                                           na.rm = TRUE)))) # 对同一亚型内的样本取中位数
  corExpMeth <- ampFreq <- delFreq <- 
    as.data.frame(matrix(NA,
                         nrow = nrow(immunomodulator),
                         ncol = n.subt, 
                         dimnames = list(rownames(immunomodulator), 
                                         unique(sinfo$subtype))))
  
  
}

## 表达谱与甲基化的相关性
for (i in rownames(immunomodulator)) {
  if(!is.element(i, rownames(expr)) | !is.element(i, rownames(meth))) { # 如果存在任意一方有缺失的基因
    corExpMeth[i,] <- NA # 则保持矩阵为NA
  } else { # 否则取出亚型样本，做表达和甲基化的相关性
    for (j in subt) {
      sam <- rownames(sinfo[which(sinfo$subtype == j),,drop = F])
      expr.subset <- as.numeric(expr[i, sam])
      meth.subset <- as.numeric(meth[i, sam])
      ct <- cor.test(expr.subset, meth.subset, method = "spearman") # 这里采用speaman相关性
      corExpMeth[i, j] <- ct$estimate
    }
  }
}
## 扩增/缺失频率
for (i in rownames(immunomodulator)) {
  if(!is.element(i, rownames(cna))) { # 同理，如果存在拷贝数中缺失某基因，则保持NA
    ampFreq[i,] <- NA 
    delFreq[i,] <- NA
  } else { # 否则
    # 计算i在总样本中的频率
    ampFreqInAll <- sum(as.numeric(cna[i,]) == 1)/ncol(cna) # 总样本中扩增的数目除以总样本数
    delFreqInAll <- sum(as.numeric(cna[i,]) == -1)/ncol(cna) # 总样本中缺失的数目除以总样本数
    for (j in subt) {
      # 计算i在亚型j中的频率
      sam <- rownames(sinfo[which(sinfo$subtype == j),,drop = F])
      cna.subset <- cna[, sam]
      ampFreqInSubt <- sum(as.numeric(cna.subset[i,]) == 1)/length(sam) # 该亚型中扩增的数目除以该亚型样本数
      delFreqInSubt <- sum(as.numeric(cna.subset[i,]) == -1)/length(sam) # 该亚型中缺失的数目除以该亚型样本数
      
      ampFreqInDiff <- ampFreqInSubt - ampFreqInAll # 根据原本，用亚型特异性扩增比例减去总扩增比例
      delFreqInDiff <- delFreqInSubt - delFreqInAll # 同理
      
      ampFreq[i, j] <- ampFreqInDiff
      delFreq[i, j] <- delFreqInDiff
    }
  }
}
write.table(expMat,"expMat.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(corExpMeth,"corExpMeth.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(ampFreq,"ampFreq.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(delFreq,"delFreq.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

##plot heatmap----

if(T){
  # 创建列注释
  annCol <- data.frame(subtype = subt,
                       row.names = subt)
  annCol <- annCol[order(annCol$subtype),,drop = F] # 按照亚型排序
  annColors <- list()
  annColors[["subtype"]] <- c("High" = "pink",
                              "Low" = "skyblue")
  top_anno <- HeatmapAnnotation(df                   = annCol,
                                col                  = annColors,
                                gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                                simple_anno_size     = unit(3.5, "mm"), # 注释高3.5毫米
                                show_legend          = F, # 不显示亚型的图例，因为一目了然
                                show_annotation_name = F, # 不显示该注释的名称
                                border               = FALSE) # 不显示注释的外边框
  
  # 创建行注释
  annRow <- immunomodulator
  annRow[which(annRow$Category == "Co-stimulator"),"Category"] <- "Co-stm" # 这里字符少一些，不会挤在一起，可以后期AI
  annRow[which(annRow$Category == "Co-inhibitor"),"Category"] <- "Co-ihb"
  annRow[which(annRow$Category == "Cell adhesion"),"Category"] <- "Cell\nadhesion" # 这里换行，不会挤在一起，可以后期AI
  annRow[which(annRow$Category == "Antigen presentation"),"Category"] <- "Antigen\npresentation"
  annRow$Category <- factor(annRow$Category, levels = c("Co-stm","Co-ihb","Ligand","Receptor","Cell\nadhesion","Antigen\npresentation","Other")) # 由于行需要按照类分割，所以需要定义因子顺序，否则按照字母表
  annRow$ICI <- factor(annRow$ICI, levels = c("Inhibitory","N/A","Stimulatory"))
  annRowColors <- list("ICI" = c("Inhibitory" = "#A77852","N/A" = "#888888","Stimulatory" = "#035F37"))
  left_anno <- HeatmapAnnotation(df                   = annRow[,"ICI",drop = F],
                                 which                = "row", # 这里是行注释（默认为列）
                                 gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                                 col                  = annRowColors,
                                 simple_anno_size     = unit(3.5, "mm"), # 注释宽3.5毫米
                                 show_annotation_name = F,
                                 border               = F)
  
  ## 绘制表达谱热图（参数下同）
  library(RColorBrewer) #RColorBrewer::display.brewer.all()
  cols=rev(brewer.pal(n=11,name="Spectral"))
  cols1=brewer.pal(n=9,name="Reds")
  cols2=brewer.pal(n=9,name="Purples")
  cols3=brewer.pal(n=9,name="Blues")
  col_expr <- colorRamp2(seq(min(na.omit(expMat)), max(na.omit(expMat)), length = 11), cols) # 创建热图颜色（将热图输入矩阵的最大最小值取5个点，分配颜色红蓝色板；注意矩阵中可能存在的NA值）
  hm.expr <- Heatmap(matrix             = as.matrix(expMat),
                     col                = col_expr,
                     border             = NA, # 无热图外边框
                     rect_gp = gpar(col = "grey80"), # 热图单元格边框为灰色
                     cluster_rows       = F, # 行不聚类
                     cluster_columns    = F, # 列不聚类
                     show_row_names     = T, # 显示行名
                     row_names_side     = "left", # 行名显示在左侧
                     row_names_gp =gpar(fontsize = 10, # 列标题大小
                                        #fontface = "bold.italic", # 列标题字体
                                        fill = "white", # 列标题背景色
                                        col = "black", # 列标题颜色
                                        border = "white"# 列标题边框色
                     ),
                     show_column_names  = F, # 不显示列名（可后期在颜色内AI使得亚型一目了然）
                     column_names_side  = "top", # 列名显示在顶部
                     row_split          = annRow$Category, # 行按照Category进行分割（因子顺序）
                     top_annotation     = top_anno, # 热图顶部注释
                     left_annotation    = left_anno, # 热图左侧注释
                     name               = "mRNA\nExpression", # 热图颜色图例的名称
                     width              = ncol(expMat) * unit(4, "mm"), # 热图单元格宽度（稍大于高度，因为所有注释都放在底部，平衡图形纵横比）
                     height             = nrow(expMat) * unit(3.5, "mm")) # 热图单元格高度
  
  col_corExprMeth <- colorRamp2(seq(min(na.omit(corExpMeth)), max(na.omit(corExpMeth)), length = 9), cols1)
  hm.corExprMeth <- Heatmap(matrix             = as.matrix(corExpMeth),
                            col                = col_corExprMeth,
                            border             = NA,
                            rect_gp = gpar(col = "grey80"),
                            cluster_rows       = F,
                            cluster_columns    = F,
                            show_row_names     = F,
                            row_names_side     = "left",
                            row_names_gp       = gpar(fontsize = 10),
                            show_column_names  = F,
                            column_names_side  = "top",
                            row_split          = annRow$Category,
                            row_title          = NULL,
                            top_annotation     = top_anno,
                            name               = "Expression\nvs. Methylation",
                            width              = ncol(expMat) * unit(4, "mm"),
                            height             = nrow(expMat) * unit(3.5, "mm"))
  col_ampFreq <- colorRamp2(seq(min(na.omit(ampFreq)), max(na.omit(ampFreq)), length = 9), cols2)
  hm.ampFreq <- Heatmap(matrix             = as.matrix(ampFreq),
                        col                = col_ampFreq,
                        border             = NA,
                        rect_gp = gpar(col = "grey80"),
                        cluster_rows       = F,
                        cluster_columns    = F,
                        show_row_names     = F,
                        row_names_side     = "left",
                        row_names_gp       = gpar(fontsize = 10),
                        show_column_names  = F,
                        column_names_side  = "top",
                        row_split          = annRow$Category,
                        row_title          = NULL,
                        top_annotation     = top_anno,
                        name               = "Amplification\nFrequency",
                        width              = ncol(expMat) * unit(4, "mm"),
                        height             = nrow(expMat) * unit(3.5, "mm"))
  col_delFreq <- colorRamp2(seq(min(na.omit(delFreq)), max(na.omit(delFreq)), length = 9), cols3)
  hm.delFreq <- Heatmap(matrix             = as.matrix(delFreq),
                        col                = col_delFreq,
                        border             = NA,
                        rect_gp = gpar(col = "grey70"),
                        cluster_rows       = F,
                        cluster_columns    = F,
                        show_row_names     = F,
                        row_names_side     = "left",
                        row_names_gp =gpar(fontsize = 10),
                        show_column_names  = F,
                        column_names_side  = "top",
                        row_split          = annRow$Category,
                        row_title          = NULL,
                        top_annotation     = top_anno,
                        name               = "Deletion\nFrequency",
                        width              = ncol(expMat) * unit(4, "mm"),
                        height             = nrow(expMat) * unit(3.5, "mm"))
  
}

# draw(hm.expr + hm.corExprMeth + hm.ampFreq + hm.delFreq, # 水平衔接各个子热图
#      heatmap_legend_side = "bottom") # 热图颜色图例显示在下方
pdf(file = "./out/pheatmap_immunomodulator-2.pdf", width = 6,height = 12)
draw(hm.expr + hm.corExprMeth + hm.ampFreq + hm.delFreq, # 水平衔接各个子热图
     heatmap_legend_side = "bottom") # 热图颜色图例显示在下方
invisible(dev.off())
#sessionInfo()
#表达-甲基化相关性Expression\nvs. Methylation
