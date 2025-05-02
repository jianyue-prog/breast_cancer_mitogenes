
rm(list = ls())
setwd("~//R//github//Mime1")
#https://mp.weixin.qq.com/s/2obBOi8Nh4PeCGr1jR6FMQ
#https://github.com/l-magnificence/Mime/tree/main

#load("~/R/my_projects/metabolism_sig/5.1-101建模/Mime1-101.RData")

#0. install----

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 
# depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' , 'ComplexHeatmap' )
# 
# for(i in 1:length(depens)){
#   depen<-depens[i]
#   if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
# }
# 
# if (!requireNamespace("CoxBoost", quietly = TRUE))
#   devtools::install_github("binderh/CoxBoost")
# 
# if (!requireNamespace("fastAdaboost", quietly = TRUE))
#   devtools::install_github("souravc83/fastAdaboost")
# 
# if (!requireNamespace("Mime", quietly = TRUE))
#   devtools::install_github("l-magnificence/Mime")#

#1.0 demo----
library(Mime1)
# load("~//R//github//Mime1/data/Example.cohort.Rdata")
# list_train_vali_Data[["Dataset1"]][1:5,1:5]
# 
# load("./data/genelist.Rdata")
# genelist
#load genelist
genelist <- readxl::read_xls("~/R/my_projects/metabolism_sig/1.scRNA/Human.MitoCarta3.0.xls",sheet = 2)
genelist <-unique(genelist$Symbol)#1136
df<-read.csv("~/R/my_projects/metabolism_sig/2.enrichment/out/差异基因.csv",row.names = 1)
#DEG

#load LUAD----
if(T){
  id <- list.files(path = "~//R//my_projects//metabolism_sig//5.确定基因//data",pattern = "csv")
  id
  cohort <- list()
  for (i in 1:length(id)) {cohort[[i]] <- fread(paste0("~//R//my_projects//metabolism_sig//5.确定基因//data//",id[i]),#id[i],
                                                header = T,data.table = F)
  cohort[[i]]  <-na.omit( column_to_rownames(cohort[[i]] ,var = colnames(cohort[[i]] )[1]))
  #print(max(cohort[[i]][,3:ncol(cohort[[i]])]))
  library(dplyr)
  cohort[[i]] <- cohort[[i]] %>%
    rownames_to_column("ID")
  colnames(cohort[[i]] )[2:3] <- c("OS.time","OS")#c("futime","fustat")
  #删除生存=0
  cohort[[i]] <- cohort[[i]][cohort[[i]]$OS.time != 0,] #futime
  }
  names(cohort) <- as.character(lapply(strsplit(id,"_"),function(i){i <- i[1]}))
  cohort <- cohort[7:1]
  names(cohort) <- toupper(names(cohort) )
} #11830 12114 gene
#names(cohort)[1] <- "train_set" #TCGA

#提取差异+线粒体基因
gene <-intersect(genelist,df[df$Group !="N.S.",]$id )#244
gene <-intersect(gene,colnames(cohort$TCGA) ) #与表达矩阵取交集 176

#重新提取表达

for (y in 1:length(cohort)) {
  cohort[[y]] <- cohort[[y]][,c("ID","OS.time","OS",gene)]
# intersect(prol$id,colnames(cohort[[y]]))
# cohort[[y]] <- cohort[[y]][cohort[[y]]$futime>0,]
} #179-3=176


#创建预后模型并选择最优模型----
#ctrl+鼠标拉小显示框
#TCGA训练集，单因素COX候选基因，后续建模
res <- Mime1::ML.Dev.Prog.Sig(train_data =cohort$TCGA,# list_train_vali_Data$Dataset1,
                              list_train_vali_Data = cohort,#list_train_vali_Data,
                              unicox.filter.for.candi = T,
                              unicox_p_cutoff = 0.05,
                              candidate_genes = genelist,
                              mode = 'all', ## 'all', 'single', and 'double'
                              nodesize =5,
                              seed = 123 
)#slow >4H;~30min
res[["Sig.genes"]] #64
#保存res
setwd("~//R//my_projects//metabolism_sig//5.1-101建模");getwd()
saveRDS(res,"./out/101_results_178input.rds")#64 unicox filtered candidate genes
#TCGA 抽样训练集-测试集？
##报错？
#Error in response.coxnet(y) :    negative event times encountered; not permitted for Cox family
#删除生存=0 成功

my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 

#可视化C-index。
pdf("./out/cindex_dis_all-1.pdf",width = 8,height = 11)
cindex_dis_all(res,
               validate_set = names(cohort)[-1],#names(list_train_vali_Data)[-1],
               order = names(cohort),#names(list_train_vali_Data),
               width = 0.15,dataset_col=my36colors,color = c("skyblue","grey","red", '#CCE0F5', "#CCC9E6") #颜色条
)
#ggsave("./out/cindex_dis_all.pdf",width = 8,height = 11) #,family="serif" #手动填上α
dev.off()

#展示指定模型的C-index。🙊
cindex_dis_select(res,
                  model="StepCox[forward] + RSF",
                  order= names(cohort))

#根据特定模型计算的风险评分，绘制不同数据集中患者的生存曲线。💪
#保存C-INDEX
my_model <-"StepCox[forward] + RSF"

survplot <- vector("list",length(cohort) ) #2
for (i in c(1:length(cohort))) { #2
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + RSF",
                              dataset = names(cohort)[i], #list_train_vali_Data
                              color=c("skyblue","pink"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="years",pval.coord=c(1,0.2))
  )
  ggsave(paste0("./out/km_",names(cohort)[i],".pdf"),family="serif",width = 4.5,height = 4.5)
  
}
aplot::plot_list(gglist=survplot,ncol=3)
#保存KM ggsave("./out/km_all.pdf",family="serif",width = 9,height = 9)


#计算模型的AUC-失败----
list_train_vali_Data <-cohort
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],#list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 0, #1 3 5报错
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],#list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],#list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")
#展示AUC
pdf("./out/auc_1year.pdf",width = 8,height = 15,family="serif") #,family="serif"
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
dev.off() #失败

#绘制不同数据集中特定模型的ROC曲线。😘
roc_vis(all.auc.1y,
        model_name = my_model,#"StepCox[forward] + plsRcox",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.55,0.15),
        year=1)
#绘制不同数据集中特定模型的1、3、5年的AUC。
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name="StepCox[forward] + plsRcox",
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))

#单变量cox回归特定模型的meta分析
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,
                                   optimal.model = "StepCox[forward] + plsRcox",
                                   type ='categorical' # 'categorical' or 'continuous'
)

metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)

meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))

#2.0 与已知模型相比较----
rs.glioma.lgg.gbm <- cal_RS_pre.prog.sig(use_your_own_collected_sig = F,
                                         type.sig = c('LGG','GBM','Glioma'), #
                                         list_input_data = list_train_vali_Data)
HR_com(rs.glioma.lgg.gbm,
       res,
       model_name=my_model,#"StepCox[forward] + plsRcox",
       dataset=names(list_train_vali_Data),
       type = "categorical")
#计算C-index
cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,
                                             type.sig = c('Glioma','LGG','GBM'),
                                             list_input_data = list_train_vali_Data)
#可视化比较C-index
cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name=my_model,#"StepCox[forward] + plsRcox",
            dataset=names(list_train_vali_Data))
#比较其他模型-c-index
ggsave("./out/cindex_comp.pdf",width = 2.5*length(list_train_vali_Data),height = 12,family="serif") #,family="serif"


# 计算AUC-失败
auc.glioma.lgg.gbm.1 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = F,
                                             type.sig = c('Glioma','LGG','GBM'),
                                             list_input_data = list_train_vali_Data,AUC_time = 0,
                                             auc_cal_method = 'KM')
auc_comp(auc.glioma.lgg.gbm.1,
         all.auc.1y,
         model_name=my_model,#"StepCox[forward] + plsRcox",
         dataset=names(list_train_vali_Data))

#3.0 Immune infiltration analysis----
#https://github.com/l-magnificence/Mime/tree/main
#remotes::install_github("icbi-lab/immunedeconv")

TME_deconvolution_all2 <- function(inputmatrix.list, deconvolution_method = c("xcell", 
                                                                              "epic", "abis", "estimate", "cibersort"), #, "cibersort_abs"
                                   microarray_names = "none", indications = NULL, tumor = TRUE, 
                                   column = "gene_symbol", rmgenes = NULL, scale_mrna = TRUE, 
                                   expected_cell_types = NULL, ...) 
{
  # 加载必要的包
  if (TRUE) {
    library(immunedeconv)
    library(dplyr)
    library(magrittr)
    library(data.table)
    library(readr)
    Sys.setenv(LANGUAGE = "en")
    options(stringsAsFactors = FALSE)
  }
  message("--- Data preprocessing ---")
  
  # 修改后的 annotate_cell_type 函数
  annotate_cell_type <- function(result_table, method) {
    cell_type_map %>% 
      filter(method_dataset == !!method) %>% 
      inner_join(result_table, by = "method_cell_type")
  }
  
  set_cibersort_binary(system.file("extdata", "CIBERSORT.R", package = "Mime1"))
  set_cibersort_mat(system.file("extdata", "LM22.txt", package = "Mime1"))
  
  TME_deconvolution <- function(gene_expression, deconvolution_method = c("xcell", 
                                                                          "epic", "abis", "estimate", "cibersort", "cibersort_abs"), 
                                tumor = TRUE, arrays = FALSE, column = "gene_symbol", 
                                rmgenes = NULL, scale_mrna = TRUE, expected_cell_types = NULL, 
                                ...) {
    deconvolution_methods <- c("quantiseq", "xcell", "epic", 
                               "abis", "mcp_counter", "estimate", "cibersort", "cibersort_abs", 
                               "timer", "consensus_tme")
    if (all(deconvolution_method %in% deconvolution_methods)) {
      if (is(gene_expression, "ExpressionSet")) {
        gene_expression <- gene_expression %>% eset_to_matrix(column)
      }
      if (!is.null(rmgenes)) {
        gene_expression <- gene_expression[!rownames(gene_expression) %in% rmgenes, ]
      }
      tme_combine <- list()
      for (method in deconvolution_method) {
        message(paste0("\n", ">>> Running ", method))
        res <- switch(method, 
                      xcell = immunedeconv::deconvolute_xcell(gene_expression, arrays = arrays, expected_cell_types = expected_cell_types, ...), 
                      epic = immunedeconv::deconvolute_epic(gene_expression, tumor = tumor, scale_mrna = scale_mrna, ...), 
                      cibersort = immunedeconv::deconvolute_cibersort(gene_expression, absolute = FALSE, arrays = arrays, ...), 
                      cibersort_abs = immunedeconv::deconvolute_cibersort(gene_expression, absolute = TRUE, arrays = arrays, ...), 
                      abis = immunedeconv::deconvolute_abis(gene_expression, arrays = arrays), 
                      estimate = immunedeconv::deconvolute_estimate(gene_expression), 
                      timer = immunedeconv::deconvolute(gene_expression, "timer", indications = indications), 
                      consensus_tme = immunedeconv::deconvolute(gene_expression, "consensus_tme", indications = indications), 
                      quantiseq = immunedeconv::deconvolute(gene_expression, "quantiseq"), 
                      mcp_counter = immunedeconv::deconvolute(gene_expression, "mcp_counter"))
        res <- res %>% as_tibble(rownames = "method_cell_type") %>% 
          annotate_cell_type(method = method)
        tme_combine[[method]] <- res
      }
      resultList <- list(tme_combine = tme_combine)
      return(resultList)
    } else {
      print("Please provide the correct parameters for deconvolution method")
    }
  }
  
  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    })
    return(x)
  })
  
  tme_decon_list <- list()
  selected_columns <- microarray_names
  
  for (i in 1:length(inputmatrix.list)) {
    train_data <- inputmatrix.list[[i]]
    # 确保列存在，然后设置行名和选择列
    if(all(c("ID", "OS.time", "OS") %in% colnames(train_data))) {
      test.matrix <- train_data %>% 
        magrittr::set_rownames(train_data$ID) %>% 
        dplyr::select(-ID, -OS.time, -OS) %>% 
        t()
    } else {
      stop("Columns 'ID', 'OS.time', and 'OS' must be present in the data.")
    }
    
    tryCatch({
      resultList <- NULL
      if (selected_columns == "none") {
        resultList <- TME_deconvolution(gene_expression = test.matrix, 
                                        deconvolution_method = deconvolution_method, 
                                        arrays = FALSE)
      } else if (all(selected_columns %in% names(inputmatrix.list))) {
        if (names(inputmatrix.list)[i] %in% selected_columns) {
          resultList <- TME_deconvolution(gene_expression = test.matrix, 
                                          deconvolution_method = deconvolution_method, 
                                          arrays = TRUE)
        } else {
          resultList <- TME_deconvolution(gene_expression = test.matrix, 
                                          deconvolution_method = deconvolution_method, 
                                          arrays = FALSE)
        }
      } else {
        cat("Invalid dataset name(s). Please try again.")
      }
      print("Success")
    }, error = function(e) {
      print(paste("An error occurred:", conditionMessage(e)))
      resultList <- NULL
    }, finally = {
      tme_decon_list[[names(inputmatrix.list)[i]]] <- resultList
    })
  }
  
  return(tme_decon_list)
}

devo <- TME_deconvolution_all2(list_train_vali_Data[1],deconvolution_method = "cibersort") #list_train_vali_Data
#slow

#devo <- TME_deconvolution_all2(cohort[[1]])
#names(devo)<-"TCGA"
#Running cibersort_abs 慢  >4H #使用训练集  

#读取之前保存的
#setwd("~//R//my_projects//metabolism_sig//5.1-101建模");getwd()
#res <-readRDS("./out/101_results.rds")

pdf("./out/immuno_heatmap.pdf",width = 16,height = 8,family="serif") #,family="serif"
immuno_heatmap(res,
               devo,col=c("skyblue","pink"),
               model_name=my_model,#"StepCox[backward] + plsRcox",
               dataset="TCGA") #Dataset1 数据集名称TCGA
dev.off() #丑

##3.1 训练集cibersort热图----
##失败，提取数据自己作图-TCGA训练集-cibersort
#处理数据
immu_cibersort <-devo$TCGA$tme_combine$cibersort[,-c(1:2)]
immu_cibersort <-as.data.frame(t(immu_cibersort ) )
colnames(immu_cibersort)<-immu_cibersort[1,];immu_cibersort<-immu_cibersort[-1,]
immu_cibersort[,1:length(immu_cibersort)]<-lapply(immu_cibersort[,1:length(immu_cibersort)],as.numeric ) #转为数值型
# gsub("\\.","-",
#      substr(colnames(immu_cibersort)[2:length(colnames(immu_cibersort))],1,12 )
#      )
immu_cibersort$ID <-gsub("\\.","-",substr(rownames(immu_cibersort),1,12) )

score_df <-res$riskscore[my_model][[1]]$TCGA 
immu_cibersort_score <-merge(immu_cibersort,score_df,by="ID") #合并riskscore
immu_cibersort_score <- immu_cibersort_score %>% tibble::column_to_rownames('ID') 
immu_cibersort_score$Group <-ifelse(immu_cibersort_score$RS >= median(immu_cibersort_score$RS),"High", "Low") #riskscore group 
table(immu_cibersort_score$Group)[[1]]
immu_cibersort_score$OS <-ifelse(immu_cibersort_score$OS==1,"Dead","Alive")
#res$Sig.genes #64

#plot pheatmap 热图
library(pheatmap);library(RColorBrewer)
#样本注释

immu_cibersort_score<- arrange(immu_cibersort_score, Group)#排序
# pheatmap(immu_cibersort_score[,1:22],
#          show_rownames = F,cluster_rows = F,cluster_cols = T,
#          scale = "column",
#          cellwidth=3,cellheight = 1,
#          annotation_row =immu_cibersort_score[,c(24:26)],
#          gaps_row= table(immu_cibersort_score$Group)[[1]] ,border_color = NA, #gaps_row= c(7), 
#          color=alpha(rev(RColorBrewer::brewer.pal(10,"PiYG") ), 0.7) #RdYlGn Spectral   #PiYG "BrBG" 褐浅蓝  "PiYG"紫绿
#          #colorRampPalette( c("skyblue", "white", "firebrick3"))(10),#
#          )

ann_colors=list(#FDR=c('[0,0.05]'='#FDDCA9','(0.05,0.5]'="#DDC9E3"),
  Group=c("Low"="skyblue","High"="pink" ))

save_pheatmap_pdf <- function(x, filename, width, height) {
  library(grid)
  x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.02)#修改聚类树线条宽度
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数

pheat<-pheatmap(t(immu_cibersort_score[,1:22]),
         show_colnames = F,cluster_rows = T,cluster_cols = F,
         scale = "column",
         #cellwidth=3,cellheight = 1,
         annotation_col =immu_cibersort_score[,c(24:26)],annotation_colors = ann_colors,
         gaps_col= table(immu_cibersort_score$Group)[[1]] ,border_color = NA, #gaps_row= c(7), 
         annotation_names_row = F, annotation_names_col = F,
         color=alpha(rev(RColorBrewer::brewer.pal(10,"PiYG") ), 0.7) #RdYlGn Spectral   #PiYG "BrBG" 褐浅蓝  "PiYG"紫绿
         #colorRampPalette( c("skyblue", "white", "firebrick3"))(10),#
)
#save_pheatmap_pdf(pheat,"./out/cibersort_pheatmap.pdf",6,3.5) 

pdf("./out/cibersort_pheatmap.pdf",6,3.5,family = "serif") 
pheat
dev.off()

#计算显著性后重新绘制热图!
#T检验显著性表格
#exp_int[,i]
High <- rownames(immu_cibersort_score[immu_cibersort_score$Group =="High",])
Low <-rownames(immu_cibersort_score[immu_cibersort_score$Group !="High",])

pvalue_interest <-data.frame() #T.test
exp_int<- immu_cibersort_score[1:22]#exprSet[rownames(exprSet) %in% gene_more,] #raw FPKM #gene_more 兴趣基因集
for (i in 1:ncol(exp_int) ){ 
  print(i)
  pwilcox <- t.test(exp_int[rownames(exp_int) %in% High,i],exp_int[rownames(exp_int) %in% Low,i]) #样本数 High vs Low
  fc <- mean(exp_int[rownames(exp_int) %in% High,i] )/mean( exp_int[rownames(exp_int) %in% Low,i] )#drug/control  手动输入？ 对照组在前！
  pvalue_interest<-rbind(pvalue_interest,data.frame(pvalue=pwilcox$p.value,Fold=fc,Symbol=colnames(exp_int)[i] ) )
}
pvalue_interest$Sign. <- ifelse(pvalue_interest$pvalue>= 0.05,"ns",
                                ifelse(pvalue_interest$pvalue < 0.001,"***",
                                       ifelse(pvalue_interest$pvalue < 0.01,"**","*")
                                ) 
)
pvalue_interest$FC <-ifelse(pvalue_interest$Fold <=0.5,"Decrease",
                            ifelse(pvalue_interest$Fold >1,"Increase","Unchange") )
rownames(pvalue_interest) <-pvalue_interest$Symbol
#图注
anno_row_p <- data.frame(
  Sign.=pvalue_interest$Sign.,FoldChange=pvalue_interest$FC) #, size = 20, replace = TRUE  adj.P.Val
rownames(anno_row_p) <- rownames(pvalue_interest)
ann_colors_p=list(FoldChange=c('Unchange'='#CCECFF','Decrease'="skyblue",'Increase'="#4a6fe3"),  #'Unchange'='#CAB2D6',
                  Sign.=c("ns"="grey","*"="pink",'**'="#D33F6A","***"="red" ) )
#plot
int_p <-pheatmap(t(exp_int),
                 show_colnames = F,cluster_rows = T,cluster_cols = F,
                 scale = "column",
                 #cellwidth=3,cellheight = 1,
                 annotation_col =immu_cibersort_score[,c(24:26)],#annotation_colors = ann_colors,
                 gaps_col= table(immu_cibersort_score$Group)[[1]] ,border_color = NA, #gaps_row= c(7), 
                 annotation_names_row = F, annotation_names_col = F,
                 annotation_row = anno_row_p,annotation_colors = ann_colors_p,
                 color=alpha(rev(RColorBrewer::brewer.pal(10,"PiYG") ), 0.7) #RdYlGn Spectral   #PiYG "BrBG" 褐浅蓝  "PiYG"紫绿
                 #colorRampPalette( c("skyblue", "white", "firebrick3"))(10),#
)
pdf("./out/cibersort_pheatmap_sign.pdf",7,4.5,family = "serif") 
int_p
dev.off()  

#save_pheatmap_pdf(int_p, "./out/cibersort_pheatmap_sign.pdf",7,4.5) 




##3.2 cibersort分组柱状图----
library(reshape2)
bar_df <-melt(immu_cibersort_score[c(1:22,26)],id.vars=c("Group") ) #exp_int, id.vars=c("Cell_type")
my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 

library(ggplot2);library(ggpubr)
p_method <-c("wilcox.test","t.test")
for (j in p_method){
  ggplot(bar_df,#
         aes(variable,value,color=Group))+  #age_grou
    geom_point(alpha=0.7,size=0.2,
               position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                             jitter.height = 0,
                                             dodge.width = 0.7))+
    geom_boxplot(alpha=1,width=0.7,fill=NA,
                 position=position_dodge(width=0.8),
                 size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    labs(x="",y="Cell proportion")+
    scale_color_manual(values =c("#EDA065","#66CCFF","#7EC7A7") )+ # #blue c("#0066CC","#66CCFF","#9DC3E6")
    theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
    theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
          axis.text.x = element_text(angle = 90,hjust = 1),
          panel.grid = element_blank(),legend.key = element_blank() ,legend.position="top" ) + #
    # stat_compare_means(#aes(group = Group) ,
    # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
    # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
    # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
    # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
    stat_compare_means(method = j,size=4, #"wilcox.test"
                       show.legend= F,label = "p.signif",#p.signif p.format
                       label.x =0.75,label.y =max(bar_df$value)-0.05)#
  ggsave(paste0("./out/cibersort_boxplot_",j,".pdf"),width = 12,height = 6,family="serif")
  
}

##3.3 单独免疫细胞boxplot----


for (i in as.character(unique(bar_df$variable) ) ){
  for (j in p_method){
    ggplot(bar_df[bar_df$variable %in% i,],#
           aes(Group,value,color=Group))+
      # geom_boxplot(alpha=1,width=0.45,fill=NA,
      #              position=position_dodge(width=0.8),
      #              size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
      #              outlier.stroke = 0.5)+
      geom_point(alpha=0.7,size=2,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.8))+
      geom_violin(alpha=0.2,width=0.9,
                  position=position_dodge(width=0.8),
                  size=0.25)+
      labs(x=i,y="Cell proportion")+
      scale_color_manual(values =c("#EDA065","#66CCFF","#7EC7A7") )+ # #blue c("#0066CC","#66CCFF","#9DC3E6")
      theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
      theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
            axis.text.x = element_text(angle = 90,hjust = 1),
            panel.grid = element_blank(),legend.key = element_blank() ,legend.position="none" ) + #
      stat_compare_means(method = j,size=8, #"wilcox.test"
                         show.legend= F,label = "p.signif",#p.signif p.format
                         label.x =1.5,label.y =max(bar_df[bar_df$variable %in% i,]$value)-0.02)#label.x =0.75,
    ggsave(paste0("./out/cibersort_boxplot_",i,"_",j,".pdf"),width = 4,height = 4,family="serif")
    
  }
}




#4.0 Construct predicting models for response----
#免疫效应预测模型
#https://github.com/l-magnificence/Mime/tree/main
load("./data/Example.ici.Rdata")
load("./data/genelist.Rdata")
res.ici <- ML.Dev.Pred.Category.Sig(train_data = list_train_vali_Data$training,
                                    list_train_vali_Data = list_train_vali_Data,
                                    candidate_genes = genelist,
                                    methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                    seed = 5201314,
                                    cores_for_parallel = 60
)
#Plot AUC of different methods among different datasets:
auc_vis_category_all(res.ici,dataset = c("training","validation"),
                     order= c("training","validation"))

#roc
plot_list<-list()
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res.ici,model_name = i,dataset = c("training","validation"),
                                   order= c("training","validation"),
                                   anno_position=c(0.4,0.25))
}
aplot::plot_list(gglist=plot_list,ncol=3)

#Compared AUC with other published models associated with immunotherapy response: !!!!!!!!!!! 比较
auc.other.pre <- cal_auc_previous_sig(list_train_vali_Data = list_train_vali_Data,seed = 5201314,
                                      train_data = list_train_vali_Data$training,
                                      cores_for_parallel = 32)
auc_category_comp(res.ici,
                  auc.other.pre,
                  model_name="svmRadialWeights",
                  dataset=names(list_train_vali_Data))

#5.0 Core feature selection----
load("./data/Example.cohort.Rdata")
load("./data/genelist.Rdata")
res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$Dataset1,
                                              candidate_genes = genelist,  
                                              mode = "all_without_SVM",nodesize =5,seed = 5201314 ) #slow
                                              #all, single, and all_without_SVM  #svm-ref This step will probably take several hours

#Upset plot of genes filtered by different methods:
core_feature_select(res.feature.all)
#The output genes are closely associated with patient outcome and higher frequence of screening means more critical.
core_feature_rank(res.feature.all, top=20)

#Here, we randomly select top two genes to analyze their correlation:
dataset_col<-c("#3182BDFF","#E6550DFF")
corplot <- list()
for (i in c(1:2)) {
  print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
                               dataset=names(list_train_vali_Data)[i],
                               color = dataset_col[i],
                               feature1="PSEN2",
                               feature2="WNT5B",
                               method="pearson"))
}
aplot::plot_list(gglist=corplot,ncol=2)

#Plot survival curve of patients according to median expression level of specific gene among different datasets:
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-core_feature_sur("PSEN2", 
                                        InputMatrix=list_train_vali_Data[[i]],
                                        dataset = names(list_train_vali_Data)[i],
                                        #color=c("blue","green"),
                                        median.line = "hv",
                                        cutoff = 0.5,
                                        conf.int = T,
                                        xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)
