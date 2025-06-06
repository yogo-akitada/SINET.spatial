
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(gridExtra)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(SingleR)
library(hdf5r)

## Place the file 'filtered_feature_bc_matrix.h5' and the image folder 'spatial' into the directory named '376A1'
Area1 <-Load10X_Spatial(
  "367A1",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",filter.matrix = TRUE
)

plot1 <- VlnPlot(Area1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Area1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

Area1@meta.data$Area <- "Area1"
table(Area1$Area)

Area2 <-Load10X_Spatial(
  "367B1",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",filter.matrix = TRUE
)

plot1 <- VlnPlot(Area2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Area2, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

Area2@meta.data$Area <- "Area2"
table(Area2$Area)

Case1 <- merge(Area1, y = Area2, add.cell.ids = c("Area1", "Area2"))
table(Case1$Area)



Case1 <- NormalizeData(Case1, normalization.method = "LogNormalize", scale.factor = 10000)
Case1 <- FindVariableFeatures(Case1, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(Case1), points = head(VariableFeatures(Case1), 25), repel = TRUE)
Case1 <- ScaleData(Case1, features = rownames(Case1))
Case1 <- RunPCA(Case1, features = VariableFeatures(object = Case1))
VizDimLoadings(Case1, dims = 1:2, reduction = "pca", balanced = T,
               nfeatures = 50)
DimPlot(Case1, reduction = "pca")

Case1 <- FindNeighbors(Case1, dims = 1:20)
Case1 <- FindClusters(Case1, resolution = 1)

table(Case1$seurat_clusters)

Case1 <- RunUMAP(Case1, dims = 1:10,n.neighbors = 200 , min.dist= 0.5)
p1 <- DimPlot(Case1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(Case1, 
                     label = TRUE, 
                     label.size = 3, 
                     pt.size.factor = 0.7, 
                     alpha = 0.5,
                     crop = FALSE)
test <-p1 + p2
test


Idents(Case1) <- "seurat_clusters"


manual.ident2 <- c(  
  "mp_outer",
  "mp_peritumoral",
  "myenteric plexus",
  "tp#1sm_&_#1mp",
  "tm",
  "mp_inner",
  "tp#2sm",
  "sm",
  "tp#1m_&_#2m",
  "m_mid_crypt",
  "m_base",
  "tp#1mp_&_#2mp",
  "tm_pni_B",
  "tm_pni_A",
  "m_superficial",
  "mes_stroma",
  "serosa",
  "mes_nerve")

names(manual.ident2) <- levels(Case1)
Case1 <- RenameIdents(Case1, manual.ident2)
levels(Case1)
Case1@meta.data$manual_ident2 <- factor(Case1@active.ident)

Idents(Case1) <- "manual_ident2"
levels(Case1)


levels(Case1) <- c(  "m_superficial",
                     "m_mid_crypt",
                     "m_base",
                     "sm",
                     "mp_inner",
                     "mp_outer",
                     "mp_peritumoral",
                     "myenteric plexus",
                     "serosa",
                     "mes_stroma",
                     "mes_nerve",
                     
                     "tp#1m_&_#2m",
                     "tp#1sm_&_#1mp",
                     "tp#2sm",
                     "tp#1mp_&_#2mp",
                     "tm_pni_A",
                     "tm_pni_B",
                     "tm")
levels(Case1)
Case1@meta.data$manual_ident2 <- factor(Case1@active.ident)


p1 <- DimPlot(Case1, reduction = "umap", label = TRUE,label.size = 2.5) 
p1
p2 <- SpatialDimPlot(Case1,
                     label = FALSE, 
                     label.size = 3, 
                     pt.size.factor = 1.5, 
                     alpha = 0.5,
                     crop = FALSE,information = NULL)+ NoLegend()
p3 <- SpatialDimPlot(Case1,
                     label = FALSE, 
                     label.size = 3, 
                     pt.size.factor = 1.5, 
                     alpha = 0,
                     crop = FALSE,information = NULL)+ NoLegend()

renameplot<-
  cowplot::plot_grid(p1,p2,p3,
                     nrow = 3,
                     ncol = 1,
                     
                     align = 'vh')


renameplot
ggsave(filename="Case1_rename.svg",width=10, height=15,plot=renameplot)



manual.ident3 <- c(  
  "m_superficial",
  "m_mid_crypt",
  "m_base",
  "sm",
  "mp_inner",
  "mp_outer",
  "mp_peritumoral",
  "myenteric plexus",
  "serosa",
  "mes_stroma",
  "mes_nerve",
  
  "tp_mucosal",
  "tp_deep",
  "tp_deep",
  "tp_deep",
  "tm",
  "tm",
  "tm")

names(manual.ident3) <- levels(Case1)
Case1 <- RenameIdents(Case1, manual.ident3)
levels(Case1)
Case1@meta.data$manual_ident3 <- factor(Case1@active.ident)

Idents(Case1) <- "manual_ident3"
levels(Case1)



### quality checking by tumor gene module scores
Idents(Case1) <- "manual_ident2"
VlnPlot(Case1 , features = "nCount_Spatial", log = FALSE, pt.size=0.1)+ NoLegend()

Idents(Case1) <- "seurat_clusters"


NT.ids <- c(  
  "normal",
  "normal",
  "normal",
  "tumor",
  "tumor",
  "normal",
  "tumor",
  "normal",
  "tumor",
  "normal",
  "normal",
  "tumor",
  "tumor",
  "tumor",
  "normal",
  "normal",
  "normal",
  "normal")
names(NT.ids) <- levels(Case1)
Case1 <- RenameIdents(Case1, NT.ids)
levels(Case1)
Case1@meta.data$NT_ident <- factor(Case1@active.ident)

Idents(Case1) <- "NT_ident"
levels(Case1)

NTall.markers <- FindMarkers(Case1, ident.1 = c("tumor"), ident.2 = c("normal"))
                             

top10 <-
  NTall.markers[order(NTall.markers$p_val_adj, decreasing = FALSE),]  %>%
  dplyr::filter(avg_log2FC > 0.3) %>%
  slice_head(n = 10) %>%
  ungroup() 
top10
rownames(top10)
[1] "CHGA"   "TTR"    "CHGB"   "PTPRN2" "TPH1"   "FEV"    "NKX2-2" "RAB3B"  "RPH3AL" "PCSK1N"

tumor_genes <- c("TTR","CHGA","CHGB","PTPRN2","TPH1","FEV","NKX2-2","RAB3B","RPH3AL","PCSK1N")

Case1 <- AddModuleScore(
  object =Case1,
  features = tumor_genes,
  name = 'tumor_genes'
)

Idents(Case1) <- "manual_ident2"
VlnPlot(Case1 , features = "tumor_genes1", log = FALSE, pt.size=0.1)+ NoLegend() 

Case1_tumor_normal <- subset(Case1,
                               tumor_genes1>3 & manual_ident2 =="tp#1m_&_#2m"|
                               tumor_genes1>3 & manual_ident2 =="tp#1sm_&_#1mp"|
                               tumor_genes1>3 & manual_ident2 =="tp#2sm"|
                               tumor_genes1>3 & manual_ident2 =="tp#1mp_&_#2mp"|
                               tumor_genes1>3 & manual_ident2 =="tm_pni_A"|
                               tumor_genes1>3 & manual_ident2 =="tm_pni_B"|
                               tumor_genes1>3 & manual_ident2 =="tm"|
                               tumor_genes1<2 & manual_ident2 =="m_superficial"|
                               tumor_genes1<2 & manual_ident2 =="m_mid_crypt"|
                               tumor_genes1<2 & manual_ident2 =="m_base"|
                               tumor_genes1<2 & manual_ident2 =="sm"|
                               tumor_genes1<2 & manual_ident2 =="mp_inner"|
                               tumor_genes1<2 & manual_ident2 =="mp_outer"|
                               tumor_genes1<2 & manual_ident2 =="mp_peritumoral"|
                               tumor_genes1<2 & manual_ident2 =="myenteric plexus"|
                               tumor_genes1<2 & manual_ident2 =="serosa"|
                               tumor_genes1<2 & manual_ident2 =="mes_stroma"|
                               tumor_genes1<2 & manual_ident2 =="mes_nerve")


Case1_tumor_normal <- subset(Case1_tumor_normal, ACTB>0 & B2M>0 &UBC>0)
save(Case1,file="Case1.RData")
save(Case1_tumor_normal,file="Case1_tumor_normal.RData")

VlnPlot(Case1_tumor_normal , features = "tumor_genes1", log = FALSE, pt.size=0.1)+ NoLegend()

SpatialDimPlot(Case1_tumor_normal, 
               label = FALSE, 
               label.size = 3, 
               pt.size.factor = 0.7, 
               alpha = 0.5,
               crop = FALSE)



Case1_tumor_normal <- RunUMAP(Case1_tumor_normal, dims = 1:10,n.neighbors = 200 , min.dist= 0.5)
UMAP_Case1_tumor_normal<-
  DimPlot(Case1_tumor_normal, reduction = "umap", label = TRUE) + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  NoLegend()
UMAP_Case1_tumor_normal
ggsave(filename="UMAP_Case1_tumor_normal.svg",width=5, height=5,plot=UMAP_Case1_tumor_normal)
UMAP_Case1_tumor_normal_w_legend<-
  DimPlot(Case1_tumor_normal, reduction = "umap", label = TRUE) + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
UMAP_Case1_tumor_normal_w_legend
ggsave(filename="UMAP_Case1_tumor_normal_w_legend.svg",width=5, height=5,plot=UMAP_Case1_tumor_normal_w_legend)











## deconvolution analysis

## https://jef.works/STdeconvolve/getting_started.html

library(tidyverse)
library(Seurat)
library(STdeconvolve)

cd <- Case1_tumor_normal@assays$Spatial@counts ## matrix of gene counts in each pixel
annot <- Case1_tumor_normal@meta.data$manual_ident3 ## annotated tissue layers assigned to each pixel

Case1_tumor_normal@images$slice1@coordinates

Idents(Case1_tumor_normal) <- "Area"
Case1_tumor_normal_area1<-subset(Case1_tumor_normal, Area == "Area1")
Case1_tumor_normal_area2<-subset(Case1_tumor_normal, Area == "Area2")

pos1<-Case1_tumor_normal_area1@images$slice1@coordinates[,4:5]
names(pos1)<-c("x","y")
pos1


pos2<-Case1_tumor_normal_area2@images$slice1.2@coordinates[,4:5]
names(pos2)<-c("x","y")
pos2

summary(pos2$x)
pos2$x<- pos2$x+10000
summary(pos2$x)

pos <- rbind(pos1,pos2)


counts <- cleanCounts(counts = cd, 
                      min.lib.size = 10,
                      min.reads = 1,
                      min.detected = 1,
                      verbose = TRUE)



odGenes <- getOverdispersedGenes(as.matrix(counts),
                                 gam.k=5,
                                 alpha=0.05,
                                 plot=FALSE,
                                 use.unadjusted.pvals=FALSE,
                                 do.par=TRUE,
                                 max.adjusted.variance=1e3,
                                 min.adjusted.variance=1e-3,
                                 verbose=FALSE, details=TRUE)


genes <- odGenes$ods


NET<-rownames(top10)
NET_overlap_fit<-rownames(counts)%in%NET
NET_overlap<-rownames(counts)[NET_overlap_fit]
gene_NET<-c(genes,NET_overlap)%>%unique()



corpus<-preprocess(t(as.matrix(counts)),
                   selected.genes = gene_NET,plot=FALSE,
                   min.reads = 1, 
                   min.lib.size = 10, 
                   min.detected = 1,
                   ODgenes = FALSE, 
                   verbose = TRUE)


ldas_case1 <- fitLDA(corpus$corpus, Ks = seq(2, 9, by = 1),
                     perc.rare.thresh = 0.05,
                     plot=TRUE,
                     verbose=TRUE)
save(ldas_case1,file="ldas_case1.RData")

## select model with minimum perplexity
optLDA <- optimalModel(models = ldas_case1, opt = "proportion")

## extract pixel cell-type proportions (theta) and cell-type gene expression profiles (beta) for the given dataset
## we can also remove cell-types from pixels that contribute less than 5% of the pixel proportion
## and scale the deconvolved transcriptional profiles by 1000 
results <- getBetaTheta(optLDA,
                        perc.filt = 0.05,
                        betaScale = 1000)

deconProp <- results$theta
deconProp 

deconGexp <- results$beta
deconProp 



deconpie_Case1_tumor_normal<-
  vizAllTopics(deconProp, pos, 
               topicCols = c("tomato","violet","springgreen","violetred2","turquoise1","springgreen4","steelblue2"),
               r=35,
               groups = annot,
               group_cols = rainbow(length(levels(annot))),)
deconpie_Case1_tumor_normal
ggsave(filename="deconpie_Case1_tumor_normal.svg",width=20, height=20,plot=deconpie_Case1_tumor_normal)



dfdeconProp <- as.data.frame(deconProp)
colnames(dfdeconProp)<- c("tumor1","tumor2","smooth_muscle","tumor3","stroma","m_base","m_superficial_crypt")
dfdeconProp$manual_ident3 <- factor(annot)
dfdeconProp 


dfdeconProp2<-
  dfdeconProp  %>%
  dplyr::group_by(manual_ident3) 


summarize(tumor1 = mean(tumor1),
          tumor2 = mean(tumor2),
          tumor3 = mean(tumor3),
          m_superficial_crypt = mean(m_superficial_crypt),
          m_base = mean(m_base),
          stroma = mean(stroma),
          smooth_muscle = mean(smooth_muscle))

dfdeconProp2

library(reshape2)

dfdeconProp3 <- melt(dfdeconProp2,id = c("manual_ident3"),measure = c("tumor1","tumor2","tumor3","m_superficial_crypt","m_base","stroma", "smooth_muscle"))
dfdeconProp3


deconbar_Case1_tumor_normal  <-
  dfdeconProp3  %>%
  ggplot()+
  aes(x= manual_ident3,
      y = value,
      fill = variable)+
  geom_bar(stat = "identity", 
           position = "fill") +
  ylab("% of cells")+  
  xlab(NULL)+  
  theme_classic()+
  scale_y_continuous(expand = c(0, 0), labels = percent)+
  scale_x_discrete(limits = c("m_superficial",
                              "m_mid_crypt",
                              "m_base",
                              "sm",
                              "mp_inner",
                              "mp_outer",
                              "mp_peritumoral",
                              "myenteric plexus",
                              "serosa",
                              "mes_stroma",
                              "mes_nerve",
                              
                              "tp_mucosal",
                              "tp_deep",
                              "tm"))+
  scale_fill_manual(values=c("tomato","violet","violetred2","steelblue2","springgreen4","turquoise1","springgreen"))+
  theme(axis.text.x=element_text(angle=75, hjust=1,  size = rel(1)),
        legend.title=element_blank())
deconbar_Case1_tumor_normal

ggsave(filename="deconbar_Case1_tumor_normal.svg",width=6, height=4,plot=deconbar_Case1_tumor_normal)




