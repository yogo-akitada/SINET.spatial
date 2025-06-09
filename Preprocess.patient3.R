
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
Area5 <-Load10X_Spatial(
  "data/356A1",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",filter.matrix = TRUE
)

plot1 <- VlnPlot(Area5, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Area5, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

Area5@meta.data$Area <- "Area5"
table(Area5$Area)

Area6 <-Load10X_Spatial(
  "data/356B1",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",filter.matrix = TRUE
)

plot1 <- VlnPlot(Area6, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Area6, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

Area6@meta.data$Area <- "Area6"
table(Area6$Area)

Case3 <- merge(Area5, y = Area6, add.cell.ids = c("Area5", "Area6"))
table(Case3$Area)



Case3 <- NormalizeData(Case3, normalization.method = "LogNormalize", scale.factor = 10000)
Case3 <- FindVariableFeatures(Case3, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(Case3), points = head(VariableFeatures(Case3), 25), repel = TRUE)
Case3 <- ScaleData(Case3, features = rownames(Case3))
Case3 <- RunPCA(Case3, features = VariableFeatures(object = Case3))
VizDimLoadings(Case3, dims = 1:2, reduction = "pca", balanced = T,
               nfeatures = 50)
DimPlot(Case3, reduction = "pca")

Case3 <- FindNeighbors(Case3, dims = 1:20)
Case3 <- FindClusters(Case3, resolution = 1.2)

table(Case3$seurat_clusters)

Case3 <- RunUMAP(Case3, dims = 1:10,n.neighbors = 200 , min.dist= 0.5)
p1 <- DimPlot(Case3, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(Case3, 
                     label = TRUE, 
                     label.size = 3, 
                     pt.size.factor = 0.7, 
                     alpha = 0.5,
                     crop = FALSE)
test <-p1 + p2
test


Idents(Case3) <- "seurat_clusters"


manual.ident2 <- c(                     
  "tp#1m_&_#2m",
  "sm_non_peritumoral",
  "m_base",
  "tp#1sm",
  "tp#1mp",
  "tm_tp#2sm_&_perithrombus",
  "mp_non_peritumoral",
  "tp#2m",
  "tp#2sm_&_thrombus",
  "sm_peritumoral",
  "tl",
  "mp_peritumoral",
  "m_lymphoid_follicle",
  "sm_deep",
  "mes_collagen",
  "mp_inner_peritumoral")

names(manual.ident2) <- levels(Case3)
Case3 <- RenameIdents(Case3, manual.ident2)
levels(Case3)
Case3@meta.data$manual_ident2 <- factor(Case3@active.ident)

Idents(Case3) <- "manual_ident2"
levels(Case3)


levels(Case3) <- c(  
  "m_base",
  "m_lymphoid_follicle",
  "sm_non_peritumoral",
  "sm_peritumoral",
  "sm_deep",
  "mp_non_peritumoral",
  "mp_peritumoral",
  "mp_inner_peritumoral",
  "mes_collagen",
  
  "tp#1m_&_#2m",
  "tp#2m",
  "tp#1sm",
  "tp#2sm_&_thrombus",
  "tp#1mp",
  "tm_tp#2sm_&_perithrombus",
  "tl")
levels(Case3)
Case3@meta.data$manual_ident2 <- factor(Case3@active.ident)


p1 <- DimPlot(Case3, reduction = "umap", label = TRUE,label.size = 2.5) 
p1
p2 <- SpatialDimPlot(Case3,
                     label = FALSE, 
                     label.size = 3, 
                     pt.size.factor = 1.5, 
                     alpha = 0.5,
                     crop = FALSE,information = NULL)+ NoLegend()
p3 <- SpatialDimPlot(Case3,
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
ggsave(filename="Case3_rename.svg",width=10, height=15,plot=renameplot)



manual.ident3 <- c(  
  "m_base",
  "m_lymphoid_follicle",
  "sm_non_peritumoral",
  "sm_peritumoral",
  "sm_deep",
  "mp_non_peritumoral",
  "mp_peritumoral",
  "mp_inner_peritumoral",
  "mes_collagen",
  
  "tp_mucosal",
  "tp_mucosal",
  "tp_other",
  "tm",
  "tp_other",
  "tm",
  "tl")

names(manual.ident3) <- levels(Case3)
Case3 <- RenameIdents(Case3, manual.ident3)
levels(Case3)
Case3@meta.data$manual_ident3 <- factor(Case3@active.ident)

Idents(Case3) <- "manual_ident3"
levels(Case3)



### quality checking by tumor gene module scores
Idents(Case3) <- "manual_ident2"
VlnPlot(Case3 , features = "nCount_Spatial", log = FALSE, pt.size=0.1)+ NoLegend()

Idents(Case3) <- "seurat_clusters"


NT.ids <- c(
  "tumor",
  "normal",
  "normal",
  "tumor",
  "tumor",
  "tumor",
  "normal",
  "tumor",
  "tumor",
  "normal",
  "tumor",
  "normal",
  "normal",
  "normal",
  "normal",
  "normal")
names(NT.ids) <- levels(Case3)
Case3 <- RenameIdents(Case3, NT.ids)
levels(Case3)
Case3@meta.data$NT_ident <- factor(Case3@active.ident)

Idents(Case3) <- "NT_ident"
levels(Case3)

NTall.markers <- FindMarkers(Case3, ident.1 = c("tumor"), ident.2 = c("normal"))
                             

top10 <-
  NTall.markers[order(NTall.markers$p_val_adj, decreasing = FALSE),]  %>%
  dplyr::filter(avg_log2FC > 0.3) %>%
  slice_head(n = 10) %>%
  ungroup() 
top10
rownames(top10)
[1] "CHGA"   "TTR"    "CHGB"   "TPH1"   "PCSK1"  "RPH3AL" "NKX2-2" "MS4A8"  "FEV"    "ECM1"  

tumor_genes <- c("TTR","CHGA","CHGB","TPH1","PCSK1","RPH3AL","NKX2-2","MS4A8","FEV","ECM1"  )

Case3 <- AddModuleScore(
  object =Case3,
  features = tumor_genes,
  name = 'tumor_genes'
)

Idents(Case3) <- "manual_ident2"
VlnPlot(Case3 , features = "tumor_genes1", log = FALSE, pt.size=0.1)+ NoLegend() 

Case3_tumor_normal <- subset(Case3,
                               tumor_genes1>2 & manual_ident2 =="tp#1m_&_#2m"|
                               tumor_genes1>2 & manual_ident2 =="tp#2m"|
                               tumor_genes1>2 & manual_ident2 =="tp#1sm"|
                               tumor_genes1>2 & manual_ident2 =="tp#2sm_&_thrombus"|
                               tumor_genes1>2 & manual_ident2 =="tp#1mp"|
                               tumor_genes1>2 & manual_ident2 =="tm_tp#2sm_&_perithrombus"|
                               tumor_genes1>2 & manual_ident2 =="tl"|
                               tumor_genes1<1.5 & manual_ident2 =="m_base"|
                               tumor_genes1<1.5 & manual_ident2 =="m_lymphoid_follicle"|
                               tumor_genes1<1.5 & manual_ident2 =="sm_non_peritumoral"|
                               tumor_genes1<1.5 & manual_ident2 =="sm_peritumoral"|
                               tumor_genes1<1.5 & manual_ident2 =="sm_deep"|
                               tumor_genes1<1.5 & manual_ident2 =="mp_non_peritumoral"|
                               tumor_genes1<1.5 & manual_ident2 =="mp_peritumoral"|
                               tumor_genes1<1.5 & manual_ident2 =="mp_inner_peritumoral"|
                               tumor_genes1<1.5 & manual_ident2 =="mes_collagen")


Case3_tumor_normal <- subset(Case3_tumor_normal, ACTB>0 & B2M>0 &UBC>0)
save(Case3,file="Case3.RData")
save(Case3_tumor_normal,file="Case3_tumor_normal.RData")

VlnPlot(Case3_tumor_normal , features = "tumor_genes1", log = FALSE, pt.size=0.1)+ NoLegend()

SpatialDimPlot(Case3_tumor_normal, 
               label = FALSE, 
               label.size = 3, 
               pt.size.factor = 0.7, 
               alpha = 0.5,
               crop = FALSE)



Case3_tumor_normal <- RunUMAP(Case3_tumor_normal, dims = 1:10,n.neighbors = 250 , min.dist= 0.5)
UMAP_Case3_tumor_normal<-
  DimPlot(Case3_tumor_normal, reduction = "umap", label = TRUE) + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  NoLegend()
UMAP_Case3_tumor_normal
ggsave(filename="UMAP_Case3_tumor_normal.svg",width=5, height=5,plot=UMAP_Case3_tumor_normal)
UMAP_Case3_tumor_normal_w_legend<-
  DimPlot(Case3_tumor_normal, reduction = "umap", label = TRUE) + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
UMAP_Case3_tumor_normal_w_legend
ggsave(filename="UMAP_Case3_tumor_normal_w_legend.svg",width=5, height=5,plot=UMAP_Case3_tumor_normal_w_legend)











## deconvolution analysis

## https://jef.works/STdeconvolve/getting_started.html

library(tidyverse)
library(Seurat)
library(STdeconvolve)

cd <- Case3_tumor_normal@assays$Spatial@counts ## matrix of gene counts in each pixel
annot <- Case3_tumor_normal@meta.data$manual_ident3 ## annotated tissue layers assigned to each pixel

Case3_tumor_normal@images$slice1@coordinates

Idents(Case3_tumor_normal) <- "Area"
Case3_tumor_normal_Area5<-subset(Case3_tumor_normal, Area == "Area5")
Case3_tumor_normal_Area6<-subset(Case3_tumor_normal, Area == "Area6")

pos1<-Case3_tumor_normal_Area5@images$slice1@coordinates[,4:5]
names(pos1)<-c("x","y")
pos1


pos2<-Case3_tumor_normal_Area6@images$slice1.2@coordinates[,4:5]
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
####



ldas_Case3 <- fitLDA(corpus$corpus, Ks = seq(2, 9, by = 1),
                     perc.rare.thresh = 0.05,
                     plot=TRUE,
                     verbose=TRUE)
save(ldas_Case3,file="ldas_Case3.RData")

## select model with minimum perplexity
optLDA <- optimalModel(models = ldas_Case3, opt = "proportion")

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



deconpie_Case3_tumor_normal<-
  vizAllTopics(deconProp, pos, 
               topicCols = c("tomato","springgreen4","violet","violetred2","turquoise1","steelblue2"),
               r=35,
               groups = annot,
               group_cols = rainbow(length(levels(annot))),)
deconpie_Case3_tumor_normal
ggsave(filename="deconpie_Case3_tumor_normal.svg",width=20, height=20,plot=deconpie_Case3_tumor_normal)

##run_me_QC(ldas_Case3,starting_k=2,ending_k=5,dir="output_Case3_tumor_normal_1/")

# the following # commented commands allow you to save the "ldas" variable (line 186 above) 
# to an RDS file
# more information about RDS files: 
# https://rstudio-education.github.io/hopr/dataio.html#saving-r-files
# remove the # (comment) to run the command:

# saveRDS(object = ldas,file = "optlDA.17p9_rep1_astrogenes.rds")

# here is the command to load our RDS file (downloaded on line 26 above) 
# remove the # (comment) to run the command:

# ldas<-readRDS(file = "optlDA.17p9_rep1_astrogenes.rds")

## exporting spatial plots and topics for the optimal model
## This function assume a gene expresison cut-off of 2.

##run_me_results(opt=18,dir = "output_Case3_tumor_normal_1/",ldas=ldas_Case3)



dfdeconProp <- as.data.frame(deconProp)
colnames(dfdeconProp)<- c("tumor1","lymphoid","tumor2","tumor3","sm_mp","mucosa")
dfdeconProp$manual_ident3 <- factor(annot)
dfdeconProp 


dfdeconProp2<-
  dfdeconProp  %>%
  dplyr::group_by(manual_ident3) 

dfdeconProp2

library(reshape2)

dfdeconProp3 <- melt(dfdeconProp2,id = c("manual_ident3"),measure = c("tumor1","tumor2","tumor3","mucosa","lymphoid","sm_mp"))
dfdeconProp3


deconbar_Case3_tumor_normal  <-
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
  scale_x_discrete(limits = c("m_base",
                              "m_lymphoid_follicle",
                              "sm_non_peritumoral",
                              "sm_peritumoral",
                              "sm_deep",
                              "mp_non_peritumoral",
                              "mp_peritumoral",
                              "mp_inner_peritumoral",
                              "mes_collagen",
                              
                              "tp_mucosal",
                              "tp_other",
                              "tm",
                              "tl"))+
  scale_fill_manual(values=c("tomato","violet","violetred2","steelblue2","springgreen4","turquoise1"))+
  theme(axis.text.x=element_text(angle=75, hjust=1,  size = rel(1)),
        legend.title=element_blank())
deconbar_Case3_tumor_normal

ggsave(filename="deconbar_Case3_tumor_normal.svg",width=6, height=4,plot=deconbar_Case3_tumor_normal)


