
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
Area3 <-Load10X_Spatial(
  "data/367C1",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",filter.matrix = TRUE
)

plot1 <- VlnPlot(Area3, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Area3, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

Area3@meta.data$Area <- "Area3"
table(Area3$Area)

Area4 <-Load10X_Spatial(
  "data/367D1",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",filter.matrix = TRUE
)

plot1 <- VlnPlot(Area4, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Area4, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

Area4@meta.data$Area <- "Area4"
table(Area4$Area)

Case2 <- merge(Area3, y = Area4, add.cell.ids = c("Area3", "Area4"))
table(Case2$Area)



Case2 <- NormalizeData(Case2, normalization.method = "LogNormalize", scale.factor = 10000)
Case2 <- FindVariableFeatures(Case2, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(Case2), points = head(VariableFeatures(Case2), 25), repel = TRUE)
Case2 <- ScaleData(Case2, features = rownames(Case2))
Case2 <- RunPCA(Case2, features = VariableFeatures(object = Case2))
VizDimLoadings(Case2, dims = 1:2, reduction = "pca", balanced = T,
               nfeatures = 50)
DimPlot(Case2, reduction = "pca")

Case2 <- FindNeighbors(Case2, dims = 1:20)
Case2 <- FindClusters(Case2, resolution = 1.2)

table(Case2$seurat_clusters)

Case2 <- RunUMAP(Case2, dims = 1:10,n.neighbors = 200 , min.dist= 0.5)
p1 <- DimPlot(Case2, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(Case2, 
                     label = TRUE, 
                     label.size = 3, 
                     pt.size.factor = 0.7, 
                     alpha = 0.5,
                     crop = FALSE)
test <-p1 + p2
test


Idents(Case2) <- "seurat_clusters"


manual.ident2 <- c(                     
"tp#2sm",
"sm",
"mp_inner",
"tp#1m_&_#2m",
"mes_dense_collagen_stroma",
"tm#1pni",
"tm#2pni",
"tp#1sm_&_#1mp",
"m_base",
"vascular_wall",
"serosa",
"mp_outer",
"tm#1_&_tp#1mp",
"m_superficial_mid_crypt",
"tm#2",
"tp#1superficialsm_&_#2superficialsm",
"mes_nerve")



names(manual.ident2) <- levels(Case2)
Case2 <- RenameIdents(Case2, manual.ident2)
levels(Case2)
Case2@meta.data$manual_ident2 <- factor(Case2@active.ident)

Idents(Case2) <- "manual_ident2"
levels(Case2)


levels(Case2) <- c(  
  "m_superficial_mid_crypt",
  "m_base",
  "sm",
  "mp_inner",
  "mp_outer",
  "serosa",
  "vascular_wall",
  "mes_dense_collagen_stroma",
  "mes_nerve",
  
  "tp#1m_&_#2m",
  "tp#1superficialsm_&_#2superficialsm",
  "tp#1sm_&_#1mp",
  "tp#2sm",
  "tm#1_&_tp#1mp",
  "tm#1pni",
  "tm#2",
  "tm#2pni")
levels(Case2)
Case2@meta.data$manual_ident2 <- factor(Case2@active.ident)


p1 <- DimPlot(Case2, reduction = "umap", label = TRUE,label.size = 2.5) 
p1
p2 <- SpatialDimPlot(Case2,
                     label = FALSE, 
                     label.size = 3, 
                     pt.size.factor = 1.5, 
                     alpha = 0.5,
                     crop = FALSE,information = NULL)+ NoLegend()
p3 <- SpatialDimPlot(Case2,
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
ggsave(filename="Case2_rename.svg",width=10, height=15,plot=renameplot)



manual.ident3 <- c(  
  "m_superficial_mid_crypt",
  "m_base",
  "sm",
  "mp_inner",
  "mp_outer",
  "serosa",
  "vascular_wall",
  "mes_dense_collagen_stroma",
  "mes_nerve",
  
  "tp_mucosal",
  "tp_deep",
  "tp_deep",
  "tp_deep",
  "tm",
  "tm",
  "tm",
  "tm")

names(manual.ident3) <- levels(Case2)
Case2 <- RenameIdents(Case2, manual.ident3)
levels(Case2)
Case2@meta.data$manual_ident3 <- factor(Case2@active.ident)

Idents(Case2) <- "manual_ident3"
levels(Case2)



### quality checking by tumor gene module scores
Idents(Case2) <- "manual_ident2"
VlnPlot(Case2 , features = "nCount_Spatial", log = FALSE, pt.size=0.1)+ NoLegend()

Idents(Case2) <- "seurat_clusters"


NT.ids <- c(
  "tumor",
  "normal",
  "normal",
  "tumor",
  "normal",
  "tumor",
  "tumor",
  "tumor",
  "normal",
  "normal",
  "normal",
  "normal",
  "tumor",
  "normal",
  "tumor",
  "tumor",
  "normal")
names(NT.ids) <- levels(Case2)
Case2 <- RenameIdents(Case2, NT.ids)
levels(Case2)
Case2@meta.data$NT_ident <- factor(Case2@active.ident)

Idents(Case2) <- "NT_ident"
levels(Case2)

NTall.markers <- FindMarkers(Case2, ident.1 = c("tumor"), ident.2 = c("normal"))
                             

top10 <-
  NTall.markers[order(NTall.markers$p_val_adj, decreasing = FALSE),]  %>%
  dplyr::filter(avg_log2FC > 0.3) %>%
  slice_head(n = 10) %>%
  ungroup() 
top10
rownames(top10)
[1] "CHGA"    "TTR"     "PCSK1"   "CHGB"    "TPH1"    "SLC18A1" "NELL2"   "PEG3"    "PTPRN2"  "KCTD12" 

tumor_genes <- c("TTR","CHGA","PCSK1","CHGB","TPH1","SLC18A1","NELL2","PEG3","PTPRN2","KCTD12" )

Case2 <- AddModuleScore(
  object =Case2,
  features = tumor_genes,
  name = 'tumor_genes'
)

Idents(Case2) <- "manual_ident2"
VlnPlot(Case2 , features = "tumor_genes1", log = FALSE, pt.size=0.1)+ NoLegend() 

Case2_tumor_normal <- subset(Case2,
                               tumor_genes1>2.5 & manual_ident2 =="tp#1m_&_#2m"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#1superficialsm_&_#2superficialsm"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#1sm_&_#1mp"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#2sm"|
                               tumor_genes1>2.5 & manual_ident2 =="tm#1_&_tp#1mp"|
                               tumor_genes1>2.5 & manual_ident2 =="tm#1pni"|
                               tumor_genes1>2.5 & manual_ident2 =="tm#2"|
                               tumor_genes1>2.5 & manual_ident2 =="tm#2pni"|
                               tumor_genes1<2 & manual_ident2 =="m_superficial_mid_crypt"|
                               tumor_genes1<2 & manual_ident2 =="m_base"|
                               tumor_genes1<2 & manual_ident2 =="sm"|
                               tumor_genes1<2 & manual_ident2 =="mp_inner"|
                               tumor_genes1<2 & manual_ident2 =="mp_outer"|
                               tumor_genes1<2 & manual_ident2 =="serosa"|
                               tumor_genes1<2 & manual_ident2 =="vascular_wall"|
                               tumor_genes1<2 & manual_ident2 =="mes_dense_collagen_stroma"|
                               tumor_genes1<2 & manual_ident2 =="mes_nerve")

Case2_tumor_normal <- subset(Case2_tumor_normal, ACTB>0 & B2M>0 &UBC>0)
save(Case2,file="Case2.RData")
save(Case2_tumor_normal,file="Case2_tumor_normal.RData")

VlnPlot(Case2_tumor_normal , features = "tumor_genes1", log = FALSE, pt.size=0.1)+ NoLegend()

SpatialDimPlot(Case2_tumor_normal, 
               label = FALSE, 
               label.size = 3, 
               pt.size.factor = 0.7, 
               alpha = 0.5,
               crop = FALSE)



Case2_tumor_normal <- RunUMAP(Case2_tumor_normal, dims = 1:10,n.neighbors = 250 , min.dist= 0.4)
UMAP_Case2_tumor_normal<-
  DimPlot(Case2_tumor_normal, reduction = "umap", label = TRUE) + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  NoLegend()
UMAP_Case2_tumor_normal
ggsave(filename="UMAP_Case2_tumor_normal.svg",width=5, height=5,plot=UMAP_Case2_tumor_normal)
UMAP_Case2_tumor_normal_w_legend<-
  DimPlot(Case2_tumor_normal, reduction = "umap", label = TRUE) + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
UMAP_Case2_tumor_normal_w_legend
ggsave(filename="UMAP_Case2_tumor_normal_w_legend.svg",width=5, height=5,plot=UMAP_Case2_tumor_normal_w_legend)











## deconvolution analysis

## https://jef.works/STdeconvolve/getting_started.html

library(tidyverse)
library(Seurat)
library(STdeconvolve)

cd <- Case2_tumor_normal@assays$Spatial@counts ## matrix of gene counts in each pixel
annot <- Case2_tumor_normal@meta.data$manual_ident3 ## annotated tissue layers assigned to each pixel

Case2_tumor_normal@images$slice1@coordinates

Idents(Case2_tumor_normal) <- "Area"
Case2_tumor_normal_Area3<-subset(Case2_tumor_normal, Area == "Area3")
Case2_tumor_normal_Area4<-subset(Case2_tumor_normal, Area == "Area4")

pos1<-Case2_tumor_normal_Area3@images$slice1@coordinates[,4:5]
names(pos1)<-c("x","y")
pos1


pos2<-Case2_tumor_normal_Area4@images$slice1.2@coordinates[,4:5]
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



ldas_Case2 <- fitLDA(corpus$corpus, Ks = seq(2, 9, by = 1),
                     perc.rare.thresh = 0.05,
                     plot=TRUE,
                     verbose=TRUE)
save(ldas_Case2,file="ldas_Case2.RData")

## select model with minimum perplexity
optLDA <- optimalModel(models = ldas_Case2, opt = "proportion")

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



deconpie_Case2_tumor_normal<-
  vizAllTopics(deconProp, pos, 
               topicCols = c("springgreen4","tomato","violet","turquoise1","steelblue2","violetred2","sienna2"),
               r=35,
               groups = annot,
               group_cols = rainbow(length(levels(annot))),)
deconpie_Case2_tumor_normal
ggsave(filename="deconpie_Case2_tumor_normal.svg",width=20, height=20,plot=deconpie_Case2_tumor_normal)

##run_me_QC(ldas_Case2,starting_k=2,ending_k=5,dir="output_Case2_tumor_normal_1/")

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

##run_me_results(opt=18,dir = "output_Case2_tumor_normal_1/",ldas=ldas_Case2)



dfdeconProp <- as.data.frame(deconProp)
colnames(dfdeconProp)<- c("stroma","tumor1","tumor2","smooth_muscle","mucosa","tumor3","tumor4")
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

dfdeconProp3 <- melt(dfdeconProp2,id = c("manual_ident3"),measure = c("tumor1","tumor2","tumor3","tumor4","mucosa","stroma","smooth_muscle"))
dfdeconProp3


deconbar_Case2_tumor_normal  <-
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
  scale_x_discrete(limits = c("m_superficial_mid_crypt",
                              "m_base",
                              "sm",
                              "mp_inner",
                              "mp_outer",
                              "serosa",
                              "vascular_wall",
                              "mes_dense_collagen_stroma",
                              "mes_nerve",
                              
                              "tp_mucosal",
                              "tp_deep",
                              "tm"))+
  scale_fill_manual(values=c("tomato","violet","violetred2","sienna2","steelblue2","springgreen4","turquoise1"))+
  theme(axis.text.x=element_text(angle=75, hjust=1,  size = rel(1)),
        legend.title=element_blank())
deconbar_Case2_tumor_normal

ggsave(filename="deconbar_Case2_tumor_normal.svg",width=6, height=4,plot=deconbar_Case2_tumor_normal)

