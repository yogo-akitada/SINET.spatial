
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
Area7 <-Load10X_Spatial(
  "data/356C1",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",filter.matrix = TRUE
)

plot1 <- VlnPlot(Area7, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Area7, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

Area7@meta.data$Area <- "Area7"
table(Area7$Area)

Area8 <-Load10X_Spatial(
  "data/356D1",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",filter.matrix = TRUE
)

plot1 <- VlnPlot(Area8, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Area8, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

Area8@meta.data$Area <- "Area8"
table(Area8$Area)

Case4 <- merge(Area7, y = Area8, add.cell.ids = c("Area7", "Area8"))
table(Case4$Area)



Case4 <- NormalizeData(Case4, normalization.method = "LogNormalize", scale.factor = 10000)
Case4 <- FindVariableFeatures(Case4, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(Case4), points = head(VariableFeatures(Case4), 25), repel = TRUE)
Case4 <- ScaleData(Case4, features = rownames(Case4))
Case4 <- RunPCA(Case4, features = VariableFeatures(object = Case4))
VizDimLoadings(Case4, dims = 1:2, reduction = "pca", balanced = T,
               nfeatures = 50)
DimPlot(Case4, reduction = "pca")

Case4 <- FindNeighbors(Case4, dims = 1:20)
Case4 <- FindClusters(Case4, resolution = 1.2)

table(Case4$seurat_clusters)

Case4 <- RunUMAP(Case4, dims = 1:10,n.neighbors = 200 , min.dist= 0.5)
p1 <- DimPlot(Case4, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(Case4, 
                     label = TRUE, 
                     label.size = 3, 
                     pt.size.factor = 0.7, 
                     alpha = 0.5,
                     crop = FALSE)
test <-p1 + p2
test


Idents(Case4) <- "seurat_clusters"


manual.ident2 <- c(                     
  "tl",
  "sm",
  "tp#1m_#2m_&_#3m",
  "tp#3sm_superficial",
  "tp#2sm",
  "unspecified_potentially_lymphatic_cell",
  "tp#1sm",
  "tp#3sm_deep",
  "m_mid_crypt",
  "tp#2mp",
  "tp#2mp_&_myenteric plexus",
  "tp#1_sm_deep_&_#2subserosa",
  "mp_outer_non_peritumoral",
  "m_base",
  "mp_inner_non_peritumoral",
  "tp#3m_&#3sm",
  "mp_peritumoral")


names(manual.ident2) <- levels(Case4)
Case4 <- RenameIdents(Case4, manual.ident2)
levels(Case4)
Case4@meta.data$manual_ident2 <- factor(Case4@active.ident)

Idents(Case4) <- "manual_ident2"
levels(Case4)


levels(Case4) <- c(  
  "m_mid_crypt",
  "m_base",
  "sm",
  "mp_inner_non_peritumoral",
  "mp_outer_non_peritumoral",
  "mp_peritumoral",
  "unspecified_potentially_lymphatic_cell",
  
  "tp#1m_#2m_&_#3m",
  "tp#3m_&#3sm",
  "tp#1sm",
  "tp#2sm",
  "tp#3sm_superficial",
  "tp#3sm_deep",
  "tp#2mp",
  "tp#2mp_&_myenteric plexus",
  "tp#1_sm_deep_&_#2subserosa",
  "tl")
levels(Case4)
Case4@meta.data$manual_ident2 <- factor(Case4@active.ident)


p1 <- DimPlot(Case4, reduction = "umap", label = TRUE,label.size = 2.5) 
p1
p2 <- SpatialDimPlot(Case4,
                     label = FALSE, 
                     label.size = 3, 
                     pt.size.factor = 1.5, 
                     alpha = 0.5,
                     crop = FALSE,information = NULL)+ NoLegend()
p3 <- SpatialDimPlot(Case4,
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
ggsave(filename="Case4_rename.svg",width=10, height=15,plot=renameplot)



manual.ident3 <- c(  
  "m_mid_crypt",
  "m_base",
  "sm",
  "mp_inner_non_peritumoral",
  "mp_outer_non_peritumoral",
  "mp_peritumoral",
  "unspecified_potentially_lymphatic_cell",
  
  "tp_mucosal",
  "tp_deep",
  "tp_deep",
  "tp_deep",
  "tp_deep",
  "tp_deep",
  "tp_deep",
  "tp_deep",
  "tp_deep",
  "tl")

names(manual.ident3) <- levels(Case4)
Case4 <- RenameIdents(Case4, manual.ident3)
levels(Case4)
Case4@meta.data$manual_ident3 <- factor(Case4@active.ident)

Idents(Case4) <- "manual_ident3"
levels(Case4)



### quality checking by tumor gene module scores
Idents(Case4) <- "manual_ident2"
VlnPlot(Case4 , features = "nCount_Spatial", log = FALSE, pt.size=0.1)+ NoLegend()

## exclude unspecified_potentially_lymphatic_cell from tumor or normal
Idents(Case4) <- "seurat_clusters"

NT.ids <- c(
  "tumor",
  "normal",
  "tumor",
  "tumor",
  "tumor",
  "unspecified_potentially_lymphatic_cell",
  "tumor",
  "tumor",
  "normal",
  "tumor",
  "tumor",
  "tumor",
  "normal",
  "normal",
  "normal",
  "tumor",
  "normal")
names(NT.ids) <- levels(Case4)
Case4 <- RenameIdents(Case4, NT.ids)
levels(Case4)
Case4@meta.data$NT_ident <- factor(Case4@active.ident)

Idents(Case4) <- "NT_ident"
levels(Case4)

NTall.markers <- FindMarkers(Case4, ident.1 = c("tumor"), ident.2 = c("normal"))
                             

top10 <-
  NTall.markers[order(NTall.markers$p_val_adj, decreasing = FALSE),]  %>%
  dplyr::filter(avg_log2FC > 0.3) %>%
  slice_head(n = 10) %>%
  ungroup() 
top10
rownames(top10)
[1] "KCTD12"  "CHGA"    "TTR"     "CHGB"    "PCSK1"   "NEUROD1" "SLC18A1" "TPH1"    "TAC1"    "RAB3B"  

tumor_genes <- c("TTR","KCTD12","CHGA","CHGB","PCSK1","NEUROD1","SLC18A1","TPH1","TAC1","RAB3B")

Case4 <- AddModuleScore(
  object =Case4,
  features = tumor_genes,
  name = 'tumor_genes'
)

Idents(Case4) <- "manual_ident2"
VlnPlot(Case4 , features = "tumor_genes1", log = FALSE, pt.size=0.1)+ NoLegend() 

Case4_tumor_normal <- subset(Case4,
                               tumor_genes1>2.5 & manual_ident2 =="tp#1m_#2m_&_#3m"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#3m_&#3sm"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#1sm"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#2sm"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#3sm_superficial"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#3sm_deep"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#2mp"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#2mp_&_myenteric plexus"|
                               tumor_genes1>2.5 & manual_ident2 =="tp#1_sm_deep_&_#2subserosa"|
                               tumor_genes1>2.5 & manual_ident2 =="tl"|
                               tumor_genes1<2 & manual_ident2 =="m_mid_crypt"|
                               tumor_genes1<2 & manual_ident2 =="m_base"|
                               tumor_genes1<2 & manual_ident2 =="sm"|
                               tumor_genes1<2 & manual_ident2 =="mp_inner_non_peritumoral"|
                               tumor_genes1<2 & manual_ident2 =="mp_outer_non_peritumoral"|
                               tumor_genes1<2 & manual_ident2 =="mp_peritumoral")


Case4_tumor_normal <- subset(Case4_tumor_normal, ACTB>0 & B2M>0 &UBC>0)
save(Case4,file="Case4.RData")
save(Case4_tumor_normal,file="Case4_tumor_normal.RData")

VlnPlot(Case4_tumor_normal , features = "tumor_genes1", log = FALSE, pt.size=0.1)+ NoLegend()

SpatialDimPlot(Case4_tumor_normal, 
               label = FALSE, 
               label.size = 3, 
               pt.size.factor = 0.7, 
               alpha = 0.5,
               crop = FALSE)



Case4_tumor_normal <- RunUMAP(Case4_tumor_normal, dims = 1:10,n.neighbors = 250 , min.dist= 0.5)
UMAP_Case4_tumor_normal<-
  DimPlot(Case4_tumor_normal, reduction = "umap", label = TRUE) + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  NoLegend()
UMAP_Case4_tumor_normal
ggsave(filename="UMAP_Case4_tumor_normal.svg",width=5, height=5,plot=UMAP_Case4_tumor_normal)
UMAP_Case4_tumor_normal_w_legend<-
  DimPlot(Case4_tumor_normal, reduction = "umap", label = TRUE) + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
UMAP_Case4_tumor_normal_w_legend
ggsave(filename="UMAP_Case4_tumor_normal_w_legend.svg",width=5, height=5,plot=UMAP_Case4_tumor_normal_w_legend)











## deconvolution analysis

## https://jef.works/STdeconvolve/getting_started.html

library(tidyverse)
library(Seurat)
library(STdeconvolve)

cd <- Case4_tumor_normal@assays$Spatial@counts ## matrix of gene counts in each pixel
annot <- Case4_tumor_normal@meta.data$manual_ident3 ## annotated tissue layers assigned to each pixel

Case4_tumor_normal@images$slice1@coordinates

Idents(Case4_tumor_normal) <- "Area"
Case4_tumor_normal_Area7<-subset(Case4_tumor_normal, Area == "Area7")
Case4_tumor_normal_Area8<-subset(Case4_tumor_normal, Area == "Area8")

pos1<-Case4_tumor_normal_Area7@images$slice1@coordinates[,4:5]
names(pos1)<-c("x","y")
pos1


pos2<-Case4_tumor_normal_Area8@images$slice1.2@coordinates[,4:5]
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



ldas_Case4 <- fitLDA(corpus$corpus, Ks = seq(2, 9, by = 1),
                     perc.rare.thresh = 0.05,
                     plot=TRUE,
                     verbose=TRUE)
save(ldas_Case4,file="ldas_Case4.RData")

## select model with minimum perplexity
optLDA <- optimalModel(models = ldas_Case4, opt = "proportion")

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


deconpie_Case4_tumor_normal<-
  vizAllTopics(deconProp, pos, 
               topicCols = c("springgreen4","tomato","turquoise1","springgreen","steelblue2"),
               r=35,
               groups = annot,
               group_cols = rainbow(length(levels(annot))),)
deconpie_Case4_tumor_normal
ggsave(filename="deconpie_Case4_tumor_normal.svg",width=20, height=20,plot=deconpie_Case4_tumor_normal)

##run_me_QC(ldas_Case4,starting_k=2,ending_k=5,dir="output_Case4_tumor_normal_1/")

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

##run_me_results(opt=18,dir = "output_Case4_tumor_normal_1/",ldas=ldas_Case4)



dfdeconProp <- as.data.frame(deconProp)
colnames(dfdeconProp)<- c("m_base","tumor","sm","smooth_muscle","m_mid_crypt")
dfdeconProp$manual_ident3 <- factor(annot)
dfdeconProp 


dfdeconProp2<-
  dfdeconProp  %>%
  dplyr::group_by(manual_ident3) 


dfdeconProp2

library(reshape2)

dfdeconProp3 <- melt(dfdeconProp2,id = c("manual_ident3"),measure = c("tumor","m_mid_crypt","m_base","sm","smooth_muscle"))
dfdeconProp3


deconbar_Case4_tumor_normal  <-
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
  scale_x_discrete(limits = c("m_mid_crypt",
                              "m_base",
                              "sm",
                              "mp_inner_non_peritumoral",
                              "mp_outer_non_peritumoral",
                              "mp_peritumoral",
                              
                              "tp_mucosal",
                              "tp_deep",
                              "tl"))+
  scale_fill_manual(values=c("tomato","steelblue2","springgreen4","turquoise1","springgreen"))+
  theme(axis.text.x=element_text(angle=75, hjust=1,  size = rel(1)),
        legend.title=element_blank())
deconbar_Case4_tumor_normal

ggsave(filename="spatial/Figs/deconbar_Case4_tumor_normal.svg",width=6, height=4,plot=deconbar_Case4_tumor_normal)
