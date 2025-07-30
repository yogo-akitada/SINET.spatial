install.packages("harmony")
library(harmony)

case.combined <- merge(Case1_tumor_normal, y = c(Case2_tumor_normal,Case3_tumor_normal,Case4_tumor_normal), add.cell.ids = c("Patient1", "Patient2","Patient3","Patient4"))
head(colnames(case.combined))
tail(colnames(case.combined))
table(case.combined$orig.ident)

case.combined$orig.ident <- ifelse(grepl("Patient1", colnames(case.combined)), "Patient1", 
                                   ifelse(grepl("Patient2", colnames(case.combined)), "Patient2", 
                                          ifelse(grepl("Patient3", colnames(case.combined)), "Patient3", 
                                                 ifelse(grepl("Patient4", colnames(case.combined)), "Patient4", "NA"))))

table(case.combined$orig.ident)

case.combined <- NormalizeData(case.combined, normalization.method = "LogNormalize", scale.factor = 10000)
case.combined <- FindVariableFeatures(case.combined, selection.method = "vst", nfeatures = 2000)
case.combined <- ScaleData(case.combined, features = rownames(case.combined))
case.combined <- RunPCA(case.combined, features = VariableFeatures(object = case.combined))
case.combined <- RunUMAP(case.combined, dims = 1:10,n.neighbors = 200 , min.dist= 0.5)

DimPlot(case.combined, reduction = "umap", label = TRUE, group.by = "orig.ident")
DefaultAssay(case.combined) <- "Spatial"


case.combined <- RunHarmony(case.combined, group.by.vars = "orig.ident")
case.combined <- RunUMAP(case.combined, reduction = "harmony", dims = 1:30)
case.combined <- FindNeighbors(case.combined, reduction = "harmony", dims = 1:30) %>% FindClusters()

UMAP_4patients <-
DimPlot(case.combined, group.by = "orig.ident")
ggsave(filename="UMAP_4patients.svg",width=7, height=5,plot=UMAP_4patients)


Idents(case.combined) <-"manual_ident3"
          
UMAP_4patients_4subtype <-
  DimPlot(case.combined, 
          cells.highlight = list(
            mucosal = WhichCells(case.combined,  idents = "tp_mucosal"),
            deep = WhichCells(case.combined,  idents = "tp_deep"),
            mesenteric = WhichCells(case.combined, idents = "tm"),
            lymphatic = WhichCells(case.combined,  idents = "tl")),
            cols.highlight = c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3" ), cols= "grey"
          )
UMAP_4patients_4subtype
ggsave(filename="UMAP_4patients_4subtype.svg",width=7, height=5,plot=UMAP_4patients_4subtype)
