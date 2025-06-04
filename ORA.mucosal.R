
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(readxl)
library(writexl)


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
library(Seurat)
library(SingleR)

library(EnhancedVolcano)
library(enrichplot)

library(MAST)

##Case1
Idents(Case1_tumor_normal) = "manual_ident2"
Case1_tumor <- subset(Case1_tumor_normal,
                      manual_ident2 == "tp#1m_&_#2m"|
                        manual_ident2 =="tp#1sm_&_#1mp"|
                        manual_ident2 =="tp#2sm"|
                        manual_ident2 =="tp#1mp_&_#2mp"|
                        manual_ident2 =="tm_pni_A"|
                        manual_ident2 == "tm_pni_B"|
                        manual_ident2 == "tm")

Case1_tumor_tp_m.markers <- FindMarkers(Case1_tumor, ident.1 = c("tp#1m_&_#2m"),
                                         ident.2 = c("tp#1sm_&_#1mp","tp#2sm","tp#1mp_&_#2mp","tm_pni_A","tm_pni_B","tm"),
                                         min.pct = 0.0, logfc.threshold = 0.0, test.use="MAST")



Case1_tumor_tp_m.markers
save(Case1_tumor_tp_m.markers,file="Case1_tumor_tp_m.markers.RData")
load(file="Case1_tumor_tp_m.markers.RData")


write.csv(Case1_tumor_tp_m.markers,"Case1_tumor_tp_m.markers.csv")

## upregulated

Case1_tumor_tp_m.markers_FC0.3 <-filter(Case1_tumor_tp_m.markers, avg_log2FC > 0.3) 
Case1_tumor_tp_m.markers_FC0.3pval<- filter(Case1_tumor_tp_m.markers_FC0.3, p_val_adj < 10e-2) 

Case1_tumor_tp_m.markers_FC0.3pval
write.csv(Case1_tumor_tp_m.markers_FC0.3pval,"Case1_tumor_tp_m.markers_FC0.3pval.csv")


HKGs <- c("ACTB","B2M","UBC")
keyvals.shape <- ifelse(rownames(Case1_tumor_tp_m.markers) %in% HKGs, 2,1)
names(keyvals.shape)[keyvals.shape == 2] <- 'Housekeeping genes'
names(keyvals.shape)[keyvals.shape == 1] <- ''



volcano_Case1_tumor_tp_m.markers_FC0.3pval<-
  EnhancedVolcano(Case1_tumor_tp_m.markers,
                  lab = rownames(Case1_tumor_tp_m.markers),
                  drawConnectors = TRUE,max.overlaps=20,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = 'Case1_tumor_tp_m',
                  pCutoff = 10e-2,
                  FCcutoff = 0.3,
                  shapeCustom = keyvals.shape,
                  pointSize = c(ifelse(rownames(Case1_tumor_tp_m.markers) %in% HKGs, 8, 1)))

volcano_Case1_tumor_tp_m.markers_FC0.3pval
ggsave(filename="volcano_Case1_tumor_tp_m.markers_FC0.3pval.svg",width=7, height=7,plot=volcano_Case1_tumor_tp_m.markers_FC0.3pval)


## download "c5.all.v2023.2.Hs.symbols.gmt" from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/
C5all<- read.gmt("c5.all.v2023.2.Hs.symbols.gmt")

ego_Case1_tumor_tp_m_c5<- enricher(gene    = rownames(Case1_tumor_tp_m.markers_FC0.3pval),
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               qvalueCutoff = 1,
               minGSSize = 10,
               maxGSSize = 500,
               TERM2GENE = C5all)

save(ego_Case1_tumor_tp_m_c5,file="ego_Case1_tumor_tp_m_c5.RData")

write_xlsx(as.data.frame(ego_Case1_tumor_tp_m_c5),"ego_Case1_tumor_tp_m_c5.xlsx")
head(as.data.frame(ego_Case1_tumor_tp_m_c5))

## download "c2.cp.v2023.2.Hs.symbols.gmt" from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/
C2all<- read.gmt("c2.cp.v2023.2.Hs.symbols.gmt")

ego_Case1_tumor_tp_m_c2<- enricher(gene    = rownames(Case1_tumor_tp_m.markers_FC0.3pval),
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               qvalueCutoff = 1,
               minGSSize = 10,
               maxGSSize = 500,
               TERM2GENE = C2all)

save(ego_Case1_tumor_tp_m_c2,file="ego_Case1_tumor_tp_m_c2.RData")

write_xlsx(as.data.frame(ego_Case1_tumor_tp_m_c2),"ego_Case1_tumor_tp_m_c2.xlsx")
head(as.data.frame(ego_Case1_tumor_tp_m_c2))


## downregulated


Case1_tumor_tp_m.markers_FCminus0.3 <-filter(Case1_tumor_tp_m.markers, avg_log2FC < -0.3) 
Case1_tumor_tp_m.markers_FCminus0.3pval<- filter(Case1_tumor_tp_m.markers_FCminus0.3, p_val_adj < 10e-2) 

Case1_tumor_tp_m.markers_FCminus0.3pval
write.csv(Case1_tumor_tp_m.markers_FCminus0.3pval,"Case1_tumor_tp_m.markers_FCminus0.3pval.csv")






ego_Case1_tumor_tp_m_c5minus<- enricher(gene    = rownames(Case1_tumor_tp_m.markers_FCminus0.3pval),
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               qvalueCutoff = 1,
               minGSSize = 10,
               maxGSSize = 500,
               TERM2GENE = C5all)
save(ego_Case1_tumor_tp_m_c5minus,file="ego_Case1_tumor_tp_m_c5minus.RData")


write_xlsx(as.data.frame(ego_Case1_tumor_tp_m_c5minus),"ego_Case1_tumor_tp_m_c5minus.xlsx")
head(as.data.frame(ego_Case1_tumor_tp_m_c5minus))




ego_Case1_tumor_tp_m_c2minus<- enricher(gene    = rownames(Case1_tumor_tp_m.markers_FCminus0.3pval),
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               qvalueCutoff = 1,
               minGSSize = 10,
               maxGSSize = 500,
               TERM2GENE = C2all)
save(ego_Case1_tumor_tp_m_c2minus,file="ego_Case1_tumor_tp_m_c2minus.RData")

write_xlsx(as.data.frame(ego_Case1_tumor_tp_m_c2minus),"ego_Case1_tumor_tp_m_c2minus.xlsx")
head(as.data.frame(ego_Case1_tumor_tp_m_c2minus))



##Case2

Case2_tumor <- subset(Case2_tumor_normal,
                      manual_ident2 == "tp#1m_&_#2m"|
                        manual_ident2 =="tp#1superficialsm_&_#2superficialsm"|
                        manual_ident2 =="tp#1sm_&_#1mp"|
                        manual_ident2 =="tp#2sm"|
                        manual_ident2 =="tm#1_&_tp#1mp"|
                        manual_ident2 == "tm#1pni"|
                        manual_ident2 =="tm#2"|
                        manual_ident2 == "tm#2pni")

Idents(Case2_tumor) = "manual_ident2"
Case2_tumor_tp_m.markers <- FindMarkers(Case2_tumor, ident.1 = c("tp#1m_&_#2m"), 
                                         ident.2 = c("tm#1_&_tp#1mp","tm#1pni","tm#2","tm#2pni","tp#1superficialsm_&_#2superficialsm","tp#1sm_&_#1mp","tp#2sm"),
                                         min.pct = 0.0, logfc.threshold = 0.0,test.use="MAST")


save(Case2_tumor_tp_m.markers,file="Case2_tumor_tp_m.markers.RData")
load(file="Case2_tumor_tp_m.markers.RData")


write.csv(Case2_tumor_tp_m.markers,"Case2_tumor_tp_m.markers.csv")

## upregulated

Case2_tumor_tp_m.markers_FC0.3 <-filter(Case2_tumor_tp_m.markers, avg_log2FC > 0.3) 
Case2_tumor_tp_m.markers_FC0.3pval<- filter(Case2_tumor_tp_m.markers_FC0.3, p_val_adj < 10e-2) 

Case2_tumor_tp_m.markers_FC0.3pval
write.csv(Case2_tumor_tp_m.markers_FC0.3pval,"Case2_tumor_tp_m.markers_FC0.3pval.csv")

HKGs <- c("ACTB","B2M","UBC")
keyvals.shape <- ifelse(rownames(Case2_tumor_tp_m.markers) %in% HKGs, 2,1)
names(keyvals.shape)[keyvals.shape == 2] <- 'Housekeeping genes'
names(keyvals.shape)[keyvals.shape == 1] <- ''


volcano_Case2_tumor_tp_m.markers_FC0.3pval<-
  EnhancedVolcano(Case2_tumor_tp_m.markers,
                  lab = rownames(Case2_tumor_tp_m.markers),
                  drawConnectors = TRUE,max.overlaps=20,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = 'Case2_tumor_tp_m',
                  pCutoff = 10e-2,
                  FCcutoff = 0.3,
                  shapeCustom = keyvals.shape,
                  pointSize = c(ifelse(rownames(Case2_tumor_tp_m.markers) %in% HKGs, 8, 1)))

volcano_Case2_tumor_tp_m.markers_FC0.3pval
ggsave(filename="volcano_Case2_tumor_tp_m.markers_FC0.3pval.svg",width=7, height=7,plot=volcano_Case2_tumor_tp_m.markers_FC0.3pval)

VlnPlot(Case2_tumor, features = HKGs, ncol = 3, y.max = NULL)


ego_Case2_tumor_tp_m_c5<- enricher(gene    = rownames(Case2_tumor_tp_m.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C5all)

save(ego_Case2_tumor_tp_m_c5,file="ego_Case2_tumor_tp_m_c5.RData")

write_xlsx(as.data.frame(ego_Case2_tumor_tp_m_c5),"ego_Case2_tumor_tp_m_c5.xlsx")
head(as.data.frame(ego_Case2_tumor_tp_m_c5))


ego_Case2_tumor_tp_m_c2<- enricher(gene    = rownames(Case2_tumor_tp_m.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C2all)

save(ego_Case2_tumor_tp_m_c2,file="ego_Case2_tumor_tp_m_c2.RData")

write_xlsx(as.data.frame(ego_Case2_tumor_tp_m_c2),"ego_Case2_tumor_tp_m_c2.xlsx")
head(as.data.frame(ego_Case2_tumor_tp_m_c2))


## downregulated


Case2_tumor_tp_m.markers_FCminus0.3 <-filter(Case2_tumor_tp_m.markers, avg_log2FC < -0.3) 
Case2_tumor_tp_m.markers_FCminus0.3pval<- filter(Case2_tumor_tp_m.markers_FCminus0.3, p_val_adj < 10e-2) 

Case2_tumor_tp_m.markers_FCminus0.3pval
write.csv(Case2_tumor_tp_m.markers_FCminus0.3pval,"Case2_tumor_tp_m.markers_FCminus0.3pval.csv")






ego_Case2_tumor_tp_m_c5minus<- enricher(gene    = rownames(Case2_tumor_tp_m.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C5all)
save(ego_Case2_tumor_tp_m_c5minus,file="ego_Case2_tumor_tp_m_c5minus.RData")


write_xlsx(as.data.frame(ego_Case2_tumor_tp_m_c5minus),"ego_Case2_tumor_tp_m_c5minus.xlsx")
head(as.data.frame(ego_Case2_tumor_tp_m_c5minus))




ego_Case2_tumor_tp_m_c2minus<- enricher(gene    = rownames(Case2_tumor_tp_m.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C2all)
save(ego_Case2_tumor_tp_m_c2minus,file="ego_Case2_tumor_tp_m_c2minus.RData")

write_xlsx(as.data.frame(ego_Case2_tumor_tp_m_c2minus),"ego_Case2_tumor_tp_m_c2minus.xlsx")
head(as.data.frame(ego_Case2_tumor_tp_m_c2minus))









##Case3

Case3_tumor <- subset(Case3_tumor_normal,
                      manual_ident2 == "tp#1m"|
                        manual_ident2 =="tp#2m"|
                        manual_ident2 =="tp#1sm"|
                        manual_ident2 =="tp#2sm_&_thrombus"|
                        manual_ident2 =="tp#1mp"|
                        manual_ident2 =="tm_tp#2sm_&_perithrombus"|
                        manual_ident2 =="tl")
Idents(Case3_tumor) = "manual_ident2"


Case3_tumor_tp_m.markers <- FindMarkers(Case3_tumor, ident.1 = c("tp#1m", "tp#2m"),
                                         ident.2 = c("tp#2sm_&_thrombus","tm_tp#2sm_&_perithrombus","tp#1sm","tp#1mp","tl"),
                                         min.pct = 0.0, logfc.threshold = 0.0 ,test.use="MAST")


Case3_tumor_tp_m.markers
save(Case3_tumor_tp_m.markers,file="Case3_tumor_tp_m.markers.RData")
load(file="Case3_tumor_tp_m.markers.RData")


write.csv(Case3_tumor_tp_m.markers,"Case3_tumor_tp_m.markers.csv")

## upregulated

Case3_tumor_tp_m.markers_FC0.3 <-filter(Case3_tumor_tp_m.markers, avg_log2FC > 0.3) 
Case3_tumor_tp_m.markers_FC0.3pval<- filter(Case3_tumor_tp_m.markers_FC0.3, p_val_adj < 10e-2) 

Case3_tumor_tp_m.markers_FC0.3pval
write.csv(Case3_tumor_tp_m.markers_FC0.3pval,"Case3_tumor_tp_m.markers_FC0.3pval.csv")


HKGs <- c("ACTB","B2M","UBC")
keyvals.shape <- ifelse(rownames(Case3_tumor_tp_m.markers) %in% HKGs, 2,1)
names(keyvals.shape)[keyvals.shape == 2] <- 'Housekeeping genes'
names(keyvals.shape)[keyvals.shape == 1] <- ''

volcano_Case3_tumor_tp_m.markers_FC0.3pval<-
  EnhancedVolcano(Case3_tumor_tp_m.markers,
                  lab = rownames(Case3_tumor_tp_m.markers),
                  drawConnectors = TRUE,max.overlaps=20,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = 'Case3_tumor_tp_m',
                  pCutoff = 10e-2,
                  FCcutoff = 0.3,
                  shapeCustom = keyvals.shape,
                  pointSize = c(ifelse(rownames(Case3_tumor_tp_m.markers) %in% HKGs, 8, 1)))


volcano_Case3_tumor_tp_m.markers_FC0.3pval
ggsave(filename="volcano_Case3_tumor_tp_m.markers_FC0.3pval.svg",width=7, height=7,plot=volcano_Case3_tumor_tp_m.markers_FC0.3pval)




ego_Case3_tumor_tp_m_c5<- enricher(gene    = rownames(Case3_tumor_tp_m.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C5all)

save(ego_Case3_tumor_tp_m_c5,file="ego_Case3_tumor_tp_m_c5.RData")

write_xlsx(as.data.frame(ego_Case3_tumor_tp_m_c5),"ego_Case3_tumor_tp_m_c5.xlsx")
head(as.data.frame(ego_Case3_tumor_tp_m_c5))


ego_Case3_tumor_tp_m_c2<- enricher(gene    = rownames(Case3_tumor_tp_m.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C2all)

save(ego_Case3_tumor_tp_m_c2,file="ego_Case3_tumor_tp_m_c2.RData")

write_xlsx(as.data.frame(ego_Case3_tumor_tp_m_c2),"ego_Case3_tumor_tp_m_c2.xlsx")
head(as.data.frame(ego_Case3_tumor_tp_m_c2))


## downregulated


Case3_tumor_tp_m.markers_FCminus0.3 <-filter(Case3_tumor_tp_m.markers, avg_log2FC < -0.3) 
Case3_tumor_tp_m.markers_FCminus0.3pval<- filter(Case3_tumor_tp_m.markers_FCminus0.3, p_val_adj < 10e-2) 

Case3_tumor_tp_m.markers_FCminus0.3pval
write.csv(Case3_tumor_tp_m.markers_FCminus0.3pval,"Case3_tumor_tp_m.markers_FCminus0.3pval.csv")






ego_Case3_tumor_tp_m_c5minus<- enricher(gene    = rownames(Case3_tumor_tp_m.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C5all)
save(ego_Case3_tumor_tp_m_c5minus,file="ego_Case3_tumor_tp_m_c5minus.RData")


write_xlsx(as.data.frame(ego_Case3_tumor_tp_m_c5minus),"ego_Case3_tumor_tp_m_c5minus.xlsx")
head(as.data.frame(ego_Case3_tumor_tp_m_c5minus))




ego_Case3_tumor_tp_m_c2minus<- enricher(gene    = rownames(Case3_tumor_tp_m.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C2all)
save(ego_Case3_tumor_tp_m_c2minus,file="ego_Case3_tumor_tp_m_c2minus.RData")

write_xlsx(as.data.frame(ego_Case3_tumor_tp_m_c2minus),"ego_Case3_tumor_tp_m_c2minus.xlsx")
head(as.data.frame(ego_Case3_tumor_tp_m_c2minus))














##Case4
Case4_tumor <- subset(Case4_tumor_normal,
                      manual_ident2=="tp#1m_#2m_&_#3m"|
                        manual_ident2=="tp#3m_&#3sm"|
                        manual_ident2=="tp#1sm"|
                        manual_ident2=="tp#2sm"|
                        manual_ident2=="tp#3sm_superficial"|
                        manual_ident2=="tp#3sm_deep"|
                        manual_ident2=="tp#2mp"|
                        manual_ident2=="tp#1_sm_deep_&_#2subserosa"|
                        manual_ident2=="tl")
Idents(Case4_tumor) = "manual_ident2"
Case4_tumor_tp_m.markers <- FindMarkers(Case4_tumor, ident.1 = c("tp#1m_#2m_&_#3m"),
                                         ident.2 = c("tp#3m_&#3sm","tp#1sm","tp#2sm","tp#3sm_superficial","tp#3sm_deep","tp#2mp","tl"),
                                         min.pct = 0.0, logfc.threshold = 0.0,test.use="MAST")


save(Case4_tumor_tp_m.markers,file="Case4_tumor_tp_m.markers.RData")
load(file="Case4_tumor_tp_m.markers.RData")


write.csv(Case4_tumor_tp_m.markers,"Case4_tumor_tp_m.markers.csv")

## upregulated

Case4_tumor_tp_m.markers_FC0.3 <-filter(Case4_tumor_tp_m.markers, avg_log2FC > 0.3) 
Case4_tumor_tp_m.markers_FC0.3pval<- filter(Case4_tumor_tp_m.markers_FC0.3, p_val_adj < 10e-2) 

Case4_tumor_tp_m.markers_FC0.3pval
write.csv(Case4_tumor_tp_m.markers_FC0.3pval,"Case4_tumor_tp_m.markers_FC0.3pval.csv")

HKGs <- c("ACTB","B2M","UBC")
keyvals.shape <- ifelse(rownames(Case4_tumor_tp_m.markers) %in% HKGs, 2,1)
names(keyvals.shape)[keyvals.shape == 2] <- 'Housekeeping genes'
names(keyvals.shape)[keyvals.shape == 1] <- ''


volcano_Case4_tumor_tp_m.markers_FC0.3pval<-
  EnhancedVolcano(Case4_tumor_tp_m.markers,
                  lab = rownames(Case4_tumor_tp_m.markers),
                  drawConnectors = TRUE,max.overlaps=20,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = 'Case4_tumor_tp_m',
                  pCutoff = 10e-2,
                  FCcutoff = 0.3,
                  shapeCustom = keyvals.shape,
                  pointSize = c(ifelse(rownames(Case4_tumor_tp_m.markers) %in% HKGs, 8, 1)))


volcano_Case4_tumor_tp_m.markers_FC0.3pval
ggsave(filename="volcano_Case4_tumor_tp_m.markers_FC0.3pval.svg",width=7, height=7,plot=volcano_Case4_tumor_tp_m.markers_FC0.3pval)



ego_Case4_tumor_tp_m_c5<- enricher(gene    = rownames(Case4_tumor_tp_m.markers_FC0.3pval),
                                 pvalueCutoff = 1,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 1,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 TERM2GENE = C5all)

save(ego_Case4_tumor_tp_m_c5,file="ego_Case4_tumor_tp_m_c5.RData")

write_xlsx(as.data.frame(ego_Case4_tumor_tp_m_c5),"ego_Case4_tumor_tp_m_c5.xlsx")
head(as.data.frame(ego_Case4_tumor_tp_m_c5))


ego_Case4_tumor_tp_m_c2<- enricher(gene    = rownames(Case4_tumor_tp_m.markers_FC0.3pval),
                                 pvalueCutoff = 1,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 1,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 TERM2GENE = C2all)

save(ego_Case4_tumor_tp_m_c2,file="ego_Case4_tumor_tp_m_c2.RData")

write_xlsx(as.data.frame(ego_Case4_tumor_tp_m_c2),"ego_Case4_tumor_tp_m_c2.xlsx")
head(as.data.frame(ego_Case4_tumor_tp_m_c2))


## downregulated


Case4_tumor_tp_m.markers_FCminus0.3 <-filter(Case4_tumor_tp_m.markers, avg_log2FC < -0.3) 
Case4_tumor_tp_m.markers_FCminus0.3pval<- filter(Case4_tumor_tp_m.markers_FCminus0.3, p_val_adj < 10e-2) 

Case4_tumor_tp_m.markers_FCminus0.3pval
write.csv(Case4_tumor_tp_m.markers_FCminus0.3pval,"Case4_tumor_tp_m.markers_FCminus0.3pval.csv")






ego_Case4_tumor_tp_m_c5minus<- enricher(gene    = rownames(Case4_tumor_tp_m.markers_FCminus0.3pval),
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 1,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      TERM2GENE = C5all)
save(ego_Case4_tumor_tp_m_c5minus,file="ego_Case4_tumor_tp_m_c5minus.RData")


write_xlsx(as.data.frame(ego_Case4_tumor_tp_m_c5minus),"ego_Case4_tumor_tp_m_c5minus.xlsx")
head(as.data.frame(ego_Case4_tumor_tp_m_c5minus))




ego_Case4_tumor_tp_m_c2minus<- enricher(gene    = rownames(Case4_tumor_tp_m.markers_FCminus0.3pval),
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 1,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      TERM2GENE = C2all)
save(ego_Case4_tumor_tp_m_c2minus,file="ego_Case4_tumor_tp_m_c2minus.RData")

write_xlsx(as.data.frame(ego_Case4_tumor_tp_m_c2minus),"ego_Case4_tumor_tp_m_c2minus.xlsx")
head(as.data.frame(ego_Case4_tumor_tp_m_c2minus))



## looking for shared gene sets



## pathway
ego_Case1_tumor_tp_m_c2_0.1<-
ego_Case1_tumor_tp_m_c2 %>%
  filter(ego_Case1_tumor_tp_m_c2@result$p.adjust<0.1 & ego_Case1_tumor_tp_m_c2@result$qvalue <0.1)
ego_Case2_tumor_tp_m_c2_0.1<-
  ego_Case2_tumor_tp_m_c2 %>%
  filter(ego_Case2_tumor_tp_m_c2@result$p.adjust<0.1 & ego_Case2_tumor_tp_m_c2@result$qvalue <0.1)
ego_Case3_tumor_tp_m_c2_0.1<-
  ego_Case3_tumor_tp_m_c2 %>%
  filter(ego_Case3_tumor_tp_m_c2@result$p.adjust<0.1 & ego_Case3_tumor_tp_m_c2@result$qvalue <0.1)
ego_Case4_tumor_tp_m_c2_0.1<-
  ego_Case4_tumor_tp_m_c2 %>%
  filter(ego_Case4_tumor_tp_m_c2@result$p.adjust<0.1 & ego_Case4_tumor_tp_m_c2@result$qvalue <0.1)

Reduce(intersect, list(ego_Case1_tumor_tp_m_c2_0.1@result$ID,
                       ego_Case2_tumor_tp_m_c2_0.1@result$ID,
                       ego_Case3_tumor_tp_m_c2_0.1@result$ID,
                       ego_Case4_tumor_tp_m_c2_0.1@result$ID))
character(0)

Reduce(intersect, list(
                       ego_Case2_tumor_tp_m_c2_0.1@result$ID,
                       ego_Case3_tumor_tp_m_c2_0.1@result$ID,
                       ego_Case4_tumor_tp_m_c2_0.1@result$ID))
character(0)

Reduce(intersect, list(ego_Case1_tumor_tp_m_c2_0.1@result$ID,
                       
                       ego_Case3_tumor_tp_m_c2_0.1@result$ID,
                       ego_Case4_tumor_tp_m_c2_0.1@result$ID))
character(0)

ego_shared_tumor_tp_m_c2_0.1<-
Reduce(intersect, list(ego_Case1_tumor_tp_m_c2_0.1@result$ID,
                       ego_Case2_tumor_tp_m_c2_0.1@result$ID,
                      
                       ego_Case4_tumor_tp_m_c2_0.1@result$ID))
save(ego_shared_tumor_tp_m_c2_0.1,file="ego_shared_tumor_tp_m_c2_0.1.RData")

[1] "REACTOME_METABOLISM_OF_FAT_SOLUBLE_VITAMINS" "KEGG_PARKINSONS_DISEASE"   

Reduce(intersect, list(ego_Case1_tumor_tp_m_c2_0.1@result$ID,
                       ego_Case2_tumor_tp_m_c2_0.1@result$ID,
                       ego_Case3_tumor_tp_m_c2_0.1@result$ID))
character(0)



## pathway minus

ego_Case1_tumor_tp_m_c2minus_0.1<-
  ego_Case1_tumor_tp_m_c2minus %>%
  filter(ego_Case1_tumor_tp_m_c2minus@result$p.adjust<0.1 & ego_Case1_tumor_tp_m_c2minus@result$qvalue <0.1)
ego_Case2_tumor_tp_m_c2minus_0.1<-
  ego_Case2_tumor_tp_m_c2minus %>%
  filter(ego_Case2_tumor_tp_m_c2minus@result$p.adjust<0.1 & ego_Case2_tumor_tp_m_c2minus@result$qvalue <0.1)
ego_Case3_tumor_tp_m_c2minus_0.1<-
  ego_Case3_tumor_tp_m_c2minus %>%
  filter(ego_Case3_tumor_tp_m_c2minus@result$p.adjust<0.1 & ego_Case3_tumor_tp_m_c2minus@result$qvalue <0.1)
ego_Case4_tumor_tp_m_c2minus_0.1<-
  ego_Case4_tumor_tp_m_c2minus %>%
  filter(ego_Case4_tumor_tp_m_c2minus@result$p.adjust<0.1 & ego_Case4_tumor_tp_m_c2minus@result$qvalue <0.1)

ego_shared_tumor_tp_m_c2minus_0.1<-
Reduce(intersect, list(ego_Case1_tumor_tp_m_c2minus_0.1@result$ID,
                       ego_Case2_tumor_tp_m_c2minus_0.1@result$ID,
                       ego_Case3_tumor_tp_m_c2minus_0.1@result$ID,
                       ego_Case4_tumor_tp_m_c2minus_0.1@result$ID))
save(ego_shared_tumor_tp_m_c2minus_0.1,file="ego_shared_tumor_tp_m_c2minus_0.1.RData")
[1] "REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA"    "WP_VEGFA_VEGFR2_SIGNALING"                                  
[3] "REACTOME_MRNA_SPLICING"                                      "REACTOME_RHO_GTPASE_CYCLE"                                  
[5] "REACTOME_SIGNALING_BY_NTRKS"                                 "REACTOME_TRANS_GOLGI_NETWORK_VESICLE_BUDDING"               
[7] "REACTOME_GOLGI_ASSOCIATED_VESICLE_BIOGENESIS"                "REACTOME_TRANSPORT_TO_THE_GOLGI_AND_SUBSEQUENT_MODIFICATION"
[9] "REACTOME_RAB_GEFS_EXCHANGE_GTP_FOR_GDP_ON_RABS"              "REACTOME_SIGNALING_BY_VEGF"                                 
[11] "REACTOME_ER_TO_GOLGI_ANTEROGRADE_TRANSPORT"      




## Biological Process


ego_Case1_tumor_tp_m_c5_0.1<-
  ego_Case1_tumor_tp_m_c5 %>%
  filter(ego_Case1_tumor_tp_m_c5@result$p.adjust<0.1 & ego_Case1_tumor_tp_m_c5@result$qvalue <0.1)
ego_Case2_tumor_tp_m_c5_0.1<-
  ego_Case2_tumor_tp_m_c5 %>%
  filter(ego_Case2_tumor_tp_m_c5@result$p.adjust<0.1 & ego_Case2_tumor_tp_m_c5@result$qvalue <0.1)
ego_Case3_tumor_tp_m_c5_0.1<-
  ego_Case3_tumor_tp_m_c5 %>%
  filter(ego_Case3_tumor_tp_m_c5@result$p.adjust<0.1 & ego_Case3_tumor_tp_m_c5@result$qvalue <0.1)
ego_Case4_tumor_tp_m_c5_0.1<-
  ego_Case4_tumor_tp_m_c5 %>%
  filter(ego_Case4_tumor_tp_m_c5@result$p.adjust<0.1 & ego_Case4_tumor_tp_m_c5@result$qvalue <0.1)

ego_shared_tumor_tp_m_c5_0.1<-
Reduce(intersect, list(ego_Case1_tumor_tp_m_c5_0.1@result$ID,
                       ego_Case2_tumor_tp_m_c5_0.1@result$ID,
                       ego_Case3_tumor_tp_m_c5_0.1@result$ID,
                       ego_Case4_tumor_tp_m_c5_0.1@result$ID))
save(ego_shared_tumor_tp_m_c5_0.1,file="ego_shared_tumor_tp_m_c5_0.1.RData")

[1] "GOBP_ACTIN_FILAMENT_ORGANIZATION"                            "GOBP_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION"       
[3] "GOBP_REGULATION_OF_ACTIN_FILAMENT_ORGANIZATION"              "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY"                 
[5] "GOBP_REGULATION_OF_PROTEIN_CONTAINING_COMPLEX_ASSEMBLY"      "GOBP_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"
[7] "GOMF_CADHERIN_BINDING"                                       "GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ORGANELLE"    
[9] "GOBP_AEROBIC_RESPIRATION"                                    "GOBP_OXIDATIVE_PHOSPHORYLATION"                             
[11] "GOBP_NEGATIVE_REGULATION_OF_PROTEIN_MODIFICATION_PROCESS"    "HP_ABNORMAL_ACTIVITY_OF_MITOCHONDRIAL_RESPIRATORY_CHAIN"    
[13] "GOBP_PROTEIN_FOLDING"   


## Biological Process minus

ego_Case1_tumor_tp_m_c5minus_0.01<-
  ego_Case1_tumor_tp_m_c5minus %>%
  filter(ego_Case1_tumor_tp_m_c5minus@result$p.adjust<0.01 & ego_Case1_tumor_tp_m_c5minus@result$qvalue <0.01)
ego_Case2_tumor_tp_m_c5minus_0.01<-
  ego_Case2_tumor_tp_m_c5minus %>%
  filter(ego_Case2_tumor_tp_m_c5minus@result$p.adjust<0.01 & ego_Case2_tumor_tp_m_c5minus@result$qvalue <0.01)
ego_Case3_tumor_tp_m_c5minus_0.01<-
  ego_Case3_tumor_tp_m_c5minus %>%
  filter(ego_Case3_tumor_tp_m_c5minus@result$p.adjust<0.01 & ego_Case3_tumor_tp_m_c5minus@result$qvalue <0.01)
ego_Case4_tumor_tp_m_c5minus_0.01<-
  ego_Case4_tumor_tp_m_c5minus %>%
  filter(ego_Case4_tumor_tp_m_c5minus@result$p.adjust<0.01 & ego_Case4_tumor_tp_m_c5minus@result$qvalue <0.01)

ego_shared_tumor_tp_m_c5minus_0.01<-
Reduce(intersect, list(ego_Case1_tumor_tp_m_c5minus_0.01@result$ID,
                       ego_Case2_tumor_tp_m_c5minus_0.01@result$ID,
                       ego_Case3_tumor_tp_m_c5minus_0.01@result$ID,
                       ego_Case4_tumor_tp_m_c5minus_0.01@result$ID))
save(ego_shared_tumor_tp_m_c5minus_0.01,file="ego_shared_tumor_tp_m_c5minus_0.01.RData")

[1] "GOCC_TRANSPORT_VESICLE"                                  "GOBP_VESICLE_LOCALIZATION"                              
[3] "GOBP_RNA_SPLICING"                                       "HP_INTERICTAL_EEG_ABNORMALITY"                          
[5] "GOCC_EXOCYTIC_VESICLE"                                   "GOCC_NUCLEAR_SPECK"                                     
[7] "HP_EEG_WITH_GENERALIZED_EPILEPTIFORM_DISCHARGES"         "GOBP_REGULATION_OF_PROTEIN_STABILITY"                   
[9] "HP_UNSTEADY_GAIT"                                        "GOBP_POSITIVE_REGULATION_OF_PROTEIN_LOCALIZATION"       
[11] "HP_MOTOR_SEIZURE"                                        "GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ORGANELLE"
[13] "GOBP_VESICLE_ORGANIZATION"                               "GOCC_VACUOLAR_MEMBRANE"                                 
[15] "GOCC_SPLICEOSOMAL_COMPLEX"                               "HP_RESTRICTED_OR_REPETITIVE_BEHAVIORS_OR_INTERESTS"     
[17] "GOBP_GOLGI_VESICLE_TRANSPORT"                            "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"          
[19] "GOBP_ENDOSOMAL_TRANSPORT" 


