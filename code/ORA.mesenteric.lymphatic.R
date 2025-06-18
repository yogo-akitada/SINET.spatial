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

Case1_tumor <- subset(Case1_tumor_normal,
                      manual_ident2 == "tp#1m_&_#2m"|
                        manual_ident2 =="tp#1sm_&_#1mp"|
                        manual_ident2 =="tp#2sm"|
                        manual_ident2 =="tp#1mp_&_#2mp"|
                        manual_ident2 =="tm_pni_A"|
                        manual_ident2 == "tm_pni_B"|
                        manual_ident2 == "tm")

Idents(Case1_tumor) = "manual_ident2"

Case1_tumor_tm.markers <- FindMarkers(Case1_tumor, ident.1 = c("tm_pni_A","tm_pni_B","tm"),
                                         ident.2 = c("tp#1sm_&_#1mp","tp#2sm","tp#1mp_&_#2mp","tp#1m_&_#2m"),
                                         min.pct = 0.0, logfc.threshold = 0.0, test.use="MAST")



Case1_tumor_tm.markers
save(Case1_tumor_tm.markers,file="Case1_tumor_tm.markers.RData")
load(file="Case1_tumor_tm.markers.RData")


write.csv(Case1_tumor_tm.markers,"Case1_tumor_tm.markers.csv")

## upregulated

Case1_tumor_tm.markers_FC0.3 <-filter(Case1_tumor_tm.markers, avg_log2FC > 0.3) 
Case1_tumor_tm.markers_FC0.3pval<- filter(Case1_tumor_tm.markers_FC0.3, p_val_adj < 10e-2) 

Case1_tumor_tm.markers_FC0.3pval
write.csv(Case1_tumor_tm.markers_FC0.3pval,"Case1_tumor_tm.markers_FC0.3pval.csv")


HKGs <- c("ACTB","B2M","UBC")
keyvals.shape <- ifelse(rownames(Case1_tumor_tm.markers) %in% HKGs, 2,1)
names(keyvals.shape)[keyvals.shape == 2] <- 'Housekeeping genes'
names(keyvals.shape)[keyvals.shape == 1] <- ''



volcano_Case1_tumor_tm.markers_FC0.3pval<-
  EnhancedVolcano(Case1_tumor_tm.markers,
                  lab = rownames(Case1_tumor_tm.markers),
                  drawConnectors = TRUE,max.overlaps=20,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = 'Case1_tumor_tm',
                  pCutoff = 10e-2,
                  FCcutoff = 0.3,
                  shapeCustom = keyvals.shape,
                  pointSize = c(ifelse(rownames(Case1_tumor_tm.markers) %in% HKGs, 8, 1)))

volcano_Case1_tumor_tm.markers_FC0.3pval
ggsave(filename="volcano_Case1_tumor_tm.markers_FC0.3pval.svg",width=7, height=7,plot=volcano_Case1_tumor_tm.markers_FC0.3pval)




ego_Case1_tumor_tm_c5<- enricher(gene    = rownames(Case1_tumor_tm.markers_FC0.3pval),
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               qvalueCutoff = 1,
               minGSSize = 10,
               maxGSSize = 500,
               TERM2GENE = C5all)

save(ego_Case1_tumor_tm_c5,file="ego_Case1_tumor_tm_c5.RData")

write_xlsx(as.data.frame(ego_Case1_tumor_tm_c5),"ego_Case1_tumor_tm_c5.xlsx")
head(as.data.frame(ego_Case1_tumor_tm_c5))


ego_Case1_tumor_tm_c2<- enricher(gene    = rownames(Case1_tumor_tm.markers_FC0.3pval),
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               qvalueCutoff = 1,
               minGSSize = 10,
               maxGSSize = 500,
               TERM2GENE = C2all)

save(ego_Case1_tumor_tm_c2,file="ego_Case1_tumor_tm_c2.RData")

write_xlsx(as.data.frame(ego_Case1_tumor_tm_c2),"ego_Case1_tumor_tm_c2.xlsx")
head(as.data.frame(ego_Case1_tumor_tm_c2))


## downregulated


Case1_tumor_tm.markers_FCminus0.3 <-filter(Case1_tumor_tm.markers, avg_log2FC < -0.3) 
Case1_tumor_tm.markers_FCminus0.3pval<- filter(Case1_tumor_tm.markers_FCminus0.3, p_val_adj < 10e-2) 

Case1_tumor_tm.markers_FCminus0.3pval
write.csv(Case1_tumor_tm.markers_FCminus0.3pval,"Case1_tumor_tm.markers_FCminus0.3pval.csv")






ego_Case1_tumor_tm_c5minus<- enricher(gene    = rownames(Case1_tumor_tm.markers_FCminus0.3pval),
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               qvalueCutoff = 1,
               minGSSize = 10,
               maxGSSize = 500,
               TERM2GENE = C5all)
save(ego_Case1_tumor_tm_c5minus,file="ego_Case1_tumor_tm_c5minus.RData")


write_xlsx(as.data.frame(ego_Case1_tumor_tm_c5minus),"ego_Case1_tumor_tm_c5minus.xlsx")
head(as.data.frame(ego_Case1_tumor_tm_c5minus))




ego_Case1_tumor_tm_c2minus<- enricher(gene    = rownames(Case1_tumor_tm.markers_FCminus0.3pval),
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               qvalueCutoff = 1,
               minGSSize = 10,
               maxGSSize = 500,
               TERM2GENE = C2all)
save(ego_Case1_tumor_tm_c2minus,file="ego_Case1_tumor_tm_c2minus.RData")

write_xlsx(as.data.frame(ego_Case1_tumor_tm_c2minus),"ego_Case1_tumor_tm_c2minus.xlsx")
head(as.data.frame(ego_Case1_tumor_tm_c2minus))



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
Case2_tumor_tm.markers <- FindMarkers(Case2_tumor, ident.1 = c("tm#1_&_tp#1mp","tm#1pni","tm#2","tm#2pni"), 
                                         ident.2 = c("tp#1m_&_#2m","tp#1superficialsm_&_#2superficialsm","tp#1sm_&_#1mp","tp#2sm"),
                                         min.pct = 0.0, logfc.threshold = 0.0,test.use="MAST")


save(Case2_tumor_tm.markers,file="Case2_tumor_tm.markers.RData")
load(file="Case2_tumor_tm.markers.RData")


write.csv(Case2_tumor_tm.markers,"Case2_tumor_tm.markers.csv")

## upregulated

Case2_tumor_tm.markers_FC0.3 <-filter(Case2_tumor_tm.markers, avg_log2FC > 0.3) 
Case2_tumor_tm.markers_FC0.3pval<- filter(Case2_tumor_tm.markers_FC0.3, p_val_adj < 10e-2) 

Case2_tumor_tm.markers_FC0.3pval
write.csv(Case2_tumor_tm.markers_FC0.3pval,"Case2_tumor_tm.markers_FC0.3pval.csv")

HKGs <- c("ACTB","B2M","UBC")
keyvals.shape <- ifelse(rownames(Case2_tumor_tm.markers) %in% HKGs, 2,1)
names(keyvals.shape)[keyvals.shape == 2] <- 'Housekeeping genes'
names(keyvals.shape)[keyvals.shape == 1] <- ''


volcano_Case2_tumor_tm.markers_FC0.3pval<-
  EnhancedVolcano(Case2_tumor_tm.markers,
                  lab = rownames(Case2_tumor_tm.markers),
                  drawConnectors = TRUE,max.overlaps=20,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = 'Case2_tumor_tm',
                  pCutoff = 10e-2,
                  FCcutoff = 0.3,
                  shapeCustom = keyvals.shape,
                  pointSize = c(ifelse(rownames(Case2_tumor_tm.markers) %in% HKGs, 8, 1)))

volcano_Case2_tumor_tm.markers_FC0.3pval
ggsave(filename="volcano_Case2_tumor_tm.markers_FC0.3pval.svg",width=7, height=7,plot=volcano_Case2_tumor_tm.markers_FC0.3pval)

VlnPlot(Case2_tumor, features = HKGs, ncol = 3, y.max = NULL)


ego_Case2_tumor_tm_c5<- enricher(gene    = rownames(Case2_tumor_tm.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C5all)

save(ego_Case2_tumor_tm_c5,file="ego_Case2_tumor_tm_c5.RData")

write_xlsx(as.data.frame(ego_Case2_tumor_tm_c5),"ego_Case2_tumor_tm_c5.xlsx")
head(as.data.frame(ego_Case2_tumor_tm_c5))


ego_Case2_tumor_tm_c2<- enricher(gene    = rownames(Case2_tumor_tm.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C2all)

save(ego_Case2_tumor_tm_c2,file="ego_Case2_tumor_tm_c2.RData")

write_xlsx(as.data.frame(ego_Case2_tumor_tm_c2),"ego_Case2_tumor_tm_c2.xlsx")
head(as.data.frame(ego_Case2_tumor_tm_c2))


## downregulated


Case2_tumor_tm.markers_FCminus0.3 <-filter(Case2_tumor_tm.markers, avg_log2FC < -0.3) 
Case2_tumor_tm.markers_FCminus0.3pval<- filter(Case2_tumor_tm.markers_FCminus0.3, p_val_adj < 10e-2) 

Case2_tumor_tm.markers_FCminus0.3pval
write.csv(Case2_tumor_tm.markers_FCminus0.3pval,"Case2_tumor_tm.markers_FCminus0.3pval.csv")






ego_Case2_tumor_tm_c5minus<- enricher(gene    = rownames(Case2_tumor_tm.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C5all)
save(ego_Case2_tumor_tm_c5minus,file="ego_Case2_tumor_tm_c5minus.RData")


write_xlsx(as.data.frame(ego_Case2_tumor_tm_c5minus),"ego_Case2_tumor_tm_c5minus.xlsx")
head(as.data.frame(ego_Case2_tumor_tm_c5minus))




ego_Case2_tumor_tm_c2minus<- enricher(gene    = rownames(Case2_tumor_tm.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C2all)
save(ego_Case2_tumor_tm_c2minus,file="ego_Case2_tumor_tm_c2minus.RData")

write_xlsx(as.data.frame(ego_Case2_tumor_tm_c2minus),"ego_Case2_tumor_tm_c2minus.xlsx")
head(as.data.frame(ego_Case2_tumor_tm_c2minus))









##Case3

Case3_tumor <- subset(Case3_tumor_normal,
                      manual_ident2 == "tp#1m_&_#2m"|
                        manual_ident2 =="tp#2m"|
                        manual_ident2 =="tp#1sm"|
                        manual_ident2 =="tp#2sm_&_thrombus"|
                        manual_ident2 =="tp#1mp"|
                        manual_ident2 =="tm_tp#2sm_&_perithrombus"|
                        manual_ident2 =="tl")
Idents(Case3_tumor) = "manual_ident2"

#tm
Case3_tumor_tm.markers <- FindMarkers(Case3_tumor, ident.1 = c("tp#2sm_&_thrombus","tm_tp#2sm_&_perithrombus"),
                                         ident.2 = c("tp#1m_&_#2m", "tp#2m","tp#1sm","tp#1mp","tl"),
                                         min.pct = 0.0, logfc.threshold = 0.0 ,test.use="MAST")


Case3_tumor_tm.markers
save(Case3_tumor_tm.markers,file="Case3_tumor_tm.markers.RData")
load(file="Case3_tumor_tm.markers.RData")


write.csv(Case3_tumor_tm.markers,"Case3_tumor_tm.markers.csv")

## upregulated

Case3_tumor_tm.markers_FC0.3 <-filter(Case3_tumor_tm.markers, avg_log2FC > 0.3) 
Case3_tumor_tm.markers_FC0.3pval<- filter(Case3_tumor_tm.markers_FC0.3, p_val_adj < 10e-2) 

Case3_tumor_tm.markers_FC0.3pval
write.csv(Case3_tumor_tm.markers_FC0.3pval,"Case3_tumor_tm.markers_FC0.3pval.csv")


HKGs <- c("ACTB","B2M","UBC")
keyvals.shape <- ifelse(rownames(Case3_tumor_tm.markers) %in% HKGs, 2,1)
names(keyvals.shape)[keyvals.shape == 2] <- 'Housekeeping genes'
names(keyvals.shape)[keyvals.shape == 1] <- ''

volcano_Case3_tumor_tm.markers_FC0.3pval<-
  EnhancedVolcano(Case3_tumor_tm.markers,
                  lab = rownames(Case3_tumor_tm.markers),
                  drawConnectors = TRUE,max.overlaps=20,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = 'Case3_tumor_tm',
                  pCutoff = 10e-2,
                  FCcutoff = 0.3,
                  shapeCustom = keyvals.shape,
                  pointSize = c(ifelse(rownames(Case3_tumor_tm.markers) %in% HKGs, 8, 1)))


volcano_Case3_tumor_tm.markers_FC0.3pval
ggsave(filename="volcano_Case3_tumor_tm.markers_FC0.3pval.svg",width=7, height=7,plot=volcano_Case3_tumor_tm.markers_FC0.3pval)




ego_Case3_tumor_tm_c5<- enricher(gene    = rownames(Case3_tumor_tm.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C5all)

save(ego_Case3_tumor_tm_c5,file="ego_Case3_tumor_tm_c5.RData")

write_xlsx(as.data.frame(ego_Case3_tumor_tm_c5),"ego_Case3_tumor_tm_c5.xlsx")
head(as.data.frame(ego_Case3_tumor_tm_c5))


ego_Case3_tumor_tm_c2<- enricher(gene    = rownames(Case3_tumor_tm.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C2all)

save(ego_Case3_tumor_tm_c2,file="ego_Case3_tumor_tm_c2.RData")

write_xlsx(as.data.frame(ego_Case3_tumor_tm_c2),"ego_Case3_tumor_tm_c2.xlsx")
head(as.data.frame(ego_Case3_tumor_tm_c2))


## downregulated


Case3_tumor_tm.markers_FCminus0.3 <-filter(Case3_tumor_tm.markers, avg_log2FC < -0.3) 
Case3_tumor_tm.markers_FCminus0.3pval<- filter(Case3_tumor_tm.markers_FCminus0.3, p_val_adj < 10e-2) 

Case3_tumor_tm.markers_FCminus0.3pval
write.csv(Case3_tumor_tm.markers_FCminus0.3pval,"Case3_tumor_tm.markers_FCminus0.3pval.csv")






ego_Case3_tumor_tm_c5minus<- enricher(gene    = rownames(Case3_tumor_tm.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C5all)
save(ego_Case3_tumor_tm_c5minus,file="ego_Case3_tumor_tm_c5minus.RData")


write_xlsx(as.data.frame(ego_Case3_tumor_tm_c5minus),"ego_Case3_tumor_tm_c5minus.xlsx")
head(as.data.frame(ego_Case3_tumor_tm_c5minus))




ego_Case3_tumor_tm_c2minus<- enricher(gene    = rownames(Case3_tumor_tm.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C2all)
save(ego_Case3_tumor_tm_c2minus,file="ego_Case3_tumor_tm_c2minus.RData")

write_xlsx(as.data.frame(ego_Case3_tumor_tm_c2minus),"ego_Case3_tumor_tm_c2minus.xlsx")
head(as.data.frame(ego_Case3_tumor_tm_c2minus))




#tl
Case3_tumor_tl.markers <- FindMarkers(Case3_tumor, ident.1 = c("tl"),
                                      ident.2 = c("tp#1m_&_#2m", "tp#2m","tp#1sm","tp#1mp","tp#2sm_&_thrombus","tm_tp#2sm_&_perithrombus"),
                                      min.pct = 0.0, logfc.threshold = 0.0 ,test.use="MAST")


Case3_tumor_tl.markers
save(Case3_tumor_tl.markers,file="Case3_tumor_tl.markers.RData")
load(file="Case3_tumor_tl.markers.RData")


write.csv(Case3_tumor_tl.markers,"Case3_tumor_tl.markers.csv")

## upregulated

Case3_tumor_tl.markers_FC0.3 <-filter(Case3_tumor_tl.markers, avg_log2FC > 0.3) 
Case3_tumor_tl.markers_FC0.3pval<- filter(Case3_tumor_tl.markers_FC0.3, p_val_adj < 10e-2) 

Case3_tumor_tl.markers_FC0.3pval
write.csv(Case3_tumor_tl.markers_FC0.3pval,"Case3_tumor_tl.markers_FC0.3pval.csv")


HKGs <- c("ACTB","B2M","UBC")
keyvals.shape <- ifelse(rownames(Case3_tumor_tl.markers) %in% HKGs, 2,1)
names(keyvals.shape)[keyvals.shape == 2] <- 'Housekeeping genes'
names(keyvals.shape)[keyvals.shape == 1] <- ''

volcano_Case3_tumor_tl.markers_FC0.3pval<-
  EnhancedVolcano(Case3_tumor_tl.markers,
                  lab = rownames(Case3_tumor_tl.markers),
                  drawConnectors = TRUE,max.overlaps=20,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = 'Case3_tumor_tl',
                  pCutoff = 10e-2,
                  FCcutoff = 0.3,
                  shapeCustom = keyvals.shape,
                  pointSize = c(ifelse(rownames(Case3_tumor_tl.markers) %in% HKGs, 8, 1)))


volcano_Case3_tumor_tl.markers_FC0.3pval
ggsave(filename="volcano_Case3_tumor_tl.markers_FC0.3pval.svg",width=7, height=7,plot=volcano_Case3_tumor_tl.markers_FC0.3pval)




ego_Case3_tumor_tl_c5<- enricher(gene    = rownames(Case3_tumor_tl.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C5all)

save(ego_Case3_tumor_tl_c5,file="ego_Case3_tumor_tl_c5.RData")

write_xlsx(as.data.frame(ego_Case3_tumor_tl_c5),"ego_Case3_tumor_tl_c5.xlsx")
head(as.data.frame(ego_Case3_tumor_tl_c5))


ego_Case3_tumor_tl_c2<- enricher(gene    = rownames(Case3_tumor_tl.markers_FC0.3pval),
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   TERM2GENE = C2all)

save(ego_Case3_tumor_tl_c2,file="ego_Case3_tumor_tl_c2.RData")

write_xlsx(as.data.frame(ego_Case3_tumor_tl_c2),"ego_Case3_tumor_tl_c2.xlsx")
head(as.data.frame(ego_Case3_tumor_tl_c2))


## downregulated


Case3_tumor_tl.markers_FCminus0.3 <-filter(Case3_tumor_tl.markers, avg_log2FC < -0.3) 
Case3_tumor_tl.markers_FCminus0.3pval<- filter(Case3_tumor_tl.markers_FCminus0.3, p_val_adj < 10e-2) 

Case3_tumor_tl.markers_FCminus0.3pval
write.csv(Case3_tumor_tl.markers_FCminus0.3pval,"Case3_tumor_tl.markers_FCminus0.3pval.csv")






ego_Case3_tumor_tl_c5minus<- enricher(gene    = rownames(Case3_tumor_tl.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C5all)
save(ego_Case3_tumor_tl_c5minus,file="ego_Case3_tumor_tl_c5minus.RData")


write_xlsx(as.data.frame(ego_Case3_tumor_tl_c5minus),"ego_Case3_tumor_tl_c5minus.xlsx")
head(as.data.frame(ego_Case3_tumor_tl_c5minus))




ego_Case3_tumor_tl_c2minus<- enricher(gene    = rownames(Case3_tumor_tl.markers_FCminus0.3pval),
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = C2all)
save(ego_Case3_tumor_tl_c2minus,file="ego_Case3_tumor_tl_c2minus.RData")

write_xlsx(as.data.frame(ego_Case3_tumor_tl_c2minus),"ego_Case3_tumor_tl_c2minus.xlsx")
head(as.data.frame(ego_Case3_tumor_tl_c2minus))














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
Case4_tumor_tl.markers <- FindMarkers(Case4_tumor, ident.1 = c("tl"),
                                         ident.2 = c("tp#3m_&#3sm","tp#1sm","tp#2sm","tp#3sm_superficial","tp#3sm_deep","tp#2mp","tp#1m_#2m_&_#3m"),
                                         min.pct = 0.0, logfc.threshold = 0.0,test.use="MAST")


save(Case4_tumor_tl.markers,file="Case4_tumor_tl.markers.RData")
load(file="Case4_tumor_tl.markers.RData")


write.csv(Case4_tumor_tl.markers,"Case4_tumor_tl.markers.csv")

## upregulated

Case4_tumor_tl.markers_FC0.3 <-filter(Case4_tumor_tl.markers, avg_log2FC > 0.3) 
Case4_tumor_tl.markers_FC0.3pval<- filter(Case4_tumor_tl.markers_FC0.3, p_val_adj < 10e-2) 

Case4_tumor_tl.markers_FC0.3pval
write.csv(Case4_tumor_tl.markers_FC0.3pval,"Case4_tumor_tl.markers_FC0.3pval.csv")

HKGs <- c("ACTB","B2M","UBC")
keyvals.shape <- ifelse(rownames(Case4_tumor_tl.markers) %in% HKGs, 2,1)
names(keyvals.shape)[keyvals.shape == 2] <- 'Housekeeping genes'
names(keyvals.shape)[keyvals.shape == 1] <- ''


volcano_Case4_tumor_tl.markers_FC0.3pval<-
  EnhancedVolcano(Case4_tumor_tl.markers,
                  lab = rownames(Case4_tumor_tl.markers),
                  drawConnectors = TRUE,max.overlaps=20,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = 'Case4_tumor_tl',
                  pCutoff = 10e-2,
                  FCcutoff = 0.3,
                  shapeCustom = keyvals.shape,
                  pointSize = c(ifelse(rownames(Case4_tumor_tl.markers) %in% HKGs, 8, 1)))


volcano_Case4_tumor_tl.markers_FC0.3pval
ggsave(filename="volcano_Case4_tumor_tl.markers_FC0.3pval.svg",width=7, height=7,plot=volcano_Case4_tumor_tl.markers_FC0.3pval)



ego_Case4_tumor_tl_c5<- enricher(gene    = rownames(Case4_tumor_tl.markers_FC0.3pval),
                                 pvalueCutoff = 1,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 1,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 TERM2GENE = C5all)

save(ego_Case4_tumor_tl_c5,file="ego_Case4_tumor_tl_c5.RData")

write_xlsx(as.data.frame(ego_Case4_tumor_tl_c5),"ego_Case4_tumor_tl_c5.xlsx")
head(as.data.frame(ego_Case4_tumor_tl_c5))


ego_Case4_tumor_tl_c2<- enricher(gene    = rownames(Case4_tumor_tl.markers_FC0.3pval),
                                 pvalueCutoff = 1,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 1,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 TERM2GENE = C2all)

save(ego_Case4_tumor_tl_c2,file="ego_Case4_tumor_tl_c2.RData")

write_xlsx(as.data.frame(ego_Case4_tumor_tl_c2),"ego_Case4_tumor_tl_c2.xlsx")
head(as.data.frame(ego_Case4_tumor_tl_c2))


## downregulated


Case4_tumor_tl.markers_FCminus0.3 <-filter(Case4_tumor_tl.markers, avg_log2FC < -0.3) 
Case4_tumor_tl.markers_FCminus0.3pval<- filter(Case4_tumor_tl.markers_FCminus0.3, p_val_adj < 10e-2) 

Case4_tumor_tl.markers_FCminus0.3pval
write.csv(Case4_tumor_tl.markers_FCminus0.3pval,"Case4_tumor_tl.markers_FCminus0.3pval.csv")






ego_Case4_tumor_tl_c5minus<- enricher(gene    = rownames(Case4_tumor_tl.markers_FCminus0.3pval),
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 1,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      TERM2GENE = C5all)
save(ego_Case4_tumor_tl_c5minus,file="ego_Case4_tumor_tl_c5minus.RData")


write_xlsx(as.data.frame(ego_Case4_tumor_tl_c5minus),"ego_Case4_tumor_tl_c5minus.xlsx")
head(as.data.frame(ego_Case4_tumor_tl_c5minus))




ego_Case4_tumor_tl_c2minus<- enricher(gene    = rownames(Case4_tumor_tl.markers_FCminus0.3pval),
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 1,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      TERM2GENE = C2all)
save(ego_Case4_tumor_tl_c2minus,file="ego_Case4_tumor_tl_c2minus.RData")

write_xlsx(as.data.frame(ego_Case4_tumor_tl_c2minus),"ego_Case4_tumor_tl_c2minus.xlsx")
head(as.data.frame(ego_Case4_tumor_tl_c2minus))



## looking for shared gene sets


#tm
## pathway
ego_Case1_tumor_tm_c2_0.1<-
ego_Case1_tumor_tm_c2 %>%
  filter(ego_Case1_tumor_tm_c2@result$p.adjust<0.1 & ego_Case1_tumor_tm_c2@result$qvalue <0.1)
ego_Case2_tumor_tm_c2_0.1<-
  ego_Case2_tumor_tm_c2 %>%
  filter(ego_Case2_tumor_tm_c2@result$p.adjust<0.1 & ego_Case2_tumor_tm_c2@result$qvalue <0.1)
ego_Case3_tumor_tm_c2_0.1<-
  ego_Case3_tumor_tm_c2 %>%
  filter(ego_Case3_tumor_tm_c2@result$p.adjust<0.1 & ego_Case3_tumor_tm_c2@result$qvalue <0.1)


Reduce(intersect, list(ego_Case1_tumor_tm_c2_0.1@result$ID,
                       ego_Case2_tumor_tm_c2_0.1@result$ID,
                       ego_Case3_tumor_tm_c2_0.1@result$ID))
character(0)

ego_shared_tumor_tm_c2_0.1 <-
Reduce(intersect, list(ego_Case1_tumor_tm_c2_0.1@result$ID,
                       
                       ego_Case3_tumor_tm_c2_0.1@result$ID))

save(ego_shared_tumor_tm_c2_0.1,file="ego_shared_tumor_tm_c2_0.1.RData")


[1] "REACTOME_NGF_STIMULATED_TRANSCRIPTION"                                        
[2] "REACTOME_NUCLEAR_EVENTS_KINASE_AND_TRANSCRIPTION_FACTOR_ACTIVATION"           
[3] "PID_AP1_PATHWAY"                                                              
[4] "WP_ROLE_OF_HYPOXIA_ANGIOGENESIS_AND_FGF_PATHWAY_IN_OA_CHONDROCYTE_HYPERTROPHY"
[5] "WP_SPINAL_CORD_INJURY"  


## pathway minus

ego_Case1_tumor_tm_c2minus_0.1<-
  ego_Case1_tumor_tm_c2minus %>%
  filter(ego_Case1_tumor_tm_c2minus@result$p.adjust<0.1 & ego_Case1_tumor_tm_c2minus@result$qvalue <0.1)
ego_Case2_tumor_tm_c2minus_0.1<-
  ego_Case2_tumor_tm_c2minus %>%
  filter(ego_Case2_tumor_tm_c2minus@result$p.adjust<0.1 & ego_Case2_tumor_tm_c2minus@result$qvalue <0.1)
ego_Case3_tumor_tm_c2minus_0.1<-
  ego_Case3_tumor_tm_c2minus %>%
  filter(ego_Case3_tumor_tm_c2minus@result$p.adjust<0.1 & ego_Case3_tumor_tm_c2minus@result$qvalue <0.1)

ego_shared_tumor_tm_c2minus_0.1<-
Reduce(intersect, list(ego_Case1_tumor_tm_c2minus_0.1@result$ID,
                       ego_Case2_tumor_tm_c2minus_0.1@result$ID,
                       ego_Case3_tumor_tm_c2minus_0.1@result$ID))
save(ego_shared_tumor_tm_c2minus_0.1,file="ego_shared_tumor_tm_c2minus_0.1.RData")

[1] "REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS"

Reduce(intersect, list(ego_Case1_tumor_tm_c2minus_0.1@result$ID,
                       ego_Case3_tumor_tm_c2minus_0.1@result$ID))


## Biological Process


ego_Case1_tumor_tm_c5_0.1<-
  ego_Case1_tumor_tm_c5 %>%
  filter(ego_Case1_tumor_tm_c5@result$p.adjust<0.1 & ego_Case1_tumor_tm_c5@result$qvalue <0.1)
ego_Case2_tumor_tm_c5_0.1<-
  ego_Case2_tumor_tm_c5 %>%
  filter(ego_Case2_tumor_tm_c5@result$p.adjust<0.1 & ego_Case2_tumor_tm_c5@result$qvalue <0.1)
ego_Case3_tumor_tm_c5_0.1<-
  ego_Case3_tumor_tm_c5 %>%
  filter(ego_Case3_tumor_tm_c5@result$p.adjust<0.1 & ego_Case3_tumor_tm_c5@result$qvalue <0.1)

ego_shared_tumor_tm_c5_0.1<-
Reduce(intersect, list(ego_Case1_tumor_tm_c5_0.1@result$ID,
                       ego_Case2_tumor_tm_c5_0.1@result$ID,
                       ego_Case3_tumor_tm_c5_0.1@result$ID))
save(ego_shared_tumor_tm_c5_0.1,file="ego_shared_tumor_tm_c5_0.1.RData")

[1] "GOBP_SYNAPSE_ORGANIZATION"


## Biological Process minus

ego_Case1_tumor_tm_c5minus_0.1<-
  ego_Case1_tumor_tm_c5minus %>%
  filter(ego_Case1_tumor_tm_c5minus@result$p.adjust<0.1 & ego_Case1_tumor_tm_c5minus@result$qvalue <0.1)
ego_Case2_tumor_tm_c5minus_0.1<-
  ego_Case2_tumor_tm_c5minus %>%
  filter(ego_Case2_tumor_tm_c5minus@result$p.adjust<0.1 & ego_Case2_tumor_tm_c5minus@result$qvalue <0.1)
ego_Case3_tumor_tm_c5minus_0.1<-
  ego_Case3_tumor_tm_c5minus %>%
  filter(ego_Case3_tumor_tm_c5minus@result$p.adjust<0.1 & ego_Case3_tumor_tm_c5minus@result$qvalue <0.1)

ego_shared_tumor_tm_c5minus_0.1<-
Reduce(intersect, list(ego_Case1_tumor_tm_c5minus_0.1@result$ID,
                       ego_Case2_tumor_tm_c5minus_0.1@result$ID,
                       ego_Case3_tumor_tm_c5minus_0.1@result$ID))
save(ego_shared_tumor_tm_c5minus_0.1,file="ego_shared_tumor_tm_c5minus_0.1.RData")

[1] "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY"




#tl
## pathway

ego_Case4_tumor_tl_c2_0.1<-
  ego_Case4_tumor_tl_c2 %>%
  filter(ego_Case4_tumor_tl_c2@result$p.adjust<0.1 & ego_Case4_tumor_tl_c2@result$qvalue <0.1)
ego_Case3_tumor_tl_c2_0.1<-
  ego_Case3_tumor_tl_c2 %>%
  filter(ego_Case3_tumor_tl_c2@result$p.adjust<0.1 & ego_Case3_tumor_tl_c2@result$qvalue <0.1)

Reduce(intersect, list(
                       ego_Case4_tumor_tl_c2_0.1@result$ID,
                       ego_Case3_tumor_tl_c2_0.1@result$ID))
character(0)



## pathway minus

ego_Case4_tumor_tl_c2minus_0.1<-
  ego_Case4_tumor_tl_c2minus %>%
  filter(ego_Case4_tumor_tl_c2minus@result$p.adjust<0.1 & ego_Case4_tumor_tl_c2minus@result$qvalue <0.1)
ego_Case3_tumor_tl_c2minus_0.1<-
  ego_Case3_tumor_tl_c2minus %>%
  filter(ego_Case3_tumor_tl_c2minus@result$p.adjust<0.1 & ego_Case3_tumor_tl_c2minus@result$qvalue <0.1)

ego_shared_tumor_tl_c2minus_0.1 <-
Reduce(intersect, list(
  ego_Case4_tumor_tl_c2minus_0.1@result$ID,
  ego_Case3_tumor_tl_c2minus_0.1@result$ID))
save(ego_shared_tumor_tl_c2minus_0.1,file="ego_shared_tumor_tl_c2minus_0.1.RData")

[1] "REACTOME_SMOOTH_MUSCLE_CONTRACTION"                                                                                              
[2] "REACTOME_SCAVENGING_BY_CLASS_A_RECEPTORS"                                                                                        
[3] "BIOCARTA_GHRELIN_PATHWAY"                                                                                                        
[4] "REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS"
[5] "WP_STRIATED_MUSCLE_CONTRACTION_PATHWAY"  


## Biological Process


ego_Case4_tumor_tl_c5_0.1<-
  ego_Case4_tumor_tl_c5 %>%
  filter(ego_Case4_tumor_tl_c5@result$p.adjust<0.1 & ego_Case4_tumor_tl_c5@result$qvalue <0.1)
ego_Case3_tumor_tl_c5_0.1<-
  ego_Case3_tumor_tl_c5 %>%
  filter(ego_Case3_tumor_tl_c5@result$p.adjust<0.1 & ego_Case3_tumor_tl_c5@result$qvalue <0.1)

Reduce(intersect, list(
  ego_Case4_tumor_tl_c5_0.1@result$ID,
  ego_Case3_tumor_tl_c5_0.1@result$ID))
character(0)

## Biological Process minus


ego_Case3_tumor_tl_c5minus_0.05<-
  ego_Case3_tumor_tl_c5minus %>%
  filter(ego_Case3_tumor_tl_c5minus@result$p.adjust<0.05 & ego_Case3_tumor_tl_c5minus@result$qvalue <0.05)
ego_Case4_tumor_tl_c5minus_0.05<-
  ego_Case4_tumor_tl_c5minus %>%
  filter(ego_Case4_tumor_tl_c5minus@result$p.adjust<0.05 & ego_Case4_tumor_tl_c5minus@result$qvalue <0.05)

ego_shared_tumor_tl_c5minus_0.05 <-
Reduce(intersect, list(
                       ego_Case3_tumor_tl_c5minus_0.05@result$ID,
                       ego_Case4_tumor_tl_c5minus_0.05@result$ID))
save(ego_shared_tumor_tl_c5minus_0.05,file="ego_shared_tumor_tl_c5minus_0.05.RData")


[1] "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX"                          
[2] "GOCC_ENDOPLASMIC_RETICULUM_LUMEN"                                       
[3] "GOBP_WOUND_HEALING"                                                     
[4] "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY"                               
[5] "GOBP_RESPONSE_TO_STEROID_HORMONE"                                       
[6] "GOBP_CELL_SUBSTRATE_ADHESION"                                           
[7] "GOCC_CONTRACTILE_FIBER"                                                 
[8] "HP_INGUINAL_HERNIA"                                                     
[9] "GOBP_POSITIVE_REGULATION_OF_SUBSTRATE_ADHESION_DEPENDENT_CELL_SPREADING"
[10] "GOCC_BLOOD_MICROPARTICLE"                                               
[11] "HP_ABNORMAL_PERISTALSIS"                                                
[12] "GOBP_MUSCLE_CONTRACTION"                                                
[13] "HP_ABNORMAL_UMBILICUS_MORPHOLOGY"                                       
[14] "HP_TALIPES_EQUINOVARUS"                                                 
[15] "GOBP_RECEPTOR_CLUSTERING"                                               
[16] "GOBP_RENAL_SYSTEM_VASCULATURE_DEVELOPMENT"                              
[17] "GOBP_MUSCLE_SYSTEM_PROCESS"                                             
[18] "GOBP_CHAPERONE_MEDIATED_AUTOPHAGY"                                      
[19] "GOBP_RESPONSE_TO_CORTICOSTEROID"                                        
[20] "GOBP_RESPONSE_TO_KETONE"                                                
[21] "HP_POSITIONAL_FOOT_DEFORMITY"                                           
[22] "HP_ABNORMAL_ATRIOVENTRICULAR_VALVE_MORPHOLOGY"                          
[23] "GOBP_RESPONSE_TO_CARBOHYDRATE"                                          
[24] "HP_ABNORMAL_TISSUE_METABOLITE_CONCENTRATION"                            
[25] "GOBP_POSITIVE_REGULATION_OF_ACTIN_FILAMENT_BUNDLE_ASSEMBLY"             
[26] "GOBP_TISSUE_MIGRATION"                                                  
[27] "HP_MEGACYSTIS"     
