## Supp Table 2
c1<-as.data.frame(table(Case1_tumor_normal$manual_ident2))
write.csv(c1,"case1spotsnumber.csv")
c2<-as.data.frame(table(Case2_tumor_normal$manual_ident2))

write.csv(c2,"case2spotsnumber.csv")
c3<-as.data.frame(table(Case3_tumor_normal$manual_ident2))

write.csv(c3,"case3spotsnumber.csv")
c4<-as.data.frame(table(Case4_tumor_normal$manual_ident2))

write.csv(c4,"case4spotsnumber.csv")


## Supp Fig. 4
##1
names(manual.ident3) <- levels(Case1_tumor_normal)
Case1_tumor_normal <- RenameIdents(Case1_tumor_normal, manual.ident3)
levels(Case1_tumor_normal)
Case1_tumor_normal@meta.data$manual_ident3 <- factor(Case1_tumor_normal@active.ident)

Idents(Case1_tumor_normal) <- "manual_ident3"
levels(Case1_tumor_normal)
Case1_tumor_normal.markers <- FindAllMarkers(Case1_tumor_normal, 
                                assay = "Spatial",
                                only.pos = TRUE,
                                logfc.threshold = 0.3,
                                test.use="MAST")
Case1_tumor_normal.markers_FC0.3pval<- filter(Case1_tumor_normal.markers, Case1_tumor_normal.markers$p_val_adj < 10e-2) 
clustered_list <- split(rownames(Case1_tumor_normal.markers_FC0.3pval), Case1_tumor_normal.markers_FC0.3pval$cluster)

c8all <- read.gmt("c8.all.v2023.2.Hs.symbols.gmt")

Case1_tumor_normalclustered_list <- compareCluster(clustered_list,
                                          fun = "enricher", TERM2GENE = c8all)

dotCase1_tumor_normalclustered_list<-
clusterProfiler::dotplot(Case1_tumor_normalclustered_list,showCategory = 2)
dotCase1_tumor_normalclustered_list
ggsave(filename="dotCase1_tumor_normalclustered_list.svg",width=20, height=8,plot=dotCase1_tumor_normalclustered_list)


##2

Idents(Case2_tumor_normal) <- "manual_ident3"
levels(Case2_tumor_normal)
Case2_tumor_normal.markers <- FindAllMarkers(Case2_tumor_normal, 
                                             assay = "Spatial",
                                             only.pos = TRUE,
                                             logfc.threshold = 0.3,
                                             test.use="MAST")
Case2_tumor_normal.markers_FC0.3pval<- filter(Case2_tumor_normal.markers, Case2_tumor_normal.markers$p_val_adj < 10e-2) 
clustered_list <- split(rownames(Case2_tumor_normal.markers_FC0.3pval), Case2_tumor_normal.markers_FC0.3pval$cluster)

c8all <- read.gmt("c8.all.v2023.2.Hs.symbols.gmt")

Case2_tumor_normalclustered_list <- compareCluster(clustered_list,
                                                   fun = "enricher", TERM2GENE = c8all)

dotCase2_tumor_normalclustered_list<-
  clusterProfiler::dotplot(Case2_tumor_normalclustered_list,showCategory = 2)
dotCase2_tumor_normalclustered_list
ggsave(filename="dotCase2_tumor_normalclustered_list.svg",width=20, height=8,plot=dotCase2_tumor_normalclustered_list)

##3

Idents(Case3_tumor_normal) <- "manual_ident3"
levels(Case3_tumor_normal)
Case3_tumor_normal.markers <- FindAllMarkers(Case3_tumor_normal, 
                                             assay = "Spatial",
                                             only.pos = TRUE,
                                             logfc.threshold = 0.3,
                                             test.use="MAST")
Case3_tumor_normal.markers_FC0.3pval<- filter(Case3_tumor_normal.markers, Case3_tumor_normal.markers$p_val_adj < 10e-2) 
clustered_list <- split(rownames(Case3_tumor_normal.markers_FC0.3pval), Case3_tumor_normal.markers_FC0.3pval$cluster)

c8all <- read.gmt("c8.all.v2023.2.Hs.symbols.gmt")

Case3_tumor_normalclustered_list <- compareCluster(clustered_list,
                                                   fun = "enricher", TERM2GENE = c8all)

dotCase3_tumor_normalclustered_list<-
  clusterProfiler::dotplot(Case3_tumor_normalclustered_list,showCategory = 2)
dotCase3_tumor_normalclustered_list
ggsave(filename="dotCase3_tumor_normalclustered_list.svg",width=20, height=8,plot=dotCase3_tumor_normalclustered_list)

##4

Idents(Case4_tumor_normal) <- "manual_ident3"
levels(Case4_tumor_normal)
Case4_tumor_normal.markers <- FindAllMarkers(Case4_tumor_normal, 
                                             assay = "Spatial",
                                             only.pos = TRUE,
                                             logfc.threshold = 0.3,
                                             test.use="MAST")
Case4_tumor_normal.markers_FC0.3pval<- filter(Case4_tumor_normal.markers, Case4_tumor_normal.markers$p_val_adj < 10e-2) 
clustered_list <- split(rownames(Case4_tumor_normal.markers_FC0.3pval), Case4_tumor_normal.markers_FC0.3pval$cluster)

c8all <- read.gmt("c8.all.v2023.2.Hs.symbols.gmt")

Case4_tumor_normalclustered_list <- compareCluster(clustered_list,
                                                   fun = "enricher", TERM2GENE = c8all)

dotCase4_tumor_normalclustered_list<-
  clusterProfiler::dotplot(Case4_tumor_normalclustered_list,showCategory = 2)
dotCase4_tumor_normalclustered_list
ggsave(filename="dotCase4_tumor_normalclustered_list.svg",width=20, height=8,plot=dotCase4_tumor_normalclustered_list)

