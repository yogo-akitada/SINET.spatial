
##creating cytoscape input data

#tp_m
ego_shared_tumor_tp_m_c5_0.1

dfego_Case1_tumor_tp_m_c5<-
  as.data.frame(ego_Case1_tumor_tp_m_c5) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_tp_m_c5_0.1)

dfego_Case1_tumor_tp_m_c5 <- data.frame(
  Name = dfego_Case1_tumor_tp_m_c5$ID,
  Description = dfego_Case1_tumor_tp_m_c5$Description,
  p.value = dfego_Case1_tumor_tp_m_c5$p.adjust,
  FDR = dfego_Case1_tumor_tp_m_c5$qvalue,
  Phenotype = "UP",
  Genes = gsub("/", ",", dfego_Case1_tumor_tp_m_c5$geneID)  # Replace "/" with "," if needed
)


dfego_Case1_tumor_tp_m_c5minus<-
  as.data.frame(ego_Case1_tumor_tp_m_c5minus) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_tp_m_c5minus_0.01)

dfego_Case1_tumor_tp_m_c5minus <- data.frame(
  Name = dfego_Case1_tumor_tp_m_c5minus$ID,
  Description = dfego_Case1_tumor_tp_m_c5minus$Description,
  p.value = dfego_Case1_tumor_tp_m_c5minus$p.adjust,
  FDR = dfego_Case1_tumor_tp_m_c5minus$qvalue,
  Phenotype = "DOWN",
  Genes = gsub("/", ",", dfego_Case1_tumor_tp_m_c5minus$geneID)  # Replace "/" with "," if needed
)

dfego_Case1_tumor_tp_m_c5updown <- rbind(dfego_Case1_tumor_tp_m_c5,dfego_Case1_tumor_tp_m_c5minus)
head(dfego_Case1_tumor_tp_m_c5updown)
write.table(dfego_Case1_tumor_tp_m_c5updown, file = "dfego_Case1_tumor_tp_m_c5updown.txt", sep = "\t", quote = FALSE, row.names = FALSE)




dfego_Case1_tumor_tp_m_c2minus<-
  as.data.frame(ego_Case1_tumor_tp_m_c2minus) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_tp_m_c2minus_0.1)

dfego_Case1_tumor_tp_m_c2minus <- data.frame(
  Name = dfego_Case1_tumor_tp_m_c2minus$ID,
  Description = dfego_Case1_tumor_tp_m_c2minus$Description,
  p.value = dfego_Case1_tumor_tp_m_c2minus$p.adjust,
  FDR = dfego_Case1_tumor_tp_m_c2minus$qvalue,
  Phenotype = "DOWN",
  Genes = gsub("/", ",", dfego_Case1_tumor_tp_m_c2minus$geneID)  # Replace "/" with "," if needed
)
dfego_Case1_tumor_tp_m_c2c5updown <- rbind(dfego_Case1_tumor_tp_m_c5updown,dfego_Case1_tumor_tp_m_c2minus)
write.table(dfego_Case1_tumor_tp_m_c2c5updown, file = "dfego_Case1_tumor_tp_m_c2c5updown.txt", sep = "\t", quote = FALSE, row.names = FALSE)




#tm
dfego_Case1_tumor_tm_c5<-
  as.data.frame(ego_Case1_tumor_tm_c5) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_tm_c5_0.1)

dfego_Case1_tumor_tm_c5 <- data.frame(
  Name = dfego_Case1_tumor_tm_c5$ID,
  Description = dfego_Case1_tumor_tm_c5$Description,
  p.value = dfego_Case1_tumor_tm_c5$p.adjust,
  FDR = dfego_Case1_tumor_tm_c5$qvalue,
  Phenotype = "UP",
  Genes = gsub("/", ",", dfego_Case1_tumor_tm_c5$geneID)  # Replace "/" with "," if needed
)


dfego_Case1_tumor_tm_c5minus<-
  as.data.frame(ego_Case1_tumor_tm_c5minus) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_tm_c5minus_0.1)

dfego_Case1_tumor_tm_c5minus <- data.frame(
  Name = dfego_Case1_tumor_tm_c5minus$ID,
  Description = dfego_Case1_tumor_tm_c5minus$Description,
  p.value = dfego_Case1_tumor_tm_c5minus$p.adjust,
  FDR = dfego_Case1_tumor_tm_c5minus$qvalue,
  Phenotype = "DOWN",
  Genes = gsub("/", ",", dfego_Case1_tumor_tm_c5minus$geneID)  # Replace "/" with "," if needed
)

dfego_Case1_tumor_tm_c5updown <- rbind(dfego_Case1_tumor_tm_c5,dfego_Case1_tumor_tm_c5minus)

dfego_Case1_tumor_tm_c2<-
  as.data.frame(ego_Case1_tumor_tm_c2) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_tm_c2_0.1)

dfego_Case1_tumor_tm_c2 <- data.frame(
  Name = dfego_Case1_tumor_tm_c2$ID,
  Description = dfego_Case1_tumor_tm_c2$Description,
  p.value = dfego_Case1_tumor_tm_c2$p.adjust,
  FDR = dfego_Case1_tumor_tm_c2$qvalue,
  Phenotype = "UP",
  Genes = gsub("/", ",", dfego_Case1_tumor_tm_c2$geneID)  # Replace "/" with "," if needed
)

dfego_Case1_tumor_tm_c2minus<-
  as.data.frame(ego_Case1_tumor_tm_c2minus) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_tm_c2minus_0.1)

dfego_Case1_tumor_tm_c2minus <- data.frame(
  Name = dfego_Case1_tumor_tm_c2minus$ID,
  Description = dfego_Case1_tumor_tm_c2minus$Description,
  p.value = dfego_Case1_tumor_tm_c2minus$p.adjust,
  FDR = dfego_Case1_tumor_tm_c2minus$qvalue,
  Phenotype = "DOWN",
  Genes = gsub("/", ",", dfego_Case1_tumor_tm_c2minus$geneID)  # Replace "/" with "," if needed
)

dfego_Case1_tumor_tm_c2updown <- rbind(dfego_Case1_tumor_tm_c2,dfego_Case1_tumor_tm_c2minus)

dfego_Case1_tumor_tm_c2c5updown <- rbind(dfego_Case1_tumor_tm_c5updown,dfego_Case1_tumor_tm_c2updown)
write.table(dfego_Case1_tumor_tm_c2c5updown, file = "dfego_Case1_tumor_tm_c2c5updown.txt", sep = "\t", quote = FALSE, row.names = FALSE)





#tl


dfego_Case4_tumor_tl_c5minus<-
  as.data.frame(ego_Case4_tumor_tl_c5minus) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_tl_c5minus_0.05)

dfego_Case4_tumor_tl_c5minus <- data.frame(
  Name = dfego_Case4_tumor_tl_c5minus$ID,
  Description = dfego_Case4_tumor_tl_c5minus$Description,
  p.value = dfego_Case4_tumor_tl_c5minus$p.adjust,
  FDR = dfego_Case4_tumor_tl_c5minus$qvalue,
  Phenotype = "DOWN",
  Genes = gsub("/", ",", dfego_Case4_tumor_tl_c5minus$geneID)  # Replace "/" with "," if needed
)





dfego_Case4_tumor_tl_c2minus<-
  as.data.frame(ego_Case4_tumor_tl_c2minus) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_tl_c2minus_0.1)

dfego_Case4_tumor_tl_c2minus <- data.frame(
  Name = dfego_Case4_tumor_tl_c2minus$ID,
  Description = dfego_Case4_tumor_tl_c2minus$Description,
  p.value = dfego_Case4_tumor_tl_c2minus$p.adjust,
  FDR = dfego_Case4_tumor_tl_c2minus$qvalue,
  Phenotype = "DOWN",
  Genes = gsub("/", ",", dfego_Case4_tumor_tl_c2minus$geneID)  # Replace "/" with "," if needed
)

dfego_Case4_tumor_tl_c2c5updown <- rbind(dfego_Case4_tumor_tl_c2minus,dfego_Case4_tumor_tl_c5minus)

write.table(dfego_Case4_tumor_tl_c2c5updown, file = "dfego_Case4_tumor_tl_c2c5updown.txt", sep = "\t", quote = FALSE, row.names = FALSE)



#deep
ego_shared_tumor_deep_c5_0.1

dfego_Case1_tumor_deep_c5<-
  as.data.frame(ego_Case1_tumor_deep_c5) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_deep_c5_0.1)

dfego_Case1_tumor_deep_c5 <- data.frame(
  Name = dfego_Case1_tumor_deep_c5$ID,
  Description = dfego_Case1_tumor_deep_c5$Description,
  p.value = dfego_Case1_tumor_deep_c5$p.adjust,
  FDR = dfego_Case1_tumor_deep_c5$qvalue,
  Phenotype = "UP",
  Genes = gsub("/", ",", dfego_Case1_tumor_deep_c5$geneID)  # Replace "/" with "," if needed
)


dfego_Case1_tumor_deep_c5minus<-
  as.data.frame(ego_Case1_tumor_deep_c5minus) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_deep_c5minus_0.1)

dfego_Case1_tumor_deep_c5minus <- data.frame(
  Name = dfego_Case1_tumor_deep_c5minus$ID,
  Description = dfego_Case1_tumor_deep_c5minus$Description,
  p.value = dfego_Case1_tumor_deep_c5minus$p.adjust,
  FDR = dfego_Case1_tumor_deep_c5minus$qvalue,
  Phenotype = "DOWN",
  Genes = gsub("/", ",", dfego_Case1_tumor_deep_c5minus$geneID)  # Replace "/" with "," if needed
)

dfego_Case1_tumor_deep_c5updown <- rbind(dfego_Case1_tumor_deep_c5,dfego_Case1_tumor_deep_c5minus)



dfego_Case1_tumor_deep_c2<-
  as.data.frame(ego_Case1_tumor_deep_c2) %>%
  select(ID,Description,p.adjust,qvalue,geneID) %>%
  filter(ID %in% ego_shared_tumor_deep_c2_0.1)

dfego_Case1_tumor_deep_c2 <- data.frame(
  Name = dfego_Case1_tumor_deep_c2$ID,
  Description = dfego_Case1_tumor_deep_c2$Description,
  p.value = dfego_Case1_tumor_deep_c2$p.adjust,
  FDR = dfego_Case1_tumor_deep_c2$qvalue,
  Phenotype = "UP",
  Genes = gsub("/", ",", dfego_Case1_tumor_deep_c2$geneID)  # Replace "/" with "," if needed
)
dfego_Case1_tumor_deep_c2c5updown <- rbind(dfego_Case1_tumor_deep_c5updown,dfego_Case1_tumor_deep_c2)
write.table(dfego_Case1_tumor_deep_c2c5updown, file = "dfego_Case1_tumor_deep_c2c5updown.txt", sep = "\t", quote = FALSE, row.names = FALSE)

