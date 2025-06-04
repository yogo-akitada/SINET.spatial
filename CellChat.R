
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

cellchat_Case1 <- createCellChat(object = Case1_tumor_normal, assay = "Spatial", group.by = "manual_ident3")
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)


CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 

cellchat_Case1@DB <- CellChatDB.use

cellchat_Case1 <- subsetData(cellchat_Case1) 

cellchat_Case1 <- identifyOverExpressedGenes(cellchat_Case1)
cellchat_Case1 <- identifyOverExpressedInteractions(cellchat_Case1)


cellchat_Case1 <- computeCommunProb(cellchat_Case1)
cellchat_Case1 <- filterCommunication(cellchat_Case1, min.cells = 10)
cellchat_Case1 <- computeCommunProbPathway(cellchat_Case1)
cellchat_Case1 <- aggregateNet(cellchat_Case1)

## same code for patient 2-4


levels(cellchat_Case1@meta$manual_ident3)

[1] "m_superficial"    "m_mid_crypt"      "m_base"           "sm"               "mp_inner"         "mp_outer"        
[7] "mp_peritumoral"   "myenteric plexus" "serosa"           "mes_stroma"       "mes_nerve"        "tp_mucosal"      
[13] "tp_deep"         "tm"    

levels(cellchat_Case2@meta$manual_ident3)

[1] "m_superficial_mid_crypt"   "m_base"                    "sm"                        "mp_inner"                 
[5] "mp_outer"                  "serosa"                    "vascular_wall"             "mes_dense_collagen_stroma"
[9] "mes_nerve"                 "tp_mucosal"                "tp_deep"                  "tm"

levels(cellchat_Case3@meta$manual_ident3)

[1] "m_base"               "m_lymphoid_follicle"  "sm_non_peritumoral"   "sm_peritumoral"       "sm_deep"             
[6] "mp_non_peritumoral"   "mp_peritumoral"       "mp_inner_peritumoral" "mes_collagen"         "tp_mucosal"          
[11] "tp_deep"             "tm"                   "tl"     

levels(cellchat_Case4@meta$manual_ident3)

[1] "m_mid_crypt"              "m_base"                   "sm"                       "mp_inner_non_peritumoral"
[5] "mp_outer_non_peritumoral" "mp_peritumoral"           "tp_mucosal"               "tp_deep"                
[9] "tl" 


pathways.show.all1 <- cellchat_Case1@netP$pathways
pathways.show.all2 <- cellchat_Case2@netP$pathways
pathways.show.all3 <- cellchat_Case3@netP$pathways
pathways.show.all4 <- cellchat_Case4@netP$pathways


## from_tumor



netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_tumor <-
  netAnalysis_contribution(cellchat_Case1, signaling = pathways.show.all1,targets.use =c(12),sources.use = c(12:14),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_tumor <-
  netAnalysis_contribution(cellchat_Case2, signaling = pathways.show.all2,targets.use =c(10),sources.use = c(10:12),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_tumor <-
  netAnalysis_contribution(cellchat_Case3, signaling = pathways.show.all3,targets.use =c(10),sources.use = c(10:13),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_tumor <-
  netAnalysis_contribution(cellchat_Case4, signaling = pathways.show.all4,targets.use =c(7),sources.use = c(7:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

shared_tp_mucosal_from_tumor<-
  Reduce(intersect, list(netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_tumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_tumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_tumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_tumor$LR.contribution$name))
[1] "MDK - NCL"

Case1_tp_mucosal_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_tumor$LR.contribution) %>%
  mutate(patient = "patient1")
Case2_tp_mucosal_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_tumor$LR.contribution) %>%
  mutate(patient = "patient2")
Case3_tp_mucosal_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_tumor$LR.contribution) %>%
  mutate(patient = "patient3")
Case4_tp_mucosal_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_tumor$LR.contribution) %>%
  mutate(patient = "patient4")

tp_mucosal_from_tumorLR <-rbind(Case1_tp_mucosal_from_tumorLR,
                                Case2_tp_mucosal_from_tumorLR,
                                Case3_tp_mucosal_from_tumorLR,
                                Case4_tp_mucosal_from_tumorLR)

tp_mucosal_from_tumorLR$contribution <- as.numeric(tp_mucosal_from_tumorLR$contribution)

tp_mucosal_from_tumorLR_sum <- data.frame(
  LR = names(tapply(tp_mucosal_from_tumorLR$contribution,
                    tp_mucosal_from_tumorLR$name,
                    sum)),
  sum = tapply(tp_mucosal_from_tumorLR$contribution,
               tp_mucosal_from_tumorLR$name,
               sum)) %>%
  filter(LR %in% shared_tp_mucosal_from_tumor) %>%
  arrange(desc(sum)) 

write.csv(tp_mucosal_from_tumorLR_sum,"tp_mucosal_from_tumorLR_sum.csv")

tp_mucosal_from_tumorLR$contribution <- as.numeric(tp_mucosal_from_tumorLR$contribution)

tp_mucosal_from_tumorLR_sum2 <- data.frame(
  LR = names(tapply(tp_mucosal_from_tumorLR$contribution,
                    tp_mucosal_from_tumorLR$name,
                    sum)),
  sum = tapply(tp_mucosal_from_tumorLR$contribution,
               tp_mucosal_from_tumorLR$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(tp_mucosal_from_tumorLR_sum2,"tp_mucosal_from_tumorLR_sum2.csv")



netAnalysis_contribution_cellchat_Case1_tm_from_tumor <-
  netAnalysis_contribution(cellchat_Case1, signaling = pathways.show.all1,targets.use =c(14),sources.use = c(12:14),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case1_tm_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case1_tm_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case1_tm_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case1_tm_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case2_tm_from_tumor <-
  netAnalysis_contribution(cellchat_Case2, signaling = pathways.show.all2,targets.use =c(12),sources.use = c(10:12),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case2_tm_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case2_tm_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case2_tm_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case2_tm_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case3_tm_from_tumor <-
  netAnalysis_contribution(cellchat_Case3, signaling = pathways.show.all3,targets.use =c(12),sources.use = c(10:13),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case3_tm_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case3_tm_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case3_tm_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case3_tm_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

shared_tm_from_tumor<-
  Reduce(intersect, list(netAnalysis_contribution_cellchat_Case1_tm_from_tumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case2_tm_from_tumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case3_tm_from_tumor$LR.contribution$name))
[1] "MDK - NCL"               "NRG1 - ERBB3"            "MDK - SDC4"              "MDK - LRP1"             
[5] "FGF9 - FGFR3"            "SEMA3B - (NRP1+PLXNA2)"  "FGF9 - FGFR1"            "GDF11 - (ACVR1B+ACVR2B)"
[9] "GDF11 - (TGFBR1+ACVR2A)" "GDF11 - (TGFBR1+ACVR2B)" "TGFB1 - (ACVR1B+TGFBR2)" "SEMA3G - (NRP2+PLXNA2)" 
[13] "SEMA3G - (NRP2+PLXNA3)" 

Case1_tm_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case1_tm_from_tumor$LR.contribution) %>%
  mutate(patient = "patient1")
Case2_tm_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case2_tm_from_tumor$LR.contribution) %>%
  mutate(patient = "patient2")
Case3_tm_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case3_tm_from_tumor$LR.contribution) %>%
  mutate(patient = "patient3")
Case1_tm_from_nontumorLR$name

tm_from_tumorLR <-rbind(Case1_tm_from_tumorLR,
                        Case2_tm_from_tumorLR,
                        Case3_tm_from_tumorLR)

tm_from_tumorLR$contribution <- as.numeric(tm_from_tumorLR$contribution)

tm_from_tumorLR_sum <- data.frame(
  LR = names(tapply(tm_from_tumorLR$contribution,
                    tm_from_tumorLR$name,
                    sum)),
  sum = tapply(tm_from_tumorLR$contribution,
               tm_from_tumorLR$name,
               sum)) %>%
  filter(LR %in% shared_tm_from_tumor)%>%
  arrange(desc(sum)) 

write.csv(tm_from_tumorLR_sum,"tm_from_tumorLR_sum.csv")

tm_from_tumorLR_sum2 <- data.frame(
  LR = names(tapply(tm_from_tumorLR$contribution,
                    tm_from_tumorLR$name,
                    sum)),
  sum = tapply(tm_from_tumorLR$contribution,
               tm_from_tumorLR$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(tm_from_tumorLR_sum2,"tm_from_tumorLR_sum2.csv")


netAnalysis_contribution_cellchat_Case3_tl_from_tumor <-
  netAnalysis_contribution(cellchat_Case3, signaling = pathways.show.all3,targets.use =c(13),sources.use = c(10:13),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case3_tl_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case3_tl_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case3_tl_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case3_tl_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case4_tl_from_tumor <-
  netAnalysis_contribution(cellchat_Case4, signaling = pathways.show.all4,targets.use =c(9),sources.use = c(7:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case4_tl_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case4_tl_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case4_tl_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case4_tl_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

shared_tl_from_tumor<-
  Reduce(intersect, list(
    netAnalysis_contribution_cellchat_Case3_tl_from_tumor$LR.contribution$name,
    netAnalysis_contribution_cellchat_Case4_tl_from_tumor$LR.contribution$name))
[1] "NRG1 - ERBB3"            "MDK - NCL"               "MDK - SDC4"              "FGF9 - FGFR3"           
[5] "FGF9 - FGFR1"            "PDGFA - PDGFRB"          "WNT4 - (FZD8+LRP6)"      "WNT4 - (FZD8+LRP5)"     
[9] "GDF11 - (ACVR1B+ACVR2B)" "GDF11 - (TGFBR1+ACVR2B)" "GDF11 - (TGFBR1+ACVR2A)" "WNT2B - (FZD8+LRP6)"    
[13] "WNT2B - (FZD8+LRP5)"   


Case3_tl_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case3_tl_from_tumor$LR.contribution) %>%
  mutate(patient = "patient3")
Case4_tl_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case4_tl_from_tumor$LR.contribution) %>%
  mutate(patient = "patient4")

tl_from_tumorLR <-rbind(
  Case3_tl_from_tumorLR,
  Case4_tl_from_tumorLR)

tl_from_tumorLR$contribution <- as.numeric(tl_from_tumorLR$contribution)

tl_from_tumorLR_sum <- data.frame(
  LR = names(tapply(tl_from_tumorLR$contribution,
                    tl_from_tumorLR$name,
                    sum)),
  sum = tapply(tl_from_tumorLR$contribution,
               tl_from_tumorLR$name,
               sum)) %>%
  filter(LR %in% shared_tl_from_tumor)%>%
  arrange(desc(sum)) 

write.csv(tl_from_tumorLR_sum,"tl_from_tumorLR_sum.csv")

tl_from_tumorLR_sum2 <- data.frame(
  LR = names(tapply(tl_from_tumorLR$contribution,
                    tl_from_tumorLR$name,
                    sum)),
  sum = tapply(tl_from_tumorLR$contribution,
               tl_from_tumorLR$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(tl_from_tumorLR_sum2,"tl_from_tumorLR_sum2.csv")






netAnalysis_contribution_cellchat_Case1_tp_other_from_tumor <-
  netAnalysis_contribution(cellchat_Case1, signaling = pathways.show.all1,targets.use =c(13),sources.use = c(12:14),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case1_tp_other_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case1_tp_other_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case1_tp_other_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case1_tp_other_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case2_tp_other_from_tumor <-
  netAnalysis_contribution(cellchat_Case2, signaling = pathways.show.all2,targets.use =c(11),sources.use = c(10:12),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case2_tp_other_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case2_tp_other_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case2_tp_other_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case2_tp_other_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case3_tp_other_from_tumor <-
  netAnalysis_contribution(cellchat_Case3, signaling = pathways.show.all3,targets.use =c(11),sources.use = c(10:13),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case3_tp_other_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case3_tp_other_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case3_tp_other_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case3_tp_other_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case4_tp_other_from_tumor <-
  netAnalysis_contribution(cellchat_Case4, signaling = pathways.show.all4,targets.use =c(8),sources.use = c(7:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case4_tp_other_from_tumor$LR.contribution,"netAnalysis_contribution_cellchat_Case4_tp_other_from_tumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case4_tp_other_from_tumor_.svg", plot=netAnalysis_contribution_cellchat_Case4_tp_other_from_tumor$gg.obj, width = 5, height = 5, units = 'in')

shared_tp_other_from_tumor<-
  Reduce(intersect, list(netAnalysis_contribution_cellchat_Case1_tp_other_from_tumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case2_tp_other_from_tumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case3_tp_other_from_tumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case4_tp_other_from_tumor$LR.contribution$name))

[1] "MDK - NCL"               "MDK - SDC4"              "NRG1 - ERBB3"            "MDK - LRP1"             
[5] "FGF9 - FGFR3"            "FGF9 - FGFR1"            "TGFB1 - (ACVR1B+TGFBR2)"


Case1_tp_other_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case1_tp_other_from_tumor$LR.contribution) %>%
  mutate(patient = "patient1")
Case2_tp_other_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case2_tp_other_from_tumor$LR.contribution) %>%
  mutate(patient = "patient2")
Case3_tp_other_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case3_tp_other_from_tumor$LR.contribution) %>%
  mutate(patient = "patient3")
Case4_tp_other_from_tumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case4_tp_other_from_tumor$LR.contribution) %>%
  mutate(patient = "patient4")

tp_other_from_tumorLR <-rbind(Case1_tp_other_from_tumorLR,
                              Case2_tp_other_from_tumorLR,
                              Case3_tp_other_from_tumorLR,
                              Case4_tp_other_from_tumorLR)

tp_other_from_tumorLR$contribution <- as.numeric(tp_other_from_tumorLR$contribution)

tp_other_from_tumorLR_sum <- data.frame(
  LR = names(tapply(tp_other_from_tumorLR$contribution,
                    tp_other_from_tumorLR$name,
                    sum)),
  sum = tapply(tp_other_from_tumorLR$contribution,
               tp_other_from_tumorLR$name,
               sum)) %>%
  filter(LR %in% shared_tp_other_from_tumor)%>%
  arrange(desc(sum)) 

write.csv(tp_other_from_tumorLR_sum,"tp_other_from_tumorLR_sum.csv")
tp_other_from_tumorLR_sum2 <- data.frame(
  LR = names(tapply(tp_other_from_tumorLR$contribution,
                    tp_other_from_tumorLR$name,
                    sum)),
  sum = tapply(tp_other_from_tumorLR$contribution,
               tp_other_from_tumorLR$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(tp_other_from_tumorLR_sum2,"tp_other_from_tumorLR_sum2.csv")


## from_nontumor



netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_nontumor <-
  netAnalysis_contribution(cellchat_Case1, signaling = pathways.show.all1,targets.use =c(12),sources.use = c(1:11),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_nontumor <-
  netAnalysis_contribution(cellchat_Case2, signaling = pathways.show.all2,targets.use =c(10),sources.use = c(1:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_nontumor <-
  netAnalysis_contribution(cellchat_Case3, signaling = pathways.show.all3,targets.use =c(10),sources.use = c(1:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_nontumor <-
  netAnalysis_contribution(cellchat_Case4, signaling = pathways.show.all4,targets.use =c(7),sources.use = c(1:6),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

shared_tp_mucosal_from_nontumor<-
  Reduce(intersect, list(netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_nontumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_nontumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_nontumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_nontumor$LR.contribution$name))
[1] "PTN - NCL"  "GCG - GIPR"

Case1_tp_mucosal_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case1_tp_mucosal_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient1")
Case2_tp_mucosal_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case2_tp_mucosal_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient2")
Case3_tp_mucosal_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case3_tp_mucosal_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient3")
Case4_tp_mucosal_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case4_tp_mucosal_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient4")

tp_mucosal_from_nontumorLR <-rbind(Case1_tp_mucosal_from_nontumorLR,
                                   Case2_tp_mucosal_from_nontumorLR,
                                   Case3_tp_mucosal_from_nontumorLR,
                                   Case4_tp_mucosal_from_nontumorLR)

tp_mucosal_from_nontumorLR$contribution <- as.numeric(tp_mucosal_from_nontumorLR$contribution)

tp_mucosal_from_nontumorLR_sum <- data.frame(
  LR = names(tapply(tp_mucosal_from_nontumorLR$contribution,
                    tp_mucosal_from_nontumorLR$name,
                    sum)),
  sum = tapply(tp_mucosal_from_nontumorLR$contribution,
               tp_mucosal_from_nontumorLR$name,
               sum)) %>%
  filter(LR %in% shared_tp_mucosal_from_nontumor)%>%
  arrange(desc(sum)) 

write.csv(tp_mucosal_from_nontumorLR_sum,"tp_mucosal_from_nontumorLR_sum.csv")
tp_mucosal_from_nontumorLR_sum2 <- data.frame(
  LR = names(tapply(tp_mucosal_from_nontumorLR$contribution,
                    tp_mucosal_from_nontumorLR$name,
                    sum)),
  sum = tapply(tp_mucosal_from_nontumorLR$contribution,
               tp_mucosal_from_nontumorLR$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(tp_mucosal_from_nontumorLR_sum2,"tp_mucosal_from_nontumorLR_sum2.csv")


netAnalysis_contribution_cellchat_Case1_tm_from_nontumor <-
  netAnalysis_contribution(cellchat_Case1, signaling = pathways.show.all1,targets.use =c(14),sources.use = c(1:11),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case1_tm_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case1_tm_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case1_tm_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case1_tm_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case2_tm_from_nontumor <-
  netAnalysis_contribution(cellchat_Case2, signaling = pathways.show.all2,targets.use =c(12),sources.use = c(1:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case2_tm_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case2_tm_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case2_tm_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case2_tm_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case3_tm_from_nontumor <-
  netAnalysis_contribution(cellchat_Case3, signaling = pathways.show.all3,targets.use =c(12),sources.use = c(1:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case3_tm_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case3_tm_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case3_tm_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case3_tm_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

shared_tm_from_nontumor<-
  Reduce(intersect, list(netAnalysis_contribution_cellchat_Case1_tm_from_nontumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case2_tm_from_nontumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case3_tm_from_nontumor$LR.contribution$name))
[1] "SEMA3C - PLXND1"         "PTN - NCL"               "SEMA3B - (NRP1+PLXNA2)"  "PTN - SDC4"             
[5] "SEMA3C - (NRP1+NRP2)"    "SEMA3B - (NRP1+PLXNA3)"  "SEMA3C - (NRP1+PLXNA2)"  "GCG - GIPR"             
[9] "SEMA3B - (NRP2+PLXNA2)"  "GRN - SORT1"             "PDGFA - PDGFRB"          "SEMA3C - (NRP1+PLXNA3)" 
[13] "GUCA2A - GUCY2C"         "SEMA3B - (NRP2+PLXNA3)"  "SEMA3C - (NRP2+PLXNA2)"  "PTN - SDC3"             
[17] "NRG1 - ERBB3"            "ADM - CALCR"             "FGF7 - FGFR1"            "SEMA3C - (NRP2+PLXNA3)" 
[21] "SEMA3G - (NRP2+PLXNA2)"  "SEMA3G - (NRP2+PLXNA3)"  "INHBA - (ACVR1B+ACVR2A)" "INHBA - (ACVR1B+ACVR2B)"
[25] "CSF1 - CSF1R"            "TGFB1 - (ACVR1B+TGFBR2)" "TGFB2 - (ACVR1B+TGFBR2)" "TGFB1 - (TGFBR1+TGFBR2)"
[29] "TGFB2 - (TGFBR1+TGFBR2)"


Case1_tm_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case1_tm_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient1")
Case2_tm_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case2_tm_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient2")
Case3_tm_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case3_tm_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient3")


tm_from_nontumorLR <-rbind(Case1_tm_from_nontumorLR,
                           Case2_tm_from_nontumorLR,
                           Case3_tm_from_nontumorLR)

tm_from_nontumorLR$contribution <- as.numeric(tm_from_nontumorLR$contribution)

tm_from_nontumorLR_sum <- data.frame(
  LR = names(tapply(tm_from_nontumorLR$contribution,
                    tm_from_nontumorLR$name,
                    sum)),
  sum = tapply(tm_from_nontumorLR$contribution,
               tm_from_nontumorLR$name,
               sum)) %>%
  filter(LR %in% shared_tm_from_nontumor)%>%
  arrange(desc(sum)) 

write.csv(tm_from_nontumorLR_sum,"tm_from_nontumorLR_sum.csv")
tm_from_nontumorLR_sum2 <- data.frame(
  LR = names(tapply(tm_from_nontumorLR$contribution,
                    tm_from_nontumorLR$name,
                    sum)),
  sum = tapply(tm_from_nontumorLR$contribution,
               tm_from_nontumorLR$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(tm_from_nontumorLR_sum2,"tm_from_nontumorLR_sum.csv")


netAnalysis_contribution_cellchat_Case3_tl_from_nontumor <-
  netAnalysis_contribution(cellchat_Case3, signaling = pathways.show.all3,targets.use =c(13),sources.use = c(1:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case3_tl_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case3_tl_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case3_tl_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case3_tl_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case4_tl_from_nontumor <-
  netAnalysis_contribution(cellchat_Case4, signaling = pathways.show.all4,targets.use =c(9),sources.use = c(1:6),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case4_tl_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case4_tl_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case4_tl_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case4_tl_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

shared_tl_from_nontumor<-
  Reduce(intersect, list(
    netAnalysis_contribution_cellchat_Case3_tl_from_nontumor$LR.contribution$name,
    netAnalysis_contribution_cellchat_Case4_tl_from_nontumor$LR.contribution$name))
[1] "PTN - NCL"               "SEMA3C - PLXND1"         "SEMA3B - (NRP1+PLXNA3)"  "CXCL12 - CXCR4"         
[5] "SEMA3B - (NRP1+PLXNA2)"  "PTN - SDC4"              "ADM - CALCR"             "GCG - GIPR"             
[9] "NAMPT - INSR"            "TGFB1 - (ACVR1B+TGFBR2)" "SEMA3C - (NRP1+PLXNA3)"  "GUCA2A - GUCY2C"        
[13] "SEMA3C - (NRP1+PLXNA2)"  "TGFB1 - (TGFBR1+TGFBR2)" "PTN - SDC3"              "TGFB2 - (ACVR1B+TGFBR2)"
[17] "PRSS3 - PARD3"           "TGFB2 - (TGFBR1+TGFBR2)"


Case3_tl_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case3_tl_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient3")
Case4_tl_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case4_tl_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient4")

tl_from_nontumorLR <-rbind(
  Case3_tl_from_nontumorLR,
  Case4_tl_from_nontumorLR)

tl_from_nontumorLR$contribution <- as.numeric(tl_from_nontumorLR$contribution)

tl_from_nontumorLR_sum <- data.frame(
  LR = names(tapply(tl_from_nontumorLR$contribution,
                    tl_from_nontumorLR$name,
                    sum)),
  sum = tapply(tl_from_nontumorLR$contribution,
               tl_from_nontumorLR$name,
               sum)) %>%
  filter(LR %in% shared_tl_from_nontumor)%>%
  arrange(desc(sum)) 

write.csv(tl_from_nontumorLR_sum,"tl_from_nontumorLR_sum.csv")
tl_from_nontumorLR_sum2 <- data.frame(
  LR = names(tapply(tl_from_nontumorLR$contribution,
                    tl_from_nontumorLR$name,
                    sum)),
  sum = tapply(tl_from_nontumorLR$contribution,
               tl_from_nontumorLR$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(tl_from_nontumorLR_sum2,"tl_from_nontumorLR_sum2.csv")


netAnalysis_contribution_cellchat_Case1_tp_other_from_nontumor <-
  netAnalysis_contribution(cellchat_Case1, signaling = pathways.show.all1,targets.use =c(13),sources.use = c(1:11),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case1_tp_other_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case1_tp_other_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case1_tp_other_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case1_tp_other_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case2_tp_other_from_nontumor <-
  netAnalysis_contribution(cellchat_Case2, signaling = pathways.show.all2,targets.use =c(11),sources.use = c(1:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case2_tp_other_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case2_tp_other_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case2_tp_other_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case2_tp_other_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case3_tp_other_from_nontumor <-
  netAnalysis_contribution(cellchat_Case3, signaling = pathways.show.all3,targets.use =c(11),sources.use = c(1:9),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case3_tp_other_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case3_tp_other_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case3_tp_other_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case3_tp_other_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case4_tp_other_from_nontumor <-
  netAnalysis_contribution(cellchat_Case4, signaling = pathways.show.all4,targets.use =c(8),sources.use = c(1:6),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case4_tp_other_from_nontumor$LR.contribution,"netAnalysis_contribution_cellchat_Case4_tp_other_from_nontumor_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case4_tp_other_from_nontumor_.svg", plot=netAnalysis_contribution_cellchat_Case4_tp_other_from_nontumor$gg.obj, width = 5, height = 5, units = 'in')

shared_tp_other_from_nontumor<-
  Reduce(intersect, list(netAnalysis_contribution_cellchat_Case1_tp_other_from_nontumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case2_tp_other_from_nontumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case3_tp_other_from_nontumor$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case4_tp_other_from_nontumor$LR.contribution$name))

[1] "PTN - NCL"               "SEMA3C - PLXND1"         "PTN - SDC4"              "SEMA3C - (NRP1+NRP2)"   
[5] "GCG - GIPR"              "GUCA2A - GUCY2C"         "ADM - CALCR"             "TGFB1 - (ACVR1B+TGFBR2)"
[9] "TGFB2 - (ACVR1B+TGFBR2)"          
[17] "ADM - CALCR"             "SEMA3G - (NRP2+PLXNA2)"  "SEMA3G - (NRP2+PLXNA3)"  "TGFB1 - (ACVR1B+TGFBR2)"
[21] "TGFB1 - (TGFBR1+TGFBR2)" "PRSS3 - F2R"             "IGF1 - IGF1R"            "TGFB2 - (ACVR1B+TGFBR2)"
[25] "TGFB2 - (TGFBR1+TGFBR2)"

Case1_tp_other_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case1_tp_other_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient1")
Case2_tp_other_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case2_tp_other_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient2")
Case3_tp_other_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case3_tp_other_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient3")
Case4_tp_other_from_nontumorLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case4_tp_other_from_nontumor$LR.contribution) %>%
  mutate(patient = "patient4")

tp_other_from_nontumorLR <-rbind(Case1_tp_other_from_nontumorLR,
                                 Case2_tp_other_from_nontumorLR,
                                 Case3_tp_other_from_nontumorLR,
                                 Case4_tp_other_from_nontumorLR)

tp_other_from_nontumorLR$contribution <- as.numeric(tp_other_from_nontumorLR$contribution)

tp_other_from_nontumorLR_sum <- data.frame(
  LR = names(tapply(tp_other_from_nontumorLR$contribution,
                    tp_other_from_nontumorLR$name,
                    sum)),
  sum = tapply(tp_other_from_nontumorLR$contribution,
               tp_other_from_nontumorLR$name,
               sum)) %>%
  filter(LR %in% shared_tp_other_from_nontumor)%>%
  arrange(desc(sum)) 

write.csv(tp_other_from_nontumorLR_sum,"tp_other_from_nontumorLR_sum.csv")
tp_other_from_nontumorLR_sum2 <- data.frame(
  LR = names(tapply(tp_other_from_nontumorLR$contribution,
                    tp_other_from_nontumorLR$name,
                    sum)),
  sum = tapply(tp_other_from_nontumorLR$contribution,
               tp_other_from_nontumorLR$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(tp_other_from_nontumorLR_sum2,"tp_other_from_nontumorLR_sum2.csv")







## m_from_m


netAnalysis_contribution_cellchat_Case1_m_from_m <-
  netAnalysis_contribution(cellchat_Case1, signaling = pathways.show.all1,targets.use =c(1:3),sources.use = c(1:3),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case1_m_from_m$LR.contribution,"netAnalysis_contribution_cellchat_Case1_m_from_m_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case1_m_from_m_.svg", plot=netAnalysis_contribution_cellchat_Case1_m_from_m$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case2_m_from_m <-
  netAnalysis_contribution(cellchat_Case2, signaling = pathways.show.all2,targets.use =c(1:2),sources.use = c(1:2),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case2_m_from_m$LR.contribution,"netAnalysis_contribution_cellchat_Case2_m_from_m_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case2_m_from_m_.svg", plot=netAnalysis_contribution_cellchat_Case2_m_from_m$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case3_m_from_m <-
  netAnalysis_contribution(cellchat_Case3, signaling = pathways.show.all3,targets.use =c(1),sources.use = c(1),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case3_m_from_m$LR.contribution,"netAnalysis_contribution_cellchat_Case3_m_from_m_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case3_m_from_m_.svg", plot=netAnalysis_contribution_cellchat_Case3_m_from_m$gg.obj, width = 5, height = 5, units = 'in')

netAnalysis_contribution_cellchat_Case4_m_from_m <-
  netAnalysis_contribution(cellchat_Case4, signaling = pathways.show.all4,targets.use =c(1:2),sources.use = c(1:2),return.data = TRUE)
write.csv(netAnalysis_contribution_cellchat_Case4_m_from_m$LR.contribution,"netAnalysis_contribution_cellchat_Case4_m_from_m_LR.contribution.csv")
ggsave(filename="netAnalysis_contribution_cellchat_Case4_m_from_m_.svg", plot=netAnalysis_contribution_cellchat_Case4_m_from_m$gg.obj, width = 5, height = 5, units = 'in')

shared_m_from_m<-
  Reduce(intersect, list(netAnalysis_contribution_cellchat_Case1_m_from_m$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case2_m_from_m$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case3_m_from_m$LR.contribution$name,
                         netAnalysis_contribution_cellchat_Case4_m_from_m$LR.contribution$name))
[1] "PRSS3 - F2R" "GCG - GIPR" 

Case1_m_from_mLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case1_m_from_m$LR.contribution) %>%
  mutate(patient = "patient1")
Case2_m_from_mLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case2_m_from_m$LR.contribution) %>%
  mutate(patient = "patient2")
Case3_m_from_mLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case3_m_from_m$LR.contribution) %>%
  mutate(patient = "patient3")
Case4_m_from_mLR <-
  as.data.frame(netAnalysis_contribution_cellchat_Case4_m_from_m$LR.contribution) %>%
  mutate(patient = "patient4")

m_from_mLR <-rbind(Case1_m_from_mLR,
                   Case2_m_from_mLR,
                   Case3_m_from_mLR,
                   Case4_m_from_mLR)

m_from_mLR$contribution <- as.numeric(m_from_mLR$contribution)

m_from_mLR_sum <- data.frame(
  LR = names(tapply(m_from_mLR$contribution,
                    m_from_mLR$name,
                    sum)),
  sum = tapply(m_from_mLR$contribution,
               m_from_mLR$name,
               sum)) %>%
  filter(LR %in% shared_m_from_m) %>%
  arrange(desc(sum)) 

write.csv(m_from_mLR_sum,"m_from_mLR_sum.csv")


m_from_mLR_sum2 <- data.frame(
  LR = names(tapply(m_from_mLR$contribution,
                    m_from_mLR$name,
                    sum)),
  sum = tapply(m_from_mLR$contribution,
               m_from_mLR$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(m_from_mLR_sum2,"m_from_mLR_sum2.csv")
