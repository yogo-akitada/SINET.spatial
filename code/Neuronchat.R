
library(NeuronChat)
library(CellChat)


Idents(Case1_tumor_normal) <- "manual_ident3"
counts.df <- Case1_tumor_normal@assays$Spatial@data %>% as.matrix %>% t %>% as.data.frame
counts.df <- tibble::rownames_to_column(counts.df, "cellnames")
clusterassignemnts <- data.frame(Case1_tumor_normal@active.ident)
clusterassignemnts <- tibble::rownames_to_column(clusterassignemnts, "cellnames")
counts.df <- merge(clusterassignemnts, counts.df, by = "cellnames")
rownames(counts.df) <- counts.df$cellnames
counts.df$cellnames <- NULL
colnames(counts.df)[1] <- c("clusters")

counts.df[2,c(1:10)]


x <- createNeuronChat(t(as.matrix(counts.df[,2:(dim(counts.df)[2]-1)])),DB='human',group.by = counts.df$clusters);


neuronchat_Case1 <- run_NeuronChat(x,M=100)

save(neuronchat_Case1,file="neuronchat_Case1.RData")


Idents(Case2_tumor_normal) <- "manual_ident3"
counts.df <- Case2_tumor_normal@assays$Spatial@data %>% as.matrix %>% t %>% as.data.frame
counts.df <- tibble::rownames_to_column(counts.df, "cellnames")
clusterassignemnts <- data.frame(Case2_tumor_normal@active.ident)
clusterassignemnts <- tibble::rownames_to_column(clusterassignemnts, "cellnames")
counts.df <- merge(clusterassignemnts, counts.df, by = "cellnames")
rownames(counts.df) <- counts.df$cellnames
counts.df$cellnames <- NULL
colnames(counts.df)[1] <- c("clusters")

counts.df[2,c(1:10)]


x <- createNeuronChat(t(as.matrix(counts.df[,2:(dim(counts.df)[2]-1)])),DB='human',group.by = counts.df$clusters);


neuronchat_Case2 <- run_NeuronChat(x,M=100)

save(neuronchat_Case2,file="neuronchat_Case2.RData")

Idents(Case3_tumor_normal) <- "manual_ident3"
counts.df <- Case3_tumor_normal@assays$Spatial@data %>% as.matrix %>% t %>% as.data.frame
counts.df <- tibble::rownames_to_column(counts.df, "cellnames")
clusterassignemnts <- data.frame(Case3_tumor_normal@active.ident)
clusterassignemnts <- tibble::rownames_to_column(clusterassignemnts, "cellnames")
counts.df <- merge(clusterassignemnts, counts.df, by = "cellnames")
rownames(counts.df) <- counts.df$cellnames
counts.df$cellnames <- NULL
colnames(counts.df)[1] <- c("clusters")

counts.df[2,c(1:10)]


x <- createNeuronChat(t(as.matrix(counts.df[,2:(dim(counts.df)[2]-1)])),DB='human',group.by = counts.df$clusters);


neuronchat_Case3 <- run_NeuronChat(x,M=100)

save(neuronchat_Case3,file="neuronchat_Case3.RData")


Idents(Case4_tumor_normal) <- "manual_ident3"
counts.df <- Case4_tumor_normal@assays$Spatial@data %>% as.matrix %>% t %>% as.data.frame
counts.df <- tibble::rownames_to_column(counts.df, "cellnames")
clusterassignemnts <- data.frame(Case4_tumor_normal@active.ident)
clusterassignemnts <- tibble::rownames_to_column(clusterassignemnts, "cellnames")
counts.df <- merge(clusterassignemnts, counts.df, by = "cellnames")
rownames(counts.df) <- counts.df$cellnames
counts.df$cellnames <- NULL
colnames(counts.df)[1] <- c("clusters")

counts.df[2,c(1:10)]


x <- createNeuronChat(t(as.matrix(counts.df[,2:(dim(counts.df)[2]-1)])),DB='human',group.by = counts.df$clusters);


neuronchat_Case4 <- run_NeuronChat(x,M=100)

save(neuronchat_Case4,file="neuronchat_Case4.RData")





library(tidyr)
library(dplyr)

OBJ <-neuronchat_Case1@net

all_dfs <- list()
obj_names <- names(OBJ)

for (i in seq_along(obj_names)) {
  name_i <- obj_names[i]
  obj_i <- OBJ[[name_i]]
  

  if (is.matrix(obj_i)) {
    df <- as.data.frame(obj_i)
  } else if (is.data.frame(obj_i)) {
    df <- obj_i
  } else {

    next
  }
  

  if (nrow(df) > 0 && ncol(df) > 0) {
    df_long <- df %>%
      tibble::rownames_to_column("rownames") %>%
      pivot_longer(cols = -rownames, names_to = "colnames", values_to = "values") %>%
      mutate(name = name_i)
    
    all_dfs[[length(all_dfs) + 1]] <- df_long
  }
}


df_combined <- bind_rows(all_dfs)


colnames(df_combined)  <- c("sources", "targets", "contribution","name")

Case1_neuronLR<-df_combined

## repeat this for other patients
OBJ <-neuronchat_Case2@net
OBJ <-neuronchat_Case3@net
OBJ <-neuronchat_Case4@net

Case2_neuronLR<-df_combined
Case3_neuronLR<-df_combined
Case4_neuronLR<-df_combined
##

Case1_neuronLR<-
Case1_neuronLR %>%
  filter(contribution>0)%>%
  mutate(patient = "patient1")

Case2_neuronLR<-
  Case2_neuronLR %>%
  filter(contribution>0)%>%
  mutate(patient = "patient2")

Case3_neuronLR<-
  Case3_neuronLR %>%
  filter(contribution>0)%>%
  mutate(patient = "patient3")

Case4_neuronLR<-
  Case4_neuronLR %>%
  filter(contribution>0)%>%
  mutate(patient = "patient4")


save(Case1_neuronLR,file="Case1_neuronLR.RData")
save(Case2_neuronLR,file="Case2_neuronLR.RData")
save(Case3_neuronLR,file="Case3_neuronLR.RData")
save(Case4_neuronLR,file="Case4_neuronLR.RData")





Case1_neuronLR_tp_mucosal_from_tumor <-
  Case1_neuronLR %>%
  filter(targets == "tp_mucosal")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))


Case2_neuronLR_tp_mucosal_from_tumor <-
  Case2_neuronLR %>%
  filter(targets == "tp_mucosal")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))


Case3_neuronLR_tp_mucosal_from_tumor <-
  Case3_neuronLR %>%
  filter(targets == "tp_mucosal")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))


Case4_neuronLR_tp_mucosal_from_tumor <-
  Case4_neuronLR %>%
  filter(targets == "tp_mucosal")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))




neuronLR_tp_mucosal_from_tumor <-rbind(Case1_neuronLR_tp_mucosal_from_tumor,
                                       Case2_neuronLR_tp_mucosal_from_tumor,
                                       Case3_neuronLR_tp_mucosal_from_tumor,
                                       Case4_neuronLR_tp_mucosal_from_tumor)


neuronLR_tp_mucosal_from_tumor_sum <- data.frame(
  LR = names(tapply(neuronLR_tp_mucosal_from_tumor$contribution,
                    neuronLR_tp_mucosal_from_tumor$name,
                    sum)),
  sum = tapply(neuronLR_tp_mucosal_from_tumor$contribution,
               neuronLR_tp_mucosal_from_tumor$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(neuronLR_tp_mucosal_from_tumor_sum,"neuronLR_tp_mucosal_from_tumor_sum_.csv")



Case1_neuronLR_tm_from_tumor <-
  Case1_neuronLR %>%
  filter(targets == "tm")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))


Case2_neuronLR_tm_from_tumor <-
  Case2_neuronLR %>%
  filter(targets == "tm")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case3_neuronLR_tm_from_tumor <-
  Case3_neuronLR %>%
  filter(targets == "tm")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))



neuronLR_tm_from_tumor <-rbind(Case1_neuronLR_tm_from_tumor,
                               Case2_neuronLR_tm_from_tumor,
                               Case3_neuronLR_tm_from_tumor)


neuronLR_tm_from_tumor_sum <- data.frame(
  LR = names(tapply(neuronLR_tm_from_tumor$contribution,
                    neuronLR_tm_from_tumor$name,
                    sum)),
  sum = tapply(neuronLR_tm_from_tumor$contribution,
               neuronLR_tm_from_tumor$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(neuronLR_tm_from_tumor_sum,"neuronLR_tm_from_tumor_sum.csv")








Case3_neuronLR_tl_from_tumor <-
  Case3_neuronLR %>%
  filter(targets == "tl")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case4_neuronLR_tl_from_tumor <-
  Case4_neuronLR %>%
  filter(targets == "tl")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))



neuronLR_tl_from_tumor <-rbind(
  Case3_neuronLR_tl_from_tumor,
  Case4_neuronLR_tl_from_tumor)


neuronLR_tl_from_tumor_sum <- data.frame(
  LR = names(tapply(neuronLR_tl_from_tumor$contribution,
                    neuronLR_tl_from_tumor$name,
                    sum)),
  sum = tapply(neuronLR_tl_from_tumor$contribution,
               neuronLR_tl_from_tumor$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(neuronLR_tl_from_tumor_sum,"neuronLR_tl_from_tumor_sum.csv")










Case1_neuronLR_tp_other_from_tumor <-
  Case1_neuronLR %>%
  filter(targets == "tp_other")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case2_neuronLR_tp_other_from_tumor <-
  Case2_neuronLR %>%
  filter(targets == "tp_other")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case3_neuronLR_tp_other_from_tumor <-
  Case3_neuronLR %>%
  filter(targets == "tp_other")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case4_neuronLR_tp_other_from_tumor <-
  Case4_neuronLR %>%
  filter(targets == "tp_other")%>%
  filter(sources %in% c("tp_mucosal","tm","tl","tp_other"))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))




neuronLR_tp_other_from_tumor <-rbind(Case1_neuronLR_tp_other_from_tumor,
                                     Case2_neuronLR_tp_other_from_tumor,
                                     Case3_neuronLR_tp_other_from_tumor,
                                     Case4_neuronLR_tp_other_from_tumor)


neuronLR_tp_other_from_tumor_sum <- data.frame(
  LR = names(tapply(neuronLR_tp_other_from_tumor$contribution,
                    neuronLR_tp_other_from_tumor$name,
                    sum)),
  sum = tapply(neuronLR_tp_other_from_tumor$contribution,
               neuronLR_tp_other_from_tumor$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(neuronLR_tp_other_from_tumor_sum,"neuronLR_tp_other_from_tumor_sum.csv")







Case1_neuronLR_tp_mucosal_from_nontumor <-
  Case1_neuronLR %>%
  filter(targets == "tp_mucosal")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case2_neuronLR_tp_mucosal_from_nontumor <-
  Case2_neuronLR %>%
  filter(targets == "tp_mucosal")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case3_neuronLR_tp_mucosal_from_nontumor <-
  Case3_neuronLR %>%
  filter(targets == "tp_mucosal")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case4_neuronLR_tp_mucosal_from_nontumor <-
  Case4_neuronLR %>%
  filter(targets == "tp_mucosal")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))



neuronLR_tp_mucosal_from_nontumor <-rbind(Case1_neuronLR_tp_mucosal_from_nontumor,
                                          Case2_neuronLR_tp_mucosal_from_nontumor,
                                          Case3_neuronLR_tp_mucosal_from_nontumor,
                                          Case4_neuronLR_tp_mucosal_from_nontumor)


neuronLR_tp_mucosal_from_nontumor_sum <- data.frame(
  LR = names(tapply(neuronLR_tp_mucosal_from_nontumor$contribution,
                    neuronLR_tp_mucosal_from_nontumor$name,
                    sum)),
  sum = tapply(neuronLR_tp_mucosal_from_nontumor$contribution,
               neuronLR_tp_mucosal_from_nontumor$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(neuronLR_tp_mucosal_from_nontumor_sum,"neuronLR_tp_mucosal_from_nontumor_sum.csv")


Case1_neuronLR_tm_from_nontumor <-
  Case1_neuronLR %>%
  filter(targets == "tm")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case2_neuronLR_tm_from_nontumor <-
  Case2_neuronLR %>%
  filter(targets == "tm")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case3_neuronLR_tm_from_nontumor <-
  Case3_neuronLR %>%
  filter(targets == "tm")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))



neuronLR_tm_from_nontumor <-rbind(Case1_neuronLR_tm_from_nontumor,
                                  Case2_neuronLR_tm_from_nontumor,
                                  Case3_neuronLR_tm_from_nontumor)


neuronLR_tm_from_nontumor_sum <- data.frame(
  LR = names(tapply(neuronLR_tm_from_nontumor$contribution,
                    neuronLR_tm_from_nontumor$name,
                    sum)),
  sum = tapply(neuronLR_tm_from_nontumor$contribution,
               neuronLR_tm_from_nontumor$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(neuronLR_tm_from_nontumor_sum,"neuronLR_tm_from_nontumor_sum.csv")








Case3_neuronLR_tl_from_nontumor <-
  Case3_neuronLR %>%
  filter(targets == "tl")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case4_neuronLR_tl_from_nontumor <-
  Case4_neuronLR %>%
  filter(targets == "tl")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))



neuronLR_tl_from_nontumor <-rbind(
  Case3_neuronLR_tl_from_nontumor,
  Case4_neuronLR_tl_from_nontumor)


neuronLR_tl_from_nontumor_sum <- data.frame(
  LR = names(tapply(neuronLR_tl_from_nontumor$contribution,
                    neuronLR_tl_from_nontumor$name,
                    sum)),
  sum = tapply(neuronLR_tl_from_nontumor$contribution,
               neuronLR_tl_from_nontumor$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(neuronLR_tl_from_nontumor_sum,"neuronLR_tl_from_nontumor_sum.csv")










Case1_neuronLR_tp_other_from_nontumor <-
  Case1_neuronLR %>%
  filter(targets == "tp_other")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case2_neuronLR_tp_other_from_nontumor <-
  Case2_neuronLR %>%
  filter(targets == "tp_other")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case3_neuronLR_tp_other_from_nontumor <-
  Case3_neuronLR %>%
  filter(targets == "tp_other")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))

Case4_neuronLR_tp_other_from_nontumor <-
  Case4_neuronLR %>%
  filter(targets == "tp_other")%>%
  filter(!(sources %in% c("tp_mucosal","tm","tl","tp_other")))%>%
  mutate(contribution2 = contribution/sum(contribution))%>%
  group_by(name) %>%
  summarise(contribution = sum(contribution2, na.rm = TRUE))



neuronLR_tp_other_from_nontumor <-rbind(Case1_neuronLR_tp_other_from_nontumor,
                                        Case2_neuronLR_tp_other_from_nontumor,
                                        Case3_neuronLR_tp_other_from_nontumor,
                                        Case4_neuronLR_tp_other_from_nontumor)


neuronLR_tp_other_from_nontumor_sum <- data.frame(
  LR = names(tapply(neuronLR_tp_other_from_nontumor$contribution,
                    neuronLR_tp_other_from_nontumor$name,
                    sum)),
  sum = tapply(neuronLR_tp_other_from_nontumor$contribution,
               neuronLR_tp_other_from_nontumor$name,
               sum)) %>%
  arrange(desc(sum)) 

write.csv(neuronLR_tp_other_from_nontumor_sum,"neuronLR_tp_other_from_nontumor_sum.csv")

