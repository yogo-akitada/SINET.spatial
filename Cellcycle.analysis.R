

s.genes <- cc.genes$s.genes
[1] "MCM5"     "PCNA"     "TYMS"     "FEN1"     "MCM2"     "MCM4"     "RRM1"     "UNG"      "GINS2"    "MCM6"     "CDCA7"    "DTL"     
[13] "PRIM1"    "UHRF1"    "MLF1IP"   "HELLS"    "RFC2"     "RPA2"     "NASP"     "RAD51AP1" "GMNN"     "WDR76"    "SLBP"     "CCNE2"   
[25] "UBR7"     "POLD3"    "MSH2"     "ATAD2"    "RAD51"    "RRM2"     "CDC45"    "CDC6"     "EXO1"     "TIPIN"    "DSCC1"    "BLM"     
[37] "CASP8AP2" "USP1"     "CLSPN"    "POLA1"    "CHAF1B"   "BRIP1"    "E2F8"    
g2m.genes <- cc.genes$g2m.genes
[1] "HMGB2"   "CDK1"    "NUSAP1"  "UBE2C"   "BIRC5"   "TPX2"    "TOP2A"   "NDC80"   "CKS2"    "NUF2"    "CKS1B"   "MKI67"   "TMPO"   
[14] "CENPF"   "TACC3"   "FAM64A"  "SMC4"    "CCNB2"   "CKAP2L"  "CKAP2"   "AURKB"   "BUB1"    "KIF11"   "ANP32E"  "TUBB4B"  "GTSE1"  
[27] "KIF20B"  "HJURP"   "CDCA3"   "HN1"     "CDC20"   "TTK"     "CDC25C"  "KIF2C"   "RANGAP1" "NCAPD2"  "DLGAP5"  "CDCA2"   "CDCA8"  
[40] "ECT2"    "KIF23"   "HMMR"    "AURKA"   "PSRC1"   "ANLN"    "LBR"     "CKAP5"   "CENPE"   "CTCF"    "NEK2"    "G2E3"    "GAS2L3" 
[53] "CBX5"    "CENPA" 
Case1_tumor_normal <- CellCycleScoring(object = Case1_tumor_normal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Case2_tumor_normal <- CellCycleScoring(object = Case2_tumor_normal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Case3_tumor_normal <- CellCycleScoring(object = Case3_tumor_normal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Case4_tumor_normal <- CellCycleScoring(object = Case4_tumor_normal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)




plot_data_raw <- Case1_tumor_normal@meta.data %>%
  group_by(manual_ident3, Phase) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  complete(manual_ident3, Phase, fill = list(n = 0))


totals <- plot_data_raw %>%
  group_by(manual_ident3) %>%
  dplyr::summarise(total_n = sum(n), .groups = "drop")

plot_data <- plot_data_raw %>%
  left_join(totals, by = "manual_ident3") %>%
  mutate(p = n / total_n)


cellcycleCase1_Ctrl_3t<-
ggplot(plot_data, aes(x = manual_ident3, y = p, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = percent(p, accuracy = 1)), 
            position = position_stack(vjust = 0.5), 
            size = 3, show.legend = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  ylab("% of cells in each subgroup") +
  theme_classic()+
  scale_x_discrete(limits = c( "m_superficial",
                               "m_base",
                               
                               "tp_mucosal",
                               "tm",
                               "tp_deep"))

cellcycleCase1_Ctrl_3t


ggsave(filename="cellcycleCase1_Ctrl_3t.svg",width=5, height=4,plot=cellcycleCase1_Ctrl_3t)


plot_data_raw <- Case2_tumor_normal@meta.data %>%
  group_by(manual_ident3, Phase) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  complete(manual_ident3, Phase, fill = list(n = 0))


totals <- plot_data_raw %>%
  group_by(manual_ident3) %>%
  dplyr::summarise(total_n = sum(n), .groups = "drop")

plot_data <- plot_data_raw %>%
  left_join(totals, by = "manual_ident3") %>%
  mutate(p = n / total_n)


cellcycleCase2_Ctrl_3t<-
  ggplot(plot_data, aes(x = manual_ident3, y = p, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = percent(p, accuracy = 1)), 
            position = position_stack(vjust = 0.5), 
            size = 3, show.legend = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  ylab("% of cells in each subgroup") +
  theme_classic()+
  scale_x_discrete(limits = c( "m_superficial_mid_crypt",
                               "m_base",
                               
                               "tp_mucosal",
                               "tm",
                               "tp_deep"))

cellcycleCase2_Ctrl_3t


ggsave(filename="cellcycleCase2_Ctrl_3t.svg",width=5, height=4,plot=cellcycleCase2_Ctrl_3t)


plot_data_raw <- Case3_tumor_normal@meta.data %>%
  group_by(manual_ident3, Phase) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  complete(manual_ident3, Phase, fill = list(n = 0))


totals <- plot_data_raw %>%
  group_by(manual_ident3) %>%
  dplyr::summarise(total_n = sum(n), .groups = "drop")

plot_data <- plot_data_raw %>%
  left_join(totals, by = "manual_ident3") %>%
  mutate(p = n / total_n)


cellcycleCase3_Ctrl_3t<-
  ggplot(plot_data, aes(x = manual_ident3, y = p, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = percent(p, accuracy = 1)), 
            position = position_stack(vjust = 0.5), 
            size = 3, show.legend = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  ylab("% of cells in each subgroup") +
  theme_classic()+
  scale_x_discrete(limits = c( "mp_inner_peritumoral",
                               "m_base",
                               
                               "tp_mucosal",
                               "tm",
                               "tl",
                               "tp_deep"))

cellcycleCase3_Ctrl_3t


ggsave(filename="cellcycleCase3_Ctrl_3t.svg",width=5, height=4,plot=cellcycleCase3_Ctrl_3t)


plot_data_raw <- Case4_tumor_normal@meta.data %>%
  group_by(manual_ident3, Phase) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  complete(manual_ident3, Phase, fill = list(n = 0))


totals <- plot_data_raw %>%
  group_by(manual_ident3) %>%
  dplyr::summarise(total_n = sum(n), .groups = "drop")

plot_data <- plot_data_raw %>%
  left_join(totals, by = "manual_ident3") %>%
  mutate(p = n / total_n)


cellcycleCase4_Ctrl_3t<-
  ggplot(plot_data, aes(x = manual_ident3, y = p, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = percent(p, accuracy = 1)), 
            position = position_stack(vjust = 0.5), 
            size = 3, show.legend = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  ylab("% of cells in each subgroup") +
  theme_classic()+
  scale_x_discrete(limits = c( "m_mid_crypt",
                               "m_base",
                               
                               "tp_mucosal",
                               "tl",
                               "tp_deep"))

cellcycleCase4_Ctrl_3t


ggsave(filename="cellcycleCase4_Ctrl_3t.svg",width=5, height=4,plot=cellcycleCase4_Ctrl_3t)
