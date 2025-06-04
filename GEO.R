
library(Biobase)
library(GEOquery)
library(limma)
library(data.table)
library(GEOquery)
library(Biobase)
library(edgeR)
library(dplyr)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(cowplot)
library(ggplot2)

urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE98894", "file=GSE98894_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&")
counts <- fread(path, header=TRUE)
counts <- as.matrix(counts[, -1, with=FALSE])
rownames(counts) <- fread(path, header=TRUE)[[1]]


gse <- getGEO("GSE98894", GSEMatrix=TRUE)[[1]]
pheno <- pData(gse)
head(pheno)  # Check available metadata
colnames(pheno)
characteristics_ch1 characteristics_ch1.1

levels(as.factor(pheno$characteristics_ch1))
"type: liver metastasis"      "type: lymph node metastasis" "type: mesenteric metastasis" "type: primary" 
levels(as.factor(pheno$characteristics_ch1.1))
"origin: pancreas"        "origin: rectum"          "origin: small intestine"


dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE)


common_samples <- intersect(colnames(logCPM), rownames(pheno))
logCPM <- logCPM[, common_samples]
pheno <- pheno[common_samples, ]
head(rownames(logCPM), 20)


entrez_ids <- rownames(logCPM)

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

logCPM <- logCPM[!is.na(gene_symbols), ]
rownames(logCPM) <- gene_symbols[rownames(logCPM)]




HTR1D_expr <- logCPM["HTR1D", ]
pheno$HTR1D_expr <- HTR1D_expr[rownames(pheno)]  
GJD2_expr <- logCPM["GJD2", ]
pheno$GJD2_expr <- GJD2_expr[rownames(pheno)]  
WNT4_expr <- logCPM["WNT4", ]
pheno$WNT4_expr <- WNT4_expr[rownames(pheno)]  
WNT2B_expr <- logCPM["WNT2B", ]
pheno$WNT2B_expr <- WNT2B_expr[rownames(pheno)]  
FZD8_expr <- logCPM["FZD8", ]
pheno$FZD8_expr <- FZD8_expr[rownames(pheno)]  
LRP5_expr <- logCPM["LRP5", ]
pheno$LRP5_expr <- LRP5_expr[rownames(pheno)]  
LRP6_expr <- logCPM["LRP6", ]
pheno$LRP6_expr <- LRP6_expr[rownames(pheno)]  
CXCR4_expr <- logCPM["CXCR4", ]
pheno$CXCR4_expr <- CXCR4_expr[rownames(pheno)]  
GDF11_expr <- logCPM["GDF11", ]
pheno$GDF11_expr <- GDF11_expr[rownames(pheno)]  
TGFBR1_expr <- logCPM["TGFBR1", ]
pheno$TGFBR1_expr <- TGFBR1_expr[rownames(pheno)]  
NRP1_expr <- logCPM["NRP1", ]
pheno$NRP1_expr <- NRP1_expr[rownames(pheno)]  
NRP2_expr <- logCPM["NRP2", ]
pheno$NRP2_expr <- NRP2_expr[rownames(pheno)]  
PLXNA2_expr <- logCPM["PLXNA2", ]
pheno$PLXNA2_expr <- PLXNA2_expr[rownames(pheno)]  
PLXNA3_expr <- logCPM["PLXNA3", ]
pheno$PLXNA3_expr <- PLXNA3_expr[rownames(pheno)]  
NCL_expr <- logCPM["NCL", ]
pheno$NCL_expr <- NCL_expr[rownames(pheno)]  
GIPR_expr <- logCPM["GIPR", ]
pheno$GIPR_expr <- GIPR_expr[rownames(pheno)]  
MDK_expr <- logCPM["MDK", ]
pheno$MDK_expr <- MDK_expr[rownames(pheno)]  
SEMA3B_expr <- logCPM["SEMA3B", ]
pheno$SEMA3B_expr <- SEMA3B_expr[rownames(pheno)]  
SEMA3G_expr <- logCPM["SEMA3G", ]
pheno$SEMA3G_expr <- SEMA3G_expr[rownames(pheno)]  
GUCY1A2_expr <- logCPM["GUCY1A2", ]
pheno$GUCY1A2_expr <- GUCY1A2_expr[rownames(pheno)]  
ERBB3_expr <- logCPM["ERBB3", ]
pheno$ERBB3_expr <- ERBB3_expr[rownames(pheno)]  
FGF9_expr <- logCPM["FGF9", ]
pheno$FGF9_expr <- FGF9_expr[rownames(pheno)]  
FGFR1_expr <- logCPM["FGFR1", ]
pheno$FGFR1_expr <- FGFR1_expr[rownames(pheno)]  
FGFR3_expr <- logCPM["FGFR3", ]
pheno$FGFR3_expr <- FGFR3_expr[rownames(pheno)]  
PDGFRB_expr <- logCPM["PDGFRB", ]
pheno$PDGFRB_expr <- PDGFRB_expr[rownames(pheno)]  





pheno$origin <- gsub("origin: ", "", pheno$characteristics_ch1.1)
pheno$type <- gsub("type: ", "", pheno$characteristics_ch1)


pheno$origin <- factor(pheno$origin)
pheno$type <- factor(pheno$type)

pheno$origin <- factor(pheno$origin,
                             levels = c("small intestine", "rectum", "pancreas"))

pheno$type <- factor(pheno$type,
                           levels = c("primary",
                                      "mesenteric metastasis",
                                      "lymph node metastasis",
                                      "liver metastasis"))

## there is only one case of mesenteric metastasis
pheno_noMM <-
  pheno %>%
  filter(type != "mesenteric metastasis")


  
  
  



# Plot with facets by origin
ggplot(pheno_noMM, aes(x = type, y = HTR1D_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "HTR1D Expression by Tumor Type and Origin",
       x = "Tumor Type",
       y = "logCPM (HTR1D)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(HTR1D_expr ~ type)$p.value)

# Plot with facets by origin
ggplot(pheno_noMM, aes(x = type, y = WNT4_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "WNT4",
       x = "Tumor Type",
       y = "logCPM (WNT4)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(WNT4_expr ~ type)$p.value)

# Plot with facets by origin
ggplot(pheno_noMM, aes(x = type, y = WNT2B_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "WNT2B",
       x = "Tumor Type",
       y = "logCPM (WNT2B)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(WNT2B_expr ~ type)$p.value)

# Plot with facets by origin
ggplot(pheno_noMM, aes(x = type, y = FZD8_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "FZD8",
       x = "Tumor Type",
       y = "logCPM (FZD8)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(FZD8_expr ~ type)$p.value)


# Plot with facets by origin
ggplot(pheno_noMM, aes(x = type, y = LRP5_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "LRP5",
       x = "Tumor Type",
       y = "logCPM (LRP5)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(LRP5_expr ~ type)$p.value)


# Plot with facets by origin
ggplot(pheno_noMM, aes(x = type, y = LRP6_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "LRP6",
       x = "Tumor Type",
       y = "logCPM (LRP6)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(LRP6_expr ~ type)$p.value)

# Plot with facets by origin

GEOCXCR4<-
ggplot(pheno_noMM, aes(x = type, y = CXCR4_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "CXCR4",
       x = "Tumor Type",
       y = "logCPM (CXCR4)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(CXCR4_expr ~ type)$p.value)
1 small intestine  0.0125
2 rectum           0.0175
3 pancreas         0.384 
ggsave(filename="GEOCXCR4.svg",width=5, height=5,plot=GEOCXCR4)
ggsave(filename="GEOCXCR4_052825.svg",width=5, height=4,plot=GEOCXCR4)

GEOWNT4<-
  ggplot(pheno_noMM, aes(x = type, y = WNT4_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "WNT4",
       x = "Tumor Type",
       y = "logCPM (WNT4)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(WNT4_expr ~ type)$p.value)
1 small intestine 0.0271 
2 rectum          0.837  
3 pancreas        0.00487
ggsave(filename="GEOWNT4.svg",width=5, height=5,plot=GEOWNT4)
ggsave(filename="GEOWNT4_052825.svg",width=5, height=4,plot=GEOWNT4)



GEONCL<-
ggplot(pheno, aes(x = origin, y = NCL_expr, fill = origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "NCL",
       x = "Tumor Origin",
       y = "logCPM (NCL)") +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

kruskal.test(NCL_expr ~ origin, data = pheno)
Kruskal-Wallis chi-squared = 3.9057, df = 2, p-value = 0.1419
ggsave(filename="GEONCL.svg",width=2, height=4,plot=GEONCL)
ggsave(filename="GEONCL_052825.svg",width=2, height=3,plot=GEONCL)


GEOMDK<-
  ggplot(pheno, aes(x = origin, y = MDK_expr, fill = origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "MDK",
       x = "Tumor Origin",
       y = "logCPM (MDK)") +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

kruskal.test(MDK_expr ~ origin, data = pheno)
Kruskal-Wallis chi-squared = 16.046, df = 2, p-value = 0.0003278
ggsave(filename="GEOMDK.svg",width=2, height=4,plot=GEOMDK)
ggsave(filename="GEOMDK_052825.svg",width=2, height=3,plot=GEOMDK)


GEOGJD2<-
  ggplot(pheno, aes(x = origin, y = GJD2_expr, fill = origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "GJD2",
       x = "Tumor Origin",
       y = "logCPM (GJD2)") +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

kruskal.test(GJD2_expr ~ origin, data = pheno)
Kruskal-Wallis chi-squared = 14.675, df = 2, p-value = 0.0006507
ggsave(filename="GEOGJD2.svg",width=2, height=4,plot=GEOGJD2)
ggsave(filename="GEOGJD2_052825.svg",width=2, height=3,plot=GEOGJD2)


GEOHTR1D<-
  ggplot(pheno, aes(x = origin, y = HTR1D_expr, fill = origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "HTR1D",
       x = "Tumor Origin",
       y = "logCPM (HTR1D)") +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

kruskal.test(HTR1D_expr ~ origin, data = pheno)
Kruskal-Wallis chi-squared = 49.227, df = 2, p-value = 2.044e-11
ggsave(filename="GEOHTR1D.svg",width=2, height=4,plot=GEOHTR1D)
ggsave(filename="GEOHTR1D_052825.svg",width=2, height=3,plot=GEOHTR1D)



ggplot(pheno_noMM, aes(x = type, y = MDK_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "MDK",
       x = "Tumor Type",
       y = "logCPM (MDK)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggplot(pheno_noMM, aes(x = type, y = GDF11_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "GDF11",
       x = "Tumor Type",
       y = "logCPM (GDF11)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(GDF11_expr ~ type)$p.value)
1 small intestine   0.282
2 rectum            0.738
3 pancreas          0.805



ggplot(pheno_noMM, aes(x = type, y = SEMA3B_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "SEMA3B",
       x = "Tumor Type",
       y = "logCPM (SEMA3B)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(SEMA3B_expr ~ type)$p.value)

ggplot(pheno_noMM, aes(x = type, y = SEMA3G_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "SEMA3G",
       x = "Tumor Type",
       y = "logCPM (SEMA3G)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(SEMA3G_expr ~ type)$p.value)

ggplot(pheno_noMM, aes(x = type, y = GUCY1A2_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "GUCY1A2",
       x = "Tumor Type",
       y = "logCPM (GUCY1A2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(GUCY1A2_expr ~ type)$p.value)




ggplot(pheno_noMM, aes(x = type, y = FGF9_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "FGF9",
       x = "Tumor Type",
       y = "logCPM (FGF9)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(FGF9_expr ~ type)$p.value)

ggplot(pheno_noMM, aes(x = type, y = FGFR1_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "FGFR1",
       x = "Tumor Type",
       y = "logCPM (FGFR1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(FGFR1_expr ~ type)$p.value)

ggplot(pheno_noMM, aes(x = type, y = FGFR3_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "FGFR3",
       x = "Tumor Type",
       y = "logCPM (FGFR3)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(FGFR3_expr ~ type)$p.value)

ggplot(pheno_noMM, aes(x = type, y = PDGFRB_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "PDGFRB",
       x = "Tumor Type",
       y = "logCPM (PDGFRB)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(PDGFRB_expr ~ type)$p.value)




ggplot(pheno_noMM, aes(x = type, y = ERBB3_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "ERBB3",
       x = "Tumor Type",
       y = "logCPM (ERBB3)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pheno_noMM %>%
  group_by(origin) %>%
  summarise(p_value = kruskal.test(ERBB3_expr ~ type)$p.value)



ggplot(pheno_noMM, aes(x = type, y = TGFBR1_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "TGFBR1",
       x = "Tumor Type",
       y = "logCPM (TGFBR1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")



ggplot(pheno_noMM, aes(x = type, y = NRP1_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "NRP1",
       x = "Tumor Type",
       y = "logCPM (NRP1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggplot(pheno_noMM, aes(x = type, y = NRP2_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "NRP2",
       x = "Tumor Type",
       y = "logCPM (NRP2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggplot(pheno_noMM, aes(x = type, y = PLXNA2_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "PLXNA2",
       x = "Tumor Type",
       y = "logCPM (PLXNA2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


ggplot(pheno_noMM, aes(x = type, y = PLXNA3_expr, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ origin) +
  labs(title = "PLXNA3",
       x = "Tumor Type",
       y = "logCPM (PLXNA3)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
