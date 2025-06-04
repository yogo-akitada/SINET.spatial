

library(data.table)
library(dplyr)
library(ggplot2)


## Load downloaded tissue-specific TPM files from GTEx Portal 
## https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
ileum <- fread("gene_tpm_v10_small_intestine_terminal_ileum.gct.gz", skip = 2)
colon <- fread("gene_tpm_v10_colon_sigmoid.gct.gz", skip = 2)
pancreas <- fread("gene_tpm_v10_pancreas.gct.gz", skip = 2)
braincortex <- fread("gene_tpm_v10_brain_cortex.gct.gz", skip = 2)
lung <- fread("gene_tpm_v10_lung.gct.gz", skip = 2)
visceral_adipose <- fread("gene_tpm_v10_adipose_visceral_omentum.gct.gz", skip = 2)
liver <- fread("gene_tpm_v10_liver.gct.gz", skip = 2)



rownames(ileum) <- ileum$Name
rownames(colon) <- colon$Name
rownames(pancreas) <- pancreas$Name
rownames(braincortex) <- braincortex$Name
rownames(lung) <- lung$Name
rownames(visceral_adipose) <- visceral_adipose$Name
rownames(liver) <- liver$Name



PTN_ileum <- as.numeric(ileum[ileum$Description == "PTN", -c(1,2)])
PTN_colon <- as.numeric(colon[colon$Description == "PTN", -c(1,2)])
PTN_pancreas <- as.numeric(pancreas[pancreas$Description == "PTN", -c(1,2)])
PTN_braincortex <- as.numeric(braincortex[braincortex$Description == "PTN", -c(1,2)])


df <- data.frame(
  Expression = c(log2(PTN_ileum + 1), log2(PTN_colon + 1), log2(PTN_pancreas + 1),log2(PTN_braincortex + 1)),
  Tissue = factor(
    c(rep("Small Intestine", length(PTN_ileum)),
      rep("Colon", length(PTN_colon)),
      rep("Pancreas", length(PTN_pancreas)),
      rep("Brain Cortex", length(PTN_braincortex))),
    levels = c("Small Intestine", "Colon", "Pancreas","Brain Cortex")
  )
)


PTN_GTEx <-
ggplot(df, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.3, alpha = 0.1, size=1) +
  labs(title = "PTN",
       y = "log2(TPM + 1)", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


ggsave(filename="PTN_GTEx.svg",width=3, height=4,plot=PTN_GTEx)


median(PTN_braincortex) /median(PTN_ileum) 
[1] 30.6785




GCG_ileum <- as.numeric(ileum[ileum$Description == "GCG", -c(1,2)])
GCG_colon <- as.numeric(colon[colon$Description == "GCG", -c(1,2)])
GCG_pancreas <- as.numeric(pancreas[pancreas$Description == "GCG", -c(1,2)])
GCG_braincortex <- as.numeric(braincortex[braincortex$Description == "GCG", -c(1,2)])


df <- data.frame(
  Expression = c(log2(GCG_ileum + 1), log2(GCG_colon + 1), log2(GCG_pancreas + 1)),
  Tissue = factor(
    c(rep("Small Intestine", length(GCG_ileum)),
      rep("Colon", length(GCG_colon)),
      rep("Pancreas", length(GCG_pancreas))),
    levels = c("Small Intestine", "Colon", "Pancreas")
  )
)
GCG_GTEx <-
  ggplot(df, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.3, alpha = 0.1, size=1) +
  labs(title = "GCG",
       y = "log2(TPM + 1)", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


ggsave(filename="GCG_GTEx.svg",width=2, height=4,plot=GCG_GTEx)





GIP_ileum <- as.numeric(ileum[ileum$Description == "GIP", -c(1,2)])
GIP_colon <- as.numeric(colon[colon$Description == "GIP", -c(1,2)])
GIP_pancreas <- as.numeric(pancreas[pancreas$Description == "GIP", -c(1,2)])
GIP_braincortex <- as.numeric(braincortex[braincortex$Description == "GIP", -c(1,2)])


df <- data.frame(
  Expression = c(log2(GIP_ileum + 1), log2(GIP_colon + 1), log2(GIP_pancreas + 1)),
  Tissue = factor(
    c(rep("Small Intestine", length(GIP_ileum)),
      rep("Colon", length(GIP_colon)),
      rep("Pancreas", length(GIP_pancreas))),
    levels = c("Small Intestine", "Colon", "Pancreas")
  )
)
GIP_GTEx <-
  ggplot(df, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.3, alpha = 0.1, size=1) +
  labs(title = "GIP",
       y = "log2(TPM + 1)", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


ggsave(filename="GIP_GTEx.svg",width=2, height=4,plot=GIP_GTEx)





FGF7_ileum <- as.numeric(ileum[ileum$Description == "FGF7", -c(1,2)])
FGF7_colon <- as.numeric(colon[colon$Description == "FGF7", -c(1,2)])
FGF7_pancreas <- as.numeric(pancreas[pancreas$Description == "FGF7", -c(1,2)])
FGF7_lung <- as.numeric(lung[lung$Description == "FGF7", -c(1,2)])


df <- data.frame(
  Expression = c(log2(FGF7_ileum + 1), log2(FGF7_colon + 1), log2(FGF7_pancreas + 1), log2(FGF7_lung + 1)),
  Tissue = factor(
    c(rep("Small Intestine", length(FGF7_ileum)),
      rep("Colon", length(FGF7_colon)),
      rep("Pancreas", length(FGF7_pancreas)),
      rep("Lung", length(FGF7_lung))),
    levels = c("Small Intestine", "Colon", "Pancreas","Lung")
  )
)
FGF7_GTEx <-
  ggplot(df, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.3, alpha = 0.1, size=1) +
  labs(title = "FGF7",
       y = "log2(TPM + 1)", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


ggsave(filename="FGF7_GTEx.svg",width=2, height=4,plot=FGF7_GTEx)







CXCL12_ileum <- as.numeric(ileum[ileum$Description == "CXCL12", -c(1,2)])
CXCL12_colon <- as.numeric(colon[colon$Description == "CXCL12", -c(1,2)])
CXCL12_pancreas <- as.numeric(pancreas[pancreas$Description == "CXCL12", -c(1,2)])
CXCL12_visceral_adipose <- as.numeric(visceral_adipose[visceral_adipose$Description == "CXCL12", -c(1,2)])


df <- data.frame(
  Expression = c(log2(CXCL12_ileum + 1), log2(CXCL12_colon + 1), log2(CXCL12_pancreas + 1), log2(CXCL12_visceral_adipose + 1)),
  Tissue = factor(
    c(rep("Small Intestine", length(CXCL12_ileum)),
      rep("Colon", length(CXCL12_colon)),
      rep("Pancreas", length(CXCL12_pancreas)),
      rep("Visceral adipose", length(CXCL12_visceral_adipose))),
    levels = c("Small Intestine", "Colon", "Pancreas","Visceral adipose")
  )
)
CXCL12_GTEx <-
  ggplot(df, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.3, alpha = 0.1, size=1) +
  labs(title = "CXCL12",
       y = "log2(TPM + 1)", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


ggsave(filename="CXCL12_GTEx.svg",width=3, height=4,plot=CXCL12_GTEx)


median(CXCL12_visceral_adipose) /median(CXCL12_ileum) 
[1] 3.385357



HTR1D_ileum <- as.numeric(ileum[ileum$Description == "HTR1D", -c(1,2)])
HTR1D_colon <- as.numeric(colon[colon$Description == "HTR1D", -c(1,2)])
HTR1D_pancreas <- as.numeric(pancreas[pancreas$Description == "HTR1D", -c(1,2)])
HTR1D_liver <- as.numeric(liver[liver$Description == "HTR1D", -c(1,2)])


df <- data.frame(
  Expression = c(log2(HTR1D_ileum + 1), log2(HTR1D_colon + 1), log2(HTR1D_pancreas + 1),log2(HTR1D_liver + 1)),
  Tissue = factor(
    c(rep("Small Intestine", length(HTR1D_ileum)),
      rep("Colon", length(HTR1D_colon)),
      rep("Pancreas", length(HTR1D_pancreas)),
      rep("Liver", length(HTR1D_liver))),
    levels = c("Small Intestine", "Colon", "Pancreas","Liver")
  )
)


HTR1D_GTEx <-
  ggplot(df, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.3, alpha = 0.1, size=1) +
  labs(title = "HTR1D",
       y = "log2(TPM + 1)", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
