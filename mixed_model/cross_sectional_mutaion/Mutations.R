# 1 Load Libraries & Import Data ---------------------------------------------------------
library(data.table)
library(tximport)
library(DESeq2)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(stringr)
library(ggpubr)
library(ggrepel)
library(splitstackshape)
library(org.Hs.eg.db)


# * 1.1 Import pheno  ------------------------

#load phenotype 
pheno <- read.table("/phenotype.txt", header = T, stringsAsFactors = F)
# Pull out control and PD individuals (0 and 1 in Status_time column)
pheno <- subset(pheno, pheno$Status_time %in% c(0, 1))

# select for last visit only
pheno <- pheno %>% group_by(PATNO) %>% mutate(last_visit = (Time == max(Time)))
pheno <- subset(pheno, pheno$last_visit == T)


# * 1.4 Read in Matrix ------------------------
countMatrix <- read.csv("/countmatrix.csv", row.names = 1, header = T)


#Removing sample with low reads
countMatrixclean<- countMatrix[(rowCounts(countMatrix[,-1]<10) < round(0.9*dim(countMatrix[,-1])[2])),]
pheno <- pheno[pheno$FILE_NAME %in% colnames(countMatrixclean),]

# Extract the sample IDs from the phenotype data
sample_ids <- pheno$FILE_NAME

# Subset the countMatrix to include only the columns in sample_ids
countMatrixfinal<- countMatrixclean[, sample_ids]
countMatrixfinal<-floor(countMatrixfinal)
#Modify pheno for mutations
pheno <- pheno %>%
  mutate(Mutation = ifelse(Status_time == "0" & is.na(Mutation), "control", Mutation),
         Mutation = ifelse(Status_time == "1" & is.na(Mutation), "non-carrier", Mutation))

#Modify count matrix for mutation
# Extract the sample IDs from the phenotype data
sample_ids <- pheno$FILE_NAME

# Subset the countMatrix to include only the columns in sample_ids
countMatrixfinal<- countMatrixclean[, sample_ids]
countMatrixfinal<-floor(countMatrixfinal)
#First we need to modify our phenotype data
pheno$phase <- ifelse(grepl("Phase1", pheno$FILE_NAME), "Phase1", ifelse(grepl("Phase2", pheno$FILE_NAME), "Phase2", NA))
pheno$phase<- factor(pheno$phase)


# 3 Differential Expression Analysis -------------------------------------------

# * 3.1 Basic DESeq ------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countMatrixfinal, colData = pheno, design = ~ factor(Gender) + AgeVisit + Mutation + factor(Status_time))
#keep <- rowSums(counts(dds) >= ) matrixStats::rowCounts(quantData$counts<10) < round(0.9*dim(quantData$counts)
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)

dds.res <- results(dds.de, tidy = T)
write.csv(dds.res, "/degs_mutation.csv")
#Make a copy of dds.res
dds.res2 <- dds.res
#Remove the dot and proceeding numbers
dds.res2$row <- sub("\\..*","",dds.res2$row)
dds.res2$row <- mapIds(org.Hs.eg.db, keys = dds.res2$row, column = "SYMBOL", keytype = "ENSEMBL")

summary(dds.res)
dds.sig <- subset(dds.res2, dds.res2$padj < 0.05)

# plotting the volcano plot
DE <- dds.res2[order(dds.res2$padj), ]
volcano.name <- "PPMI Linear RNA Expression Last Visit ca/co"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
plotDE <- DE %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up",
                                              log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns"))
plotDE <- plotDE[order(plotDE$pvalue), ]
#Get the 9 circ degs
result_9circ<- read.table("/9circ.txt", header = T) 
result_9circ$linear <- sub("\\..*","",result_9circ$linear)
result_9circ$linear <- mapIds(org.Hs.eg.db, keys = result_9circ$linear, column = "SYMBOL", keytype = "ENSEMBL")
circ_names<-result_9circ$linear

nine_circ<- plotDE[plotDE$row %in% circ_names, ]
cols <- densCols(plotDE$log2FoldChange, plotDE$padj)
cols[plotDE$gene_type=='up']<-"#E69F00"
cols[plotDE$gene_type=='down']<-"#56B4E9"
cols[plotDE$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- plotDE[plotDE$pvalue<=plim & plotDE$gene_type!='ns',]

top_hits_num <- 10
sig <- sig[1:top_hits_num, ]



ggplot(plotDE, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 15))+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = row), force = 2, nudge_y = 1) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")


