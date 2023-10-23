library(data.table)
library(tximport)
library(DESeq2)
library(stringr)
library(biomaRt)
library("lmerTest")
library(reshape2)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(org.Hs.eg.db)
#load phenotype 
pheno <- read.table('phenotype.file.txt', header = T, stringsAsFactors = F)

# Base directory where subdirectories with quant.genes.sf files are located
base_directory <- "Base/directory"

# List all subdirectories
subdirectories <- list.dirs(base_directory, full.names = TRUE, recursive = TRUE)

# Initialize empty list to store data frames and row names
data_frames <- list()
row_names <- NULL

# Loop through each subdirectory
for (subdir in subdirectories) {
  quant_file <- file.path(subdir, "quant.genes.sf")
  
  if (file.exists(quant_file)) {
    individual_name <- basename(subdir)   # Get the name of the individual
    
    # Read the quant.genes.sf file
    quant_data <- fread(quant_file)
    
    # Extract the NumReads column and add it to the data frame list
    data_frames[[individual_name]] <- quant_data$NumReads
    
  }
}

# Combine data frames into a single matrix
countMatrix <- do.call(cbind, data_frames)

# Set column names to individual names without .salmon suffix and change - to .
colnames(countMatrix) <- sub("\\.salmon$", "", names(data_frames))
colnames(countMatrix) <- gsub("-", ".", colnames(countMatrix))


# Set row names to Names column from quant_data
row.names(countMatrix)<-quant_data$Name

# Print the count matrix
print(countMatrix)

#save the count matrix for easier loading
write.csv(countMatrix,"save/count/matrix.csv")

#Load count matrix
count_matrix<-read.csv("load/count/matrix.csv", row.names = 1)

#Removing sample with low reads from coun matrix and phenotype file.
countMatrixclean<- count_matrix[(rowCounts(count_matrix[,-1]<10) < round(0.9*dim(count_matrix[,-1])[2])),]
pheno <- pheno[pheno$sample_id %in% colnames(countMatrixclean),]

# Extract the sample IDs from the phenotype data
sample_ids <- pheno$sample_id

# Subset the countMatrix to include only the columns in sample_ids
countMatrixfinal<- countMatrixclean[, sample_ids]
countMatrixfinal<-floor(countMatrixfinal)

#save the cleaned count matrix for easier loading
write.csv(countMatrixfinal,"save/count/matrix.csv")
count_matrix<-read.csv("load/count/matrix.csv", row.names = 1)

#count normalization
rnaDDS <- DESeqDataSetFromMatrix(countData =countMatrixfinal , colData = pheno, design = ~ 1)
vstDDS<-vst(rnaDDS)
vst_countmatrix<-assay(vstDDS)

#save the normalized count matrix for easier loading
write.table(vst_countmatrix,"save/count/matrix.csv")
vstDDS<-read.table("load/count/matrix.csv", header = T)



# plotPCA(vst.all, intgroup="sex")
pcaData <- plotPCA(vstDDS, intgroup="Status",ntop=500 ,returnData = T)
write.csv(pcaData, "save/PCA.csv", row.names = T)

rnaPercentVar <- round(100 * attr(pcaData, "percentVar"))
pc1.sd <- sd(pcaData$PC1)
pc1.mean <- mean(pcaData$PC1)
pc2.sd <- sd(pcaData$PC2)
pc2.mean <- mean(pcaData$PC2)
dist = 3 # choose of 3 sd or 2 sd
range.pc1 <- c(pc1.mean - dist*pc1.sd, pc1.mean + dist*pc1.sd)
range.pc2 <- c(pc2.mean - dist*pc2.sd, pc2.mean + dist*pc2.sd)
# check PCA plots with SD borders
panel <- c("#E69F00", "#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plotStatus <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) + geom_point() +
  xlab(paste0("PC1: ", rnaPercentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", rnaPercentVar[2], "% variance")) + 
  geom_hline(yintercept = range.pc2, linetype = 'dashed') +
  geom_vline(xintercept = range.pc1, linetype = 'dashed') +
  scale_color_manual(values=panel) + theme_bw()



#Generation of Longtable

# load normalized count matrix
counts<-read.table("load/count/matrix.csv",header = T, stringsAsFactors = F)
counts<-as.matrix(counts)
# load cleaned phenotype data
pheno <- read.table('load/phenotype.csv', header = T, stringsAsFactors = F)

#subset phenotype file 
pheno<- pheno[colnames(counts) %in% pheno$sample_id, ]
# check time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12

# select covariates in pheno that are needed for mixed models
mixed_model_pheno <- pheno %>% select(sample_id, sex, status = Pt.Category, age_at_baseline, time_in_years, PATNO = participant_id)
mixed_model_pheno$PATNO <- str_replace_all(mixed_model_pheno$PATNO, "PD-", "")
mixed_model_pheno <- mixed_model_pheno[order(mixed_model_pheno$PATNO, mixed_model_pheno$time_in_years), ]
tmp <- mixed_model_pheno %>% group_by(PATNO) %>% mutate(visit_counts = length(time_in_years))
table(tmp$visit_counts)

mixed_model_count_matirx <- data.frame(t(counts))
mixed_model_count_matirx$sample_id <- rownames(mixed_model_count_matirx)

mixed_model_longtable <- merge(mixed_model_pheno, mixed_model_count_matirx, by = "sample_id")

# create the long table
mixed_model_longtable <- melt(mixed_model_longtable, id.vars = c('sample_id', "PATNO", 'sex', 'age_at_baseline', 'status', "time_in_years"))

colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "variable"] <- "RNA"
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "value"] <- "counts"
str(mixed_model_longtable)
mixed_model_longtable$RNA <- as.character(mixed_model_longtable$RNA)

mixed_model_longtable<- mixed_model_longtable[order(mixed_model_longtable$PATNO), ]

# remove 0 counts
mixed_model_longtable[mixed_model_longtable$counts == 0, 'counts'] <- NA
mixed_model_longtable <- na.omit(mixed_model_longtable)

# add baseline counts 
mixed_model_longtable$counts_at_baseline <- -9

mixed_model_longtable <- mixed_model_longtable %>% group_by(PATNO, RNA) %>% mutate(counts_at_baseline = counts[which(time_in_years == min(time_in_years))])

any(mixed_model_longtable$counts_at_baseline == -9) # should be FALSE

# checking
## mixed_model_longtable %>% filter(PATNO == "PDAA503EF5", circRNA == "circASXL1")

patno <- unique(mixed_model_longtable$PATNO)
linear_names <- unique(mixed_model_longtable$RNA)

rand_linear <- linear_names[sample(1:length(linear_names), 1)]
rand_sample <- patno[sample(1:length(patno), 1)]
mixed_model_longtable %>% filter(PATNO == rand_sample, RNA == rand_linear)

write.table(mixed_model_longtable, 'save/longtable.txt', row.names = F, col.names = T)


#MIXED MODEL WITH LONGTABLE
# load in the Longtable
data<- read.table("load/longtable.txt", header = T, stringsAsFactors = F)
data <- data[order(data$PATNO, data$RNA),]
data <- group_by(data, PATNO, RNA)
data <- mutate(data, visit_counts = length(time_in_years))



# filtered with as least two visits
atleast_2_visit <- subset(data, data$visit_counts > 1)

loop_linear<- as.character(unique(atleast_2_visit$RNA))

full_result_table <- NULL
for (one_linear in loop_linear){
  if (which(loop_linear == one_linear) %% 100 == 0){print(paste("Progress:", one_linear, which(loop_linear == one_linear)))}
  tmp <- subset(atleast_2_visit, atleast_2_visit$RNA == one_linear)
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + sex + status + time_in_years +  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', 'conditioncontrol', 'years', 'conditioncontrol:years')
  res$covar <- rownames(res)
  res$linear <- one_linear
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}

# save result table
full_result_table %>% write.table("save/reuslts.txt", quote = F, col.names = T, row.names = F, sep = '\t')


#Read in the results table
full_result_table<- read.table("load/results.txt", header=T) 
full_result_table$beta <-as.numeric(full_result_table$beta)
full_result_table$p <-as.numeric(full_result_table$p)

#Calculate padj value
full_result_table <- full_result_table %>%
  group_by(covar) %>%
  mutate(padj = p.adjust(p, method = "fdr"))

#Create a table for interaction (conditioncontrol:years)
full_result_table<- subset(full_result_table, full_result_table$covar == 'conditioncontrol:years')


#Creating an specific table
selected_columns <- full_result_table[, c("linear","beta", "p", "padj")]

# Create a new data frame with the selected columns
result_table<- as.data.frame(selected_columns)

#write.table(result_table,"save/results.txt", quote = F, col.names = T, row.names = F, sep = '\t')


#Ploting for gene symbols intea of ENSEMBL ids
DE<-full_result_table
#extracting ensembl-ids
ensembl_ids <- DE$linear
# Remove everything after the dot (including the dot)
primary_gene_id <- sub("\\..*", "", ensembl_ids)
# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = primary_gene_id, column = "SYMBOL", keytype = "ENSEMBL")
# Add the gene symbols as a new column in your DEG data frame
DE$gene_symbol <- gene_symbols


#Reorganizing the df
#Creating an specific table
selected_columns <- DE[, c("gene_symbol","beta", "p","padj", "linear")]

# Create a new data frame with the selected columns
DE<- as.data.frame(selected_columns)
#Omit unmapped genes
DE<-na.omit(DE)

DE <- DE[order(DE$padj), ]
volcano.name <- "PDBP  Longitudinal differential expression"
l2f.lim <- 0
plim <- 0.05
plotDE <- DE %>% mutate(gene_type = case_when(beta >= l2f.lim & padj <= plim ~ "up",
                                              beta <= -l2f.lim & padj <= plim ~ "down", TRUE ~ "ns"))


#Get the 9 circ degs
result_9circ<- read.table("nince/circs.txt", header = T)  
circ_names<-result_9circ$linear



#Subset the result_table based on circRNA names to create nine_circ
nine_circ<- plotDE[plotDE$linear %in% circ_names, ]

#Volcano plot for degs
cols <- densCols(plotDE$beta, plotDE$padj)
cols[plotDE$gene_type=='up']<-"#E69F00"
cols[plotDE$gene_type=='down']<-"#56B4E9"
cols[plotDE$gene_type=='ns']<-"#999999"
sizes <- c("up" = 2, "down" = 2, "ns" = 1)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)


ggplot(plotDE, aes(x=beta, y=-log10(padj), size=gene_type, alpha = gene_type))+
  geom_point(col= cols)+
  scale_y_continuous(expand = c(0,0), limits = c(0, 6))+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = gene_symbol), force = 2, nudge_y = 1, size = 3) + #labeling our 9 circs
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: beta") +
  ylab("-log10(p)")+
  theme_bw()+
  theme(legend.position = "none")
