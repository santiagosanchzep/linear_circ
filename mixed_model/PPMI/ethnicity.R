library(reshape2)
library(dplyr)
library(stringr)
library(matrixStats)
library(org.Hs.eg.db)
Load PPMI count matrix
countMatrix<-read.csv("/countmatrix.csv", row.names = 1)
#Load PPMI pheno
pheno<- read_excel("/phenotype.xlsx", "PPMI_FINAL" )
#Subsetting for African Americans
pheno<-pheno[pheno$race=="Black or African American",]
#removing Prodromal cases
pheno <- subset(pheno, pheno$`Pt-Category` != "Prodromal")
#Modify phenotype in order to match it with colnames in count matrix.
pheno$time<-pheno$visit_month
pheno$time <- ifelse(pheno$time == "0", "BL", pheno$time)
pheno$time <- ifelse(pheno$time == "6", "V02", pheno$time)
pheno$time <- ifelse(pheno$time == "12", "V04", pheno$time)
pheno$time <- ifelse(pheno$time == "24", "V06", pheno$time)
pheno$time <- ifelse(pheno$time == "36", "V08", pheno$time)
pheno$ID<- paste(pheno$participant_id, ".",pheno$time, sep = "")
pheno$ID <- gsub("PP-", "", pheno$ID)


#Vector of all sample IDs
sample_ids <- ID
# Data frame to store matching colnames and IDs
matching_data <- data.frame(ID = character(), Colnames = character(), stringsAsFactors = FALSE)

# Loop through each sample ID and find matching colnames
for (sample_id in sample_ids) {
  matching_colnames <- grep(sample_id, colnames(countMatrix), value = TRUE)
  # If matching colnames are found, add them to the data frame
  if (length(matching_colnames) > 0) {
    matching_data <- rbind(matching_data, data.frame(ID = sample_id, Colnames = matching_colnames))
  }
}

#Remove problematic sample
pheno_PPMI<-pheno[pheno$ID !="sample",]
#Merge
pheno_PPMI<-merge(pheno_PPMI, matching_data, by="ID")
# Now 'pheno_PPMI' contains a data frame with ID and matching colnames


# Subset the countMatrix based on matching_colnames
countMatrix<- countMatrix[, pheno_PPMI$Colnames]

#Removing sample with low reads
countMatrixclean<- countMatrix[(rowCounts(countMatrix[,-1]<10) < round(0.9*dim(countMatrix[,-1])[2])),]
countMatrix_PPMI<-countMatrixclean



#Load PDBP count matrix
countMatrix<-read.csv("/countmatrix.csv", row.names = 1)
colnames(countMatrix) <- gsub("\\.", "-", colnames(countMatrix))
#Load PDBO pheno
pheno<- read_excel("/phenotype.xlsx", "PDBP")
#Subset for African American
pheno<-pheno[pheno$race=="Black or African American",]

#Removing sample with low reads
countMatrixclean<- countMatrix[(rowCounts(countMatrix[,-1]<10) < round(0.9*dim(countMatrix[,-1])[2])),]
pheno_PDBP<-pheno[pheno$sample_id %in% colnames(countMatrixclean),]

# Extract the sample IDs from the phenotype data
sample_ids <- pheno_PDBP$sample_id

# Subset the countMatrix to include only the columns in sample_ids
countMatrixfinal<- countMatrixclean[, sample_ids]
countMatrix_PDBP<-floor(countMatrixfinal)



# Combine both matrixes PDBP and PPMI
shared_genes<-countMatrix_PDBP[rownames(countMatrix_PDBP) %in% rownames(countMatrix_PPMI),]
countMatrix_PDBP$rownames_col <- rownames(countMatrix_PDBP)
countMatrix_PPMI$rownames_col <- rownames(countMatrix_PPMI)
countMatrix_combined <- merge(countMatrix_PDBP, countMatrix_PPMI, by= "rownames_col")
rownames(countMatrix_combined) <- countMatrix_combined$rownames_col

# Remove the "rowname_col" column from the data frame
countMatrix_combined <- countMatrix_combined[,-1]


# Get the common column names between pheno_PPMI and pheno_PDBP to create a combined phenotype
common_columns <- intersect(colnames(pheno_PPMI), colnames(pheno_PDBP))

# Subset pheno_PPMI to include only the common columns
pheno_PPMI_cleaned <- pheno_PPMI[common_columns]

# Combine the cleaned pheno_PPMI and pheno_PDBP data frames
combined_pheno <- rbind(pheno_PPMI_cleaned, pheno_PDBP)

write.table(combined_pheno, "/combined_phenotype.txt" )
write.table(countMatrix_combined, "/combined_countmatrix.txt" )

#count normalization
countMatrix_combined<- round(countMatrix_combined)
combined_pheno$visit_month<-as.factor(combined_pheno$visit_month)
rnaDDS <- DESeqDataSetFromMatrix(countData = countMatrix_combined, colData = combined_pheno, design = ~ 1)


vstDDS<-vst(rnaDDS)
vstDDS<-varianceStabilizingTransformation(rnaDDS)
vst_countmatrix<-assay(vstDDS)
write.table(vst_countmatrix,"/norm_countmatrix.txt")

#Principal Components analysis

plotPCA(vst.all, intgroup="sex")
pcaData <- plotPCA(vstDDS, intgroup="visit_month",ntop=500 ,returnData = T)
#write.csv(pcaData, "/pcaData.csv", row.names = T)

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
plot(plotStatus)


# load normalized count matrix
data<-read.csv("/norm_countmatrix.csv", row.names = 1)
data<- data[(rowCounts(data[,-1]<10) < round(0.9*dim(data[,-1])[2])),]
#load combined phenotype data
pheno <- read.table("/combined_pheno.txt", header = T, stringsAsFactors = F)


#In mixed model we need to merge count matrix and phenotype, so wee need to assign same values  

# Create a column of sample_id_modified similar to what we do at the begining to merge PPMi phenotype and count matrix.
pheno$sample_id_modified <- ifelse(pheno$Cohort == "PDBP",
                                  as.character(pheno$sample_id),
                                  paste0(sub("^PP-", "", as.character(pheno$participant_id)),
                                         ifelse(pheno$visit_month == 0, ".BL",
                                                ifelse(pheno$visit_month == 6, ".V02",
                                                       ifelse(pheno$visit_month == 12, ".V04",
                                                              ifelse(pheno$visit_month == 24, ".V06", ".V08"))))))


#assign the new sample id  as  colnames in count matrix, only for PPMI samples
sample_names <- pheno$sample_id_modified[1:67]
colnames(data)[17:83] <- sample_names



pheno$time_in_years <- pheno$visit_month / 12
pheno$Pt-Category <- ifelse(pheno$Pt-Category == "Case", 1,
                              ifelse(pheno$Pt-Category == "Control", 0, pheno$Pt-Category))

mixed_model_pheno <- pheno %>% dplyr::select(sample_id=sample_id_modified, sex, status = `Pt-Category`, age_at_baseline, time_in_years, PATNO = participant_id)

mixed_model_pheno <- mixed_model_pheno %>% group_by(PATNO) %>% mutate(age_at_baseline)
# merge the matrix and phenotype
mixed_model_count_matirx <- data.frame(t(data))
mixed_model_count_matirx$sample_id<- rownames(mixed_model_count_matirx)

mixed_model_longtable <- merge(mixed_model_pheno, mixed_model_count_matirx, by = "sample_id")
mixed_model_longtable <- melt(mixed_model_longtable, id.vars = c( "sample_id", "PATNO", 'sex', "age_at_baseline", 'status', 'time_in_years'))

colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "variable"] <- "circRNA"
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "value"] <- "counts"
str(mixed_model_longtable)

mixed_model_longtable$circRNA <- as.character(mixed_model_longtable$circRNA)
mixed_model_longtable <- mixed_model_longtable[order(mixed_model_longtable$PATNO, mixed_model_longtable$time_in_years), ]

mixed_model_longtable[mixed_model_longtable$counts == 0, 'counts'] <- NA
mixed_model_longtable <- na.omit(mixed_model_longtable)



# create count at baseline
mixed_model_longtable$counts_at_baseline <- -9

mixed_model_longtable <- mixed_model_longtable %>% group_by(PATNO, circRNA) %>% mutate(counts_at_baseline = counts[which(time_in_years == min(time_in_years))])

any(mixed_model_longtable$counts_at_baseline == -9) # should be FALSE

# checking
patno <- unique(mixed_model_longtable$PATNO)
circ_names <- unique(mixed_model_longtable$circRNA)

rand_circ <- circ_names[sample(1:length(circ_names), 1)]
rand_sample <- patno[sample(1:length(patno), 1)]
mixed_model_longtable %>% filter(PATNO == rand_sample, circRNA == rand_circ)

write.table(mixed_model_longtable, '/combined_longtable_AA.csv', quote= F, sep = '\t', row.names = F, col.names = T)


data<-read.table('/combined_longtable_AA.csv', header = T)

# filtered with as least two visits
data<- data %>% group_by(PATNO, circRNA) %>% mutate(visit_counts = length(time_in_years))
atleast_2_visit <- subset(data, data$visit_counts > 1)
str(atleast_2_visit)
atleast_2_visit <- atleast_2_visit %>%
  filter(circRNA != "ENSG00000279809.1") #Problematic gene, only present in PDBP samples

library("lmerTest")

loop_circs <- as.character(unique(atleast_2_visit$circRNA))

full_result_table <- NULL
for (one_circ in loop_circs){
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(atleast_2_visit, atleast_2_visit$circRNA == one_circ)
  print(paste("Number of rows in tmp for circRNA", one_circ, ":", nrow(tmp)))
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + status + time_in_years +  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}


# save result table
write.table(full_result_table, "/results.AA.txt", quote = F, col.names = T, row.names = F, sep = '\t')



#Read table
full_result_table<- read.table("/results.AA.txt", header=T)
full_result_table$beta <-as.numeric(full_result_table$beta)
full_result_table$p <-as.numeric(full_result_table$p)

full_result_table <- full_result_table %>%
  group_by(covar) %>%
  mutate(padj = p.adjust(p, method = "fdr"))


#Create a table for  interaction
full_result_table<- subset(full_result_table, full_result_table$covar == 'status:time_in_years')

#Creating an specific table
selected_columns <- full_result_table[, c("linear","beta", "p", "padj")]

# Create a new data frame with the selected columns
result_table<- as.data.frame(selected_columns)


#Obtaining  the  gene symbols
DE<-full_result_table
#extracting degs
ensembl_ids <- DE$linear
# Remove everything after the dot (including the dot)
primary_gene_id <- sub("\\..*", "", ensembl_ids)
# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = primary_gene_id, column = "SYMBOL", keytype = "ENSEMBL")
# Add the gene symbols as a new column in your DEG data frame
DE$gene_symbol <- gene_symbols



#Reorganizing the df
#Creating an specific table
selected_columns <- DE[, c("gene_symbol","beta", "p", "padj", "linear")]

# Create a new data frame with the selected columns
DE<- as.data.frame(selected_columns)
# DE<-na.omit(DE)



#VOLCANO PLOT
DE <- DE[order(DE$padj), ]
volcano.name <- "PPMI & PDBP  Longitudinal differential expression African Americans"
l2f.lim <- 0
plim <- 0.05
plotDE <- DE %>% mutate(gene_type = case_when(beta >= l2f.lim & padj <= plim ~ "up",
                                              beta <= -l2f.lim & padj <= plim ~ "down", TRUE ~ "ns"))

#Get the 9 circ degs
result_9circ<- read.table("/9circs.txt", header = T)

circ_names<-result_9circ$linear

#Subset the result_table based on circRNA names
nine_circ<- plotDE[plotDE$linear %in% circ_names, ]

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
  geom_label_repel(data = nine_circ, aes(label = gene_symbol), force = 2, nudge_y = 1, size = 3) + #Labelling 9 circs
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: beta") +
  ylab("-log10(p)")+
  theme_bw()+
  theme(legend.position = "none")

