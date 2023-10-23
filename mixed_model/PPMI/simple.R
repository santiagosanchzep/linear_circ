library(data.table)
library(tximport)
library(DESeq2)
library(readxl)
library(ggplot2)
# pheno<- read.table("load/pheno.txt", header = T, stringsAsFactors = F)
# pheno <- subset(pheno, pheno$Status_time != 2)
# # Base directory where subdirectories with quant.genes.sf files are located
# base_directory <- "base/directory"

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

countMatrix<-read.csv("load/count/matrix.csv", row.names = 1)

#Removing sample with low reads
countMatrixclean<- countMatrix[(rowCounts(countMatrix[,-1]<10) < round(0.9*dim(countMatrix[,-1])[2])),]
pheno <- pheno[pheno$FILE_NAME %in% colnames(countMatrixclean),]

# Extract the sample IDs from the phenotype data
sample_ids <- pheno$FILE_NAME

# Subset the countMatrix to include only the columns in sample_ids
countMatrixfinal<- countMatrixclean[, sample_ids]
countMatrixfinal<-floor(countMatrixfinal)


#First create a phase column using the File_NAME
pheno$phase <- ifelse(grepl("Phase1", pheno$FILE_NAME), "Phase1", ifelse(grepl("Phase2", pheno$FILE_NAME), "Phase2", NA))
pheno$phase<- factor(pheno$phase)



#count normalization
rnaDDS <- DESeqDataSetFromMatrix(countData =countMatrixfinal , colData = pheno, design = ~ 1)


vstDDS<-vst(rnaDDS)
# vstDDS<-varianceStabilizingTransformation(rnaDDS)
vst_countmatrix<-assay(vstDDS)

# #save the normalized count matrix for easier loading
# write.table(vst_countmatrix,"save/count/matrix.txt")
# vstDDS<-read.table("load/count/matrix.txt", header = T)



#Principal Components Analyses

# plotPCA(vst.all, intgroup="sex")
pcaData <- plotPCA(vstDDS, intgroup="phase",ntop=22218 ,returnData = T)
write.csv(pcaData, "save/pcaData.txt", row.names = T)

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

#LONGTABLE

# load normalized count matrix
data <- read.table("load/count/matrix.csv", header = T, stringsAsFactors = F)

# load cleaned phenotype data
pheno <- read.table('read/pheno.txt', header = T, stringsAsFactors = F)
pheno <- subset(pheno, pheno$Status_time != 2)
data <- data[ , colnames(data) %in% pheno$FILE_NAME]
pheno<- pheno[colnames(data) %in% pheno$FILE_NAME,]
mixed_model_pheno <- pheno %>% dplyr::select(sample_id = FILE_NAME, sex = Gender, status = Status_time, AgeVisit, time_in_years, PATNO)

# set age at baseline to their age at first visit
mixed_model_pheno <- mixed_model_pheno %>% group_by(PATNO) %>% mutate(age_at_baseline = min(AgeVisit))

# merge the matrix and phenotype
mixed_model_count_matirx <- data.frame(t(data))
mixed_model_count_matirx$sample_id <- rownames(mixed_model_count_matirx)

mixed_model_longtable <- merge(mixed_model_pheno, mixed_model_count_matirx, by = "sample_id")
mixed_model_longtable <- melt(mixed_model_longtable, id.vars = c('sample_id', "PATNO", 'sex', "AgeVisit", "age_at_baseline", 'status', 'time_in_years'))

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

write.table(mixed_model_longtable, 'save/longtable.txt', quote= F, sep = '\t', row.names = F, col.names = T)

