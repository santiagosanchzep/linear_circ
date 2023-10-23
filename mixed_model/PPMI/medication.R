library(reshape2)
library(dplyr)
library(DESeq2)
library(stringr)


# load normalized count matrix
data <- read.table("load/count/matrix.txt", header = T, stringsAsFactors = F)

# load cleaned phenotype data
pheno <- read.table('load/phenotype.txt', header = T, stringsAsFactors = F)



pheno <- subset(pheno, pheno$Status_time != 2)
data <- data[ , colnames(data) %in% pheno$FILE_NAME]

mixed_model_pheno <- pheno %>% dplyr::select(sample_id = FILE_NAME, sex = Gender, status = Status_time, AgeVisit, time_in_years, PATNO, PDMed)

# set age at baseline to their age at first visit
mixed_model_pheno <- mixed_model_pheno %>% group_by(PATNO) %>% mutate(age_at_baseline = min(AgeVisit))

# merge the matrix and phenotype
mixed_model_count_matirx <- data.frame(t(data))
mixed_model_count_matirx$sample_id <- rownames(mixed_model_count_matirx)

mixed_model_longtable <- merge(mixed_model_pheno, mixed_model_count_matirx, by = "sample_id")
mixed_model_longtable <- melt(mixed_model_longtable, id.vars = c('sample_id', "PATNO", 'sex', "AgeVisit", "PDMed", "age_at_baseline", 'status', 'time_in_years' ))

colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "variable"] <- "RNA"
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "value"] <- "counts"
str(mixed_model_longtable)

mixed_model_longtable$RNA <- as.character(mixed_model_longtable$RNA)
mixed_model_longtable <- mixed_model_longtable[order(mixed_model_longtable$PATNO, mixed_model_longtable$time_in_years), ]

mixed_model_longtable[mixed_model_longtable$counts == 0, 'counts'] <- NA
mixed_model_longtable<- na.omit(mixed_model_longtable)

# create count at baseline
mixed_model_longtable$counts_at_baseline <- -9

mixed_model_longtable <- mixed_model_longtable %>% group_by(PATNO, RNA) %>% mutate(counts_at_baseline = counts[which(time_in_years == min(time_in_years))])

any(mixed_model_longtable$counts_at_baseline == -9) # should be FALSE

# checking
patno <- unique(mixed_model_longtable$PATNO)
li_names <- unique(mixed_model_longtable$RNA)

rand_li <- li_names[sample(1:length(li_names), 1)]
rand_sample <- patno[sample(1:length(patno), 1)]
mixed_model_longtable %>% filter(PATNO == rand_sample, RNA == rand_li)

write.table(mixed_model_longtable, 'save/longtable.txt', quote= F, sep = '\t', row.names = F, col.names = T)



#RUN MIXED MODEL
data <- read.table("load/longtable.txt", header = T, stringsAsFactors = F)
data <- data %>% group_by(PATNO, RNA) %>% mutate(visit_counts = length(time_in_years))
# filtered with as least two visits
atleast_2_visit <- subset(data, data$visit_counts > 1)

library("lmerTest")

loop_circs <- as.character(unique(atleast_2_visit$RNA))

full_result_table <- NULL
for (one_circ in loop_circs){
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(atleast_2_visit, atleast_2_visit$RNA == one_circ)
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + sex + PDMed  + status + time_in_years +  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', "PDMed", 'conditioncontrol', 'years', 'conditioncontrol:years')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}

# save result table
write.table(full_result_table, "save/results.txt", quote = F, col.names = T, row.names = F, sep = '\t')

#Read table
full_result_table<- read.table("load/results.txt", header =T)
full_result_table$beta <-as.numeric(full_result_table$beta)
full_result_table$p <-as.numeric(full_result_table$p)

full_result_table <- full_result_table %>%
  group_by(covar) %>%
  mutate(padj = p.adjust(p, method = "fdr"))




#Create a table for  interaction
full_result_table<- subset(full_result_table, full_result_table$covar == 'conditioncontrol:years')


#Creating an specific table
selected_columns <- full_result_table[, c("linear","beta", "p", "padj")]

# Create a new data frame with the selected columns
result_table<- as.data.frame(selected_columns)


#Getting  the  gene symbols
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
DE <- DE[order(DE$padj), ]


#VOLCANO PLOT
volcano.name <- "PPMI  Longitudinal differential expression + Medication"
l2f.lim <- 0
plim <- 0.05
plotDE <- DE %>% mutate(gene_type = case_when(beta >= l2f.lim & padj <= plim ~ "up",
                                              beta <= -l2f.lim & padj <= plim ~ "down", TRUE ~ "ns"))

#Get the 9 circ degs
result_9circ<- read.table("load/9circ.txt", header = T)
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
  geom_label_repel(data = nine_circ, aes(label = gene_symbol), force = 2, nudge_y = 1, size = 3) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: beta") +
  ylab("-log10(p)")+
  theme_bw()+
  theme(legend.position = "none")

