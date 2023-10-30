llll# load normalized count matrix
counts<-read.table("load/pdbp_normalized_matrix.txt",header = T, stringsAsFactors = F)
counts<-as.matrix(counts)
# load cleaned phenotype data
pheno <- read.table('load/pdbp_phenotype.txt', header = T, stringsAsFactors = F)
# load normalized cell counts
cell_counts<-read.csv("load/pdbp_cellcounts.csv")

#merge cell counts to pheno
pheno <- merge(cell_counts, pheno, by = "sample_id")

# check time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12

# select covariates in pheno that are needed for mixed models
mixed_model_pheno <- pheno %>% dplyr::select(sample_id, sex, status = Pt.Category, age_at_baseline, time_in_years, PATNO = participant_id,B_cell, T_cell_CD4, T_cell_CD8, Monocyte, NK_cell, Neutrophil)
mixed_model_pheno$PATNO <- str_replace_all(mixed_model_pheno$PATNO, "PD-", "")
mixed_model_pheno <- mixed_model_pheno[order(mixed_model_pheno$PATNO, mixed_model_pheno$time_in_years), ]
tmp <- mixed_model_pheno %>% group_by(PATNO) %>% mutate(visit_counts = length(time_in_years))
table(tmp$visit_counts)

mixed_model_count_matirx <- data.frame(t(counts))
mixed_model_count_matirx$sample_id <- rownames(mixed_model_count_matirx)

mixed_model_longtable <- merge(mixed_model_pheno, mixed_model_count_matirx, by = "sample_id")

# create the long table
mixed_model_longtable <- melt(mixed_model_longtable, id.vars = c('sample_id', "PATNO", 'sex', 'age_at_baseline', 'status', "time_in_years", "B_cell", "T_cell_CD4", "T_cell_CD8", "Monocyte", "NK_cell",  "Neutrophil"))

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

write.table(mixed_model_longtable, 'save/cell_counts_pdbp_longtable.txt', quote= F, sep = '\t', row.names = F, col.names = T)

#MIXED MODEL WITH LONGTABLE

# load normalized count matrix
data<- read.table("read/cell_counts_pdbp_longtable.txt", header = T, stringsAsFactors = F)
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
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + B_cell + T_cell_CD4 + T_cell_CD8 + Monocyte + NK_cell + Neutrophil + sex + status + time_in_years +  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', "B_cell", "T_cell_CD4", "T_cell_CD8", "Monocyte", "NK_cell",  "Neutrophil", 'sex', 'conditioncontrol', 'years', 'conditioncontrol:years')
  res$covar <- rownames(res)
  res$linear <- one_linear
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}


# save result table
full_result_table %>% write.table("save/cell_counts_pdbp_results.txt", quote = F, col.names = T, row.names = F, sep = '\t')


#Read table
full_result_table<- read.table("load/cell_counts_pdbp_results.txt", header =T) 
full_result_table$beta <-as.numeric(full_result_table$beta)
full_result_table$p <-as.numeric(full_result_table$p)

full_result_table <- full_result_table %>%
  group_by(covar) %>%
  mutate(padj = p.adjust(p, method = "fdr"))

#Create a table for intercept and interaction
full_result_table<- subset(full_result_table, full_result_table$covar == 'conditioncontrol:years')


#Creating an specific table
selected_columns <- full_result_table[, c("linear","beta", "p", "padj")]

# Create a new data frame with the selected columns
result_table<- as.data.frame(selected_columns)

#write.table(result_table,"save/cell_counts_results.txt", quote = F, col.names = T, row.names = F, sep = '\t')

#Ploting for gene symbols instead of ENSEMBL ids
DE<-result_table
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

DE <- DE[order(DE$padj), ]
volcano.name <- "PDBP  Longitudinal differential expression + Cell counts"
l2f.lim <- 0
plim <- 0.05
plotDE <- DE %>% mutate(gene_type = case_when(beta >= l2f.lim & padj <= plim ~ "up",
                                              beta <= -l2f.lim & padj <= plim ~ "down", TRUE ~ "ns"))
#Get the 9 circ degs
result_9circ<- read.table("load/9circs.txt", header = T)  

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
  geom_label_repel(data = nine_circ, aes(label = gene_symbol), force = 2, nudge_y = 1, size = 3) + #label the nine circ
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: beta") +
  ylab("-log10(p)")+
  theme_bw()+
  theme(legend.position = "none")



