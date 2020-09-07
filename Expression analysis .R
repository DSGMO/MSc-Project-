library(sva)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggpubr)

#Expression analysis 
exp_seq.PRAD.CA <- read.delim("~/Downloads/exp_seq.PRAD-CA.tsv")
load("/Users/katefodder/Desktop/MSC_PROJECT/exp_array.RData")
exp_seq <- exp_seq.PRAD.CA

length(unique(exp_seq.PRAD.CA$icgc_donor_id)) #144
length(unique(exp_array.PRAD.CA$icgc_donor_id)) #213



exp_array$samples <- paste(exp_array$icgc_donor_id, exp_array$analysis_id, sep = "-")
array_mat1 <- exp_array[,c("samples","Gene.name","normalized_expression_value")]
array_mat2 <- melt(array_mat1)
array_mat3 <- dcast(array_mat2, samples~Gene.name, value.var="value", fun.aggregate = sum)

exp_seq$samples <- paste(exp_seq$icgc_donor_id, exp_seq$analysis_id, seq = "-")
seq_mat1 <- exp_seq[,c("samples","gene_id","normalized_read_count")]
seq_mat2 <- melt(seq_mat1)
seq_mat3 <- dcast(seq_mat2, samples~gene_id, value.var = "value", fun.aggregate = sum)
 
seq_samples <- as.data.frame(seq_mat3$samples)
array_samples <- as.data.frame(array_mat3$samples)

seq_samples$label <- paste("RNASeq")
array_samples$label <- paste("Array")

all_data <- rbind(array_samples, seq_samples)

setnames(seq_samples, old = c("seq_mat3$samples"), new = c("Sample"))
setnames(array_samples, old = c("array_mat3$samples"), new = c("Sample"))
cancer_type1$samples <- paste(seq_mat3$samples )
cancer_type2$samples <- array_mat3$samples

array_genes <- colnames(array_mat3)
seq_genes <- colnames(seq_mat3)

seq_genes <- as.data.frame(seq_genes)
array_genes <- as.data.frame(array_genes)
setnames(seq_genes, old= c("seq_genes"), new=c("genes"))
setnames(array_genes, old = c("array_genes"), new=c("genes"))

all_genes <- merge(seq_genes, array_genes, by = c("genes"), all.x=F, all.y=F)

all_genes = as.data.table(all_genes)
setkey(all_genes, genes)
all_genes <- all_genes[!"samples"]
all_genes_list <- as.vector(all_genes$genes)

#adjusting datasets to matrices 
setnames(seq_samples, old = c("seq_mat3$samples"), new = c("Sample"))
setnames(array_samples, old = c("array_mat3$samples"), new = c("Sample"))
cancer_type1$samples <- paste(seq_mat3$samples )
cancer_type2$samples <- array_mat3$samples

array_genes <- colnames(array_mat3)
seq_genes <- colnames(seq_mat3)

a_mat <- as.data.frame(array_mat3)
s_mat <- as.data.frame(seq_mat3)

#Finding common genes between two matrices and binding two datasets together 

a_mat1 <- a_mat[,all_genes_list]
s_mat1 <- s_mat[,all_genes_list]

combined_data <- rbind(a_mat1, s_mat1)


combined_data2 <- as.data.frame(combined_data[,-17272])
rownames(combined_data) <- combined_samples$Sample

combined_samples <- rbind(array_samples, seq_samples)


#Adding gene names to array data 
library(biomaRt)
df_symbol<-getBM(filters="refseq_dna", attributes="external_gene_id", values=refseq, mart=)
ensembl <- useMart("ensembl")
useMart("ensembl", exp_array.PRAD.CA,
        archive=FALSE, ensemblRedirect = NULL, version, verbose = FALSE)
mart= useMart("ensembl")


mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

array_values <- as.data.frame(exp_array.PRAD.CA$gene_id)
df_symbol <- getBM(filters = "refseq_mrna", attributes="external_gene_name", values= array_values, mart = mart)


array_values = array_values[!is.na(array_values)]
write.table(array_values, "array_val.txt", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

attributes = listAttributes(mart)
filters = listFilters(mart)

all$gene_id <- all$RefSeq.mRNA.ID
all <- all[-which(all$gene_id == ""), ]
exp_array <- merge(exp_array.PRAD.CA, all, by = c("gene_id","gene_id"), all.x=T, all.y=F)

save(exp_array, file="exp_array.RData")

exp_seq <- exp_seq.PRAD.CA

#Combat and batch 
platform<-c(rep(1,n_samples_seq_tcga),rep(0,n_samples_seq_cgga))
input_ge_tcga_cgga<-merge(log(matrix_to_plot_hs_nocombat2+1,2),log(CGGA_table2[,-1]+1,2),by="row.names")
n_samples_seq_tcga<-ncol(matrix_to_plot_hs_nocombat2)
rownames(CGGA_table2)<-CGGA_table2[,1]


n_samples_seq <- nrow(exp_seq)
n_samples_array <- nrow(exp_array)
platform <- c(rep(1, n_samples_array), rep(0, n_samples_seq))

exp_array_rid <- log(exp_array[,2+8+9])
exp_seq_rid <- log(exp_seq[,1+2+7+9])
input_exp <- merge(exp_array, exp_seq, by="icgc_donor_id", all.x=F, all.y=F)



rownames(exp_array) <- exp_array[,2]

input_ge_tcga_cgga<-merge(log(matrix_to_plot_hs_nocombat2+1,2),log(CGGA_table2[,-1]+1,2),by="row.names")

combined_data2 = as.data.frame(combined_data2)
setwd("/Users/katefodder/Desktop/MSC_PROJECT")

write.table(combined_data, file="combined_data.txt")

#Remove columns/genes with 0 variance
zeroVar <- function(data, useNA = 'ifany') {
  out <- apply(data, 2, function(x) {length(table(x, useNA = useNA))})
  which(out==1)
}
combined_data2 <- combined_data[,-zeroVar(combined_data)]


#Combat - Anna 


combined_data <- as.matrix(combined_data)
combined_data <- t(combined_data)
modcombat = model.matrix(~1, data=all_data)
batch = all_data$label

scaled_data = ComBat(dat=combined_data, batch=batch, mod=modcombat, par.prior=TRUE)
scaled_data <- data.frame(t(scaled_data))

input_combat <- ComBat(combined_data2, )

#Combat - Guidantonio 
combined_data2 <- data.matrix(combined_data2)
combined_data2 <- t(combined_data2)

n_samples_seq <- nrow(exp_seq)
n_samples_array <- nrow(exp_array)
platform <- c(rep(1, n_samples_seq), rep(0, n_samples_array))

input_combat2<-ComBat(combined_data2, platform)


library(ggfortify)

print("PCA")

#Pre combat
input_pca<-data.frame(t(combined_data))
pca_res <- prcomp(input_pca, scale. = TRUE)
input_pca$assay<-as.factor(all_data$label)

print(dim(input_pca))
#Post combat 

input_pca_combat<-data.frame(scaled_data)
pca_res_combat <- prcomp(input_pca_combat, scale. = TRUE)
input_pca_combat$assay<-as.factor(all_data$label)

output_pca<-paste(paste(list_directories[i],'method.fc',pval,lfc,sep='.'),'.PCA_TCGA_CGGA.pdf',sep='')

#PCA plots pre and post combat 
pdf(output_pca)
autoplot(pca_res, data = input_pca, colour = 'assay',main="No Combat")
autoplot(pca_res_combat,data=input_pca_combat,colour = 'assay',main="Combat")


#Saving Data 

write.table(scaled_data, "scaled_expression_data.txt", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)




#Analysing differeing expression between mutated and non-mutated genes 

expression_data <- as.data.frame(scaled_data)

expression_data$icgc_donor_id <- rownames(expression_data)
expression_data$icgc_donor_id <- sub("-PRAD_CA-RNA_ARRAY","", expression_data$icgc_donor_id)
expression_data$icgc_donor_id <- sub("PRAD_CA-WTS -", "", expression_data$icgc_donor_id)

rownames(expression_data) <- expression_data$icgc_donor_id
exp_data <- expression_data %>% distinct(icgc_donor_id, A1BG ,keep_all=TRUE)

expression_data2 <- distinct(expression_data, rownames(expression_data))
#If rownames the same - remove the 
require(dplyr)


#donor_ids which are present in both datasets 

exp_array_id_list <- exp_array$icgc_donor_id  
exp_seq_id_list <- exp_seq$icgc_donor_id 

duplicated_donors <- intersect(exp_array_id_list, exp_seq_id_list)

exp_data <- expression_data %>% group_by(icgc_donor_id) 


# Genes of interest 
IGHD3_10 <- expression_data$IGHD3_10
IGHD3_10 <- as.data.frame(IGHD3_10)
IGHD3_10$icgc_donor_id <- paste(expression_data$icgc_donor_id)
IGHD3_10_muts <- motif_df[motif_df$Gene_name == "IGHD3_10",]
IGHD3_10_muts <- IGHD3_10_muts[,c("icgc_donor_id")]
IGHD3_10$mutation <-  IGHD3_10$icgc_donor_id%in%IGHD3_10_muts 
IGHD3_10$mutation <- sub("TRUE","Mutated", IGHD3_10$mutation)
IGHD3_10$mutation <- sub("FALSE","Not Mutated", IGHD3_10$mutation)

mut_IGHD3_10 <- subset(IGHD3_10, mutation=="Mutated")
nmut_IGHD3_10 <- subset(IGHD3_10, mutation=="Not Mutated")

setnames(IGHD3_10, old = c("IGHD3_10"), new = c("Val"))
ggplot(IGHD3_10, aes(y=Val, x=mutation, fill=mutation)) + 
  geom_jitter() + 
  geom_boxplot(alpha = 0.7) + 
  ggtitle("IGHD3_10") + 
  ylab("Expression Value") + 
  xlab("Motif alteration present/not present") + 
  stat_compare_means(method = "t.test") +
  scale_fill_aaas() + theme_bw()


colnames(scaled_data) <- colnames(combined_data)
rownames(scaled_data) <- rownames(combined_data)
  
