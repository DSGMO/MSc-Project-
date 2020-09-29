#Filtering for protein coding genes and other relevant genes from list_sigs file 
library(ggrepel)
library("ggsci")
library("ggplot2")
library("gridExtra")
library(ggplot2)
library(wesanderson)
library(tidyverse)
library(gapminder)
library(forcats)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(tidyr)
setwd('/Users/katefodder/Desktop/MSC_PROJECT')

Sig_pro_overlap$Gene_name <- sub("-","_", Sig_pro_overlap$Gene_name)
Sig_pro_overlap$Gene_name <- sub("[.]","_", Sig_pro_overlap$Gene_name)

pc_genes <- Sig_pro_overlap[Sig_pro_overlap$Gene_type == "protein_coding",]
IG_C_genes <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="IG_C_gene",]
IG_V_genes <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="IG_V_gene",]
IG_J_genes <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="IG_J_gene",]
ID_D_gene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="ID_D_gene",]
rRNA <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="rRNA",]
miRNA <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="miRNA",]
lncRNA <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="lncRNA",]
TR_C_gene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="TR_C_gene",]
TR_D_gene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="TR_D_gene",]
TR_J_gene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="TR_J_gene",]
misc_RNA <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="misc_RNA",]


IG_V_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="IG_V_pseudogene",]
IG_C_pseudogene   <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="IG_C_pseudogene  ",]
IG_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="IG_pseudogene",]
polymorphic_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="polymorphic_pseudogene",]
processed_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="processed_pseudogene",]
processed_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="processed_pseudogene",]
rRNA_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="rRNA_pseudogene",]
transcribed_processed_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="transcribed_processed_pseudogene",]
transcribed_unprocessed_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="transcribed_unprocessed_pseudogene",]
unprocessed_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="unprocessed_pseudogene",]
unitary_pseudogene <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="unitary_pseudogene",]

pseudogenes <- rbind(IG_V_pseudogene, IG_C_pseudogene, IG_pseudogene, polymorphic_pseudogene, processed_pseudogene, unprocessed_pseudogene, rRNA_pseudogene, transcribed_processed_pseudogene, transcribed_unprocessed_pseudogene, unprocessed_pseudogene, unitary_pseudogene )




geneList <- unique(rel_genes$Gene_name)
write.table(rel_genes, "rel_genes.RData")

# code so only one mutated gene per sample 
miRNA_unique <- miRNA[,-c("icgc_donor_id")]
miRNA_unique <- unique(miRNA_unique)
miRNA_unique_res <- miRNA_unique %>% group_by(Gene_name) %>% summarise(Freq=n())


rel_genes_unique <- rel_genes %>% distinct(Gene_name, icgc_donor_id, .keep_all = TRUE) ##### line so that only one mutated gene per patient 

rel_genes_unique_res <- rel_genes_unique %>% group_by(Gene_name) %>% summarise(freq=n())




rel_genes_unique_res <- rel_genes_unique %>% group_by(Gene_name) %>% summarise(Freq=n())
#Plots with filtered genes 
rel_genes_res <- rel_genes %>% group_by(Gene_name, Gene_type, sig, Signature_type) %>% summarise(freq=n())
rel_genes_rid <- subset(rel_genes_unique_res, freq>)

ggplot(rel_genes_rid, aes(y=freq, x=Gene_name, fill=sig)) + geom_col() + scale_fill_npg()





############################HEATMAP OF ALL GENES AND ALL SIGNATURES #################
Sig_pro_overlap_unique <- Sig_pro_overlap %>% distinct(Gene_name, icgc_donor_id, sig, .keep_all = T)

genes_distinct <- Sig_pro_overlap %>% distinct(Gene_name, icgc_donor_id)
genes_distinct_five <- genes_distinct %>% group_by(Gene_name) %>% summarise(Freq=n())
genes_five <- subset(genes_distinct_five, Freq>4)

genes_list <- as.list(genes_five$Gene_name)

genes_res <- Sig_pro_overlap_unique %>% group_by(sig, Gene_name) %>% summarise(Freq=n())
genes_res_subset <- subset(genes_res, Gene_name %in% genes_list)

for(i in 1:nrow(genes_res_subset)){
  genes_res_subset$log_freq <- log10((genes_res_subset$Freq)+1)
}

genes_res_mat <- acast(genes_res_subset, Gene_name~sig, value.var="log_freq")
genes_res_mat[is.na(genes_res_mat)] <- 0
genes_res_mat_scale <- scale(genes_res_mat)
pheatmap(genes_res_mat, cluster_cols = T, rownames=F, cluster_rows = T, show_rownames = F, annotation_col = list_sigs)

list_sigs <- as.data.frame(colnames(genes_res_mat))
setnames(list_sigs, old = c('colnames(genes_res_mat)'), new = c('sig'))
list_sigs$Signature_type <- paste("f")
for(i in 1:nrow(list_sigs)){
  if(list_sigs$sig[i]=="SBS3"||list_sigs$sig[i]== "SBS30"||
     list_sigs$sig[i]== "SBS44"||list_sigs$sig[i]== "SBS6"||
     list_sigs$sig[i]== "SBS15"||list_sigs$sig[i]=="SBS36"){
    list_sigs$Signature_type[i] <- paste("DNA Repair Related")
  }
  else if(list_sigs$sig[i] == "SBS1"|| list_sigs$sig[i]== "SBS5"||list_sigs$sig[i]=="SBS40"){
    list_sigs$Signature_type[i] <- paste("Ageing Related")
  }
  else if(list_sigs$sig[i] == "SBS2"){
    list_sigs$Signature_type[i] <- paste("APOBEC Related")
  }
  else if(list_sigs$sig[i] == "SBS7a"||list_sigs$sig[i]=="SBS7b"){
    list_sigs$Signature_type[i] <- paste("UV light Related")
  }
  else if(list_sigs$sig[i]=="SBS32"||list_sigs$sig[i]=="SBS31"||
          list_sigs$sig[i]=="SBS11"||list_sigs$sig[i]=="SBS25"||
          list_sigs$sig[i]=="SBS35"){
    list_sigs$Signature_type[i] <- paste("Treatment Related")
  }
  else if(list_sigs$sig[i] == "SBS42"||list_sigs$sig[i]=="SBS22"){
    list_sigs$Signature_type[i] <- paste("Endogenous Exposure Related")
  }
  else if(list_sigs$sig[i] == "SBS37"||list_sigs$sig[i]=="SBS9"||
          list_sigs$sig[i]=="SBS16"||list_sigs$sig[i]=="SBS19"||
          list_sigs$sig[i]=="SBS23"||list_sigs$sig[i]=="SBS28"||
          list_sigs$sig[i]=="SBS41"||list_sigs$sig[i]=="SBS10b"||
          list_sigs$sig[i]=="SBS17b"||list_sigs$sig[i]=="SBS33"||
          list_sigs$sig[i]=="SBS39"||list_sigs$sig[i]=="SBS8"){
    list_sigs$Signature_type[i] <- paste("Other/Unknown")
  }
  else if(list_sigs$sig[i]=="SBS84"||list_sigs$sig[i]=="SBS85"){
    list_sigs$Signature_type[i] <- paste("Cytidine deaminase Related")
  }
}

list_sigs <- list_sigs %>% remove_rownames %>% column_to_rownames(var="sig")









rel_genes_res <- rel_genes %>% group_by(Gene_name, sig, Gene_type) %>% summarise(Freq=n())
rel_genes_all_res <- rel_genes %>% group_by(Gene_name, Gene_type) %>% summarise(Freq=n())
rel_genes_filt <- subset(rel_genes_all_res, Freq>10)

genes_df <- merge(rel_genes_res, rel_genes_filt, by = c("Gene_name","Gene_name"), all.x=F, all.y=T)


#Adding annotation rows 
gene_types <- genes_df[,c("Gene_name","Gene_type.x")]
gene_types <- unique(gene_types)
gene_types <- gene_types %>% remove_rownames %>% column_to_rownames(var="Gene_name")

#Adding annotation columns 

list_sigs <- as.data.frame(colnames(genes_mat))
setnames(list_sigs, old = c('colnames(genes_mat)'), new = c('sig'))
all_SIGS <- as.
for(i in 1:nrow(all_SIGS)){
  if(all_SIGS$sig[i]=="SBS3"||all_SIGS$sig[i]== "SBS30"||all_SIGS$sig[i]== "SBS44"||all_SIGS$sig[i]== "SBS6"){
    all_SIGS$Signature_type[i] <- paste("DNA Repair Related")
  }
  else if(all_SIGS$sig[i] == "SBS1"|| all_SIGS$sig[i]== "SBS5"||all_SIGS$sig[i]=="SBS40"){
    all_SIGS$Signature_type[i] <- paste("Ageing Related")
  }
  else if(all_SIGS$sig[i] == "SBS2"){
    all_SIGS$Signature_type[i] <- paste("APOBEC Related")
  }
  else if(all_SIGS$sig[i] == "SBS7a"||all_SIGS$sig[i]=="SBS7b"){
    all_SIGS$Signature_type[i] <- paste("UV light Related")
  }
  else if(all_SIGS$sig[i]=="SBS32"){
    all_SIGS$Signature_type[i] <- paste("Treatment Related")
  }
  else if(all_SIGS$sig[i] == "SBS42"){
    all_SIGS$Signature_type[i] <- paste("Endogenous Exposure Related")
  }
  else if(all_SIGS$sig[i] == "SBS37"||all_SIGS$sig[i]=="SBS9"||all_SIGS$sig[i]=="SBS16"||all_SIGS$sig[i]=="SBS19"||all_SIGS$sig[i]=="SBS23"||all_SIGS$sig[i]=="SBS28"||all_SIGS$sig[i]=="SBS41"){
    all_SIGS$Signature_type[i] <- paste("Other/Unknown")
  }
  else if(all_SIGS$sig[i]=="SBS84"||all_SIGS$sig[i]=="SBS85")
    all_SIGS$Signature_type[i] <- paste("Cytidine deaminase Related")
  }
}

list_sigs <- list_sigs %>% remove_rownames %>% column_to_rownames(var="sig")

genes_mat <- acast(genes_df, Gene_name~sig, value.var="Freq.x")
genes_mat[is.na(genes_mat)] <- 0
mypal = pal_npg
mypal = list(category = mypal)
genes_mat_scale <- scale(genes_mat)
pheatmap(genes_mat,cluster_cols = T, cluster_rows = T, annotation_row = gene_types, main = "Mutated Genes and Signatures", annotation_col = list_sigs, annotation_colors = scale_fill_npg())

geneList <- as.list(unique(rel_genes$Gene_name))
me <- `mart_export.(6)`
gene_list <- as.data.frame(me$EntrezGene.ID)
library(ReactomePA)

data(gene_list)




#Final HEATMAP 

 ##### line so that only one mutated gene per patient 
rel_genes_unique_res <- rel_genes_unique %>% group_by(Gene_name, sig) %>% summarise(freq=n())
rel_genes_res <- rel_genes_unique %>% group_by(Gene_name) %>% summarise(freq=n())
rel_genes_rid <- subset(rel_genes_res, freq>5)


genes_to_keep <- as.list(rel_genes_rid$Gene_name)

rel_genes <- subset(rel_genes_unique_res, Gene_name %in% genes_to_keep)

genes_mat <- acast(rel_genes, Gene_name~sig, value.var = "freq")
genes_mat[is.na(genes_mat)] <- 0
genes_mat_scale <- scale(genes_mat)

#annotation cols 
list_sigs <- as.data.frame(colnames(genes_mat))
setnames(list_sigs, old = c('colnames(genes_mat)'), new = c('sig'))

for(i in 1:nrow(list_sigs)){
  if(list_sigs$sig[i]=="SBS3"||list_sigs$sig[i]== "SBS30"||list_sigs$sig[i]== "SBS44"||list_sigs$sig[i]== "SBS6"){
    list_sigs$Signature_type[i] <- paste("DNA Repair Related")
  }
  else if(list_sigs$sig[i] == "SBS1"|| list_sigs$sig[i]== "SBS5"||list_sigs$sig[i]=="SBS40"){
    list_sigs$Signature_type[i] <- paste("Ageing Related")
  }
  else if(list_sigs$sig[i] == "SBS2"){
    list_sigs$Signature_type[i] <- paste("APOBEC Related")
  }
  else if(list_sigs$sig[i] == "SBS7a"||list_sigs$sig[i]=="SBS7b"){
    list_sigs$Signature_type[i] <- paste("UV light Related")
  }
  else if(list_sigs$sig[i]=="SBS32"){
    list_sigs$Signature_type[i] <- paste("Treatment Related")
  }
  else if(list_sigs$sig[i] == "SBS42"){
    list_sigs$Signature_type[i] <- paste("Endogenous Exposure Related")
  }
  else if(list_sigs$sig[i] == "SBS37"||list_sigs$sig[i]=="SBS9"||list_sigs$sig[i]=="SBS16"||list_sigs$sig[i]=="SBS19"||list_sigs$sig[i]=="SBS23"||list_sigs$sig[i]=="SBS28"||list_sigs$sig[i]=="SBS41"){
    list_sigs$Signature_type[i] <- paste("Other/Unknown")
  }
  else if(list_sigs$sig[i]=="SBS84"||list_sigs$sig[i]=="SBS85")
    list_sigs$Signature_type[i] <- paste("Cytidine deaminase Related")
}


list_sigs <- list_sigs %>% remove_rownames %>% column_to_rownames(var="sig")
#Annotation rows 

gene_types1 <- subset(rel_genes_unique, Gene_name %in% genes_to_keep)
gene_types2 <- gene_types1[,c("Gene_name","Gene_type")]
gene_types <- gene_types2 %>% distinct(Gene_name, .keep_all = TRUE)
gene_types <- gene_types %>% remove_rownames %>% column_to_rownames(var="Gene_name")

#HEATMAP 
pheatmap(genes_mat_scale, annotation_col = list_sigs, annotation_row = gene_types, fontsize_row = 2.4, cluster_rows = T, show_rownames = F)
#GGPLOT OF FREQUENTLY MUTATED GENES 

SDHB <- Sig_pro_overlap[Sig_pro_overlap$Gene_name== "SDHB",]
SAMD9 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name== "SAMD9",]
SCGB1C1 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name== "SCGB1C1",]
ZNF595 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name== "ZNF595",]
CD68 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name== "CD68",]
OR4N5 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name== "OR4N5",]
RARS2 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name== "RARS2",]
ARHGEF11 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="ARHGEF11",]
TPTE <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="TPTE",]
ORC3 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="ORC3",]
ODF3 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="ODF3",]
SYT3 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="SYT3",]
GGT8P <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="GGT8P",]
IGKV1OR2_1 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="IGKV1OR2_1",]
BX284668_5 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="BX284668_5",]
AF0254981_1 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="AF0254981_1",]
SLC25A15P4 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="SLC25A15P4",]
HERC2P3 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="HERC2P3",]
AC027612_5 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="AC027612_5",]
BAGE2 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="BAGE2",]
RNU1_2 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="RNU1_2",]
IGHD3_10 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="IGHD3_10",]
IGHD2_8 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="IGHD2_8",]
IGHD3_9 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="IGHD3_9",]
IGHD4_11 <- Sig_pro_overlap[Sig_pro_overlap$Gene_name=="IGHD4_11",]

freq_genes <- NULL
freq_genes <- rbind(SDHB, SAMD9, SCGB1C1, ZNF595, CD68, OR4N5, RARS2, ARHGEF11, TPTE, ORC3, ODF3, SYT3,
                    GGT8P, IGKV1OR2_1, BX284668_5, AF254981_1, SLC25A15P4, HERC2P3, ACO27612_5, 
                    BAGE2, RNU1_2, IGHD3_10, IGHD2_8, IGHD3_9, IGHD4_11)

freq_genes <- freq_genes %>% distinct(Gene_name, icgc_donor_id, sig, .keep_all = T)
freq_genes_res <- freq_genes %>% group_by(sig, Signature_type, Gene_name) %>% summarise(freq=n())

ggplot(freq_genes_res, aes(y=freq, x=Gene_name, fill=sig)) + geom_col() + 
   ggtitle("Frequently mutated genes and Signatures associated with Mutations") + 
  xlab("Genes") + ylab("Number of Mutations") + 
  theme(axis.text.x=element_text(angle=90, hjust=1), panel.background = element_rect(fill="white")) + labs(fill="Signature Aetiology") 

freq_genes_res$Gene_name <- factor(freq_genes_res$Gene_name, levels = c("IGHD4_11","IGHD3_10","IGHD3_9","IGHD2_8","ODF3","ARHGEF11","SYT3","ORC3","TPTE","OR4N5","CD68","RARS2","SAMD9","ZNF595","SCGB1C1","SDHB","RNU1_2","BAGE2","HERC2P3","SLC25A15P4", "AF0254981_1","BX284668_5","GGT8P","IGKV1OR2_1"))

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

col_vector16 <- col_vector[2:17]


