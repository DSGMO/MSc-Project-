#Script for analysing signatures associated with most affected transcription factors 
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
library(reshape2)
library(tidyr)
setwd('/Users/katefodder/Desktop/MSC_PROJECT')
Homo_sapiens_TF <- read.delim("~/Downloads/Homo_sapiens_TF.txt")
load("/Users/katefodder/Desktop/MSC_PROJECT/Sig_pro_overlap_annotated.RDATA")
load("/Users/katefodder/Desktop/MSC_PROJECT/PRAD_CA_hocomoco_nopvalue.RData")


load("/Users/katefodder/Desktop/MSC_PROJECT/Sig_pro_overlap_annotated.RDATA")
Sig_pro_overlap$chr <- paste(Sig_pro_overlap$Promoter_chromosome)
Sig_pro_overlap$chr <- sub("","chr", Sig_pro_overlap$chr)
Sig_pro_overlap$chr <- gsub("chr23", "chrX", Sig_pro_overlap$chr)
Sig_pro_overlap$chr <- gsub("chr24", "chrY", Sig_pro_overlap$chr)

Sig_pro_overlap$names <- paste(Sig_pro_overlap$chr, Sig_pro_overlap$mutation_end,
                               Sig_pro_overlap$ref_sig, Sig_pro_overlap$alt_sig,
                               Sig_pro_overlap$icgc_donor_id, Sig_pro_overlap$sig, sep = ":")

motif_df2 <- NULL
motif_df2 <- merge(Sig_pro_overlap, ALL_RESULTS, by = c("names","names"), all.x=F, all.y=T, allow.cartesian=T, fill=TRUE)

motif_df <- motif_df2[, c("names","Gene_name","Gene_type","ref_sig","alt_sig",
                          "icgc_donor_id","sig",
                          "snpPos","geneSymbol","seqMatch","effect","Signature_type"),]

setnames(motif_df, old = c('geneSymbol'), 
         new = c('TF'), skip_absent = T)

setnames(Homo_sapiens_TF, old = c("Symbol"), new = ("TF"), skip_absent = T)
motif_df <- merge(motif_df, Homo_sapiens_TF, by=c("TF","TF"), all.x=T, all.y=F)


motif_df_strong <- subset(motif_df, effect=="strong")
motif_df_weak <- subset(motif_df, effect =="weak")

motif_df[is.na(motif_df)] <- "FALSE"
motif_df <- motif_df[!(motif_df$Gene_name=="FALSE"),]

motif_df_unique <- motif_df %>% distinct(TF, icgc_donor_id, .keep_all = T)

sig_tf_unique <- motif_df_unique %>% group_by(TF, sig, Signature_type) %>% summarise(Freq=n())
sig_unique <- motif_df_unique %>% group_by(TF) %>% summarise(Freq=n())

sig_tf <- motif_df %>% group_by(TF, sig, Signature_type) %>% summarise(Freq=n())

for(i in 1:nrow(sig_tf)){
  sig_tf$log_freq <- log10((sig_tf$Freq)+1)
}

sig_tf2 <- subset(sig_tf, Freq>50)
sig_tf_mat <- acast(sig_tf, TF~sig, value.var="log_freq")
sig_tf_mat[is.na(sig_tf_mat)] <- 0
pheatmap(sig_tf_mat, cluster_cols = T, rownames=F, cluster_rows = T)

#removing ageing signatures 

sig_tf_rd <- sig_tf[!(sig_tf$sig == "SBS40"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS5"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS1"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS37"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS7a"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS7b"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS7c"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS7d"),]
 # Heatmap 

sig_tf_rd2 <- subset(sig_tf_rd, Freq>50)

sig_tf_rd_mat <- acast(sig_tf_rd2, TF~sig, value.var = "Freq")
scale_sig_tf_rd_mat <- scale(sig_tf_rd_mat)
scale_sig_tf_rd_mat[is.na(scale_sig_tf_rd_mat)] <- 0

tfs <- as.data.frame(rownames(sig_tf_rd_mat))
setnames(tfs, old = c("rownames(sig_tf_rd_mat)"), new = c("TF"))
tfs <- merge(tfs, Homo_sapiens_TF, by=c("TF","TF"), all.x=T, all.y=F)
tfs <- tfs[,c("TF","Family")]
tfs <- tfs %>% remove_rownames %>% column_to_rownames(var="TF")


pheatmap(scale_sig_tf_rd_mat, annotation_row = tfs, annotation_col = list_sigs, fontsize_row = 3)
pheatmap(scale_sig_tf_rd_mat, annotation_row = tfs, 
         annotation_col = list_sigs, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         fontsize_row = 3, )

#Adding annotation rows and columns to heatmap 

tfs <- as.data.frame(rownames(sig_tf_rd_mat))
setnames(tfs, old = c("rownames(sig_tf_rd_mat)"), new = c("TF"))
tfs <- merge(tfs, Homo_sapiens_TF, by=c("TF","TF"), all.x=T, all.y=F)


for(i in 1:nrow(tfs)){
  print(i)
  if(tfs$TF[i] == "BPTF"||tfs$TF[i] == "BRCA1"||tfs$TF[i] == "ENO1"||tfs$TF[i] == "FUBP1"||
     tfs$TF[i] == "HLTF"||tfs$TF[i] == "T"){
    tfs$Family[i] <- paste("Miscellaneous")
  }
  else if(tfs$TF[i] =="ATF2+ATF4"){
    tfs$Family[i] <- paste("TF_bZIP")
  }
  else if(tfs$TF[i] =="GABPB1+GABPB2"){
    tfs$Family[i] <- paste("ETS")
  }
  else if(tfs$TF[i] =="NFIA+NFIB+NFIC+NFIX"){
    tfs$Family[i] <- paste("Nuclear Factor")
  }
  else if(tfs$TF[i] =="TBP"){
    tfs$Family[i] <- paste("TBF_related")
  }
}

tfs <- tfs[,c("TF","Family")]

list_sigs <- as.data.frame(colnames(sig_tf_rd_mat))
setnames(list_sigs, old = c('colnames(sig_tf_rd_mat)'), new = c('sig'))




tfs <- tfs %>% remove_rownames %>% column_to_rownames(var="TF")

list_sigs <- list_sigs %>% remove_rownames %>% column_to_rownames(var="sig")

sig_tf_mat_scale <- scale(sig_tf_mat)
pheatmap(sig_tf_mat_scale, annotation_row = tfs, annotation_col = list_sigs,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         fontsize_row = 3
)
#new df with 'relevant' genes


pc_genes <- motif_df[motif_df$Gene_type == "protein_coding",]
IG_C_genes <- motif_df[motif_df$Gene_type=="IG_C_gene",]
IG_V_genes <- motif_df[motif_df$Gene_type=="IG_V_gene",]
IG_J_genes <- motif_df[motif_df$Gene_type=="IG_J_gene",]
ID_D_gene <- motif_df[motif_df$Gene_type=="ID_D_gene",]
rRNA <- motif_df[motif_df$Gene_type=="rRNA",]
miRNA <- motif_df[motif_df$Gene_type=="miRNA",]
lncRNA <- motif_df[motif_df$Gene_type=="lncRNA",]
TR_C_gene <- motif_df[motif_df$Gene_type=="TR_C_gene",]
TR_D_gene <- motif_df[motif_df$Gene_type=="TR_D_gene",]
TR_J_gene <- motif_df[motif_df$Gene_type=="TR_J_gene",]
misc_RNA <- motif_df[motif_df$Gene_type=="misc_RNA",]


IG_V_pseudogene <- motif_df[motif_df$Gene_type=="IG_V_pseudogene",]
IG_C_pseudogene   <- motif_df[motif_df$Gene_type=="IG_C_pseudogene  ",]
IG_pseudogene <- motif_df[motif_df$Gene_type=="IG_pseudogene",]
polymorphic_pseudogene <- motif_df[motif_df$Gene_type=="polymorphic_pseudogene",]
processed_pseudogene <- motif_df[motif_df$Gene_type=="processed_pseudogene",]
processed_pseudogene <- motif_df[motif_df$Gene_type=="processed_pseudogene",]
rRNA_pseudogene <- motif_df[motif_df$Gene_type=="rRNA_pseudogene",]
transcribed_processed_pseudogene <- motif_df[motif_df$Gene_type=="transcribed_processed_pseudogene",]
transcribed_unprocessed_pseudogene <- motif_df[motif_df$Gene_type=="transcribed_unprocessed_pseudogene",]
unprocessed_pseudogene <- motif_df[motif_df$Gene_type=="unprocessed_pseudogene",]
unitary_pseudogene <- motif_df[motif_df$Gene_type=="unitary_pseudogene",]

pseudogenes <- rbind(IG_V_pseudogene, IG_C_pseudogene, IG_pseudogene, polymorphic_pseudogene, processed_pseudogene, unprocessed_pseudogene, rRNA_pseudogene, transcribed_processed_pseudogene, transcribed_unprocessed_pseudogene, unprocessed_pseudogene, unitary_pseudogene )
rel_genes <- rbind(pc_genes, ID_D_gene, IG_C_genes, rRNA, miRNA, lncRNA, TR_C_gene, TR_D_gene, TR_J_gene, misc_RNA)



pc_genes <- motif_df[motif_df$Gene_type == "protein_coding",]
IG_C_genes <- motif_df[motif_df$Gene_type=="IG_C_gene",]
ID_D_gene <- motif_df[motif_df$Gene_type=="ID_D_gene",]
rRNA <- motif_df[motif_df$Gene_type=="rRNA",]
miRNA <- motif_df[motif_df$Gene_type=="miRNA",]
lncRNA <- motif_df[motif_df$Gene_type=="lncRNA",]
TR_C_gene <- motif_df[motif_df$Gene_type=="TR_C_gene",]
TR_D_gene <- motif_df[motif_df$Gene_type=="TR_D_gene",]
TR_J_gene <- motif_df[motif_df$Gene_type=="TR_J_gene",]
misc_RNA <- motif_df[motif_df$Gene_type=="misc_RNA",]
un_pseudo <- motif_df[motif_df$Gene_type == "unprocessed_pseudogene",]
pseudo <- motif_df[motif_df$Gene_type == "processed_pseudogene",]


motif_df3 <- rbind(pc_genes, ID_D_gene, IG_C_genes, rRNA, miRNA, lncRNA, TR_C_gene, TR_D_gene, TR_J_gene, misc_RNA)

I_genes <- rbind(IG_C_genes, ID_D_gene)
rel_map <- motif_df3 %>% group_by(TF, sig, Signature_type) %>% summarise(Freq=n())

rel_map <- acast(rel_map, TF~sig, value.var="Freq")
rel_map[is.na(rel_map)] <- 0
pheatmap(rel_map)




#New heatmap with removed ageing signatures 



#New heatmap with removed ageing signatures and only relevant genes 

#Adding annotation rows and columns to heatmap 

tfs <- as.data.frame(rownames(rel_map))
setnames(tfs, old = c("rownames(rel_map)"), new = c("TF"))
tfs <- merge(tfs, Homo_sapiens_TF, by=c("TF","TF"), all.x=T, all.y=F)


for(i in 1:nrow(tfs)){
  print(i)
  if(tfs$TF[i] == "BPTF"||tfs$TF[i] == "BRCA1"||tfs$TF[i] == "ENO1"||tfs$TF[i] == "FUBP1"||
     tfs$TF[i] == "HLTF"||tfs$TF[i] == "T"){
    tfs$Family[i] <- paste("Miscellaneous")
  }
  else if(tfs$TF[i] =="ATF2+ATF4"){
    tfs$Family[i] <- paste("TF_bZIP")
  }
  else if(tfs$TF[i] =="GABPB1+GABPB2"){
    tfs$Family[i] <- paste("ETS")
  }
  else if(tfs$TF[i] =="NFIA+NFIB+NFIC+NFIX"){
    tfs$Family[i] <- paste("Nuclear Factor")
  }
  else if(tfs$TF[i] =="TBP"){
    tfs$Family[i] <- paste("TBF_related")
  }
}

tfs <- tfs[,c("TF","Family")]

list_sigs <- as.data.frame(colnames(rel_map))
setnames(list_sigs, old = c('colnames(rel_map)'), new = c('sig'))

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
  else if(list_sigs$sig[i]=="SBS84"||list_sigs$sig[i]=="SBS85"){
    list_sigs$Signature_type[i] <- paste("Cytidine deaminase Related")
  }
}




list_sigs <- list_sigs %>% remove_rownames %>% column_to_rownames(var="sig")

rel_map_scale <- scale(rel_map)
pheatmap(rel_map_scale, annotation_row = tfs, annotation_col = list_sigs, fontsize_row = 2)

#Barplot of transcription factors most affected 

tf_freqs <- motif_df %>% group_by(TF) %>% summarise(Freq=n())

#380 Transcription factors affected 
#Top 20 TFs
#TF
#Freq
#FOXJ3	ARID3A IRF4 IRF5	FOXJ2	FUBP1 ONECUT2 FOXF1	
#HOXC6	KLF15	FOXQ1 FOXM1 IRF3	MEF2A	SP1	MEF2D	FOXO1	FOXD3
# STAT2	IRF2	IRF1#

#Adding transcription factor families to data 

#TF Families - ARID, FOX (Forkhead Box), IRF, ONECUT, HOX, KLF/SP, MEF2, SOX


promoter_mut_genes <- as.vector(unique(Sig_pro_overlap$Gene_name))

motif_dis_genes <- as.vector(unique(motif_df$Gene_name))
Names <- as.list("promoter_mut_genes", "motif_dis_genes")
venn.diagram(x = Names,filename = '#14_venn_diagramm.png',
             output=TRUE )
venn.diagram(
  x = list(motif_dis_genes, promoter_mut_genes),
  category.names = c("Set 1" , "Set 2 "),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)

#Analysis of frequency of tfs etc 

freq_tfs <- motif_df %>% group_by(Signature_type, TF, sig) %>% summarise(Freq=n())
freq_tfs_rid <- subset[freq_tfs, ]



#FINAL SCRIPTS FOR PLOTS 
#1. All - all signatures and all data 

sig_tf <- motif_df %>% group_by(TF, sig, Signature_type) %>% summarise(Freq=n())
#######adding log scale ##########
for(i in 1:nrow(sig_tf)){
  sig_tf$log_freq <- log10((sig_tf$Freq)+1)
}

sig_tf_mat <- acast(sig_tf, TF~sig, value.var="log_freq")
sig_tf_mat[is.na(sig_tf_mat)] <- 0

tfs <- as.data.frame(rownames(sig_tf_mat))
setnames(tfs, old = c("rownames(sig_tf_mat)"), new = c("TF"))
tfs <- merge(tfs, Homo_sapiens_TF, by=c("TF","TF"), all.x=T, all.y=F)
tfs <- tfs[,c("TF","Family")]
for(i in 1:nrow(tfs)){
  print(i)
  if(tfs$TF[i] == "BPTF"||tfs$TF[i] == "BRCA1"||tfs$TF[i] == "ENO1"||tfs$TF[i] == "FUBP1"||
     tfs$TF[i] == "HLTF"||tfs$TF[i] == "T"){
    tfs$Family[i] <- paste("Miscellaneous")
  }
  else if(tfs$TF[i] =="ATF2+ATF4"){
    tfs$Family[i] <- paste("TF_bZIP")
  }
  else if(tfs$TF[i] =="GABPB1+GABPB2"){
    tfs$Family[i] <- paste("ETS")
  }
  else if(tfs$TF[i] =="NFIA+NFIB+NFIC+NFIX"){
    tfs$Family[i] <- paste("Nuclear Factor")
  }
  else if(tfs$TF[i] =="TBP"){
    tfs$Family[i] <- paste("TBF_related")
  }
}

tfs <- tfs %>% remove_rownames %>% column_to_rownames(var="TF")

list_sigs <- as.data.frame(colnames(sig_tf_mat))
setnames(list_sigs, old = c('colnames(sig_tf_mat)'), new = c('sig'))
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

mypal = pal_npg("nrc")(7)
mypal = list(mypal)
list_sigs <- list_sigs %>% remove_rownames %>% column_to_rownames(var="sig")
sig_tf_mat_scale <- scale(sig_tf_mat)
pheatmap(sig_tf_mat, annotation_col = list_sigs, show_rownames = F, fontsize_col = 9, annotation_colors = mypal)
pheatmap(sig_tf_mat, annotation_row = tfs, annotation_col = list_sigs, show_rownames = F, fontsize_col = 9)


#2. Ageing Signature and SBS7a/b/c/d removed 
TFs <- motif_df %>% group_by(TF) %>% summarise(Freq=n())
TFs1 <- subset(sigs, Freq>2500)
TFs_to_keep <- as.character(TFs1$TF)

sig_tf <- NULL
sig_tf <-subset(motif_df, TF %in% TFs_to_keep)

sig_tf <- sig_tf %>% group_by(TF, sig, Signature_type) %>% summarise(Freq=n())


sig_tf_rd <- sig_tf[!(sig_tf$sig == "SBS40"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS5"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS1"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS37"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS7a"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS7b"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS7c"),]
sig_tf_rd <- sig_tf_rd[!(sig_tf_rd$sig =="SBS7d"),]



for(i in 1:nrow(sig_tf_rd)){
  sig_tf_rd$log_freq <- log10((sig_tf_rd$Freq)+1)
}

sig_tf_mat <- acast(sig_tf_rd, TF~sig, value.var="log_freq")
sig_tf_mat[is.na(sig_tf_mat)] <- 0

tfs <- as.data.frame(rownames(sig_tf_mat))
setnames(tfs, old = c("rownames(sig_tf_mat)"), new = c("TF"))
tfs <- merge(tfs, Homo_sapiens_TF, by=c("TF","TF"), all.x=T, all.y=F)
tfs <- tfs[,c("TF","Family")]
for(i in 1:nrow(tfs)){
  print(i)
  if(tfs$TF[i] == "BPTF"||tfs$TF[i] == "BRCA1"||tfs$TF[i] == "ENO1"||tfs$TF[i] == "FUBP1"||
     tfs$TF[i] == "HLTF"||tfs$TF[i] == "T"||tfs$TF[i]=="CXXC1"){
    tfs$Family[i] <- paste("Miscellaneous")
  }
  else if(tfs$TF[i] =="ATF2+ATF4"){
    tfs$Family[i] <- paste("TF_bZIP")
  }
  else if(tfs$TF[i] =="GABPB1+GABPB2"){
    tfs$Family[i] <- paste("ETS")
  }
  else if(tfs$TF[i] =="NFIA+NFIB+NFIC+NFIX"){
    tfs$Family[i] <- paste("Nuclear Factor")
  }
  else if(tfs$TF[i] =="TBP"){
    tfs$Family[i] <- paste("TBF_related")
  }
}
for(i in 1:nrow(tfs)){
  print(i)
  if(tfs$TF[i] == "CXXC1"){
    tfs$Family[i] <- paste("Miscellaneous")
  }
}

tfs <- tfs %>% remove_rownames %>% column_to_rownames(var="TF")

list_sigs <- as.data.frame(colnames(sig_tf_mat))
setnames(list_sigs, old = c('colnames(sig_tf_mat)'), new = c('sig'))
list_sigs$Signature_type <- paste("f")

for(i in 1:nrow(list_sigs)){
  if(list_sigs$sig[i]=="SBS3"||list_sigs$sig[i]== "SBS30"||
     list_sigs$sig[i]== "SBS10a"||list_sigs$sig[i]== "SBS10b"||
     list_sigs$sig[i]== "SBS36"||list_sigs$sig[i]== "SBS15"||list_sigs$sig[i]== "SBS44"){
    list_sigs$Signature_type[i] <- paste("DNA Repair Related")
  }
  else if(list_sigs$sig[i] == "SBS2"){
    list_sigs$Signature_type[i] <- paste("APOBEC Related")
  }
  else if(list_sigs$sig[i]=="SBS32"||list_sigs$sig[i] =="SBS35"||
          list_sigs$sig[i]== "SBS25"||
          list_sigs$sig[i]== "SBS11"||list_sigs$sig[i]== "SBS31"){
    list_sigs$Signature_type[i] <- paste("Treatment Related")
  }
  else if(list_sigs$sig[i] == "SBS42"||list_sigs$sig[i] == "SBS22"){
    list_sigs$Signature_type[i] <- paste("Endogenous Exposure Related")
  }
  else if(list_sigs$sig[i]=="SBS9"||list_sigs$sig[i]=="SBS16"
          ||list_sigs$sig[i]=="SBS19"||list_sigs$sig[i]=="SBS23"||list_sigs$sig[i]=="SBS28"
          ||list_sigs$sig[i]=="SBS41"||list_sigs$sig[i]=="SBS17b"
          ||list_sigs$sig[i]=="SBS38"||list_sigs$sig[i]== "SBS8"||list_sigs$sig[i]== "SBS39"
          ||list_sigs$sig[i]== "SBS33"){
    list_sigs$Signature_type[i] <- paste("Other/Unknown")
  }
  else if(list_sigs$sig[i]=="SBS84"||list_sigs$sig[i]=="SBS85"){
    list_sigs$Signature_type[i] <- paste("Cytidine deaminase Related")
  }
}
for(i in 1:nrow(list_sigs)){
  if(list_sigs$sig[i]=="SBS6"){
    list_sigs$Signature_type[i] <- paste("DNA Repair Related")
  }
} 
list_sigs <- list_sigs %>% remove_rownames %>% column_to_rownames(var="sig")
sig_tf_mat_scale <- scale(sig_tf_mat)
pheatmap(sig_tf_mat, annotation_row = tfs, 
         annotation_col = list_sigs, fontsize_col = 10, 
         fontsize_row = 8, 
         cluster_cols = T, 
         cluster_rows = T)

#3. "Relevant" genes only 

pc_genes <- motif_df[motif_df$Gene_type == "protein_coding",]
IG_C_genes <- motif_df[motif_df$Gene_type=="IG_C_gene",]
ID_D_gene <- motif_df[motif_df$Gene_type=="ID_D_gene",]
rRNA <- motif_df[motif_df$Gene_type=="rRNA",]
miRNA <- motif_df[motif_df$Gene_type=="miRNA",]
lncRNA <- motif_df[motif_df$Gene_type=="lncRNA",]
TR_C_gene <- motif_df[motif_df$Gene_type=="TR_C_gene",]
TR_D_gene <- motif_df[motif_df$Gene_type=="TR_D_gene",]
TR_J_gene <- motif_df[motif_df$Gene_type=="TR_J_gene",]
misc_RNA <- motif_df[motif_df$Gene_type=="misc_RNA",]

AR <- motif_df[motif_df$TF =="AR",]


motif_df3 <- rbind(pc_genes, ID_D_gene, IG_C_genes, rRNA, miRNA, lncRNA, TR_C_gene, TR_D_gene, TR_J_gene, misc_RNA)

un_pseudo <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="unprocessed_pseudogene",]
pseudo <- Sig_pro_overlap[Sig_pro_overlap$Gene_type=="processed_pseudogene",]
all_pseudo <- rbind(un_pseudo, pseudo)

rel_map <- motif_df3 %>% group_by(TF, sig, Signature_type) %>% summarise(Freq=n())
rel_map <- sig_tf[!(sig_tf$sig == "SBS40"),]
rel_map <- rel_map[!(rel_map$sig =="SBS5"),]
rel_map <- rel_map[!(rel_map$sig =="SBS1"),]
rel_map <- rel_map[!(rel_map$sig =="SBS37"),]
rel_map <- rel_map[!(rel_map$sig =="SBS7a"),]
rel_map <- rel_map[!(rel_map$sig =="SBS7b"),]
rel_map <- rel_map[!(rel_map$sig =="SBS7c"),]
rel_map <- rel_map[!(rel_map$sig =="SBS7d"),]
rel_map <- subset(rel_map, Freq>10)
sig_tf_mat <- acast(rel_map, TF~sig, value.var="Freq")
sig_tf_mat[is.na(sig_tf_mat)] <- 0

tfs <- as.data.frame(rownames(sig_tf_mat))
setnames(tfs, old = c("rownames(sig_tf_mat)"), new = c("TF"))
tfs <- merge(tfs, Homo_sapiens_TF, by=c("TF","TF"), all.x=T, all.y=F)
tfs <- tfs[,c("TF","Family")]
for(i in 1:nrow(tfs)){
  print(i)
  if(tfs$TF[i] == "BPTF"||tfs$TF[i] == "BRCA1"||tfs$TF[i] == "ENO1"||tfs$TF[i] == "FUBP1"||
     tfs$TF[i] == "HLTF"||tfs$TF[i] == "T"){
    tfs$Family[i] <- paste("Miscellaneous")
  }
  else if(tfs$TF[i] =="ATF2+ATF4"){
    tfs$Family[i] <- paste("TF_bZIP")
  }
  else if(tfs$TF[i] =="GABPB1+GABPB2"){
    tfs$Family[i] <- paste("ETS")
  }
  else if(tfs$TF[i] =="NFIA+NFIB+NFIC+NFIX"){
    tfs$Family[i] <- paste("Nuclear Factor")
  }
  else if(tfs$TF[i] =="TBP"){
    tfs$Family[i] <- paste("TBF_related")
  }
}

tfs <- tfs %>% remove_rownames %>% column_to_rownames(var="TF")

list_sigs <- as.data.frame(colnames(sig_tf_mat))
setnames(list_sigs, old = c('colnames(sig_tf_mat)'), new = c('sig'))
list_sigs$Signature_type <- paste("f")
for(i in 1:nrow(list_sigs)){
  if(list_sigs$sig[i]=="SBS44"){
    list_sigs$Signature_type[i] <- paste("DNA Repair Related")
  }
}

for(i in 1:nrow(list_sigs)){
  if(list_sigs$sig[i]=="SBS3"||list_sigs$sig[i]== "SBS30"||
     list_sigs$sig[i]== "SBS6"||list_sigs$sig[i]== "SBS10a"||list_sigs$sig[i]== "SBS10b"||
     list_sigs$sig[i]== "SBS36"||list_sigs$sig[i]== "SBS15"){
    list_sigs$Signature_type[i] <- paste("DNA Repair Related")
  }
  else if(list_sigs$sig[i] == "SBS1"|| list_sigs$sig[i]== "SBS5"||list_sigs$sig[i]=="SBS40"){
    list_sigs$Signature_type[i] <- paste("Ageing Related")
  }
  else if(list_sigs$sig[i] == "SBS2"){
    list_sigs$Signature_type[i] <- paste("APOBEC Related")
  }
  else if(list_sigs$sig[i]=="SBS32"||list_sigs$sig[i] =="SBS35"||
          list_sigs$sig[i]== "SBS32"||list_sigs$sig[i]== "SBS25"||
          list_sigs$sig[i]== "SBS11"||list_sigs$sig[i]== "SBS31"){
    list_sigs$Signature_type[i] <- paste("Treatment Related")
  }
  else if(list_sigs$sig[i] == "SBS42"||list_sigs$sig[i] == "SBS22"){
    list_sigs$Signature_type[i] <- paste("Endogenous Exposure Related")
  }
  else if(list_sigs$sig[i] == "SBS37"||list_sigs$sig[i]=="SBS9"||list_sigs$sig[i]=="SBS16"
          ||list_sigs$sig[i]=="SBS19"||list_sigs$sig[i]=="SBS23"||list_sigs$sig[i]=="SBS28"
          ||list_sigs$sig[i]=="SBS41"||list_sigs$sig[i]=="SBS17b"
          ||list_sigs$sig[i]=="SBS38"||list_sigs$sig[i]== "SBS8"||list_sigs$sig[i]== "SBS39"
          ||list_sigs$sig[i]== "SBS33"){
    list_sigs$Signature_type[i] <- paste("Other/Unknown")
  }
  else if(list_sigs$sig[i]=="SBS84"||list_sigs$sig[i]=="SBS85"){
    list_sigs$Signature_type[i] <- paste("Cytidine deaminase Related")
  }
}

list_sigs <- list_sigs %>% remove_rownames %>% column_to_rownames(var="sig")
sig_tf_mat_scale <- scale(sig_tf_mat)
pheatmap(sig_tf_mat_scale, annotation_row = tfs, annotation_col = list_sigs, fontsize = 6, fontsize_col = 10)

#4. Only "Strong" alterations 

motif_df_strong <- subset(motif_df, effect=="strong")

sig_tf <- motif_df_strong %>% group_by(TF, sig, Signature_type) %>% summarise(Freq=n())

sig_tf <- subset(sig_tf, Freq>100)
sig_tf_mat <- acast(sig_tf, TF~sig, value.var="Freq")
sig_tf_mat[is.na(sig_tf_mat)] <- 0

tfs <- as.data.frame(rownames(sig_tf_mat))
setnames(tfs, old = c("rownames(sig_tf_mat)"), new = c("TF"))
tfs <- merge(tfs, Homo_sapiens_TF, by=c("TF","TF"), all.x=T, all.y=F)
tfs <- tfs[,c("TF","Family")]
for(i in 1:nrow(tfs)){
  print(i)
  if(tfs$TF[i] == "BPTF"||tfs$TF[i] == "BRCA1"||tfs$TF[i] == "ENO1"||tfs$TF[i] == "FUBP1"||
     tfs$TF[i] == "HLTF"||tfs$TF[i] == "T"){
    tfs$Family[i] <- paste("Miscellaneous")
  }
  else if(tfs$TF[i] =="ATF2+ATF4"){
    tfs$Family[i] <- paste("TF_bZIP")
  }
  else if(tfs$TF[i] =="GABPB1+GABPB2"){
    tfs$Family[i] <- paste("ETS")
  }
  else if(tfs$TF[i] =="NFIA+NFIB+NFIC+NFIX"){
    tfs$Family[i] <- paste("Nuclear Factor")
  }
  else if(tfs$TF[i] =="TBP"){
    tfs$Family[i] <- paste("TBF_related")
  }
}

tfs <- tfs %>% remove_rownames %>% column_to_rownames(var="TF")

list_sigs <- as.data.frame(colnames(sig_tf_mat))
setnames(list_sigs, old = c('colnames(sig_tf_mat)'), new = c('sig'))
list_sigs$Signature_type <- paste("f")
for(i in 1:nrow(list_sigs)){
  if(list_sigs$sig[i]=="SBS3"||list_sigs$sig[i]== "SBS30"||
     list_sigs$sig[i]== "SBS44"||list_sigs$sig[i]== "SBS6"||
     list_sigs$sig[i]== "SBS15"){
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
  else if(list_sigs$sig[i]=="SBS84"||list_sigs$sig[i]=="SBS85"){
    list_sigs$Signature_type[i] <- paste("Cytidine deaminase Related")
  }
}


list_sigs <- list_sigs %>% remove_rownames %>% column_to_rownames(var="sig")
sig_tf_mat_scale <- scale(sig_tf_mat)
pheatmap(sig_tf_mat_scale, annotation_row = tfs, annotation_col = list_sigs, fontsize = 5)

#Looking at which transcription factors more affected 

transcription_factors <- motif_df_strong
transcription_factors <- transcription_factors %>% distinct(icgc_donor_id, TF, sig, .keep_all=T)

trans_factors <- transcription_factors %>% group_by(TF, sig) %>% summarise(Freq=n())


trans_subset <- subset(trans_factors, Freq>100)

ggplot(trans_subset, aes(y=Freq, x=reorder(TF,Freq), fill=sig)) + geom_col() + scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        # Remove panel border
        panel.border = element_blank(),  
        # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank())


transfacts1 <- motif_df %>% group_by(sig, Signature_type, icgc_donor_id) %>% summarise(Freq=n())

transfacts2 <- motif_df %>% group_by(sig, Signature_type) %>% summarise(Freq=n())

ggplot(transfacts2, aes(y=Freq, x=reorder(sig, Freq), fill=Signature_type)) + geom_col() + 
  scale_fill_brewer(palette = "Paired") + theme(axis.text.x=element_text(angle=90, hjust=1),
                                                # Remove panel border
                                                panel.border = element_blank(),  
                                                # Remove panel grid lines
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                # Remove panel background
                                                panel.background = element_blank()) + 
  ylab("Number of times signature associated with a motif alteration") + xlab("Signature")

ggplot(transfacts, aes(y=Freq, x=reorder(sig, Freq), fill=Signature_type)) + geom_boxplot() + 
   theme(axis.text.x=element_text(angle=90, hjust=1),
                                                           # Remove panel border
                                                           panel.border = element_blank(),  
                                                           # Remove panel grid lines
                                                           panel.grid.major = element_blank(),
                                                           panel.grid.minor = element_blank(),
                                                           # Remove panel background
                                                           panel.background = element_blank()) + 
          scale_fill_brewer(palette = "Paired") + coord_trans((y = "log1p"))
  
protein_coding <- 


