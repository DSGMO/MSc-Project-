library(MutationalPatterns)
library(data.table)
library(GenomicRanges)
library(BSgenome)
setwd('/Users/katefodder/Desktop/MSC_PROJECT')
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"

vcf<-as.data.frame(fread('simple_somatic_mutation.open.PRAD-CA.tsv'))

vcf<-vcf[vcf$mutation_type=='single base substitution',]

vcf_reduce<-vcf[,c('icgc_donor_id','chromosome','chromosome_start','reference_genome_allele','mutated_to_allele')]

list_samples<-NULL

for(i in unique(vcf_reduce$icgc_donor_id)){
  
  print(i)
  
  myvcf4<-vcf_reduce[vcf_reduce[,1]==i,-1]
  
  myvcf4<-myvcf4[myvcf4$reference_genome_allele!='-',]
  myvcf4<-myvcf4[myvcf4$mutated_to_allele!='-',]
  
  if(length(which(myvcf4$reference_genome_allele == myvcf4$mutated_to_allele))!=0)
  {	
    myvcf4<-myvcf4[-which(myvcf4$reference_genome_allele == myvcf4$mutated_to_allele),]
  }else{
    myvcf4<-myvcf4
  }
  myvcf4_parse<-data.frame(CHROM=myvcf4[,'chromosome'],
                           POS=myvcf4[,'chromosome_start'],
                           ID=paste(paste(myvcf4[,'chromosome'],myvcf4[,'chromosome_start'],sep=':'),paste(myvcf4[,'reference_genome_allele'],myvcf4[,'mutated_to_allele'],sep='/'),sep='_'),
                           REF=myvcf4[,'reference_genome_allele'],
                           ALT=myvcf4[,'mutated_to_allele'],
                           QUAL=rep('.',nrow(myvcf4)),
                           FILTER=rep('.',nrow(myvcf4)),
                           INFO=rep('.',nrow(myvcf4)),
                           FORMAT=rep('PASS',nrow(myvcf4)),
                           samples=rep('.',nrow(myvcf4)))
  
  colnames(myvcf4_parse)[ncol(myvcf4_parse)]<-i
  colnames(myvcf4_parse)[1]<-'#CHROM'
  
  setwd('/Users/katefodder/Desktop/MSC_PROJECT/vcf_files')
  writeLines('##fileformat=VCFv4.1',paste(i,'.prad-ca.vcf',sep=''))
  write.table(myvcf4_parse,file=paste(i,'.prad-ca.vcf',sep=''),sep='\t',row.names=F,quote=F,col.names=T,append=T)
  
  list_samples<-c(list_samples,i)
}
setwd('/Users/katefodder/Desktop/MSC_PROJECT/vcf_files')

write.table(list_samples,file='PRAD-CA_samples_list.txt',sep='\t',row.names=F,quote=F,col.names=F)


library(dplyr)
library(magrittr)

library(MutationalPatterns)
library(data.table)
library(GenomicRanges)
library(BSgenome)
library('BSgenome.Hsapiens.UCSC.hg19')
library(VariantAnnotation)

setwd('/Users/katefodder/Desktop/MSC_PROJECT')
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
sample_names<-as.character(read.delim(file='PRAD-CA_samples_list.txt',header=F)[,1])

setwd('/Users/katefodder/Desktop/MSC_PROJECT/vcf_files')
vcf_files <- NULL
vcf_files<-grep(dir(),pattern='vcf',value=T)
vcfs <- NULL
vcfs <- read_vcfs_as_granges(vcf_files[grep(vcf_files,pattern="PRAD",invert=T)], sample_names, ref_genome)

#Create context file for each samples

final_table_context_patients<-data.frame()

for(i in names(vcfs)){

print(i)

print(i)

context = mut_context(vcfs[[i]], ref_genome)

start_context<-substr(context,1,1)
end_context<-substr(context,4,4)

REF<-unlist(lapply(mcols(vcfs[[i]])$REF,as.character))
ALT<-unlist(lapply(mcols(vcfs[[i]])$ALT,as.character))

MUTATIONS<-paste('[',paste(REF,ALT,sep='>'),']',sep='')
VARIANT_CLASS<-paste(start_context,MUTATIONS,end_context,sep='')

table_context_patient<-data.frame(
	   CHROM=data.frame(vcfs[[i]])[,1],
	   POS=start(vcfs[[i]]),
	   REF=REF,
	   ALT=ALT,
	   VARIANT_CLASS=VARIANT_CLASS,
	   SAMPLE=rep(i,length(REF)))

final_table_context_patients<-rbind(final_table_context_patients,table_context_patient)

}

write.table(final_table_context_patients,file='PRAD-CA.snvs.vcf',sep='\t',row.names=F,quote=F)

setwd('/Users/katefodder/Desktop/MSC_PROJECT')

sigProfiler_SBS_signatures_2019_05_22 <- read.csv("~/Desktop/MSC_PROJECT/sigProfiler_SBS_signatures_2019_05_22.txt", header=FALSE)
cosmic_table <- sigProfiler_SBS_signatures_2019_05_22
colnames(cosmic_table) <- cosmic_table[1,]
cosmic_table <- cosmic_table[-1, ] 


cosmic_table2<-as.data.frame(cosmic_table[,-c(1:2)])


start_context<-substr(cosmic_table$SubType,1,1)
end_context<-substr(cosmic_table$SubType,2,2)

MUTATIONS<-paste('[',cosmic_table$Type,']',sep='')
VARIANT_CLASS<-paste(start_context,MUTATIONS,end_context,sep='')
rownames(cosmic_table2)<-VARIANT_CLASS

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

cosmic_table2 <- cosmic_table2 %>% mutate_if(is.character,as.numeric)
cosmic_table2 <- as.matrix(cosmic_table2)

cosmic_table2 <- as.data.frame(cosmic_table2)
rownames(cosmic_table2) <- rownames(cosmic_table)
#General 

type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
p1 <- plot_spectrum(type_occurrences)
p4 <- plot_spectrum(type_occurrences, CT = TRUE)
#next stage after mut_mat

pdf('test.pdf')
plot_96_profile(cosmic_table4[,1:4], condensed = TRUE, ymax = 0.4)
dev.off()

mut_mat <- mut_mat + 0.0001
library("NMF")
estimate <- nmf(mut_mat, rank=4, method="brunet", nrun=10, seed=144456)

plot_96_profile(mut_mat[,c(1,2)], condensed = T)

estimate <- nmf(mut_mat, rank=4:6, method="brunet", nrun=10, seed=144456)
plot(estimate)
nmf_res5 <- extract_signatures(mut_mat, rank = 5, nrun = 10)
rank=4

estimate4 <- nmf(mut_mat, rank=4, method="brunet")

setwd('/home/guidantoniomt/kate/PRAD.snvs')

# get matrix W
w <- basis(estimate4)
dim(w)
colnames(w)<-paste('x',1:ncol(w))
write.table(w,file=paste('processes_',rank,sep=''),sep='\t',row.names=T,quote=F,col.names=T)

wavg<-t(data.frame(apply(w,4,mean)))
colnames(wavg)<-paste('x',1:ncol(wavg))
write.table(wavg,file=paste('processesStabAvg_',rank,sep=''),sep='\t',row.names=F,quote=F,col.names=T)

# get matrix H
h <- coef(estimate4)
dim(h)
write.table(h,file=paste('exposures_',rank,sep=''),sep='\t',row.names=F,quote=F,col.names=T)
write.table(h,file=paste('exposures_fitting_',rank,sep=''),sep='\t',row.names=F,quote=F,col.names=T)




sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
                + "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
cosmic_table2 <- as.data.frame(cosmic_table2)


#Chose rank 5 

nmf_res <- extract_signatures(mut_mat, rank = 5, nrun = 10)
plot_contribution(nmf_res$contribution[,1:22], nmf_res$signature, mode="relative") + 
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8))

fit_res <- fit_to_signatures(mut_mat, cosmic_table4)
fit_res <- fit_to_signatures(mut_mat, cosmic_table4)
#Results with mutational Patterns
sigs_cosmic <- sigs_cosmic %>% mutate_if(is.character,as.numeric)
sigs_cosmic <- as.matrix(sigs_cosmic)
fit_res <- fit_to_signatures(mut_mat, sigs_cosmic)

mp_sigs <- as.data.frame(t(fit_res$contribution))

signatures.cosmic <- as.data.frame(signatures.cosmic)
sigs_cosmic <- as.data.frame(t(signatures.cosmic))


#Genomic Distribution 

chromosomes <- seqnames(get(ref_genome))[1:22]
plot_rainfall(vcfs[[1]], title = names(vcfs[5]),
                  chromosomes = chromosomes, cex = 1.5, ylim = 1e+09)

#Fit res code 
sigs_cosmic <- as.matrix(t(signatures.cosmic))
cosmic_table2 <- as.matrix(cosmic_table4)
fit_res <- fit_to_signatures(mut_mat, sigs_cosmic)
select <- which(rowSums(fit_res$contribution) > 10)
plot_contribution(fit_res$contribution[select,1:10], cosmic_table2[select,],coord_flip = FALSE,mode = "absolute") + theme(axis.text.x=element_text(angle=90, hjust=1))

plot_contribution_heatmap(fit_res$contribution,
                          cluster_samples = TRUE,
                          method = "complete") + theme(axis.text.x=element_text(angle=90, hjust=1, size=0.5), axis.text.y=element_text(size=0.5))


mp_signatures <- as.data.frame(fit_res$contribution)


# pie charts 

first_don <- as.data.frame(mp_signatures$DO229393)
first_don$signature <- rownames(mp_signatures)
first_don$sig <- first_don$signature
setnames(first_don, old=c("mp_signatures$DO229393"), new=c("weight"))
ggplot(first_don, aes(x="", y=weight, fill=sig)) + geom_bar(stat = "identity", width=1) + 
  coord_polar("y", start=0) + scale_fill_aaas() + theme_void() + ggtitle("D0229393")

for(i in 1:nrow(first_don)){ 
  if(first_don$weight[i]<100){
    first_don$sig[i] <- paste("Other")
  }
}
    
first_don = as.data.table(first_don)
setkey(first_don, sig)
first_don <- first_don[!"Other"]

ggplot(first_don, aes(x="", y=weight, fill=sig)) + geom_bar(stat = "identity", width=1) + 
  coord_polar("y", start=0)

### Second Donor #####

second_don <- as.data.frame(mp_signatures$DO229392)
second_don$signature <- rownames(mp_signatures)
second_don$sig <- second_don$signature
setnames(second_don, old=c("mp_signatures$DO229392"), new=c("weight"))
ggplot(second_don, aes(x="", y=weight, fill=sig)) + geom_bar(stat = "identity", width=1) + 
  coord_polar("y", start=0) + scale_fill_lancet() + theme_void() + ggtitle("D0229392")

for(i in 1:nrow(second_don)){ 
  if(second_don$weight[i]<100){
    second_don$sig[i] <- paste("Other")
  }
}

second_don = as.data.table(second_don)
setkey(second_don, sig)
first_don <- first_don[!"Other"]

ggplot(first_don, aes(x="", y=weight, fill=sig)) + geom_bar(stat = "identity", width=1) + 
  coord_polar("y", start=0)

####### third donor ######
third_don <- as.data.frame(mp_signatures$DO229394)
third_don$signature <- rownames(mp_signatures)
third_don$sig <- third_don$signature
setnames(third_don, old=c("mp_signatures$DO229394"), new=c("weight"))
ggplot(third_don, aes(x="", y=weight, fill=sig)) + geom_bar(stat = "identity", width=1) + 
  coord_polar("y", start=0) + scale_fill_lancet() + theme_void() + ggtitle("D0229394")

for(i in 1:nrow(third_don)){ 
  if(third_don$weight[i]<400){
    third_don$sig[i] <- paste("Other")
  }
}



mp <- as.data.frame(t(mp_signatures))
ds <- allsignatures
install.packages("xlsx")
write.xlsx(x, file, sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
library(xlsx)
write.xlsx(mp, "mp.xlsx", col.names=T, row.names =F)


mp_df <- tibble::rownames_to_column(mp, "Donor")
ds_df <- tibble::rownames_to_column(ds, "Donor")
library(janitor)
mp_sum <- mp_df %>%
  adorn_totals("row")
ds_sum <- ds_df %>% adorn_totals("row")

mp_sum <- as.data.frame(t(mp_sum))
ds_sum <- as.data.frame(t(ds_sum))

mp_sum <- mp_sum[,c(307)]
ds_sum <- ds_sum[,c(307)]

mp_sum <- as.data.frame(mp_sum)
ds_sum <- as.data.frame(ds_sum)

mp_sum <- mp_sum[-1,]
ds_sum <- ds_sum[-1,]

mp_sum <- as.data.frame(mp_sum)
ds_sum <- as.data.frame(ds_sum)

rownames(mp_sum) <- rownames(signatures.cosmic)

combined_sigs <- cbind(ds_sum, mp_sum)
setnames(combined_sigs, old=c("mp_sum","ds_sum"), new=c("MutPat", "DecSigs"))

combined_sigs$MutPat <- as.numeric(combined_sigs$MutPat)
combined_sigs$DecSigs <- as.numeric(combined_sigs$DecSigs)
combined_sigs <- tibble::rownames_to_column(combined_sigs, "Signature")

ggplot(combined_sigs, aes(y=MutPat, x=DecSigs, fill=Signature)) + 
  geom_point(aes(fill=Signature)) + 
  geom_smooth(method = "auto") 

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

library(ggraph)
library(ggpubr)
ggscatter(combined_sigs, x="MutPat", y="DecSigs", fill ="Signature", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) + ylab("deconstructSigs Weights") + 
  xlab("MutationalPatterns Weights") 
  

#removing the most prevalent signatures to see whats happening 

combined_sigs2 <- combined_sigs

combined_sigs2 = as.data.table(combined_sigs2)
setkey(combined_sigs2, Signature)
combined_sigs2 <- combined_sigs2[!"Signature.1"]
combined_sigs2 <- combined_sigs2[!"Signature.3"]
combined_sigs2 <- combined_sigs2[!"Signature.5"]
combined_sigs2 <- combined_sigs2[!"Signature.8"]
combined_sigs2 <- combined_sigs2[!"Signature.9"]
combined_sigs2 <- combined_sigs2[!"Signature.16"]

ggscatter(combined_sigs2, x="MutPat", y="DecSigs", fill ="Signature", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) + ylab("deconstructSigs Weights") + 
  xlab("MutationalPatterns Weights") 
  




ds_sig1 <- as.data.frame(ds_df$Signature.1)
mp_sig1 <- as.data.frame(mp_df$Signature.1)

sig1 <- cbind(mp_sig1, ds_sig1)

setnames(sig1, old=c("ds_df$Signature.1", "mp_df$Signature.1"), new=c("DeconstructSigs", "MutationalPatterns"))
#SIG5
ds_sig8 <- as.data.frame(ds_df$Signature.8)
mp_sig8 <- as.data.frame(mp_df$Signature.8)

sig8 <- cbind(mp_sig8, ds_sig8)

setnames(sig8, old=c("ds_df$Signature.8", "mp_df$Signature.8"), new=c("DeconstructSigs", "MutationalPatterns"))
ggplot(sig8, aes(y=MutationalPatterns, x=DeconstructSigs)) + geom_jitter()

#Strand Bias 
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
strand = mut_strand(vcfs[[1]], genes_hg19)
mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19)

strand_counts <- strand_occurrences(mut_mat_s)
strand_bias <- strand_bias_test(strand_counts)
                                
ps1 <- plot_strand(strand_counts, mode = "relative")                                

#Fit res with 40 Signatures
sigs_cosmic <- as.matrix(t(signatures.cosmic))
fit_res4 <- fit_to_signatures(mut_mat, sigs_cosmic)
mp_allsignatures <- as.data.frame(fit_res4$contribution)

