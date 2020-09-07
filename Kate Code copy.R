#Genomic Ranges Code

library("GenomicRanges")
#Creating range for Promoter regions 
promoters <- makeGRangesFromDataFrame(Genes, keep.extra.columns = TRUE,
                                      ignore.strand = TRUE, seqinfo = NULL, 
                                      seqnames.field = c("Chromosome"), 
                                      start.field = "Promoter_start", 
                                      end.field = c("Promoter_end", "stop"), 
                                      starts.in.df.are.0based = TRUE
)

#Sorting X and Y Chromosomes 
simple_somatic_mutation.open.PRAD.CA$chromosome <- gsub(simple_somatic_mutation.open.PRAD.CA$chromosome,pattern="X",replacement=23)
simple_somatic_mutation.open.PRAD.CA$chromosome <- gsub(simple_somatic_mutation.open.PRAD.CA$chromosome,pattern="Y",replacement=24)

#Creating mutation positions 
mutations <- makeGRangesFromDataFrame(simple_somatic_mutation.open.PRAD.CA, keep.extra.columns = TRUE, 
                                      ignore.strand = TRUE, seqinfo = NULL, 
                                      seqnames.field = c("chromosome"), 
                                      start.field = "chromosome_start",
                                      end.field = c("chromosome_start","stop"),
                                      starts.in.df.are.0based = TRUE)

# Finding overlap 
findOverlaps(promoters, mutations)

overlap <- findOverlaps(promoters, mutations)

indices_subject <- subjectHits(overlap)
indices_query <- queryHits(overlap)

indices_subject[1:10]
promoters_df<-as.data.frame(promoters)
mutations_df<-as.data.frame(mutations)

promoters_ovp<-promoters_df[indices_query,]
mutations_ovp<-mutations_df[indices_subject,]

res_overlap<-cbind(promoters_ovp,mutations_ovp)
res_overlap$position <- res_overlap$start

#adding column with mutation to dataframe 
res_overlap$alteration <- paste(res_overlap$mutated_from_allele, res_overlap$mutated_to_allele, sep = '>')



head(res_overlap)
str(res_overlap)

res_overlap[,c(1:10)]
head(res_overlap[1:10,c(1:10)])

colnames(res_overlap)
head(res_overlap[1:10,c(1:30)])

mutations[mutations$icgc_donor_id=="DO51070" & mutations$icgc_mutation_id=="MU31677070",]

nrow(res_overlap)
res_overlap_no_syn <- res_overlap[res_overlap$consequence_type!="synonymous_variant",]

write.table(res_overlap, file = "output_annotation_syn.txt", sep="\t", row.names = F, quote = F)
freq_genes<-data.frame(table(res_overlap$Gene_name))
head(freq_genes)
freq_genes2<- freq_genes[order(freq_genes[,2], decreasing = T),]
write.table(res_overlap_no_syn,file="output_annotation.txt",sep="\t",row.names=F,quote=F)

#Assigning Mut-sigs with Signatures 

setwd('/Users/katefodder/Desktop/MSC_PROJECT')

sig_mut_files <- list.files(path = "~/Desktop/MSC_PROJECT/PRAD-CA-ICGC_assigment_mutsig_with_muts")
sig_mut_files <-gsub(".tsv", "", sig_mut_files)

library(data.table)
for(i in sig_mut_files){
  print(i)
  filepath <- file.path("~/Desktop/MSC_PROJECT/PRAD-CA-ICGC_assigment_mutsig_with_muts",paste(i,".tsv",sep=""))
  assign(i, read.delim(filepath, header = FALSE, col.names = c("chromosome", "position","ref", "alt", "icgc_donor_id"),
                       sep = "\t"))
  
}


SBS1$Signature <- paste('SBS1')
SBS2$Signature <- paste('SBS2')
SBS3$Signature <- paste('SBS3')
SBS5$Signature <- paste('SBS5')
SBS6$Signature <- paste('SBS6')
SBS7a$Signature <- paste('SBS7a')
SBS7b$Signature <- paste('SBS7b')
SBS7c$Signature <- paste('SBS7c')
SBS7d$Signature <- paste('SBS7d')
SBS84$Signature <- paste('SBS84')
SBS85$Signature <- paste('SBS85')
SBS9$Signature <- paste('SBS9')
SBS10b$Signature <- paste('SBS10b')
SBS11$Signature <- paste('SBS11')
SBS15$Signature <- paste('SBS15')
SBS16$Signature <- paste('SBS16')
SBS17a$Signature <- paste('SBS17a')
SBS17b$Signature <- paste('SB17b')
SBS18$Signature <- paste('SBS18')
SBS19$Signature <- paste('SBS19')
SBS20$Signature <- paste('SBS20')
SBS22$Signature <- paste('SBS22')
SBS23$Signature <- paste('SBS23')
SBS24$Signature <- paste('SBS24')
SBS25$Signature <- paste('SBS25')
SBS26$Signature <- paste('SBS26')
SBS28$Signature <- paste('SBS28')
SBS30$Signature <- paste('SBS30')
SBS31$Signature <- paste('SBS31')
SBS32$Signature <- paste('SBS32')
SBS35$Signature <- paste('SBS35')
SBS36$Signature <- paste('SBS36')
SBS37$Signature <- paste('SBS37')
SBS38$Signature <- paste('SBS38')
SBS39$Signature <- paste('SBS39')
SBS40$Signature <- paste('SBS40')
SBS41$Signature <- paste('SBS41')
SBS42$Signature <- paste('SBS42')
SBS43$Signature <- paste('SBS43')
SBS44$Signature <- paste('SBS44')
SBS46$Signature <- paste('SBS46')
SBS47$Signature <- paste('SBS47')
SBS50$Signature <- paste('SBS50')
SBS51$Signature <- paste('SBS51')
SBS52$Signature <- paste('SBS52')
SBS53$Signature <- paste('SBS53')
SBS54$Signature <- paste('SBS54')
SBS55$Signature <- paste('SBS55')
SBS57$Signature <- paste('SBS57')
SBS58$Signature <- paste('SBS58')
SBS59$Signature <- paste('SBS59')
SBS60$Signature <- paste('SBS60')

all_sig_mut_files = NULL
all_sig_mut_files <- rbind(SBS1, SBS10b, SBS11, SBS15, SBS16, SBS17a, SBS17b, SBS18, SBS19, SBS2, SBS20, SBS22, SBS23, SBS24, SBS25, SBS26, SBS28, SBS3, SBS30, SBS31, SBS32, SBS35, SBS36, SBS37, SBS39, SBS38, SBS40, SBS41, SBS41, SBS42, SBS43, SBS44, SBS46, SBS47, SBS50, SBS51, SBS52, SBS53, SBS54, SBS55, SBS57, SBS58, SBS59, SBS60, SBS6, SBS60, SBS7a, SBS7b, SBS7c, SBS7d, SBS84, SBS85, SBS9)
setwd("/Users/katefodder/Desktop/MSC_PROJECT")
save(all_sig_mut_files, file="all_sig_mut_files.RData")


#Merging promoter data with signature data 
setwd('/Users/katefodder/Desktop/MSC_PROJECT')
merged_sigs_promoters2 <- merge(res_overlap, all_sig_mut_files, by.x=c("icgc_donor_id","position","mutated_from_allele","mutated_to_allele"),by.y=c("icgc_donor_id","position","ref","alt"),all.x=FALSE, all.y=FALSE)

save(merged_sigs_promoters, file="merged_sigs_promoters.RData")




