library(GenomicRanges)
setwd('/Users/katefodder/Desktop/MSC_PROJECT')
`mart_export.(6)` <- read.delim("~/Desktop/MSC_PROJECT/mart_export (6).txt")
Genes <- `mart_export.(6)`
Genes <- Genes[,c(2,3,4,5,6,7)]
colnames(Genes)[1:6]<-c("Chromosome", "Gene_start","Gene_end", "Strand", "Gene_name","Gene_type")

Genes$Strand <- gsub(Genes$Strand, pattern = "-1", replacement = "negative")
Genes$Strand <- gsub(Genes$Strand, pattern = "1", replacement = "positive")



for(i in 1:nrow(Genes)){
  if(Genes$Strand[i]=="negative"){
    Genes$Promoter_end[i] <- paste((Genes$Gene_end[i])+3000)
    Genes$Promoter_start[i] <- paste((Genes$Gene_end)[i]-2000)
  }
  if(Genes$Strand[i]=="positive"){
    Genes$Promoter_start[i] <- paste((Genes$Gene_start)[i]-3000)
    Genes$Promoter_end[i]<- paste((Genes$Gene_start)[i]+2000)
  }
}

Genes$Chromosome <- sub(pattern = "X", "23", Genes$Chromosome)
Genes$Chromosome <- sub(pattern = "Y", "24", Genes$Chromosome)

#creating Bed file to check on promoter regions 
# in write up - mention check with IGV (promoter regions) - 'to check the validity of this - promoter regions were checked against gene start sites using the IGV track annotation tool -then reference' 

Genes$Promoter_end = as.numeric(Genes$Promoter_end)
Genes$Promoter_start = as.numeric(Genes$Promoter_start)
Genes_parse<-Genes[,c("Chromosome","Promoter_start","Promoter_end","Gene_name")]
Genes_parse2 <- Genes_parse
Genes_parse2$Chromosome <- sub("", "chr", Genes_parse2$Chromosome)
Genes$Chromosome <- sub(pattern = "23", "X", Genes$Chromosome)
Genes$Chromosome <- sub(pattern = "24", "Y", Genes$Chromosome)
Genes_parse2 <- Genes_parse2[,c("Chromosome", "Promoter_start", "Promoter_end")]

rownames(Genes_parse2) = NULL
colnames(Genes_parse2) <- NULL
Genes_parse2 = as.data.frame(Genes_parse2)

GuidPromoters <- as.data.frame(read.table("promoters.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))


write.table(Genes_parse2, "Promoter_check.txt", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#Assigning Mut-sigs with Signatures 
setwd("~/Desktop/MSC_PROJECT/PRAD-CA-ICGC_assigment_mutsig_with_muts")

sig_mut_files<-grep(grep(dir(),pattern=".tsv",value=T),pattern="SBS",value=T)

all_SIGS<-data.frame()

for(i in sig_mut_files){
  
  current_sig<-gsub(i,pattern=".tsv",replacement="")
  
  current_sig_tab<-read.delim(i, header = FALSE, col.names = c("chromosome", "position","ref", "alt", "icgc_donor_id"),
                              sep = "\t",stringsAsFactors = F)
  
  if(nrow(current_sig_tab)!=0){
    
    current_sig_tab2<-data.frame(current_sig_tab,sig=current_sig)
    
    
    all_SIGS<-rbind(all_SIGS,current_sig_tab2)
  }
  
}

setwd("/Users/katefodder/Desktop/MSC_PROJECT")
save(all_SIGS, file="all_sig_mut_files.RData")

load("/Users/katefodder/Desktop/MSC_PROJECT/all_sig_mut_files.RData")
colnames(all_SIGS)[1:4]<-c("chromosome_sig","start_sig","ref_sig","alt_sig")
all_SIGS[,1]<-gsub(gsub(gsub(all_SIGS[,1],pattern="chr",replacement=""),pattern="X",replacement=23),pattern="Y",replacement=24)

#Genomic Ranges packages 
Genes_parse<-Genes[,c("Chromosome","Promoter_start","Promoter_end","Gene_name")]
Genes_parse$Chromosome <- gsub(Genes_parse$Chromosome, pattern='X',replacement=23)
Genes_parse$Chromosome <- gsub(Genes_parse$Chromosome, pattern='Y',replacement=24)
promoters <- makeGRangesFromDataFrame(Genes, keep.extra.columns = TRUE,
                                      ignore.strand = TRUE, seqinfo = NULL, 
                                      seqnames.field = c("Chromosome"), 
                                      start.field = "Promoter_start", 
                                      end.field = c("Promoter_end"), 
                                      starts.in.df.are.0based = TRUE)

Genes$Chromosome <- gsub(Genes$Chromosome, pattern='X', replacement=23)
Genes$Chromosome <- gsub(Genes$Chromosome, pattern='Y', replacement=24)

#Sorting X and Y Chromosomes 
simple_somatic_mutation.open.PRAD.CA$chromosome <- gsub(simple_somatic_mutation.open.PRAD.CA$chromosome,pattern="X",replacement=23)
simple_somatic_mutation.open.PRAD.CA$chromosome <- gsub(simple_somatic_mutation.open.PRAD.CA$chromosome,pattern="Y",replacement=24)

mut_rid<-simple_somatic_mutation.open.PRAD.CA[,c("chromosome","chromosome_start","chromosome_end","mutated_from_allele","mutated_to_allele","icgc_donor_id")]

#Creating Sigs positions 

all_SIGS$chromosome <- gsub("chr","", all_SIGS$chromosome)
all_SIGS$chromosome <- gsub("X", "23", all_SIGS$chromosome)
all_SIGS$chromosome <- gsub("Y", "24", all_SIGS$chromosome)


sigs <- makeGRangesFromDataFrame(all_SIGS, keep.extra.columns = TRUE, 
                                 ignore.strand = TRUE, seqinfo = NULL, 
                                 seqnames.field = c("chromosome_sig"),
                                 start.field = c("start_sig"),
                                 end.field = c("start_sig"),
                                 starts.in.df.are.0based = TRUE)

#Finding overlap of promoters and Signatures 
findOverlaps(promoters, sigs)
overlap_sigs_proms <- findOverlaps(sigs, promoters)

SIGSindices_subject <- subjectHits(overlap_sigs_proms)
SIGSindices_query <- queryHits(overlap_sigs_proms)

Sig_pro_overlap <- data.frame(subject = promoters[SIGSindices_subject], query = sigs[SIGSindices_query])

save(Sig_pro_overlap, file = "Sig_pro_overlap.RData")

head(Sig_pro_overlap)

#cleaning up Sig_pro_overlap File t 
library(data.table)


setnames(Sig_pro_overlap, old = c('subject.seqnames', 'subject.start', 'subject.end', 'subject.width', 'subject.strand','subject.Gene_Stable_ID','subject.Gene_start','subject.Gene_end', 'subject.Strand', 'subject.Gene_name', 'subject.Gene_type', 'query.seqnames','query.start','query.end','query.width','query.strand','query.ref_sig','query.alt_sig','query.icgc_donor_id','query.sig','query.chromosome'), new = c('Promoter_chromosome','Promoter_start','Promoter_end','Promoter_width','Promoter_strandxxx','Gene_Stable_ID', 'Gene_start','Gene_end', 'Gene_strand', 'Gene_name','Gene_type','mutation_chromosome','mutation_start','mutation_end','mutation_width','mutation_strand','ref_sig','alt_sig','icgc_donor_id','sig','chromosome'), skip_absent = T)
head(Sig_pro_overlap)

Sig_pro_overlap <- Sig_pro_overlap[,-5,8]

Sig_pro_overlap$Promoter_strand <- gsub("negative","-", Sig_pro_overlap$Promoter_strand)
Sig_pro_overlap$Promoter_strand <- gsub("positive","+", Sig_pro_overlap$Promoter_strand)

save(Sig_pro_overlap, file = "Sig_pro_overlap.RData")
Sig_pro_overlap <- Sig_pro_overlap[,c("Promoter_chromosome","Gene_name","Gene_type","mutation_start","mutation_end","ref_sig","alt_sig","icgc_donor_id","sig","alteration","signature_type")]


Mut_sig_pro_overlap <- Sig_pro_overlap[,-1,2]

#Looking at frequency of mutations in genes 
freq_genes<-data.frame(table(Sig_pro_overlap$subject.Gene_name))

head(freq_genes)

freq_genes2<- freq_genes[order(freq_genes[,2], decreasing = T),]
table(freq_genes2)

#Analysis of signatures present 

head(Sig_pro_overlap)
freq_sigs <- data.frame(table(Sig_pro_overlap$query.sig))
head(freq_sigs)
freq_sigs2<- freq_sigs[order(freq_sigs[,2])]
table(freq_sigs2)
