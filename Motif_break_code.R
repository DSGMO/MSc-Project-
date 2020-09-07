#Motifbreak code 
install.packages("motifbreakR")
library("motifbreakR")
library("BSgenome")

#Creating bed file to input into package 

mutation_input <- Sig_pro_overlap
mutation_parse <- mutation_input[,c("Promoter_chromosome","Gene_name","Promoter_strand","mutation_start","mutation_end","ref_sig","alt_sig","icgc_donor_id","sig")]
mutation_parse$Promoter_chromosome <- sub("","chr", mutation_parse$Promoter_chromosome)
mutation_parse$name <- paste(mutation_parse$Promoter_chromosome, mutation_parse$mutation_end, mutation_parse$ref_sig, mutation_parse$alt_sig, sep = ':')
mutation_bed <- mutation_parse[,c("Promoter_chromosome","mutation_start","mutation_end","name","Promoter_strand")]
mutation_bed$score <- paste("0")
mutation_bed <- mutation_bed[,c("Promoter_chromosome","mutation_end","mutation_start","name","score","Promoter_strand")]
save(mutation_bed, file = "mutation_bed.txt", header = F)

mutation_bed <- mutation_bed[,c(1,2,3,4,6,5)]
head(mutation_bed)

#motif break code 

library(motifbreakR)
library(data.table)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37) # dbSNP142 in hg19
library(BSgenome.Hsapiens.UCSC.hg19)     # hg19 genome
library(MotifDb)

#
# GMT: the input is the matrix with the promoters (and mutations), associated with mutational signatures)
#
setwd("/Users/katefodder/Desktop/MSC_PROJECT")
read.table("simple_somatic_mutation.open.PRAD_CA.tsv")
mut_sig<-fread("simple_somatic_mutation.open.PRAD-CA.tsv",data.table=F)
mut_sig <- Sig_pro_overlap
mut_sig$Promoter_chromosome<-paste("chr",mut_sig$Promoter_chromosome,sep="")

mut_sig<-mut_sig[grep(mut_sig$Promoter_chromosome,pattern="chrM",invert=T),]

mut_sig<-mut_sig[mut_sig$ref_sig%in%c("A","C","G","T") & mut_sig$alt_sig%in%c("A","C","G","T"),]

#
# GMT: replace the name of the columns
#
id_var<-paste(mut_sig[,"chromosome"],
              mut_sig[,""],
              mut_sig[,"ref_sig"],
              mut_sig[,"alt_sig"],
              sep=":")

mut_sig$id_var<-id_var

#
# GMT: create the input for snps.from.file from motifbreakR
#
mut_rid<-data.frame(mut_sig[,c("chromosome",
                               "chromosome_start",
                               "chromosome_end",
                               "id_var")],
                    score=0,
                    strand="+")

mut_rid2<-mut_rid


#snps.from.file import the end column, 
#in this way the position of the variants does not change and remain the same
#see object snps.from.file, the ranges are the same of the start values of mut_rid
mut_rid2$chromosome_start<-mut_rid2$chromosome_start-1 # start = start -1 (not real start coordinate)
mut_rid2$chromosome_end<-mut_rid2$chromosome_start+1   # end = stat+1 (become the original start)

#
# GMT: only 20000 variants to test the code
#
write.table(mut_rid2[1:20000,],file="input_for_motif_breaker.bed",sep="\t",row.names=F,quote=F,col.names=F)

setwd("/Users/katefodder/Desktop/MSC_PROJECT")

snps.mb.frombed <- snps.from.file(file = "input_for_motif_breaker.bed",
                                  search.genome = BSgenome.Hsapiens.UCSC.hg19,
                                  format = "bed")
snps.mb.frombed <- snps.from.file(file = "mutation_bed_for_motif_breaker.bed", 
                                 search.genome = BSgenome.Hsapiens.UCSC.hg19, 
                                 format = "bed")
#
# GMT: see if PWM matrices from Encode are the best solution
#
data(encodemotif)

results <- motifbreakR(snpList = snps.mb.frombed[1:5], filterp = TRUE,
                       pwmList = encodemotif, # check what is the best database
                       threshold = 0.85,
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())

# Katherine: create an iterative process to compute the p-values for each variant

rs1006140 <- results[names(results) %in% "rs1006140"]
calculatePvalue(rs1006140)

# Katherine: add the p-values to the data.frame with the results
ranges_df<-as.data.frame(results@ranges)
motif_df<-as.data.frame(results@elementMetadata)

res_tot<-cbind(ranges_df,motif_df)
write.table(res_tot,file="motif_breaker_test.txt",row.names=F,quote=F,sep="\t")





# Katherine: create plots as in the vignette of the most relevant events (most distrupted motif)



