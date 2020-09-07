#Stage 3 - Motif_break Code 

library(motifbreakR)
library(data.table)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)    
library(MotifDb)


mut_sig <- Sig_pro_overlap[,c("Promoter_chromosome", "Promoter_start","Promoter_end","Gene_name", "Gene_strand", "mutation_start","mutation_end","ref_sig","alt_sig","icgc_donor_id","sig")]

mut_sig <- setNames(mut_sig, c("chromosome","pstart","pend", "Gene_name", "strand" ,"chromosome_start","chromosome_end","mutated_from_allele","mutated_to_allele","icgc_donor_id","sig"))
mut_sig <- mut_sig[,c("chromosome","strand","chromosome_start","chromosome_end","mutated_from_allele","mutated_to_allele")]
mut_sig$chromosome<-paste("chr",mut_sig$chromosome,sep="")
mut_sig<-mut_sig[grep(mut_sig$chromosome,pattern="chrM",invert=T),]

mut_sig$chromosome_start<-mut_sig$chromosome_start-2

mut_sig<-mut_sig[mut_sig$mutated_from_allele%in%c("A","C","G","T") & mut_sig$mutated_to_allele%in%c("A","C","G","T"),]
id_var<-paste(mut_sig[,"chromosome"],
              mut_sig[,"chromosome_end"],
              mut_sig[,"mutated_from_allele"],
              mut_sig[,"mutated_to_allele"],
              sep=":")

mut_sig$id_var<-id_var

mut_sig$strand <- sub("negative", "-", mut_sig$strand)
mut_sig$strand <- sub("positive", "+", mut_sig$strand)

mut_rid<-data.frame(mut_sig[,c("chromosome",
                               "chromosome_start",
                               "chromosome_end",
                               "id_var")],
                    score=0,
                    strand="+")

mut_rid$chromosome_end = as.numeric(mut_rid$chromosome_end)
mut_rid$chromosome_start = as.numeric(mut_rid$chromosome_start)
mut_rid2 <- mut_rid
 

mut_rid2$chromosome <- sub("chr23", "chrX", mut_rid2$chromosome)
mut_rid2$id_var <- sub("chr23", "chrX", mut_rid2$id_var)

mut_rid2$chromosome <- sub("chr24", "chrY", mut_rid2$chromosome)
mut_rid2$id_var <- sub("chr24", "chrY", mut_rid2$id_var)


setwd("/Users/katefodder/Desktop/MSC_PROJECT")
write.table(mut_rid2[1:2,],file="input_for_motif_breaker.bed",sep="\t",row.names=F,quote=F,col.names=F)

snps.mb.frombed <- snps.from.file(file = "input_for_motif_breaker.bed",
                                  search.genome = BSgenome.Hsapiens.UCSC.hg19,
                                  format = "bed")
data(encodemotif)

results <- motifbreakR(snpList = snps.mb.frombed[1:5], filterp = TRUE,
                       pwmList = encodemotif,
                       threshold = 0.85,
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())


calculatePvalue(results)

ranges_df<-as.data.frame(results@ranges)
motif_df<-as.data.frame(results@elementMetadata)

res_tot<-cbind(ranges_df,motif_df)


write.table(res_tot,file="motif_breaker_results.txt",row.names=F,quote=F,sep="\t")



plotMB(results = ALL_RESULTS, rsid = "chr1_4409259_G_A_DO229393_SBS1", effect = "strong")
ALL_RESULTS$rsid <- ALL_RESULTS$names
ALL_RESULTS$rsid <- gsub(":","_", ALL_RESULTS$rsid)

ALL_RESULTS$rsid=="chr1_4409259_G_A_DO229393_SBS1" <- sub("chr1_4409259_G_A_DO229393_SBS1","x", ALL_RESULTS$rsid=="chr1_4409259_G_A_DO229393_SBS1")
plotMB(ALL_RESULTS, "chr1:4409259:G:A:DO229393:SBS1", reverseMotif = TRUE, effect = c("strong", "weak"))

plotMB(results = results, rsid = "rs1006140", effect = "strong")
