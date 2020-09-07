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

# Creating new data frame so only number of times gene mutated rather than total number of mutations in a gene 
unique_genes <- unique(res_overlap[,c("icgc_donor_id","Gene_name")])

freq_genes_bysample <-data.frame(table(unique_genes$Gene_name))
freq_genes_bysample2 <- freq_genes_bysample[order(freq_genes_bysample[,2], decreasing = T),]

head(freq_genes_bysample)

head(freq_genes)






dim(res_overlap)
dim(res_overlap_no_syn)

freq_genes<-data.frame(table(res_overlap_no_syn$Gene_name))
head(freq_genes)

freq_genes2<-freq_genes[order(freq_genes[,2],decreasing=T),]

head(freq_genes2)



table(freq_genes2)
table(res_overlap_no_syn$Genename)

table(res_overlap_no_syn$consequence_type)
nrow(res_overlap_no_syn)
Consequences <- res_overlap_no_syn$consequence_type
head(res_overlap_no_syn$consequence_type)
table(Consequences)

Consequences = as.data.frame(Consequences)

head(Consequences)

table(Consequences)

save(res_overlap_no_syn,file="results_annotation_promoters_with_mutations_PRAD.RData")
table(res_overlap_no_syn$consequence_type)
table(res_overlap_no_syn$consequence_type)
test<-res_overlap_no_syn[res_ovelap_no_syn$consequence_type="intron_variant",]

test<-res_ovelap_no_syn[which(res_ovelap_no_syn$consequence_type%in%"intron_variant"),]

head(test[1:10,c(1:10)])


