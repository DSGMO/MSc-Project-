library(deconstructSigs)

setwd()
sigs.input <- mut.to.sigs.input(mut.ref = simple_somatic_mutation.open.PRAD.CA, 
                                sample.id = "icgc_donor_id", 
                                chr = "chromosome", 
                                pos = "chromosome_start", 
                                ref = "mutated_from_allele", 
                                alt = "mutated_to_allele")

samples<-rownames(sigs.input)
head(signatures.cosmic)
head(sigProfiler_SBS_signatures_2019_05_22)
allsignatures <- NULL
for(i in samples[1:328]){
  print(i)
  
  x <-whichSignatures(tumor.ref = sigs.input, signatures.ref = cosmic_table3, i, contexts.needed = TRUE)
  
  #x=as.data.frame(x)
  allsignatures <- rbind(allsignatures,x$weights)

 
}


save(allsignatures, file="allsignatures.RData")


dt2 <- allsignatures%>%
  +     rownames_to_column(allsignatures) %>%
  +     gather(colname, value, -rowname) 

head(dt2)

ggplot(dt2, aes(x = colname, y = rowname, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "black") + ggtitle("Signature Weights Across All Samples") + ylab("Sample") + xlab("Signature")



library(deconstructSigs)
DO229395 <- whichSignatures(tumor.ref = sigs.input, "DO229395", signatures.ref = cosmic_table3, 
                            contexts.needed = TRUE, tri.counts.method = "default")
plotSignatures(DO229395)
makePie(DO229395)

head(cosmic_table2)


setwd("/Users/katefodder/Desktop/MSC_PROJECT")

write.table(allsignatures, file="all.signatures.txt")


sappl(ex, function(x) )

write.table(allsignatures, file="....txt", quote=FALSE, sep="\t")

pdf('output.pdf')


dev.off()
sapply(1:90, function(x) whichSignatures(tumor.ref = sigs.input,
                                         signatures.ref = signatures.cosmic ,
                                         sample.id = x,
                                         contexts.needed = TRUE,
                                         tri.counts.method = 'default'))


table <- cbind(t2[!sapply(t2, is.list)], 
               +                (t(apply(t2[sapply(t2, is.list)], 1, unlist))))