#Merging promoter data with signature data 
setwd('/Users/katefodder/Desktop/MSC_PROJECT')
merged_sigs_promoters2 <- merge(res_overlap, all_sig_mut_files, by.x=c("icgc_donor_id","position","mutated_from_allele","mutated_to_allele"),by.y=c("icgc_donor_id","position","ref","alt"),all.x=FALSE, all.y=FALSE)

save(merged_sigs_promoters, file="merged_sigs_promoters.RData")
