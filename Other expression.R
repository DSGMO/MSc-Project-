# For genes that are only in sequencing data not array 

Sig_pro_overlap$Gene_name <- gsub("-","_", Sig_pro_overlap$Gene_name)
Sig_pro_overlap$Gene_name <- gsub(".","_", Sig_pro_overlap$Gene_name)

motif_df1 <- motif_df
motif_df1$Gene_name <- gsub("-","_", motif_df1$Gene_name)
setnames(s_mat, old=c("IGHD3-10"), new = c("IGHD2_8"))
names(s_mat) <- gsub("\\.", "", names(s_mat))

colnames(s_mat) <- gsub(".","_", colnames(s_mat))
head(s_mat)

IGHD2_8  <- s_mat$IGHD2_8 
IGHD2_8 <- as.data.frame(IGHD2_8)
IGHD2_8$icgc_donor_id <- paste(s_mat$samples)
IGHD2_8$icgc_donor_id <- sub(" PRAD_CA-WTS -","",IGHD2_8$icgc_donor_id)

IGHD2_8muts <- motif_df1[motif_df1$Gene_name == "IGHD2_8",]
IGHD2_8muts <-  IGHD2_8muts[,c("icgc_donor_id")]
IGHD2_8$mutation <-   IGHD2_8$icgc_donor_id%in%IGHD2_8muts 
IGHD2_8$mutation <- sub("TRUE","Mutated",IGHD2_8$mutation)
IGHD2_8$mutation <- sub("FALSE","Not Mutated", IGHD2_8$mutation)

setnames(IGHD2_8, old = c("IGHD2_8"), new = c("Val"))
ggplot(IGHD2_8, aes(y=Val, x=mutation, fill=mutation)) + 
  geom_jitter() + 
  geom_boxplot(alpha = 0.7) + 
  ggtitle("IGHD2_8") + 
  ylab("Expression Value") + 
  xlab("Mutation present/not present") + 
  stat_compare_means(method = "t.test") + 
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 2, type = "discrete")) + theme_bw()


genes <- Sig_pro_overlap %>% group_by(Gene_name) %>% summarise(Freq=n())
