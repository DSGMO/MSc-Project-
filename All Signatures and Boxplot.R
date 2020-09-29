#Code for barplot with signatures and prevalence 

for(i in 1:nrow(Sig_pro_overlap)){
  if(Sig_pro_overlap$sig[i]=="SBS3"||Sig_pro_overlap$sig[i]== "SBS30"||Sig_pro_overlap$sig[i]== "SBS44"||Sig_pro_overlap$sig[i]== "SBS6"){
    Sig_pro_overlap$Signature_type[i] <- paste("DNA Repair Related")
  }
  else if(Sig_pro_overlap$sig[i] == "SBS1"|| Sig_pro_overlap$sig[i]== "SBS5"||Sig_pro_overlap$sig[i]=="SBS40"){
    Sig_pro_overlap$Signature_type[i] <- paste("Ageing Related")
  }
  else if(Sig_pro_overlap$sig[i] == "SBS2"){
    Sig_pro_overlap$Signature_type[i] <- paste("APOBEC Related")
  }
  else if(Sig_pro_overlap$sig[i] == "SBS7a"||Sig_pro_overlap$sig[i]=="SBS7b"){
    Sig_pro_overlap$Signature_type[i] <- paste("UV light Related")
  }
  else if(Sig_pro_overlap$sig[i]=="SBS32"){
    Sig_pro_overlap$Signature_type[i] <- paste("Treatment Related")
  }
  else if(Sig_pro_overlap$sig[i] == "SBS42"){
    Sig_pro_overlap$Signature_type[i] <- paste("Endogenous Exposure Related")
  }
  else if(Sig_pro_overlap$sig[i] == "SBS37"||Sig_pro_overlap$sig[i]=="SBS9"||Sig_pro_overlap$sig[i]=="SBS16"||Sig_pro_overlap$sig[i]=="SBS19"||Sig_pro_overlap$sig[i]=="SBS23"||Sig_pro_overlap$sig[i]=="SBS28"||Sig_pro_overlap$sig[i]=="SBS41"){
    Sig_pro_overlap$Signature_type[i] <- paste("Other/Unknown")
  }
  else if(Sig_pro_overlap$sig[i]=="SBS84"||Sig_pro_overlap$sig[i]=="SBS85")
    Sig_pro_overlap$Signature_type[i] <- paste("Cytidine deaminase Related")
  }
}


Sigs_res <- Sig_pro_overlap %>% group_by(sig, Signature_type, icgc_donor_id) %>%  summarise(Freq=n())

ggplot(Sigs_res, aes(x=reorder(sig, Freq), y=Freq, fill=Signature_type)) + 
  geom_boxplot() + 
  scale_fill_npg() + 
  xlab("Mutational Signatures") + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  ggtitle("All promoter mutations associated with Signatures") +
  coord_trans(y = "log1p") + 
  labs(fill="Signature Type") +
  ylab("Mutations/Sample (Log1p)") + 
  xlab("Mutational Signature")





