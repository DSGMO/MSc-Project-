Sig_pro_overlap = as.data.table(Sig_pro_overlap)
setkey(Sig_pro_overlap, sig)
Sig_pro_overlap <- Sig_pro_overlap[!"SBS43"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS27"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS45"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS48"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS49"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS52"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS56"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS59"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS46"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS47"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS50"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS51"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS53"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS54"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS55"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS57"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS58"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS60"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS7a"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS7b"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS7c"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS7d"]
Sig_pro_overlap <- Sig_pro_overlap[!"SBS38"]

all_SIGS = as.data.frame(all_SIGS)

for(i in 1:nrow(all_SIGS)){
  print(i)
  if(all_SIGS$sig[i]=="SBS3"||all_SIGS$sig[i]== "SBS15"||
     all_SIGS$sig[i]== "SBS21"||all_SIGS$sig[i]== "SBS26"||
     all_SIGS$sig[i]== "SBS30"||all_SIGS$sig[i]== "SBS36"||
     all_SIGS$sig[i]== "SBS44"||all_SIGS$sig[i]== "SBS6"||
     all_SIGS$sig[i]== "SBS14"){
    all_SIGS$Signature_type[i] <- paste("DNA Repair Related")
  }
  else if(all_SIGS$sig[i] == "SBS1"|| all_SIGS$sig[i]== "SBS5"||
          all_SIGS$sig[i]=="SBS40"||all_SIGS$sig[i]=="SBS20"){
    all_SIGS$Signature_type[i] <- paste("Ageing Related")
  }
  else if(all_SIGS$sig[i] == "SBS2"||all_SIGS$sig[i]== "SBS13"){
    all_SIGS$Signature_type[i] <- paste("APOBEC Related")
  }
  else if(all_SIGS$sig[i] == "SBS7a"||all_SIGS$sig[i]=="SBS7b"||
          all_SIGS$sig[i]=="SBS7c"||all_SIGS$sig[i]=="SBS7d"||
          all_SIGS$sig[i]=="SBS38"){
    all_SIGS$Signature_type[i] <- paste("UV light Related")
  }
  else if(all_SIGS$sig[i] == "SBS11"||all_SIGS$sig[i]=="SBS25"||
          all_SIGS$sig[i]=="SBS31"||all_SIGS$sig[i]=="SBS32"||
          all_SIGS$sig[i]=="SBS35"){
    all_SIGS$Signature_type[i] <- paste("Treatment Related")
  }
  else if(all_SIGS$sig[i] == "SBS42"||all_SIGS$sig[i]=="SBS24"||
          all_SIGS$sig[i]=="SBS22"){
    all_SIGS$Signature_type[i] <- paste("Endogenous Exposure Related")
  }
  else if(all_SIGS$sig[i] == "SBS4"||all_SIGS$sig[i]=="SBS29"){
    all_SIGS$Signature_type[i] <- paste("Tobacco Related")
  }
  else if(all_SIGS$sig[i] == "SBS8"||all_SIGS$sig[i]=="SBS53"|all_SIGS$sig[i]=="SBS47"||
          all_SIGS$sig[i]=="SBS46"||all_SIGS$sig[i]=="SBS9"||all_SIGS$sig[i]=="SBS57"||
          all_SIGS$sig[i]=="SBS58"||all_SIGS$sig[i]=="SBS60"||
          all_SIGS$sig[i]=="SBS54"||all_SIGS$sig[i]=="SBS10a"||
          all_SIGS$sig[i]=="SBS10b"||all_SIGS$sig[i]=="SBS12"|
          all_SIGS$sig[i]=="SBS16"||all_SIGS$sig[i]=="SBS17a"||
          all_SIGS$sig[i]=="SBS17b"||all_SIGS$sig[i]=="SBS18"||
          all_SIGS$sig[i]=="SBS19"||all_SIGS$sig[i]=="SBS23"||
          all_SIGS$sig[i]=="SBS28"||all_SIGS$sig[i]=="SBS41"||all_SIGS$sig[i]=="SBS33"||
          all_SIGS$sig[i]=="SBS34"||all_SIGS$sig[i]=="SBS37"||all_SIGS$sig[i]=="SBS39"){
    all_SIGS$Signature_type[i] <- paste("Other/Unknown")
  }
  else if(all_SIGS$sig[i] =="SBS85"||all_SIGS$sig[i] =="SBS84"){
    all_SIGS$Signature_type[i] <- paste("Cytidine Deaminase Related")
  }
  else if(all_SIGS$sig[i] == "SBS43"||all_SIGS$sig[i] == "SBS46"||all_SIGS$sig[i] == "SBS47"
          ||all_SIGS$sig[i] == "SBS50"||all_SIGS$sig[i] == "SBS51"||all_SIGS$sig[i] == "SBS53"
          ||all_SIGS$sig[i] == "SBS54"||all_SIGS$sig[i] == "SBS55"||all_SIGS$sig[i] == "SBS57"
          ||all_SIGS$sig[i] == "SBS58"||all_SIGS$sig[i] == "SBS60"||all_SIGS$sig[i] == "SBS27"
          ||all_SIGS$sig[i]=="SBS45"||all_SIGS$sig[i]=="SBS48"||all_SIGS$sig[i]=="SBS49"
          ||all_SIGS$sig[i]=="SBS52"||all_SIGS$sig[i]=="SBS56"||all_SIGS$sig[i]=="SBS59"){
    all_SIGS$Signature_type[i] <- paste("Possible Sequencing Artefacts")
  }
}


save(all_SIGS, file = "all_SIGS_annotated.RDATA")
geneList = as.data.frame(rel_genes$Gene_name)

save(geneList, file = "geneList.txt", sep = "", quote = F, header = F, row.names = F)
write.table(geneList, "GeneList.txt", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

