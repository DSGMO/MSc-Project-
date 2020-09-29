#Alteration Types code 
#Preparing frequency of alterations tables for promoters and all samples to effectively compare
all_SIGS$alteration <- paste(all_SIGS$ref, all_SIGS$alt, sep = ">")
all_alteration_res <- all_SIGS %>% group_by(alteration) %>% summarise(Freq=n())
all_alteration_res <- all_alteration_res[-(13),]

Sig_pro_overlap$alteration <- paste(Sig_pro_overlap$ref_sig, Sig_pro_overlap$alt_sig, sep = ">")
pro_alteration_res <- Sig_pro_overlap %>% group_by(alteration) %>% summarise(Freq=n())
pro_alteration_res <- pro_alteration_res[-(13),]


#Code for stats 

sum(all_alteration_res$Freq)
sum(pro_alteration_res$Freq)

A_C <- matrix(c(2435, 26852, 56832, 636662), nrow = 2)
A_G <- matrix(c(5448, 63619 , 53819, 599895), nrow = 2)
A_T <- matrix(c(3146, 38762, 56121 , 624752), nrow = 2)
C_A <- matrix(c(4858, 56363, 54409 , 607151), nrow = 2)
C_G <- matrix(c(2399, 25630, 56868 , 637884), nrow = 2)
C_T <- matrix(c(10951, 117779, 48316, 545735), nrow = 2)
G_A <- matrix(c(11460, 122124, 47807, 541390), nrow = 2)
G_C <- matrix(c(2462, 25598, 56805, 637916), nrow = 2)
G_T<- matrix(c(4736, 56354, 54531, 607160), nrow = 2)
T_A <- matrix(c(3269, 39211, 55998 ,624303), nrow = 2)
T_C <- matrix(c(5618, 63636, 53649, 657896), nrow = 2)
T_G <- matrix(c(2485, 27585, 56782 , 635929), nrow = 2)

alt_p_values <- NULL
alt_MATS <- list(A_C, A_G, A_T, C_A, C_G, C_T, G_A, G_C, G_T, T_A, T_C, T_G)

for(i in alt_MATS){
  print(i)
  x <- fisher.test(i)
  alt_p_values <- rbind(alt_p_values,x$p.value)
  
}

alt_p_values = as.data.frame(alt_p_values)
alt_ORs <- NULL
alt_ORs <- as.data.frame(alt_ORs)


for(i in alt_MATS){
  print(i)
  x <- fisher.test(i)
  alt_ORs <- rbind(alt_ORs,x$estimate)
}
alt_ORs <- as.data.frame(alt_ORs)

alt_conf_i <- NULL
for(i in alt_MATS){
  print(i)
  x <- fisher.test(i)
  alt_conf_i <- rbind(alt_conf_i,x$conf.int)
  
}
alt_conf_i <- as.data.frame(alt_conf_i)

alt_stats <- NULL
setnames(alt_conf_i, old = c("V1","V2"), new=c("lower","upper"))
setnames(alt_ORs, old = c("X1.01587207453593"), new = c("OR"))
setnames(alt_p_values, old=c("V1"), new=c("P_value"))

alt_stats <- cbind(alt_conf_i, alt_ORs, alt_p_values)
alt_stats$alt <- c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G")

#Forest Plot for Alteration Types 
alt_fp <- ggplot(data=alt_stats, aes(x=reorder(alt, OR), y=OR, ymin=lower, ymax=upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Alteration") + ylab("Odds Ratio (95% CI)") +
  theme_bw()  # use a white background

alt_fp

alt_stats <- alt_stats %>% mutate_if(is.numeric, ~round(., 4))
for(i in 1:nrow(alt_stats)){
  if(sig_stats$p_value[i] >0.05){
    alt_stats$label[i] <- paste("P > 0.05")
  }
  else if(alt_stats$p_value[i]<0.05){
    alt_stats$label[i] <- paste("P < 0.05")
  }
}

alt_p_values <- alt_p_values %>% mutate_if(is.numeric, ~round(., 4))

for(i in 1:nrow(alt_stats)){
  alt_stats$conf[i] <- paste(alt_stats$lower[i], alt_stats$upper[i], sep = "-")
}

#Chart for comparing all sigs and promoter sigs 

all_SIGS = as.data.table(all_SIGS)
setkey(all_SIGS, sig)
all_SIGS <- all_SIGS[!"SBS43"]
all_SIGS <- all_SIGS[!"SBS46"]
all_SIGS <- all_SIGS[!"SBS47"]
all_SIGS <- all_SIGS[!"SBS50"]
all_SIGS <- all_SIGS[!"SBS51"]
all_SIGS <- all_SIGS[!"SBS53"]
all_SIGS <- all_SIGS[!"SBS54"]
all_SIGS <- all_SIGS[!"SBS55"]
all_SIGS <- all_SIGS[!"SBS57"]
all_SIGS <- all_SIGS[!"SBS58"]
all_SIGS <- all_SIGS[!"SBS60"]

