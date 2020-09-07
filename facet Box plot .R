#Facet plot with boxplots of altnatures comparing all samples and Promoter samples 

library(ggrepel)
library("ggsci")
library("ggplot2")
library("gridExtra")
library(ggplot2)
library(wesanderson)
library(tidyverse)
library(gapminder)
library(forcats)
library(dplyr)
library(data.table)
library(tidyverse)
library(ggpubr)
library(rstatix)

alt_pro_res <- alt_pro_overlap %>% group_by(alt, altnature_type) %>% summarise(Freq=n())
alt_all_res <- all_altS %>% group_by(alt, altnature_type) %>% summarise(Freq=n())




sum(alt_all_res$Freq) #663518
sum(alt_pro_res$Freq) #59267


alt_pro_res$label <- paste("Promoter Mutations")
alt_all_res$label <- paste("All Mutations")

alt_pro_res$adj_f <- paste((alt_pro_res$Freq)/sum(alt_pro_res$Freq))
alt_all_res$adj_f <- paste((alt_all_res$Freq)/sum(alt_all_res$Freq))

pro_all_res <- rbind(alt_pro_res, alt_all_res)

ggplot(pro_all_res, aes(x=label, y=adj_f, fill=altnature_type)) + 
  geom_col() + 
  scale_fill_npg() + 
  facet_wrap(~alt) +
  theme(axis.text.x=element_text(angle=90, hjust=1), axis.text.y = element_blank()) + 
  ggtitle("Mutations associated with altnatures/Sample") + 
  ylab("Mutations") + 
  xlab("altnature")

p_values$p_val <- p_values$V1
p_values$V1 <- NULL
p_values$group1 <- paste("Promoter Mutations")
p_values$group2 <- paste("All Mutations")

p_values$y <- paste("adj_f")
stat.test <- p_values %>%
  group_by() %>%
  t_test(len ~ supp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_altnificance()

stat.test <- p_values

#Adding ORs and P-values 


alts <- as.character(unique(alt_all_res$alt))
SBS1 <- matrix(c(11396, 115926, 47871, 547592), nrow = 2)
SBS10b <- matrix(c(109, 1332, 59158, 662186),nrow = 2)
SBS11 <- matrix(c(138, 1370, 59129, 662148),nrow = 2)
SBS15 <- matrix(c(1198, 12707, 58069, 650811),nrow = 2)
SBS16 <- matrix(c(34, 337, 59233, 663181),nrow = 2)
SBS17b <- matrix(c(5, 97, 59262, 663421),nrow = 2)
SBS19 <- matrix(c(1579, 18866, 57688, 644652),nrow = 2)
SBS2 <- matrix(c(1496, 17479, 57771, 646039),nrow = 2)
SSBS20 <- matrix(c(0, 16, 59267, 663502),nrow = 2)
SBS22 <- matrix(c(86, 1025, 59181, 662493),nrow = 2)
SBS23 <- matrix(c(129, 1470, 59138, 662048),nrow = 2)
SBS25 <- matrix(c(51, 574, 59216, 662944),nrow = 2)
SBS26 <- matrix(c(0, 5, 59267, 663513),nrow = 2)
SBS28 <- matrix(c(300, 3692, 58967, 659826),nrow = 2)
SBS3 <- matrix(c(207, 2091, 59060, 661427),nrow = 2)
SBS30 <- matrix(c(2091, 22173, 57176, 641345),nrow = 2)
SBS31 <- matrix(c(294, 2997, 58973, 660521),nrow = 2)
SBS32 <- matrix(c(2141, 25678, 57126, 637840),nrow = 2)
SBS33 <- matrix(c(21, 210, 59246, 663308),nrow = 2)
SBS35 <- matrix(c(4, 81, 59263, 663437),nrow = 2)
SBS36 <- matrix(c(5, 33, 59262, 663485),nrow = 2)
SBS37 <- matrix(c(4682, 55861, 54585, 607657),nrow = 2)
SBS39 <- matrix(c(61, 859, 59206, 662659),nrow = 2)
SBS40 <- matrix(c(19855, 227545, 39412, 435973),nrow = 2)
SBS41 <- matrix(c(996, 12389, 58271, 651129),nrow = 2)
SBS42 <- matrix(c(197, 2129, 59070, 661389),nrow = 2)
SBS44 <- matrix(c(196, 2213, 59071, 661305),nrow = 2)
SBS5 <- matrix(c(9789, 110199, 49478, 553319),nrow = 2)
SBS6 <- matrix(c(808, 8166, 58459, 655352), nrow=2)
SBS8 <- matrix(c(22, 277, 59245, 663241), nrow=2)
SBS84 <- matrix(c(356, 3766, 58911, 659752), nrow=2)
SBS85 <- matrix(c(130, 1514, 59137, 662004), nrow=2)
SBS9 <- matrix(c(891, 10402, 58376, 653116), nrow=2)

p_values <- NULL
MATS <- list(SBS1, SBS10b, SBS11, SBS15, SBS16, SBS17b, 
             SBS19, SBS2, SBS22, SBS23, SBS25, SBS28, SBS3, 
             SBS30, SBS31, SBS32, SBS33, SBS35, SBS36, 
             SBS37, SBS39, SBS40, SBS41, SBS42,
             SBS44, SBS5, SBS6, SBS84, SBS85)

for(i in MATS){
  print(i)
  x <- fisher.test(i)
  p_values <- rbind(p_values,x$p.value)

}
ORs <- NULL
ORs <- as.data.frame(ORs)

for(i in MATS){
  print(i)
  x <- fisher.test(i)
  ORs <- rbind(ORs,x$estimate)
  
}



conf_i <- NULL
for(i in MATS){
  print(i)
  x <- fisher.test(i)
  conf_i <- rbind(conf_i,x$conf.int)
  
}

test <- fisher.test(SBS1)

p_values = as.data.frame(p_values)
ORs = as.data.frame(ORs)
conf_i = as.data.frame(conf_i)

p_values$alt <- c("SBS1", "SBS10b", "SBS11", "SBS15", "SBS16", "SBS17b", 
                "SBS19", "SBS2", "SBS22", "SBS23", "SBS25","SBS28", "SBS3", 
                "SBS30", "SBS31", "SBS32","SBS33", "SBS35", "SBS36", 
                "SBS37", "SBS39", "SBS40", "SBS41", "SBS42",
                "SBS44", "SBS5", "SBS6", "SBS84", "SBS85")
ORs$alt <- c("SBS1", "SBS10b", "SBS11", "SBS15", "SBS16", "SBS17b", 
                  "SBS19", "SBS2", "SBS22", "SBS23", "SBS25","SBS28", "SBS3", 
                  "SBS30", "SBS31", "SBS32","SBS33", "SBS35", "SBS36", 
                  "SBS37", "SBS39", "SBS40", "SBS41", "SBS42",
                  "SBS44", "SBS5", "SBS6", "SBS84", "SBS85")

conf_i$alt <- c("SBS1", "SBS10b", "SBS11", "SBS15", "SBS16", "SBS17b", 
                "SBS19", "SBS2", "SBS22", "SBS23", "SBS25","SBS28", "SBS3", 
                "SBS30", "SBS31", "SBS32","SBS33", "SBS35", "SBS36", 
                "SBS37", "SBS39", "SBS40", "SBS41", "SBS42",
                "SBS44", "SBS5", "SBS6", "SBS84", "SBS85")

p_values <- p_values %>% mutate_if(is.numeric, ~round(., 4))
setnames(ORs, old = c("X1.12451616845076"), new = c("OR"))


alt_stats <- cbind(ORs, p_values)
setnames(conf_i, old  = c("V1","V2"), new=c("lower","upper"))

alt_stats <- cbind(alt_stats, conf_i)
setnames(alt_stats, old=c("V1"), new = c("p_value"))
alt_stats <- alt_stats[,-2]

alt_stats <- alt_stats %>% mutate_if(is.numeric, ~round(., 4))
for(i in 1:nrow(alt_stats)){
  alt_stats$conf[i] <- paste(alt_stats$lower[i], alt_stats$upper[i], sep = "-")
}
#Forest plots for stats

ggplot(data=alt_stats, aes(x=reorder(alt, OR), y=OR, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(colour=label)) + 
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("altnature") + ylab("Odds Ratio (95% CI)") +
  theme_bw()  # use a white background

for(i in 1:nrow(alt_stats)){
  if(alt_stats$p_value[i]>0.05){
    alt_stats$label[i] <- paste("P > 0.05")
  }
  else if(alt_stats$p_value[i]<0.05){
    alt_stats$label[i] <- paste("P < 0.05")
  }
}




