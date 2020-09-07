#Graph plots 
install.packages("ggplot2")
library(ggplot2)

Conseq <- ggplot(Consequences, aes(x="", y=Frequency, fill=Consequence))+     geom_bar(width = 1, stat = "identity")
Conseq + coord_polar("y", start=0) + ylab(NULL) + xlab(NULL) + ggtitle("Mutation Type Consequence")

Chrom <- res_ovelap_no_syn$seqnames
table(res_ovelap_no_syn$seqnames)

Chrom <- ggplot(Chroms, aes(x="", y=Frequency, fill=Chromosome)) + 
  geom_bar(width = 1, stat = "identity")
Chrom + coord_polar("y", start=0)
library(wesanderson)
Substitutions <- ggplot(AlleleAlterations, aes(x="", y=Frequency, fill=Alteration))+     geom_bar(width = 1, stat = "identity")

Substitutions + coord_polar("y", start=0) + ylab(NULL) + xlab(NULL) + ggtitle("Substitutions") + scale_fill_manual(values = wes_palette(13, name = "Darjeeling1", type = "continuous")) + theme(axis.text.x=element_blank()) + geom_text(aes(label = paste0(AlleleAlterations$Perc, "%")))


