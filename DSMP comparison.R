

ds_sigs <- rownames_to_column(allsignatures, var = "donor")
ds_sigs2 <- melt(ds_sigs, id.vars = c("donor"), measure.vars = signatures )

allsigs2 <- as.data.frame(t(allsignatures))



signatures <- as.vector(colnames(allsignatures))

mp_sigs <- rownames_to_column(mp_sigs, var = "donor")
mp_sigs2 <- melt(mp_sigs, id.vars=c("donor"), measure.vars = signatures)

setnames(mp_sigs2, old=c("value"), new=c("mutpat_val"))
setnames(ds_sigs2, old=c("value"), new=c("desigs_val"))

comb_sigs <- cbind(ds_sigs2, mp_sigs2)
comb_sigs <- comb_sigs[,-4]
comb_sigs <- comb_sigs[,-4]

pc_genes <- Sig_pro_overlap[Sig_pro_overlap$Gene_type == "protein_coding",]
sig_1 <- comb_sigs[comb_sigs$variable =="Signature.1",]
ggplot(sig_1, aes(y=mutpat_val, x=desigs_val, fill=variable)) + geom_jitter()
