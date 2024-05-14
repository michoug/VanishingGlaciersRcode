library(tidyverse)

## From EcolUtils package (https://github.com/GuillemSalazar/EcolUtils)
spec.gen_log<-function(comm.tab,niche.width.method="levins",perm.method="quasiswap",n=1000,probs=c(0.025,0.975)){
  require(spaa)
  require(vegan)
  occurrence<-function(x){apply(ceiling(x/max(x)),2,sum)}
  n<-n
  if (niche.width.method=="occurrence") levin.index.real<-occurrence(comm.tab) else levin.index.real<-as.numeric(niche.width(comm.tab,method=niche.width.method))
  names(levin.index.real)<-colnames(comm.tab)
  
  levin.index.simul<-matrix(NA,ncol=dim(comm.tab)[2],nrow=n)
  for (i in 1:n){
    print(i)
    if (niche.width.method=="occurrence") levin.index.simul[i,]<-occurrence(permatswap(comm.tab,perm.method,times=1)$perm[[1]]) else levin.index.simul[i,]<-as.numeric(niche.width(permatswap(comm.tab,perm.method,times=1)$perm,method=niche.width.method))
  }
  colnames(levin.index.simul)<-colnames(comm.tab)
  levin.index.simul<-as.data.frame(levin.index.simul)
  media<-apply(levin.index.simul,2,mean)
  ci<-apply(levin.index.simul,2,quantile,probs=probs)
  resultats<-data.frame(observed=levin.index.real,mean.simulated=media,lowCI=ci[1,],uppCI=ci[2,],sign=NA)
  for (j in 1:dim(resultats)[1]){
    if (resultats$observed[j]>resultats$uppCI[j]) resultats$sign[j]<-"GENERALIST"
    if (resultats$observed[j]<resultats$lowCI[j]) resultats$sign[j]<-"SPECIALIST"
    if (resultats$observed[j]>=resultats$lowCI[j] & resultats$observed[j]<=resultats$uppCI[j]) resultats$sign[j]<-"NON SIGNIFICANT"
  }
  resultats$sign<-as.factor(resultats$sign)
  resultats}


mags_tab <- read.table("data/pMAGs_cov_norm.txt.gz",sep="\t",row.names=1,header=TRUE,comment.char="@")
meta <- read_tsv("data/metadata_NOMIS_sed_chem.txt")

colnames(mags_tab) = gsub('DownB', 'DN', colnames(mags_tab))
colnames(mags_tab) = gsub('UpB', 'UP', colnames(mags_tab))
colnames(mags_tab) = gsub('Down', 'DN', colnames(mags_tab))
colnames(mags_tab) = gsub('Up', 'UP', colnames(mags_tab))
colnames(mags_tab) = gsub('(GL140)\\d$', '\\1_UP', colnames(mags_tab))
colnames(mags_tab) <- gsub("_\\d$", "", colnames(mags_tab))
colnames(mags_tab) <- gsub("X\\d+R$", "GLR140_UP", colnames(mags_tab))


mags_tab[mags_tab < 100] = 0

# write.table(mags_tab, "../MAGs_cov_norm_100.txt", row.names = TRUE, quote = F, sep = "\t")

spec_tab = mags_tab

colnames(spec_tab) = as.character(vapply(colnames(spec_tab), function(x) strsplit(x, '_')[[1]][1], FUN.VALUE = character(1)))

spec_tab = as.data.frame(do.call(cbind, by(t(spec_tab),INDICES=names(spec_tab),FUN=colMeans)))

spec_tab = spec_tab %>% reframe(across(where(is.numeric), as.integer))

rownames(spec_tab) = rownames(mags_tab)
spec_tab = spec_tab[rowSums(spec_tab) > 0,]

Sys.time()
spec_gen = spec.gen_log(t(spec_tab), niche.width.method = 'levins', perm.method = 'quasiswap', n=1000)
Sys.time()

rownames(spec_gen) = rownames(spec_tab)

write_tsv(spec_gen, "spec_gen.tsv")

