library(lme4)
library(dplyr)
library(ggplot2)
library(evir)
library(miqtl)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("At least two argument must be supplied", call.=FALSE)
} else if (length(args)==2) {
  FROM=args[1]
  END=args[2]
}

Expdf.full <- read.table(file="../data/CCexpsum.full.adj.txt",sep="\t",header=TRUE)

CClist <- sort(unique(intersect(Expdf.full$CC[Expdf.full$Study=='TAK'],Expdf.full$CC[Expdf.full$Study=='TOV'])))

Expdf <- Expdf.full[Expdf.full$CC %in% CClist,]

info <- Expdf[,c(1:7)]

info$Study <- as.character(info$Study)
info$Study2 <- info$Study
info$Study[which(info$Study=='TAK')] <- 0.5
info$Study[which(info$Study=='TOV')] <- -0.5
info$Study <- as.numeric(info$Study)

probelist <- colnames(Expdf)[8:dim(Expdf)[2]]

anno4 <- read.table(file="../data/TOV_anno4.txt",sep="\t",header=TRUE)

all(probelist %in% paste0('X',anno4$Probe.Set.ID))
anno4$ProbeX <- paste0('X',anno4$Probe.Set.ID)

# make a table of LR and LRlist
thisprobelist <- probelist
outpd <- data.frame()

info -> info.ori
for (pb in thisprobelist[FROM:END]){
  info <- info.ori
  
  message(pb)
  info$y <- as.numeric(Expdf[,pb])
  probechr <- anno4$Achr[which(anno4$ProbeX==pb)] 
  probep1 <- anno4$Apos[which(anno4$ProbeX==pb)]
  probep2 <- anno4$Apos2[which(anno4$ProbeX==pb)]

  info$y.adj[which(info$Study==0.5)] <- info$y[which(info$Study==0.5)] - mean(info$y[which(info$Study==0.5)])
  info$y.adj[which(info$Study== -0.5)] <- info$y[which(info$Study== -0.5)] - mean(info$y[which(info$Study== -0.5)])
  # mega peak height
  info <- info %>% group_by(CC,Study) %>% summarize(weight=length(y),y=mean(y),y.adj=mean(y.adj))
  info$An <- paste0(info$CC,info$Study) 

  myweight <- info$weight;names(myweight) <- info$An;  

  TAKscan <- scan.h2lmm(genomecache="~/IDSS/segments_happy_format_mm10/",return.allele.effects = TRUE,weights=myweight,
                        data=info[which(info$Study==0.5),], formula= rint(y)~ 1,geno.id="CC",pheno.id = 'An',use.multi.impute=FALSE,print.locus.fit=FALSE)
  TOVscan <- scan.h2lmm(genomecache="~/IDSS/segments_happy_format_mm10/",return.allele.effects = TRUE,weights=myweight,
                        data=info[which(info$Study==-0.5),], formula= rint(y)~ 1,geno.id="CC",pheno.id = 'An',use.multi.impute=FALSE,print.locus.fit=FALSE)
  Studyscan <- scan.h2lmm(genomecache="~/IDSS/segments_happy_format_mm10/",return.allele.effects = TRUE,weights=myweight,
                          data=info, formula= rint(y)~ 1 + Study,geno.id="CC",pheno.id = 'An',use.multi.impute=FALSE,print.locus.fit=FALSE)
  Studyscan2 <- scan.h2lmm(genomecache="~/IDSS/segments_happy_format_mm10/",return.allele.effects = TRUE,weights=myweight,
                          data=info, formula= rint(y.adj)~ 1 + Study,geno.id="CC",pheno.id = 'An',use.multi.impute=FALSE,print.locus.fit=FALSE)
# TAKTOVMega
  TAKout <- as.data.frame( t( round(-log10(TAKscan$p.value),4) ) )
  write.table(TAKout,file=paste0("../result/eqtlcistrans_allp/TAKscan_",FROM,"_",END,".txt"),sep="\t",col.names=F,row.names=F,quote=F,append=T)
  
  TOVout <- as.data.frame( t( round(-log10(TOVscan$p.value),4) ) )
  write.table(TOVout,file=paste0("../result/eqtlcistrans_allp/TOVscan_",FROM,"_",END,".txt"),sep="\t",col.names=F,row.names=F,quote=F,append=T)

  Megaout <- as.data.frame( t( round(-log10(Studyscan$p.value),4) ) )
  write.table(Megaout,file=paste0("../result/eqtlcistrans_allp/Megascan_",FROM,"_",END,".txt"),sep="\t",col.names=F,row.names=F,quote=F,append=T)

  Megaout2 <- as.data.frame( t( round(-log10(Studyscan2$p.value),4) ) )
  write.table(Megaout2,file=paste0("../result/eqtlcistrans_allp/Megascanadj_",FROM,"_",END,".txt"),sep="\t",col.names=F,row.names=F,quote=F,append=T)
}

 
