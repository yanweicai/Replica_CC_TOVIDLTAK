library(lme4)
library(lmerTest)
library(qvalue)
library(varComp)
library(miqtl)
library(ggplot2)
library(tidyverse)

setwd("~/Dropbox/ValdarLab/IDSScross_git/src/")
# read in data table
info <- read.table(file="../data/info_3study.txt",header=TRUE,sep = "\t")
info$CCbyStudy <- paste0(info$CC,'_',info$Study)
options(contrasts = rep("contr.sum", 2)) # Set contrasts -- see below

### Section 1. Variance components of phenotype models.
phone_list=c("D.BW","BWratio","LWratio","AST","ALT","TBIL")

outdf <- data.frame() # output data frame
for (ph in phone_list){
  message(ph)
  info$y <- as.numeric(info[,ph])

  # Adjust for Batch effects 
  fit0 <- lmer( y ~ 1 + Study + (1|CCline) + (1|CCbyStudy) + (1|Dosing.Date), data=info,REML=FALSE)
  batef <- ranef(fit0)$Dosing.Date;batch2num <- batef[[1]]; names(batch2num) <- rownames(batef)
  info$y <- info$y - batch2num[info$Dosing.Date]
  
  ### single study linear mix model by lmer 
  fit1 <- lmer( y ~ 1 +(1|CCline), data=(info[which(info$Study=='TOV'),]),REML=FALSE)
  varout1 <- as.data.frame(VarCorr(fit1))
  h2.tov <- round(varout1$vcov[[1]]/sum(varout1$vcov),3)
  h2.tov.p <- ranova(fit1)$"Pr(>Chisq)"[2]
    
  fit2 <- lmer( y ~ 1 + (1|CCline), data=info[which(info$Study=='TAK'),],REML=FALSE)
  varout2 <- as.data.frame(VarCorr(fit2))
  h2.tak <- round(varout2$vcov[[1]]/sum(varout2$vcov),3)
  h2.tak.p <- ranova(fit2)$"Pr(>Chisq)"[2]

  fit3 <- lmer( y ~ 1 + (1|CCline), data=info[which(info$Study=='GLD'),],REML=FALSE)
  varout3<- as.data.frame(VarCorr(fit3))
  h2.gld <- round(varout3$vcov[[1]]/sum(varout3$vcov),3)
  h2.gld.p <- ranova(fit3)$"Pr(>Chisq)"[2]
  

  # phenotype analysis for three studies
  fit5 <- lmer( y ~ 1 + Study + (1|CCline) + (1|CCbyStudy), data=info,REML=FALSE)
  fit7 <- lmer( y ~ 1 +  (1|CCline) + (1|CCbyStudy), data=info,REML=FALSE)
    
  h2.SxS.p <- ranova(fit5)$"Pr(>Chisq)"[3]
  h2.C.p <- ranova(fit5)$"Pr(>Chisq)"[2]
  h2.S.p <- anova(fit5,fit7)$"Pr(>Chisq)"[2]
    
  varout5 <- as.data.frame(VarCorr(fit5))
  v.CC <- varout5$vcov[which(varout5$grp=="CCline")]
  v.R <- varout5$vcov[which(varout5$grp=="Residual")]
  v.CCS <- varout5$vcov[which(varout5$grp=='CCbyStudy')] # var(X)+var(Y)+2Cov(X+Y)

  v.All <- var(info$y,use='na.or.complete')
  v.Study <- v.All-v.CC-v.R-v.CCS

  h2.mega.CS <- round(v.CCS/v.All,3)
  h2.mega.CC <- round(v.CC/v.All,3)
  h2.mega.Study <- round(v.Study/v.All,3)

  p2p <- function(plist){
    plisto <- rep('',length(plist))
    plisto[which(plist>=0.1)]<-''
    plisto[which(plist<=0.1)]<-'.'
    plisto[which(plist<=0.05)]<-'*'
    plisto[which(plist<=0.01)]<-'**'
    plisto[which(plist<=0.001)]<-'***'
    plisto
  }

  thisdf <- data.frame(ph=ph,h2.tov=paste0(h2.tov,p2p(h2.tov.p)),h2.tak=paste0(h2.tak,p2p(h2.tak.p)),h2.gld=paste0(h2.gld,p2p(h2.gld.p)),
                       CC.123=paste0(h2.mega.CC,p2p(h2.C.p)),
                       Study.123=paste0(h2.mega.Study,p2p(h2.S.p)),
                       SxS.123=paste0(h2.mega.CS,p2p(h2.SxS.p)),
                       h2.tov.p=h2.tov.p,h2.tak.p=h2.tak.p,h2.gld.p=h2.gld.p,
                       h2.C.p=h2.C.p,h2.S.p=h2.S.p,h2.SxS.p=h2.SxS.p)
  
  #thisdf2 <- data.frame(ph=ph,h2.tov=h2.tov,h2.tov.p=h2.tov.p,h2.tak=h2.tak,h2.tak.p=h2.tak.p,h2.gld=h2.gld,h2.gld.p=h2.gld.p,
  #                     Study.123=mega123[1],Study.123.p=mega123[2],CC.123=mega123[3],CC.123.p=mega123[4],SxS.123=mega123[5],SxS.123.p=mega123[6])
                       
  outdf <- rbind(outdf,thisdf)
}
# outdf --> Variance components of phenotype models 


### Section 2. shrinkage estimation of variance components (Body weight)
library(limma)

BWsum <- info %>% filter(!is.na(D.BW)) %>% group_by(Study,CCline) %>% dplyr::summarize(mean = mean(D.BW),var = var(D.BW), obs=length(D.BW))

BWsum$x <- 0;BWsum$x[which(BWsum$Study=='TAK')]<-1;BWsum$x[which(BWsum$Study=='GLD')]<-2;
BWsum$var.shink <- BWsum$var
BWsum$var.shink <- squeezeVar(BWsum$var,(BWsum$obs-1))$var.post
BWsum$sd <- sqrt(BWsum$var.shink)
BWsum$cv <- BWsum$sd/BWsum$mean

BWcount <- as.data.frame.matrix(table(BWsum[,c('CCline','x')]))
strainfull <- rownames(BWcount)[which(BWcount[[1]]==1 & BWcount[[2]]==1 & BWcount[[3]]==1)]
BWsum$col <- 0
BWsum$col[which(BWsum$CCline %in% strainfull)] <- 1

BWsumsub <- BWsum[which(BWsum$x==0),]
BWsumsub$x <- 3
BWsum <- rbind(BWsum,BWsumsub)

# plot for the within-strain differences figure
#pdf(file="Result/BWsd.pdf",width=4.5,height=4.5)
p<- ggplot(BWsum,aes(x=x,y=sd,group=CCline))
for (xi in 0:3){p <- p + geom_vline(xintercept = xi, color = "gray", size=1)}
p <- p + geom_point(aes(color=CCline),show.legend = F) + 
  geom_line(aes(color=CCline),show.legend = F,subset(BWsum,col==1)) +
  geom_line(aes(),color='gray',show.legend = F,subset(BWsum,x %in% c(0,1) & col==0)) +
  geom_line(aes(),color='gray',show.legend = F,subset(BWsum,x %in% c(1,2) & col==0)) +
  geom_line(aes(),color='gray',show.legend = F,subset(BWsum,x %in% c(2,3) & col==0))
  
p+theme(text=element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_blank())
#dev.off()

### Section 3. QTL mapping for ALT/AST phenotype
for (ph in c('ALT','AST')){
  message(ph)
  info.ph <- cbind(info[,1:7],data.frame(ph=info[,ph]))
  info.ph$ph.adj <- info.ph$ph 
  
  # debatch and take strain mean within each 
  Tov_df <- info.ph[which(info.ph$Study=='TOV'),]
  batcheff <- ranef( ( lmer(ph~ (1|CCline) + (1|Dosing.Date),data=Tov_df) ))$Dosing.Date
  batcheffv <- batcheff[[1]];names(batcheffv) <- rownames(batcheff)
  Tov_df$ph.adj <- Tov_df$ph - batcheffv[as.character(Tov_df$Dosing.Date)]
  Tov_df2 <- merge( aggregate(ph.adj ~ CCline,Tov_df,length),aggregate(ph.adj ~ CCline,Tov_df,mean),by="CCline")
  colnames(Tov_df2) <- c('CC','weight','ph')
  Tov_df2$Study <- 'TOV'
  
  Tak_df <- info.ph[which(info.ph$Study=='TAK'),]
  batcheff <- ranef( ( lmer(ph~ (1|CCline) + (1|Dosing.Date),data=Tak_df) ))$Dosing.Date
  batcheffv <- batcheff[[1]];names(batcheffv) <- rownames(batcheff)
  Tak_df$ph.adj <- Tak_df$ph - batcheffv[as.character(Tak_df$Dosing.Date)]
  Tak_df2 <- merge( aggregate(ph.adj ~ CCline,Tak_df,length),aggregate(ph.adj ~ CCline,Tak_df,mean),by="CCline")
  colnames(Tak_df2) <- c('CC','weight','ph')
  Tak_df2$Study <- 'TAK'
  
  Gld_df <- info.ph[which(info.ph$Study=='GLD'),]
  batcheff <- ranef( ( lmer(ph~ (1|CCline) + (1|Dosing.Date),data=Gld_df) ))$Dosing.Date
  batcheffv <- batcheff[[1]];names(batcheffv) <- rownames(batcheff)
  Gld_df$ph.adj <- Gld_df$ph - batcheffv[as.character(Gld_df$Dosing.Date)]
  Gld_df2 <- merge( aggregate(ph.adj ~ CCline,Gld_df,length),aggregate(ph.adj ~ CCline,Gld_df,mean),by="CCline")
  colnames(Gld_df2) <- c('CC','weight','ph')
  Gld_df2$Study <- 'GLD'
  
  info.adj <- rbind(Tov_df2,Tak_df2,Gld_df2)
  info.adj$Study <- as.factor(info.adj$Study)
  ph.df <- info.adj
  ph.df$pheno.id <- paste(ph.df$Study,ph.df$CC,sep=".")
  ph.df$SUBJECT.NAME <- as.factor(ph.df$CC)
  ph.df[-which(ph.df$CC=='CC078'),] -> ph.df
  
  Tov_df2 <- ph.df[which(ph.df$Study=='TOV'),]
  Tak_df2 <- ph.df[which(ph.df$Study=='TAK'),]
  Gld_df2 <- ph.df[which(ph.df$Study=='GLD'),]

  # mapping
  MI=TRUE
    
  myweights <- Tov_df2$weight; names(myweights) <- Tov_df2$CC 
  Scan_Tov <- scan.h2lmm(genomecache="~/Dropbox/ValdarLab/Takeda_copy/segments_happy_format_mm10/", 
                           data=Tov_df2, formula= ph ~ 1,geno.id="CC",pheno.id = "CC",weights = myweights,use.multi.impute=MI,num.imp=15,print.locus.fit=FALSE)
    
  myweights <- Tak_df2$weight; names(myweights) <- Tak_df2$CC 
  Scan_Tak <- scan.h2lmm(genomecache="~/Dropbox/ValdarLab/Takeda_copy/segments_happy_format_mm10/", 
                           data=Tak_df2, formula= ph ~ 1,geno.id="CC",pheno.id = "CC",weights = myweights,use.multi.impute=MI,num.imp=15,print.locus.fit=FALSE)
    
  myweights <- Gld_df2$weight; names(myweights) <- Gld_df2$CC 
  Scan_Gld <- scan.h2lmm(genomecache="~/Dropbox/ValdarLab/Takeda_copy/segments_happy_format_mm10/",
                           data=Gld_df2, formula= ph ~ 1,weights = myweights,geno.id="CC",pheno.id = "CC",use.multi.impute=MI,num.imp=15,print.locus.fit=FALSE)
   
  myweights <- ph.df$weight;names(myweights) <- ph.df$pheno.id
  Scan_mega <- scan.h2lmm(genomecache="~/Dropbox/ValdarLab/Takeda_copy/segments_happy_format_mm10/",
                            data=ph.df, formula= ph ~ 1+Study,weights = myweights,pheno.id="pheno.id",use.multi.impute=MI,num.imp=20,print.locus.fit=FALSE)
  saveRDS(Scan_mega,file=paste0("~/Dropbox/ValdarLab/IDSScross/Result/PhenoAuto/",ph,"scan.RDS"))

  #HPeak <- ceiling(max( -log10(Scan_Tov$p.value),-log10(Scan_Tak$p.value), -log10(Scan_Gld$p.value),-log10(Scan_mega$p.value)))
  pdf(paste0("~/Dropbox/ValdarLab/IDSScross/Result/Pheno/",ph,".mapping.mm10.pdf"),width=10,height=14)
  par(mfrow=c(4,1))
  genome.plotter.whole(scan.list=list(MI=Scan_Tov),main="TOV")
  genome.plotter.whole(scan.list=list(MI=Scan_Tak),main="TAK")
  genome.plotter.whole(scan.list=list(MI=Scan_Gld),main="GLD")
  genome.plotter.whole(scan.list=list(MI=Scan_mega),main="Mega")
  dev.off()
}


# threshold for ALT: 9.265661   (gev:-0.2303773  1.2742609  6.5247375)
# threshold for AST: 6.516393   (gev:-0.01232742  0.77762467  4.24846875)

pdf(paste0("~/Dropbox/ValdarLab/IDSScross/Result/Pheno/",ph,".mapping.mm10.pdf"),width=8.5,height=7)
par(mfrow=c(2,1))
genome.plotter.whole(scan.list=list(TOV=Scan_Tov,TAK=Scan_Tak,IDL=Scan_Gld),main="",y.max.manual=6,
                     main.colors = c("blue", "orange", "gray48"))
genome.plotter.whole(scan.list=list(Mega=Scan_mega),main="",y.max.manual=HPeak)
dev.off()

pdf(paste0("~/Dropbox/ValdarLab/IDSScross/Result/Pheno/ALTAST.mapping.mm10.pdf"),width=8.5,height=7)
par(mfrow=c(2,1))
genome.plotter.whole(scan.list=list(Mega=ALT.scan),main="",hard.thresholds=9.265661)
genome.plotter.whole(scan.list=list(Mega=AST.scan),main="",hard.thresholds=6.516393)
dev.off()

