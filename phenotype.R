library(lme4)
library(lmerTest)
library(qvalue)
library(varComp)
library(miqtl)
library(ggplot2)
library(tidyverse)

# read in data table
info <- read.table(file="../data/info_3study.txt",header=TRUE,sep = "\t")

options(contrasts = rep("contr.sum", 2)) # Set contrasts -- see below

### Section 1. Variance components of phenotype models.
phone_list=c("D.BW","BWratio","LWratio","AST","ALT","TBIL")

outdf <- data.frame() # output data frame
for (ph in phone_list){
  info$y <- as.numeric(info[,ph])

  # Adjust for Batch effects 
  fit0 <- lmer( y ~ 1 + Study + (1|CCline) + (0+Study|CCline) + (1|Dosing.Date), data=info,REML=FALSE)
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
  
  # pairwise comparison model 
  options(contrasts = c("contr.sum","contr.poly"))
  info <- info[which(info$Study!='GLD'),]
  fit4 <- lmer( y ~ 1 + StudyN + (1|CCline) + StudyN*(1|CCline), data=info,REML=FALSE)
  fit4.2 <- lmer( y ~ 1 + StudyN + StudyN*(1|CCline), data=info,REML=FALSE)
  anova(fit4,fit4.2)
  ranova(fit4)
  as.data.frame(VarCorr(fit4))
  
  # phenotype analysis for three studies
  fit5 <- lmer( y ~ 1 + Study + (1+Study|CCline), data=info,REML=FALSE)
  fit7 <- lmer( y ~ 1 +  (1|CCline) + (1+Study|CCline), data=info,REML=FALSE)
    
  h2.SxS.p <- ranova(fit5)$"Pr(>Chisq)"[3]
  h2.C.p <- ranova(fit5)$"Pr(>Chisq)"[2]
  h2.S.p <- anova(fit5,fit7)$"Pr(>Chisq)"[2]
    
  varout5 <- as.data.frame(VarCorr(fit5))
  v.CC <- varout5$vcov[which(varout5$grp=="CCline" & varout5$var1=="(Intercept)" & is.na(varout5$var2))]
  v.R <- varout5$vcov[which(varout5$grp=="Residual")]
  v.CCS <- sum(varout5$vcov[which(varout5$grp=='CCline.1')][c(1,2,3,4,4,5,5,6,6)]) # var(X)+var(Y)+2Cov(X+Y)

  v.Study <- var(predict(lm(y~Study,data=info)))
  v.sum <- v.CC + v.R + v.CCS + v.Study
  
  h2.mega.CS <- round(v.CCS/v.sum,3)
  h2.mega.CC <- round(v.CC/v.sum,3)
  h2.mega.Study <- round(v.Study/v.sum,3)

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
                       SxS.123=paste0(h2.mega.CS,p2p(h2.SxS.p)))
  
  #thisdf2 <- data.frame(ph=ph,h2.tov=h2.tov,h2.tov.p=h2.tov.p,h2.tak=h2.tak,h2.tak.p=h2.tak.p,h2.gld=h2.gld,h2.gld.p=h2.gld.p,
  #                     Study.123=mega123[1],Study.123.p=mega123[2],CC.123=mega123[3],CC.123.p=mega123[4],SxS.123=mega123[5],SxS.123.p=mega123[6])
                       
  outdf <- rbind(outdf,thisdf)
}

# shrinkage estimation of 
library(limma)

BWsum <- info %>% filter(!is.na(D.BW)) %>% group_by(Study,CCline) %>% dplyr::summarize(mean = mean(D.BW),var = var(D.BW), obs=length(D.BW))

BWsum$x <- 0;BWsum$x[which(BWsum$Study=='TAK')]<-1;BWsum$x[which(BWsum$Study=='GLD')]<-2;
BWsum$var.shink <- BWsum$var
BWsum$var.shink <- squeezeVar(BWsum$var,(BWsum$obs-1))$var.post
#for (x in 0:2){
#  ix <- which(BWsum$x==x)
#  BWsum$var.shink[ix] <- squeezeVar(BWsum$var[ix],(BWsum$obs[ix]-1))$var.post
#}
BWsum$sd <- sqrt(BWsum$var.shink)
BWsum$cv <- BWsum$sd/BWsum$mean

BWcount <- as.data.frame.matrix(table(BWsum[,c('CCline','x')]))
strainfull <- rownames(BWcount)[which(BWcount[[1]]==1 & BWcount[[2]]==1 & BWcount[[3]]==1)]
BWsum$col <- 0
BWsum$col[which(BWsum$CCline %in% strainfull)] <- 1

BWsumsub <- BWsum[which(BWsum$x==0),]
BWsumsub$x <- 3
BWsum <- rbind(BWsum,BWsumsub)

pdf(file="Result/BWsd.pdf",width=4.5,height=4.5)
p<- ggplot(BWsum,aes(x=x,y=sd,group=CCline))
for (xi in 0:3){p <- p + geom_vline(xintercept = xi, color = "gray", size=1)}
p <- p + geom_point(aes(color=CCline),show.legend = F) + 
    geom_line(aes(),color='gray',show.legend = F,subset(BWsum,x %in% c(0,1) & col==0)) +
  geom_line(aes(),color='gray',show.legend = F,subset(BWsum,x %in% c(1,2) & col==0)) +
  geom_line(aes(),color='gray',show.legend = F,subset(BWsum,x %in% c(2,3) & col==0))
  
p+theme(text=element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_blank())
dev.off()

BW.13 <- merge(BWsum[BWsum$Study=='GLD',],BWsum[BWsum$Study=='TOV',],by="CCline")
thisdf <- BW.13
plot(thisdf$mean.x,thisdf$mean.y,xlab='Body Weight GLD',ylab='Body Weight TOV',
     xlim=c(14,34),ylim=c(14,34))
arrows(thisdf$mean.x-thisdf$sd.x,thisdf$mean.y,thisdf$mean.x+thisdf$sd.x,thisdf$mean.y,code=3,angle = 90,length=0.05)
arrows(thisdf$mean.x,thisdf$mean.y-thisdf$sd.y,thisdf$mean.x,thisdf$mean.y+thisdf$sd.y,code=3,angle = 90,length=0.05)

t.test(thisdf$sd.x/thisdf$mean.x,thisdf$sd.y/thisdf$mean.y,paired = TRUE)
# no difference of standard error

for (ph in phone_list){
  message(ph)
  info.ph <- cbind(info[,1:5],data.frame(ph=info[,ph]))
  # foreach Study, remove Batch effect
  info.ph$ph.adj <- info.ph # add pesudo ph.adj
  
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
  
  #info.adj.raw <- rbind(Tov_df,Tak_df,Gld_df)
  #info.adj.raw$ph <- info.adj.raw$ph.adj
  #info.adj.raw <- info.adj.raw[,c("Animal_Number","Study","CCline","ph")]
  #write.table(info.adj.raw,file=paste0("~/Dropbox/ValdarLab/IDSScross/Result/PhenoAuto/",ph,"phraw.txt"),quote=F,col.names = T,row.names = F)
  
  info.adj <- rbind(Tov_df2,Tak_df2,Gld_df2)
  info.adj$Study <- as.factor(info.adj$Study)
  
  df.list[[ph]] <- info.adj
}
write.table(df.list$AST,file="~/Dropbox/ValdarLab/IDSScross/Result/PhenoAuto/ASTph.txt",quote=F,row.names = F,col.names = T)
write.table(df.list$ALT,file="~/Dropbox/ValdarLab/IDSScross/Result/PhenoAuto/ALTph.txt",quote=F,row.names = F,col.names = T)


for (thistudy in c('TOV','TAK','GLD')){
  tmplist = list()
  for (ph in phone_list){ 
    tmplist[[ph]] <- df.list[[ph]][which(df.list[[ph]]$Study==thistudy),c('CC','ph')] 
    colnames(tmplist[[ph]])[2] <- ph
  }
  thisdf <- tmplist %>% reduce(left_join,by='CC')
  thisdf$Study <- thistudy
  df.list[[thistudy]] <- thisdf
}

#### Calculate h2
ph.df <- read.table(file="~/Dropbox/ValdarLab/IDSScross/Result/PhenoAuto/ALTph.txt",header = T)

### Here is for each phenotype the plot to make
for (ph in phone_list){
  ph.df <- df.list[[ph]]
  ph.df$pheno.id <- paste(ph.df$Study,ph.df$CC,sep=".")
  ph.df$SUBJECT.NAME <- as.factor(ph.df$CC)
  ph.df[-which(ph.df$CC=='CC078'),] -> ph.df
  
  Tov_df2 <- ph.df[which(ph.df$Study=='TOV'),]
  Tak_df2 <- ph.df[which(ph.df$Study=='TAK'),]
  Gld_df2 <- ph.df[which(ph.df$Study=='GLD'),]

  #### Pairwise comparison
  pdf(file=paste0("~/Dropbox/ValdarLab/IDSScross/Result/Pheno/",ph,".comparewisecompare.pdf"),width=10,height=10)
  par(mfrow=c(2,2))
  mergedf <- merge(Tov_df2,Tak_df2,by='CC')
  plot(mergedf$ph.x,mergedf$ph.y,xlab='TOV',ylab='TAK')
  mergedf <- merge(Tov_df2,Gld_df2,by='CC')
  plot(mergedf$ph.x,mergedf$ph.y,xlab='TOV',ylab='GLD')
  mergedf <- merge(Tak_df2,Gld_df2,by='CC')
  plot(mergedf$ph.x,mergedf$ph.y,xlab='TAK',ylab='GLD')
  
  plot(density(c(Tov_df2$ph,Tak_df2$ph,Gld_df2$ph)),ylim=c(0,max(density(Gld_df2$ph)$y,density(Tak_df2$ph)$y,density(Tov_df2$ph)$y)),
               type='n',main='n',xlab='ph')
  lines(density(Tov_df2$ph),col='red',lwd=2)
  lines(density(Tak_df2$ph),col='blue',lwd=2)
  lines(density(Gld_df2$ph),col='green',lwd=2)
  dev.off()
  
  #### Mapping plots
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
  # 
  #  HPeak <- ceiling(max( -log10(Scan_Tov$p.value),-log10(Scan_Tak$p.value), -log10(Scan_Gld$p.value),-log10(Scan_mega$p.value)))
  # 
    pdf(paste0("~/Dropbox/ValdarLab/IDSScross/Result/Pheno/",ph,".mapping.mm10.pdf"),width=10,height=14)
    par(mfrow=c(4,1))
    genome.plotter.whole(scan.list=list(MI=Scan_Tov),main="TOV",y.max.manual=HPeak)
    genome.plotter.whole(scan.list=list(MI=Scan_Tak),main="TAK",y.max.manual=HPeak)
    genome.plotter.whole(scan.list=list(MI=Scan_Gld),main="GLD",y.max.manual=HPeak)
    genome.plotter.whole(scan.list=list(MI=Scan_mega),main="Mega",y.max.manual=HPeak)
  #  dev.off()
}

pdf(paste0("~/Dropbox/ValdarLab/IDSScross/Result/Pheno/",ph,".mapping.mm10.pdf"),width=8.5,height=7)
par(mfrow=c(2,1))
genome.plotter.whole(scan.list=list(TOV=Scan_Tov,TAK=Scan_Tak,IDL=Scan_Gld),main="",y.max.manual=6,
                     main.colors = c("blue", "orange", "gray48"))
genome.plotter.whole(scan.list=list(Mega=Scan_mega),main="",y.max.manual=HPeak)
dev.off()

### Variance component analysis
phone_list=c("D.BW","ALT","AST","TBIL","miR122")
outdf <- data.frame()
for (ph in phone_list){
  info$y <- as.numeric(info[,ph])
  vartotal <- var(info$y)
  
  fit4 <- lmer( y ~ 1 + Study + (1|CCline) + (1|Dosing.Date), data=info,REML=FALSE)
  fit4.1 <- lmer( y ~ 1 + Study + (1+Study|CCline) + (1|Dosing.Date), data=info,REML=FALSE)
  
  myanova<-anova(fit4,fit4.1,refit=FALSE)
  p1 <- -log10(myanova$`Pr(>Chisq)`[2])
  
  varout4 <- as.data.frame(VarCorr(fit4))
  v.CC <- varout4$vcov[which(varout4$grp=="CCline" & varout4$var1=="(Intercept)" & is.na(varout4$var2))]
  v.DD <- varout4$vcov[which(varout4$grp=="Dosing.Date")]
  v.R <- varout4$vcov[which(varout4$grp=="Residual")]
  v.CCS <- varout4$vcov[which(varout4$grp=="CCline" & varout4$var1=="StudyN" & is.na(varout4$var2))]
  varfix4 <- var(predict(lm(y~Study,data=info)))
  v.list <- c(varfix4,v.DD,v.CC,v.CCS,v.R)
  v.list.p <- v.list/sum(v.list)
  
  thisdf <- data.frame(ph=ph,v.S=v.list.p[1],v.Batch=v.list.p[2],v.CC=v.list.p[3],v.CCS=v.list.p[4],v.R=v.list.p[5])
  outdf <- rbind(outdf,thisdf)
}


#ph       v.S     v.Batch       v.CC      v.CCS        v.R
#1   D.BW 0.1633564 0.029044059 0.50794591 0.03433468 0.26531896
#2    ALT 0.2627233 0.047154487 0.27178998 0.02356175 0.39477047
#3    AST 0.2597767 0.020439184 0.23426211 0.02742971 0.45809233
#4   TBIL 0.1975941 0.016268062 0.06546853 0.13847044 0.58219883
#5 miR122 0.7944671 0.004561951 0.04382017 0.07803026 0.07912053
# Variance components of phenotype models analysis
## CC MTY 
#CCmty <- read.table(file="UNCcsbioCCMTY.txt",header=TRUE)
#CCmty$CCline <- substring(CCmty$Strain,1,5)
#CCmty$MTraw <-CCmty$MT
#info <- merge(info,CCmty[,c("CCline","MT","MTraw","Y")],by='CCline')
#info <- info[order(info$SAN),]

fit4 <- lmer( ALT ~ 1 + Study + status + (1|CCline) + (1|Dosing.Date), data=info,REML=FALSE)
summary(fit4)

Var_Random_effect <- as.numeric(VarCorr(fit4))
Var_Residual <- attr(VarCorr(fit4), "sc")^2
Var_Study <- var(predict(lm(ALT~Study,data=info)))
Var_status <- var(predict(lm(ALT~status,data=info)))

Var_status/sum(c(Var_Random_effect,Var_Residual,Var_Study,Var_status))

var(info$ALT)


