library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(evir)
library(miqtl)
Expdf.full <- read.table(file="../data/CCexpsum.full.adj.txt",sep="\t",header=TRUE)

CClist <- sort(unique(intersect(Expdf.full$CC[Expdf.full$Study=='TAK'],Expdf.full$CC[Expdf.full$Study=='TOV'])))

Expdf <- Expdf.full[Expdf.full$CC %in% CClist,]
rm(Expdf.full)
#Expdf <- read.table(file="CCexpsum.txt",sep="\t",header=TRUE)

info <- Expdf[,c(1:7)]
info$CC <- as.factor(info$CC)


info$Study <- as.character(info$Study)
info$Study2 <- info$Study
info$Study[which(info$Study=='TAK')] <- 0.5
info$Study[which(info$Study=='TOV')] <- -0.5
info$Study <- as.numeric(info$Study)

probelist <- colnames(Expdf)[8:dim(Expdf)[2]]
anno4 <- read.table(file="../data/TOV_anno4.txt",sep="\t",header=TRUE)
anno4 <- anno4[which(paste0('X',anno4$Probe.Set.ID) %in% probelist),]
all( paste0('X',anno4$Probe.Set.ID) == probelist)

################################################
#### Analysis of lmer models
outdf <- data.frame()
j=0
for (pb in probelist){
  j<- j+1; if (j %% 1000 ==0) {message(j); write.table(outdf,file="../result/Res2.1.expstrainlevel.txt",quote=TRUE,sep="\t",row.names = FALSE,col.names = TRUE)
}
  info$y <- as.numeric(Expdf[,pb])
  
  fit1 <- lmer( y ~ (1|CC),data=info[which(info$Study==0.5),],REML = F)
  fit1.v <- as.data.frame(summary(fit1)$varcor)[,c('grp','vcov')]
  h2.1 <- fit1.v$vcov[1]/sum(fit1.v$vcov)
  h2.1.p <- rand(fit1)$'Pr(>Chisq)'[2]
  
  fit2 <- lmer( y ~ (1|CC),data=info[which(info$Study==-0.5),],REML = F)
  fit2.v <- as.data.frame(summary(fit2)$varcor)[,c('grp','vcov')]
  h2.2 <- fit2.v$vcov[1]/sum(fit2.v$vcov)
  h2.2.p <- rand(fit1)$'Pr(>Chisq)'[2]
  
  fit3 <- lmer( y ~ 1 + Study + (1|CC),data=info,REML = F)
  fit3.v <- as.data.frame(summary(fit3)$varcor)[,c('grp','vcov')]
  h2.3 <- fit3.v$vcov[1]/sum(fit3.v$vcov)
  h2.3.p <- rand(fit3)$'Pr(>Chisq)'[2]
  
  fit4 <- lmer( y ~ 1 + Study + (1+Study|CC), data=info,REML=FALSE)
  myanova<-anova(fit4,fit3,refit=FALSE)
  logPSxS <- -log10(myanova$`Pr(>Chisq)`[2])
  
  ## 
  corf1f2df <- cbind(ranef(fit1)$CC,ranef(fit2)$CC)
  corf1f2df <- corf1f2df[complete.cases(corf1f2df),]
  thiscor <- cor(corf1f2df[[1]],corf1f2df[[2]])
  
  fit5 <- fit4
  varout5 <- as.data.frame(VarCorr(fit5))
  v.CC <- varout5$vcov[1]
  v.CCS <- varout5$vcov[2]
  v.CCScor <- varout5$vcov[3]
  v.R <- varout5$vcov[4]

  v.Study <- var(predict(lm(y~Study,data=info)))

  var.y <- var(info$y,na.rm=TRUE)
  
  thisdf <- data.frame(pb=pb,
                       h2.1=h2.1,h2.1.p= -log10(h2.1.p),
                       h2.2=h2.2,h2.2.p= -log10(h2.2.p),
                       h2.3=h2.3,h2.3.p= -log10(h2.3.p),
                       logPSxS=logPSxS,pcc=thiscor,
                       v.CC=v.CC,v.CCS=v.CCS,v.CCScor=v.CCScor,v.R=v.R,v.Study=v.Study,var.y=var.y)
  outdf <- rbind(outdf,thisdf)
}
dim(outdf)
outdf$CCtotal <- outdf$v.CC + outdf$v.CCS + 2*outdf$v.CCScor
tmp <- outdf[,c('v.CC','v.CCS','v.CCScor','v.Study','v.R')]
tmp$CCtotal <- tmp$v.CC + tmp$v.CCS + 2*tmp$v.CCScor
boxplot(tmp)

write.table(outdf,file="../result/Res2.1.expstrainlevel.txt",quote=TRUE,sep="\t",row.names = FALSE,col.names = TRUE)


outdf <- read.table(file="Result/Res1.2.1.strainlevel.txt",sep="\t",header=TRUE)
fitmp<- glm(pcc ~ h2.1 + h2.1.p  + StrainStudy + h2.1:StrainStudy,family = gaussian,data=outdf)

plot(outdf$h2.mega,outdf.strain$StrainStudy,xlab='Single Study Heritability (TAK)',ylab='Two Study Strain-by-study effect',
     cex=0.5,pch=19,xlim=c(0,1),ylim=c(0,1))
lines(c(-1,2),c(-1,2),lty=2,col='blue')

plot(outdf$h2.3,outdf$StrainStudy,xlab='Single Study Heritability (TAK)',ylab='Two Study Strain-by-study effect',
     cex=0.5,pch=19,xlim=c(0,1),ylim=c(0,1))
lines(c(-1,2),c(-1,2),lty=2,col='blue')

plot(outdf$h2.1, (outdf$h2.mega+outdf$StrainStudy),xlim=c(0,1),ylim=c(0,1),
     xlab='Single Study Heritability (TAK)',ylab='Two Study: Strain + Strain-by-study effect',cex=0.5,pch=19)
points(c(-1,2),c(-1,2),lty=2,type='l',col='blue')

### plot 1.1.1 Variance explained by each parts, histgram of SxS effect
library(reshape2)
outdf2 <- outdf[,c('h2.mega','StudyEff','StrainStudy','Residual')]
colnames(outdf2) <- c("Strain","Study","Strain-by-Study","Residual")
outdf2 <- melt(outdf2,value.name='Variance')
ggplot(outdf2,aes(variable,Variance)) + 
  geom_violin(scale = 'width') + 
  stat_summary(geom = 'crossbar', fun.data = median_hilow, width = 0.1, fun.args = list(conf.int = 0.5)) +
  ggtitle("y ~ 1 + Study + (1 + Study|CC)")
ggsave(filename = "modelvarsum/boxplot4.png", width = 7, height = 5)

## figure 2.1.1 h2 h2 plot
plot(outdf$h2.1,outdf$h2.2,xlab='h2 TAK Study',ylab='h2 TOV Study',pch=19,cex=0.4)
#ix <- outdf$logPSxS>15
#points(outdf$h2.1[ix],outdf$h2.2[ix],pch=19,cex=0.5,col='blue') 

## figure 2.1.2 BLUP cor hitgram
plot(outdf$h2.1,outdf$pcc,xlab='h2 TAK Study',ylab='Correlation of Strain Effect (BLUP pcc)',pch=19,cex=0.4)

# BLUPcor correlated with 
plot(outdf$StrainStudy,outdf$pcc,xlab='Strain-by-Study effect (Variance percentage)',ylab='Correlation of Strain Effect (BLUP pcc)',pch=19,cex=0.4)


## figure 2.1.2 
plot(predict(fitmp,data.frame(h2.1=0.1,h2.1.p=0.01,StrainStudy=seq(0.01,0.1,0.01))),
     xlab='Gene-by-Study Effect (Percentage of variance)',ylab='Reproducibility of strain effect (BLUPpcc)',
     xlim=c(0,75),ylim=c(0,1),type='l')
for (ix in 2:9){points(predict(fitmp,data.frame(h2.1=0.1*ix,h2.1.p=0.01,StrainStudy=seq(0.01,0.1*ix,0.01))),type='l',lty=ix)}



# peak a probe list of height probe for each gene
anno4$mean <- apply(as.matrix(Expdf[which(Expdf$Study=='TAK'),8:dim(Expdf)[2]]),2,mean)
maxpb<-anno4 %>% group_by(Gene.Symbol) %>% summarize(max.pb = Probe.Set.ID[which.max(mean)])
maxpb <- paste0('X',as.character(maxpb$max.pb))

TAK.cor <- cor(as.matrix(Expdf[which(Expdf$Study=='TAK'),8:dim(Expdf)[2]]),method='spearman')
TOV.cor <- cor(as.matrix(Expdf[which(Expdf$Study=='TOV'),8:dim(Expdf)[2]]),method = 'spearman')

library('diffcoexp')
exprs.1 <- t(as.matrix(Expdf[which(Expdf$Study=='TAK'),maxpb]))
exprs.2 <- t(as.matrix(Expdf[which(Expdf$Study=='TOV'),maxpb]))
tmp <- coexpr(exprs.1, exprs.2,rth = 0.5, qth = 0.1)
str(tmp)

diffix <- which(tmp$q.diffcor<0.1)
cor(tmp$cor.1[diffix],tmp$cor.2[diffix])

cor1.df <- melt(cor1)[1:dim(melt(cor1))[1]/2,]
cor2.df <- melt(cor2)[1:dim(melt(cor1))[1]/2,]

library('DiffCorr')

comp.2.cc.fdr(output.file="Result/DiffCorr.test.txt",
              exprs.1,exprs.2,
              method="pearson",
              p.adjust.methods='BH',
              threshold=1)

DFc.obs <- read.table(file="Result/DiffCorr.test.txt",header=TRUE,sep="\t")

sig.ix <- which(abs(DFc.obs$r1)>0.5 | abs(DFc.obs$r2)>0.5)

cor(DFc.obs$r1[sig.ix],DFc.obs$r2[sig.ix])

library('cocor')

data("aptitude")
cocor(~logic + intelligence.a | logic + intelligence.a, aptitude)



# counld get a DiffCorr count
set.seed(12345)
for (j in 1:100){
  message(j)
comp.2.cc.fdr(output.file="Result/DiffCorr.test.rep.txt",
              TAKin[1:8000,],TOVin[sample(1:8000),],
              method="pearson",
              threshold=0.05)
  DFc.rep <- dim(read.table(file="Result/DiffCorr.test.rep.txt",header=TRUE,sep="\t"))[1]
}


Fobs <- sqrt(sum(  (TAK.cor-TOV.cor)^2 ))

# Frobeneous norm

set.seed(1234)
outF <- c()
for (i in 1:1000){
  message(i)
  thiso <- sample(1:dim(TOV.cor)[1])
  thisF <- sqrt(sum(  (TAK.cor-TOV.cor[thiso,thiso])^2 ))  
  outF <- c(outF,thisF)
}

################ Supplimentory Code ########################
## CC MTY 
## dont include this section if not want CC MTY involved
CCmty <- read.table(file="UNCcsbioCCMTY.txt",header=TRUE)
CCmty$CCline <- substring(CCmty$Strain,1,5)
CCmty$MTraw <-CCmty$MT
CCmty$MT <- substring(CCmty$MTraw,1,1)
CCmty <- CCmty[,c("CCline","MT","MTraw","Y")]
colnames(CCmty)[1] <- "CC"
info <- left_join(info,CCmty,by="CC")

##### Test effects of MT genome and whether to include it or not. 
### The standard MT genome mapping
outdf <- data.frame()
for (pb in probelist){
  message(pb)
  info$y[info$Study==0.5] <- rint((Expdf[info$Study==0.5,pb]))
  info$y[info$Study==-0.5] <- rint((Expdf[info$Study==-0.5,pb]))
  
  TAK.df <- merge( aggregate(y ~ CC,info[info$Study==0.5,],length),aggregate(y ~ CC,info[info$Study==0.5,],mean),by="CC")
  TAK.df$Study <- 0.5
  TOV.df <- merge( aggregate(y ~ CC,info[info$Study==-0.5,],length),aggregate(y ~ CC,info[info$Study==-0.5,],mean),by="CC")
  TOV.df$Study <- -0.5
  thisdf <- rbind(TAK.df,TOV.df)
  colnames(thisdf)[2:3] <- c("weight","y")
  thisdf <- merge(thisdf, CCmty[,c("CC","MT")],by="CC")
  
  fit.mega <- lm( y ~ Study + MT,data=thisdf,weights=weight)
  fit.mega.0 <- lm( y ~ Study,data=thisdf,weights=weight)
  F.mega <- -log10(anova(fit.mega,fit.mega.0)$"Pr(>F)"[2])
  
  fit.tak <- lm( y ~ MT,data=thisdf[thisdf$Study==0.5,],weights=weight)
  fit.tak.0 <- lm( y ~ 1,data=thisdf[thisdf$Study==0.5,],weights=weight)
  F.tak <- -log10(anova(fit.tak,fit.tak.0)$"Pr(>F)"[2])
  
  fit.tov <- lm( y ~ MT,data=thisdf[thisdf$Study==-0.5,],weights=weight)
  fit.tov.0 <- lm( y ~ 1,data=thisdf[thisdf$Study==-0.5,],weights=weight)
  F.tov <- -log10(anova(fit.tov,fit.tov.0)$"Pr(>F)"[2])
  
  thisdf <- data.frame(pb=pb,F.mega=F.mega,F.tak=F.tak,F.tov=F.tov)
  outdf <- rbind(outdf,thisdf)
}
write.table(outdf,file="Result/R2.2.MTtablerint.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)

MTtable <- read.table(file="Result/R2.2.MTtablerint.txt",header=TRUE)

MTtable$F.mega.p <- apply(as.matrix(MTtable$F.mega),1,function(x) {length(which(x<MTrintsample$F.mega))})/dim(MTtable)[1]
MTtable$F.tak.p <- apply(as.matrix(MTtable$F.tak),1,function(x) {length(which(x<MTrintsample$F.tak))})/dim(MTtable)[1]
MTtable$F.tov.p <- apply(as.matrix(MTtable$F.tov),1,function(x) {length(which(x<MTrintsample$F.tov))})/dim(MTtable)[1]

megasig <- which(p.adjust(MTtable$F.mega.p,method='BH')<0.05)
tak.sig <- which(p.adjust(MTtable$F.tak.p,method='BH')<0.05)
tov.sig <- which(p.adjust(MTtable$F.tov.p,method='BH')<0.05)
#### analysis of significant MT rint probes ####
outdf <- data.frame()
for (pb in probelist){
  message(pb)
  info$y[info$Study==0.5] <- sample(rint((Expdf[info$Study==0.5,pb])))
  info$y[info$Study==-0.5] <- sample(rint((Expdf[info$Study==-0.5,pb])))
  
  TAK.df <- merge( aggregate(y ~ CC,info[info$Study==0.5,],length),aggregate(y ~ CC,info[info$Study==0.5,],mean),by="CC")
  TAK.df$Study <- 0.5
  TOV.df <- merge( aggregate(y ~ CC,info[info$Study==-0.5,],length),aggregate(y ~ CC,info[info$Study==-0.5,],mean),by="CC")
  TOV.df$Study <- -0.5
  thisdf <- rbind(TAK.df,TOV.df)
  colnames(thisdf)[2:3] <- c("weight","y")
  thisdf <- merge(thisdf, CCmty[,c("CC","MT")],by="CC")
  
  fit.mega <- lm( y ~ Study + MT,data=thisdf,weights=weight)
  fit.mega.0 <- lm( y ~ Study,data=thisdf,weights=weight)
  F.mega <- -log10(anova(fit.mega,fit.mega.0)$"Pr(>F)"[2])
  
  fit.tak <- lm( y ~ MT,data=thisdf[thisdf$Study==0.5,],weights=weight)
  fit.tak.0 <- lm( y ~ 1,data=thisdf[thisdf$Study==0.5,],weights=weight)
  F.tak <- -log10(anova(fit.tak,fit.tak.0)$"Pr(>F)"[2])
  
  fit.tov <- lm( y ~ MT,data=thisdf[thisdf$Study==-0.5,],weights=weight)
  fit.tov.0 <- lm( y ~ 1,data=thisdf[thisdf$Study==-0.5,],weights=weight)
  F.tov <- -log10(anova(fit.tov,fit.tov.0)$"Pr(>F)"[2])
  
  thisdf <- data.frame(pb=pb,F.mega=F.mega,F.tak=F.tak,F.tov=F.tov)
  outdf <- rbind(outdf,thisdf)
}
outdf -> MTrintsample


### tmp for BLUPcor 
# pb, h2.1, h2.2, h2.3,h2.4,cor12,cor13,cor23
set.seed(12345)
nzszlist <- c(1,2,3,4)
outm <- matrix(0,nrow=length(probelist),ncol=length(nzszlist))
row.names(outm) <- probelist
for (k in 1:length(nzszlist)){
  nzsz <- nzszlist[k]
  w4 <- c(rep(1,nzsz),rep(0,(4-nzsz)))
  info$weight <- 0
  for (ws in c(0.5,-0.5)){
  for (wcc in unique(info$CC)){
    ix <- which(info$Study==ws & info$CC==wcc)
    if (length(ix)==4){info$weight[ix]<- sample(w4)}
    if (length(ix)==3){info$weight[ix]<- sample(w4[1:3])}
  }
  }
thispbc <- c()
for (pb in probelist){
  info$y <- as.numeric(Expdf[,pb])
  
  fit1 <- lmer( y ~ (1|CC),data=info[which(info$Study==0.5),],weights = weight,REML = F)

  fit2 <- lmer( y ~ (1|CC),data=info[which(info$Study==-0.5),],weights = weight,REML = F)

  corf1f2df <- cbind(ranef(fit1)$CC,ranef(fit2)$CC)
  corf1f2df <- corf1f2df[complete.cases(corf1f2df),]
  thiscor <- as.numeric(cor(corf1f2df[[1]],corf1f2df[[2]]))
  thispbc <- c(thispbc,thiscor)
}
outm[,k] <- thispbc
save(outm,file="outm_sizenoisetmp.Rdata")
}
plot(density(outm[,4],na.rm =T),lwd=2,ylim=c(0,1.8),xlab='Strain Mean Correlation',main='')
lines(density(outm[,3]),lwd=2,col='red')
lines(density(outm[,2]),lwd=2,col='orange')
lines(density(outm[,1]),lwd=2,col='blue')
abline(v=0,lty=3)
legend('topright',c('Rep 4','Rep 3','Rep 2','Rep 1'),col=c('black','red','orange','blue'),lwd=2,)

outdf <- data.frame()

### Check info of y ###
library(shrink)
outdf <- data.frame()

shinkvar <- function(thisinfo){
  fit.tmp <- lm(y~CC,data=thisinfo,x=TRUE,y=TRUE)
  thissk <-shrink(fit.tmp,postfit = T)
  thisinfo$yp <- predict(thissk,type='response')
  outsum <- thisinfo %>% group_by(CC) %>% dplyr::summarise(sum=sum((y-yp)^2),size=length(y))
  outsum$var <- outsum$sum/max((outsum$size-1),1)
  outsum[,c("CC",'var')]
}
  
vdf <- data.frame(CC=rownames(ssm))

for (pb in probelist){
  message(dim(vdf))
  info$y <- as.numeric(Expdf[,pb])
  
  fit4 <- lmer( y ~ 1 + Study + (1+Study|CC), data=info,REML=FALSE)
  info$y <- info$y-info$Study*fixef(fit4)[2]
  
  ss <- merge(shinkvar(info[info$Study==0.5,]),shinkvar(info[info$Study==-0.5,]),by='CC')
  ssm <- ss[,2:3]
  row.names(ssm) <- ss[[1]]
  colnames(ssm) <- paste0(pb,c('.tak','.tov'))
  
  vdf <- cbind(vdf,ssm)
}
  
load(file = "Result/vdf.tmp.Rdata")

vdfm <- as.matrix(vdf[,-1])
vdfm.tak <- vdfm[,2*c(1:16615)-1]
vdfm.tov <- vdfm[,2*c(1:16615)]
vdfm.tak.rank <- apply(vdfm.tak,2,rank)
vdfm.tov.rank <- apply(vdfm.tov,2,rank)
plot(apply(vdfm.tak.rank,1,mean),apply(vdfm.tov.rank,1,mean),type='n',xlab='Rank mean TAK',ylab='Rank mean TOV')
text(apply(vdfm.tak.rank,1,mean),apply(vdfm.tov.rank,1,mean),substr(vdf$CC,4,5))

# get entropy score across strains?
entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}
Chisq.my <- function(target){
  freq <- table(target)/length(target)
  vec <- as.data.frame(freq)[,2]
  E <- rep(1/length(vec),length(vec))
  sum((vec-E)^2/E)
}
ie.tak <- apply(vdfm.tak.rank,1,entropy)
ie.tov <- apply(vdfm.tov.rank,1,entropy)

chisq.tak <- apply(vdfm.tak.rank,1,Chisq.my)
chisq.tov <- apply(vdfm.tov.rank,1,Chisq.my)

plot(ie.tak,ie.tov,type='n',xlab='Shannon Entropy TAK',ylab='Entropy TOV')
text(ie.tak,ie.tov,substr(vdf$CC,4,5))
plot(chisq.tak,chisq.tov,type='n',xlab='Chisq TAK',ylab='Chisq TOV')
text(chisq.tak,chisq.tov,substr(vdf$CC,4,5))
# correlation table of each CC
vdf2 <- vdf[,1:2]
vdf2$rank.cor <- 0;vdf2$var.cor <- 0;
vdf2 <- vdf2[,-2]
for (j in 1:dim(vdf2)[1]){
  vdf2[j,2] <- cor(vdfm.tak.rank[j,],vdfm.tov.rank[j,])
  vdf2[j,3] <- cor(vdfm.tak[j,],vdfm.tov[j,])
}
ggplot(vdf2) + 
  geom_col(aes(x = CC, y = rank.cor), size = 1, color = "darkblue", fill = "white") +
  geom_line(aes(x = CC, y = var.cor), size = 1.5, color="red", group = 1) +
  theme(axis.text.x = element_text(angle = 60))


Rank.CC055 <- data.frame(tak=vdfm.tak.rank[31,],tov=vdfm.tov.rank[31,])
plot(table(Rank.CC055))

Rank.CC052 <- data.frame(tak=vdfm.tak.rank[29,],tov=vdfm.tov.rank[29,])
plot(table(Rank.CC052))

Rank.CC039 <- data.frame(tak=vdfm.tak.rank[24,],tov=vdfm.tov.rank[24,])
plot(table(Rank.CC039))

Rank.CC055$mean.TAK <- apply(Expdf[53:56,-(1:7)],2,mean)
ix <- which(Rank.CC055$tak>32 & Rank.CC055$tov>32)
Rank.CC055$mean.TAK[ix]
probemean <- apply(Expdf[,-c(1:7)],2,mean)

plot(probemean,Rank.CC055$mean.TAK,xlab='Probe Mean TAK',ylab='CC055 mean TAK',main='High variance of CC005 isnot due to high mean
     (Red: rank>32)  ')
points(c(-1,1000),c(-1,1000),type='l',lty=2)
points(probemean[ix],Rank.CC055$mean.TAK[ix],col='red')





boxplot(t(vdfm.tak.rank),horizontal=TRUE)

  sddf <- info %>% group_by(Study,CC) %>% dplyr::summarise(mean=mean(y),sd=sd(y))
  tt <- merge(sddf[sddf$Study=='0.5',],sddf[sddf$Study=='-0.5',],by='CC')
  tt <- tt[complete.cases(tt),]
  thisdf <- data.frame(pb=pb,mean.tak=mean(tt$mean.x),sd.tak=mean(tt$sd.x),mean.tov=mean(tt$mean.y),sd.tov=mean(tt$sd.y),
                       cor.sd=cor(tt$sd.x,tt$sd.y),cor.var.shink=cor(ss$var.x,ss$var.y,use="complete.obs"))
  outdf <- rbind(outdf,thisdf)
}
write.table(outdf,file="Result/Res1.2.1.2.sdvartable.txt",sep="\t",col.names = T,row.names = FALSE,quote = F)

plot(outdf$sd.tak,outdf$sd.tov,xlim=c(0,1.83),ylim=c(0,1.83),pch=19,cex=0.5,xlab='SD TAK',ylab='SD TOV')
points(c(-10,10),c(-10,10),lty=2,type='l')

outdf -> sddf
outdf <- read.table(file="Result/R2.2.MTtable.txt",sep="\t",header=TRUE)
outdf <- merge(outdf, sddf,by='pb')

