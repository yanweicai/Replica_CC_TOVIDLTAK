library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(evir)
library(miqtl)

chrlist <- c(1:19,'X')
chrsize <-c(195.471971,182.113224,160.03968,156.508116,151.834684,149.736546,145.441459,129.401213,124.59511,130.694993,122.082543,120.129022,120.421639,124.902244,104.043685,98.207768,94.987271,90.702639,61.431566,171.031299)

chrtototal <- function(chr,pos){
  input <- data.frame(chr=chr,pos=pos)
  input$out <- 0
  for (i in 1:dim(input)[1]){
    if (input$chr[i] %in% chrlist){
      input$out[i] <- input$pos[i] + sum(chrsize[0:(which(chrlist==input$chr[i])-1)])
    }
  }
  return(as.numeric(input$out))
}

setwd("~/Dropbox/ValdarLab/IDSScross/")
Expdf.full <- read.table(file="CCexpsum.full.txt",sep="\t",header=TRUE,stringsAsFactors = F)

CClist <- sort(unique(intersect(Expdf.full$CC[Expdf.full$Study=='TAK'],Expdf.full$CC[Expdf.full$Study=='TOV'])))

info.full <- Expdf.full[,c(1:7)]
Expdf <- Expdf.full[Expdf.full$CC %in% CClist,]

info <- Expdf[,c(1:7)]

info$Study <- as.character(info$Study)
info$Study2 <- info$Study
info$Study[which(info$Study=='TAK')] <- 0.5
info$Study[which(info$Study=='TOV')] <- -0.5
info$Study <- as.numeric(info$Study)

probelist <- colnames(Expdf)[8:dim(Expdf)[2]]

### section 0. plots
outdf.strain <- read.table(file="Result/Res1.2.1.strainlevel.txt",sep="\t",header=TRUE)
outdf.strain$quan <- '50-100'
outdf.strain$quan[which(outdf.strain$h2.mega>=quantile(outdf.strain$h2.mega,0.5))] <- '25-50%'
outdf.strain$quan[which(outdf.strain$h2.mega>=quantile(outdf.strain$h2.mega,0.75))] <- '10-25%'
outdf.strain$quan[which(outdf.strain$h2.mega>=quantile(outdf.strain$h2.mega,0.9))] <- '05-10%'
outdf.strain$quan[which(outdf.strain$h2.mega>=quantile(outdf.strain$h2.mega,0.95))] <- '0-5%'

# --> Violin plot of strainmean correlation
ggplot(outdf.strain,aes(x=quan,y=pcc))+
geom_violin(width=0.8,draw_quantiles = c(0.25, 0.5, 0.75))

outdf.gglasso <- read.table(file="Result/ggLASSO_RES_outsample_LRE.txt",sep="\t")[,1:2]
colnames(outdf.gglasso) <- c('pb','outsR')
outdf.strain <- merge(outdf.strain,outdf.gglasso,by='pb')
outdf.gglasso <- read.table(file="Result/ggLASSO_insample_LRE.txt",sep="\t")[,1:2]
colnames(outdf.gglasso) <- c('pb','insR')
outdf.strain <- merge(outdf.strain,outdf.gglasso,by='pb')

# --> Violin plot of R-square of out-of-sample andin-sample group lasso prediction
ggplotdf <- outdf.strain[which(outdf.strain$quan!='50-100'),c('quan','outsR','insR')]
library(reshape2)
ggplotdf <- melt(ggplotdf)
ggplotdf <- ggplotdf[c(which(ggplotdf$variable=='insR'),which(ggplotdf$variable=='outsR')),]
ggplot(ggplotdf,aes(x=quan,y=value,fill=variable,order=variable))+
  geom_violin(width=0.8,draw_quantiles = c(0.25, 0.5, 0.75))

# compare group-Lasso, Lasso, Ridge, and Elastic net (0.5)
outdf.lassos <- read.table(file="Result/ggLASSO_RES_outsample_LRE.txt",sep="\t")
colnames(outdf.lassos) <- c('pb','o.gLss','o.Lss','o.Rdg','o.Els')
outdf.strain <- merge(outdf.strain,outdf.lassos,by='pb')

ggplotdf <- outdf.strain[which(outdf.strain$quan %in% c('0-5%')),c('o.gLss','o.Lss','o.Rdg','o.Els')]
boxplot(ggplotdf)

### Section 1. eQTL mapping ###
anno4 <- read.table(file="TOV_anno4.txt",sep="\t",header=TRUE)
anno4$ProbeX <- paste0('X',anno4$Probe.Set.ID)

anno4 <- anno4[which(anno4$ProbeX %in% probelist),]
all(probelist==anno4$ProbeX)

anno4$plot.pb <- chrtototal(anno4$Achr,(anno4$Apos+anno4$Apos2)/2)

thisprobelist <- probelist
outpd <- data.frame()

outdf <- data.frame()

chrlist <- c(1:19,'X')

# Annotation of markers
load(file="CC.scan.Rdata")
mk2chr <- CC.mean.scan$chr
names(mk2chr) <- as.character(names(CC.mean.scan$p.value))
mk2pos <- CC.mean.scan$pos$Mb 
names(mk2pos) <- as.character(names(CC.mean.scan$p.value))

tablepos <- data.frame(mk=CC.mean.scan$loci,chr=CC.mean.scan$chr,pos=CC.mean.scan$pos$Mb)
tablepos$plotpos <- chrtototal(tablepos$chr,tablepos$pos)

outdf.strain <- read.table(file="Result/Res1.2.1.strainlevel.txt",sep="\t",header=TRUE)

# GEV parameters in order of: shape parameter xi=,scale parameter sigma=,location parameter mu=
TAK.gev <- c(-0.03154331,0.64054052,4.09265555) 
TOV.gev <- c(-0.04022183,0.62374341,4.03302492) 
Mega.gev <- c(-0.08838203,1.34577122,3.55855865) 
Mega.adj.gev <- c(-0.09964554,1.29155727,3.52161113)

# QTL analysis from full Matrix of p-values
for (Stdy in c('TAK','TOV','Mega')){
  message(Stdy)
  thispm <- read.table(file=paste0("~/Dropbox (ValdarLab)/YanweiData/IDSScross/",Stdy,"scan.rint.pmatrix.txt"))
  if (Stdy=='TAK'){this.gev <- c(-0.03154331,0.64054052,4.09265555)}
  if (Stdy=='TOV'){this.gev <- c(-0.04022183,0.62374341,4.03302492)}
  if (Stdy=='Mega'){this.gev <- c(-0.08838203,1.34577122,3.55855865)}
  if (Stdy=='Mega.agj'){this.gev <- c(-0.09964554,1.29155727,3.52161113)}

  this.heights <- apply(thispm,1,max) 

  pgev <- evir::pgev(this.heights,xi=this.gev[1],mu=this.gev[3],sigma=this.gev[2])
  if (length(is.na(pgev))>0){pgev[is.na(pgev)] <- 1-1e-16}

  this.maxdf <- data.frame(height=this.heights,gev.p=1-pgev)
  this.maxdf$q <- p.adjust(this.maxdf$gev.p,method='BH')
  write.table(this.maxdf,file=paste0("Result/eQTLmapping_plot/",Stdy,".maxdf.txt"),col.names = T,row.names = F,quote=F,sep="\t")
  
  # function to get p and q cutoff for 
  getafromb <- function(anumber,aname,bname,this.maxdf){
    tmpdf <- this.maxdf
    tmpdf$delta <- abs(anumber-tmpdf[,aname])
    tmpdf <- tmpdf[order(tmpdf$delta),]
    lmdf <- tmpdf[which(tmpdf$delta %in% unique(tmpdf$delta)[1:3]),]
    fitc <- as.numeric(coefficients(lm(lmdf[,bname]~lmdf[,aname])))
    anumber*fitc[2]+fitc[1]
  }
  
  plotcut1 <- getafromb(0.1,'gev.p','height',this.maxdf)
  plotcut2 <- getafromb(0.2,'q','height',this.maxdf)
  plotcut3 <-getafromb(0.1,'q','height',this.maxdf)

  message(paste0(plotcut1,"\t",plotcut2,"\t",plotcut3,"\n"))

  InMx <- as.matrix(thispm)
  outdf <- data.frame()
  
  for (j in 1:dim(InMx)[1]){
    message(j)
    ix <- which(InMx[j,]>=plotcut1)
    if (length(ix)==0){next}
    thisdf <- data.frame(pbix=j,mkix=ix,height=InMx[j,ix])
    # get alleic effects
    info$y <- as.numeric(Expdf[,(7+j)])
    if (Stdy %in% c('TAK','TOV')){
      Thisscan <- scan.h2lmm(genomecache="~/Dropbox/ValdarLab/Takeda_copy/segments_happy_format_mm10/",
                        data=info[which(info$Study2==Stdy),], just.these.loci= tablepos$mk[ix],
                        formula= rint(y)~ 1, return.allele.effects=TRUE,geno.id="CC",pheno.id = 'An',use.multi.impute=FALSE,print.locus.fit=FALSE)
    }else{
      Thisscan <- scan.h2lmm(genomecache="~/Dropbox/ValdarLab/Takeda_copy/segments_happy_format_mm10/",
                           data=info, just.these.loci= tablepos$mk[ix],
                           formula= rint(y)~ 1 + Study, return.allele.effects=TRUE,geno.id="CC",pheno.id = 'An',use.multi.impute=FALSE,print.locus.fit=FALSE)
    }
    thisdf$afv <- apply(Thisscan$allele.effects,2,function(x) {paste0(x,collapse = ',')})
    outdf <- rbind(outdf,thisdf)
  }
  pgev <- evir::pgev(outdf$height,xi=this.gev[1],mu=this.gev[3],sigma=this.gev[2])
  if (length(is.na(pgev))>0){pgev[is.na(pgev)] <- 1-1e-16}
  outdf$gev.p <- 1-pgev
  for (j in 1:dim(outdf)[1]){ 
    message(j)
    outdf$q[j] <- getafromb(outdf$height[j],'height','q',this.maxdf)
  }
  write.table(outdf,file=paste0("Result/eQTLmapping_plot/",Stdy,".sigpair.txt"),col.names = T,row.names = F,quote=F,sep="\t")

  Indf <- outdf[which(outdf$height>=plotcut1),]
  Indf$plot.pb <- anno4$plot.pb[Indf$pbix]
  Indf$plot.m <- tablepos$plotpos[Indf$mkix]

  # --> eQTL mapping plot
  pdf(paste0("Result/eQTLmapping_plot/",Stdy,".matrixplot.pdf"),height=9,width=9)
  plot(Indf$plot.m,Indf$plot.pb,type='n',xaxt='n',yaxt='n',bty='n',xlab="QTL location",ylab="Transcript Location",xlim=c(0,sum(chrsize)),ylim=c(0,sum(chrsize)))
  abline(h=chrtototal(chrlist,0),lwd=1,col='gray',lty=2)
  abline(v=chrtototal(chrlist,0),lwd=1,col='gray',lty=2)
  axis(1, at=chrtototal(chrlist,chrsize/2), las=1,,lwd=0,labels =chrlist)
  axis(2, at=chrtototal(chrlist,chrsize/2), las=2,lwd=0,labels =chrlist )
  # add color bar
  col=c('#E7C1A5','#E1B357','#83C7ED','#C9DC91','#7AE0D5','#6FE09E','#D843D4','#998F83','#D67DD2','#806CD9','#C9E5CE','#DFA3C6','#8A3DEA','#8F99DE','#DAE053','#D7D3E3','#7FE751','#D15185','#DD694D','#61A0B2')
  segments(chrtototal(chrlist,0),-4,c(chrtototal(chrlist[2:20],0),sum(chrsize)),-4,lwd=4,col=col)
  segments(-4,chrtototal(chrlist,0),-4,c(chrtototal(chrlist[2:20],0),sum(chrsize)),lwd=4,col=col)

  segments(chrtototal(chrlist,0),sum(chrsize)+4,c(chrtototal(chrlist[2:20],0),sum(chrsize)),sum(chrsize)+4,lwd=4,col=col)
  segments(sum(chrsize)+4,chrtototal(chrlist,0),sum(chrsize)+4,c(chrtototal(chrlist[2:20],0),sum(chrsize)),lwd=4,col=col)

  #ix1 <- which(tablepos$q.eqtl<0.1)#
  ix1 <- which(Indf$height>=plotcut1)
  points(Indf$plot.m[ix1],Indf$plot.pb[ix1],pch=19,col='gray',cex=0.5)
  #ix2 <- which(tablepos$q.eqtl<0.01)#
  ix2 <- which(Indf$height>=plotcut2)
  points(Indf$plot.m[ix2],Indf$plot.pb[ix2],pch=19,col='blue',cex=0.7)
  #ix3 <- which(tablepos$q.eqtl<0.001)#
  ix3 <- which(Indf$height>=plotcut3)
  points(Indf$plot.m[ix3],Indf$plot.pb[ix3],pch=19,col='red',cex=0.7)
  dev.off()
}

# Section 2. Compare Replicability of eQTLs
sigdftopeak <- function(sigdf){
  outdf <- data.frame()
  for (thispb in unique(sigdf$pbix)){
    thisdf <- sigdf[sigdf$pbix==thispb,]
    while(dim(thisdf)[1]>0){
      maxix <- which.max(thisdf$height)
      outdf <- rbind(outdf,thisdf[maxix,])
      # find connect peaks with max
      ixleft <- maxix
      ixright <- maxix
      if (maxix==1){
        ixleft <- maxix
      }else{
        for (it in 1:1000){
          if (maxix-it==0){break}
          if ( thisdf$mkix[(maxix-it)] != (thisdf$mkix[maxix]-it) ){break}
          ixleft <- maxix-it
        }}
      if (maxix==dim(thisdf)[1]){
        ixright <- maxix
      }else{
        for (it in 1:1000){
          if(maxix+it>dim(thisdf)[1]){break}
          if ( thisdf$mkix[(maxix+it)] != (thisdf$mkix[maxix]+it) ){break}
          ixright <- maxix+it
        }}
      thisdf <- thisdf[-(ixleft:ixright),]
    }
  }
  outdf
}

sigdflcldstl <- function(TAKsig){
  TAKsig$lcldstl <- 'local'
  for (i in 1:dim(TAKsig)[1]){
    if (tablepos$chr[TAKsig$mkix[i]] != anno4$Achr[TAKsig$pbix[i]]){
      TAKsig$lcldstl[i]='distal'
    }else{
      if ( abs(tablepos$pos[TAKsig$mkix[i]]-anno4$Apos[TAKsig$pbix[i]])>=15){
        TAKsig$lcldstl[i]='distal'
      }
    }
  }
  TAKsig
}

TAKsig <- read.table(file="Result/eQTLmapping_plot/TAK.sigpair.txt",header=TRUE,stringsAsFactors = FALSE)
TOVsig <- read.table(file="Result/eQTLmapping_plot/TOV.sigpair.txt",header=TRUE,stringsAsFactors = FALSE)

for (j in 1:dim(tablepos)[1]){
  tablepos$TAKc[j] <- length(which(TAKsig$mkix==j))
  tablepos$TOVc[j] <- length(which(TOVsig$mkix==j))
}

# get pairs of significant peaks
TAKsig <- read.table(file="Result/eQTLmapping_plot/TAK.sigpair.txt",header=TRUE,stringsAsFactors = FALSE)
TAKsig <- TAKsig[which(TAKsig$height>5.48413662157165),]
TAKsig <- sigdflcldstl(TAKsig)
TAKsig <- sigdftopeak(TAKsig) # 37distal 475local

Megasig <- read.table(file="Result/eQTLmapping_plot/Mega.sigpair.txt",header=TRUE,stringsAsFactors = FALSE)
Megasig <- Megasig[which(Megasig$height>7.99231178568917),]
Megasig <- sigdflcldstl(Megasig)
Megasig <- sigdftopeak(Megasig) 

comdf <- merge(TAKsig[,1:6],Megasig[1:6],by=c('pbix','mkix'),all=TRUE)
colnames(comdf)[c(3,7)] <- c('TAK.height','Mega.height')
# add TAKdf and Megadf significance
this.maxdf <- read.table(file="Result/eQTLmapping_plot/TAK.maxdf.txt",header=TRUE)
for (j in which(is.na(comdf$TAK.height))){
  comdf$TAK.height[j] <- TAKpm[comdf$pbix[j],comdf$mkix[j]]
  comdf$q.x[j] <- getafromb(comdf$TAK.height[j],'height','q',this.maxdf)
}
this.maxdf <- read.table(file="Result/eQTLmapping_plot/Mega.maxdf.txt",header=TRUE)
for (j in which(is.na(comdf$Mega.height))){
  comdf$Mega.height[j] <- Megapm[comdf$pbix[j],comdf$mkix[j]]
  comdf$q.y[j] <- getafromb(comdf$Mega.height[j],'height','q',this.maxdf)
}
plot(-log10(comdf[,c('q.x','q.y')]),xlab='-log10(q.TAK)',ylab='-log10(q.Mega)')
points(c(-10,10),c(-10,10),type='l',lty=2)
# order: xi sigma mu

corcutoff <- 0.7
checkdf1df2 <- function(sigdf1,sigdf2){
  sigdf1$REP <- 'FALSE'
  localix <- which(sigdf1$lcldstl=='local')
  for (j in localix){
    x1 <- as.numeric(strsplit(as.character(sigdf1$afv[j]),split = ',')[[1]])
    afvs <- as.character(sigdf2$afv[which(sigdf2$pbix==sigdf1$pbix[j] & sigdf2$lcldstl=='local')])
    if (length(afvs)==0){next}
    for(x2 in afvs){
      thisx2 <- as.numeric(strsplit(x2,split = ',')[[1]]) # are NAs
      if ( cor(x1,thisx2,use = "complete.obs",method='pearson')>corcutoff){
         sigdf1$REP[j] <- 'TRUE';break;}
    }
  }
  
  distalix <- which(sigdf1$lcldstl=='distal')
  for (j in distalix){
    x2df <- sigdf2[which(sigdf2$pbix==sigdf1$pbix[j]),]
    if (dim(x2df)[1]==0){next;}
    # filter by postion
    thismdf <- tablepos[sigdf1$mkix[j],]
    othermdf <- tablepos[x2df$mkix,]
    x2dfix <- which(othermdf$chr==thismdf$chr & abs(othermdf$pos-thismdf$pos)<=30)
    x2df <- x2df[x2dfix,]
    if (dim(x2df)[1]==0){next;}
    afvs <- x2df$afv
    x1 <- as.numeric(strsplit(as.character(sigdf1$afv[j]),split = ',')[[1]])
    for(x2 in afvs){
      thisx2 <- as.numeric(strsplit(x2,split = ',')[[1]]) # are NAs
      if ( cor(x1,thisx2,use = "complete.obs",method='pearson')>corcutoff){
        sigdf1$REP[j] <- 'TRUE';break;}
    }
  }
  sigdf1
}
TAKsig <- checkdf1df2(TAKsig,Megasig)
table(TAKsig[,c('lcldstl','REP')])

thisout<- data.frame(threshold=c(5.48413662157165,5.908781,7.25814538725819,7.87659507878913),peakhv=0,cishv=0,transhv=0)
for (i in 1:dim(thisout)[1]){
  thisdf <- TAKsig[which(TAKsig$height>thisout$threshold[i]),]
  thisout$peakhv[i] <- sum(thisdf$REP=='TRUE')/dim(thisdf)[1]
  thisout$cishv[i] <- sum(thisdf$REP[which(thisdf$lcldstl=='local')]=='TRUE')/sum(thisdf$lcldstl=='local')
  thisout$transhv[i] <- sum(thisdf$REP[which(thisdf$lcldstl=='distal')]=='TRUE')/sum(thisdf$lcldstl=='distal')
}

row.names(thisout) <- c('p<0.1','p<0.05','q<0.2','q<0.1')

# make replicability plot
pdf(file="Result/eQTLmapping_plot//TAKmegacompare.pdf",height = 5,width=5)
plot(thisout$peakhv,type='b',xlab='',ylab='Replicate Percentage',ylim=c(0,1),pch=1,xaxt='n',yaxt='n',lwd=2.5,main='TAK v.s. Mega(q<0.2)')
axis(side=1,at=1:4,labels=rownames(thisout),las=1)
axis(side=2,las=1)
lines(thisout$cishv, type = "b", col = "blue",pch=3,lwd=2.5)
lines(thisout$transhv, type = "b", col = "orange",pch=2,lwd=2.5)
legend('bottomright',col=c("black","blue","orange"),bg="transparent",legend = c('All eQTLs','local eQTLs','distal eQTLs'),lwd=2,pch=c(1,3,2),cex=1)
dev.off()

### Section 3. QTL hotpots analysis 
TAKpm <- read.table(file="~/Dropbox (ValdarLab)/YanweiData/IDSScross/TAKscan.rint.pmatrix.txt")
TAKqm <- apply(t(TAKpm),1,function(x) {
  p.adjust(10^-x,method='BH')
})
out1 <- apply(TAKqm,2,function(x){  sum(-log10(p.adjust(x,method='BH'))) })

TOVpm <- read.table(file="~/Dropbox (ValdarLab)/YanweiData/IDSScross/TOVscan.rint.pmatrix.txt")

tablepos$TAKc <- apply(TAKpm,2,function(x){  length(which(p.adjust(10^-x,method='BH')<0.05)) })
tablepos$TOVc <- apply(TOVpm,2,function(x){  length(which(p.adjust(10^-x,method='BH')<0.05)) })

# plots of Figure 5B
pdf("Result/QTPhotspot_mega.pdf",height=4,width=9)
plot(tablepos$plotpos,tablepos$Megac,ylim=c(1,max(tablepos$TAKc,tablepos$TOVc)),xlim=c(0,max(tablepos$plotpos)),type='n',xaxt='n',ylab='Counts of adjust.p<0.05',xlab='',bty="l")
for (chr in c(1:19,'X')){
  thistable <- tablepos[which(tablepos$chr==chr),]
  for (j in 1:dim(thistable)[1]){
    points(thistable$plotpos[j],thistable$TOVc[j],type='p',col='blue',pch=1)
    points(thistable$plotpos[j],thistable$TAKc[j],type='p',col='orange',pch=2)
  }
}
axis(1, at=chrtototal(chrlist,0), las=1,lwd=1,labels =rep("",20))
axis(1, at=chrtototal(chrlist,chrsize/2), las=1,lwd=0,labels =chrlist )
abline(v=chrtototal(chrlist,0),lwd=1,col='gray',lty=2)
legend('topleft',c('TOV','TAK'),pch=c(1,2),col=c('blue','orange'),cex=0.9)
dev.off()


Megapm <- read.table(file="~/Dropbox (ValdarLab)/YanweiData/IDSScross/Megascan.rint.pmatrix.txt")

# Section 4. 
