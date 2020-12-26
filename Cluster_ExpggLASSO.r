library(lme4)
library(dplyr)
library(ggplot2)
library(evir)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("At least two argument must be supplied", call.=FALSE)
} else if (length(args)==2) {
  FROM=args[1]
  END=args[2]
}

Expdf.full <- read.table(file="../data/CCexpsum.full.adj.txt",sep="\t",header=TRUE)

CClist <- sort(unique(intersect(Expdf.full$CC[Expdf.full$Study=='TAK'],Expdf.full$CC[Expdf.full$Study=='TOV'])))

info <- Expdf.full[,c(1:7)]
probelist <- colnames(Expdf.full)[8:dim(Expdf.full)[2]]

QTLccm <- as.matrix(t(read.table(file="../data/CC54eqltable.txt",header = T,sep="\t")))
QTLccm <- QTLccm[,- 8*c(1:7641)]

info$Study <- as.character(info$Study)
info$Study2 <- info$Study
info$Study[which(info$Study=='TAK')] <- 0.5
info$Study[which(info$Study=='TOV')] <- -0.5
info$Study <- as.numeric(info$Study)

thisgroup <- c(rep(c(1:7641),each=7))

test.ix <- which(info$Study== 0.5 & !(info$CC %in% CClist))
train.ix <-  which(!(info$Study== 0.5 & !(info$CC %in% CClist)))

library(gglasso)
library(glmnet)

for(pb in probelist[FROM:END]){
  message(pb)
  info$y <- as.numeric(Expdf.full[,pb])
  fit3 <- lmer( y ~ 1 + Study + (1|CC), data=info,REML=FALSE)
  Y <- info$y - info$Study*fixef(fit3)['Study']
  X <- QTLccm[info$CC,]

  # group Lasso
  cv_fit <- cv.gglasso(X[train.ix,], Y[train.ix], thisgroup)
  fit <- gglasso(X[train.ix,], Y[train.ix], thisgroup,lambda=cv_fit$lambda.min)
  R.gl <- cor(predict(fit, s = cv_fit$lambda.min, newx = X[test.ix,]),Y[test.ix])^2
  # lasso
  cv_fit <- cv.glmnet(X[train.ix,], Y[train.ix], type.measure="mse", family="gaussian", alpha=1)
  fit <- glmnet(X[train.ix,], Y[train.ix], family="gaussian", alpha=1 ,lambda=cv_fit$lambda.min)
  R.ls <- cor(predict(fit, newx = X[test.ix,]),Y[test.ix])^2
  # Ridge
  cv_fit <- cv.glmnet(X[train.ix,], Y[train.ix], type.measure="mse", family="gaussian", alpha=0)
  fit <- glmnet(X[train.ix,], Y[train.ix], family="gaussian", alpha=0 ,lambda=cv_fit$lambda.min)
  R.rd <- cor(predict(fit, newx = X[test.ix,]),Y[test.ix])^2
  # EN
  cv_fit <- cv.glmnet(X[train.ix,], Y[train.ix], type.measure="mse", family="gaussian", alpha=0.5)
  fit <- glmnet(X[train.ix,], Y[train.ix], family="gaussian", alpha=0.5 ,lambda=cv_fit$lambda.min)
  R.en <- cor(predict(fit, newx = X[test.ix,]),Y[test.ix])^2
  
  thisdf <- data.frame(pb=pb,R.gl=R.gl,R.ls=R.ls,R.rd=R.rd,R.en=R.en)

  write.table(thisdf,file=paste0("result/ggLASSO/ggLASSO_Study_LRE",FROM,"_",END,".txt"),sep="\t",col.names=F,row.names=F,quote=F,append=T)
}


 
