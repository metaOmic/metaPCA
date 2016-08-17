## By SungHwan
## A: Tuning parameter is chosen by a rule proposed in the manuscript. I attach an following example code to follow at ease (see Figure S3).

####################################################################################################
## Searching the optimal tuning parameter based on the proportion of increased explained variance
####################################################################################################
library(doMC)
sequence = seq(1,10,1)
var.tmp <-foreach(i= sequence,.combine=rbind) %do% {
  res <- meta.pca(DList=Spellman, method="SSC", Meta.Dim=2, is.auto.Dim = TRUE, is.sparse=TRUE, Lambda=i)
  non.zero <- sum(res$v != 0)
  c(sum(foreach(ii = 1: length(res$coord),.combine=c) %do% {
    sum(diag(var(res$coord[[ii]])))
  }), non.zero)
}
num.nonzero <- var.tmp[,2]
var.tmp <- var.tmp[,1]
# plot(var.tmp ~ num.nonzero, type='o',ylab="Explained Variance",xlab="The number of non-penalized features",cex.lab=2,lwd=2,cex.axis=1.5)
scree <-foreach(i= 1:10,.combine=c) %do% {
  (var.tmp[i+1] - var.tmp[i]) / var.tmp[i+1]
}
par(mfrow=c(1,1), mar=c(5, 5, 2, 4))
plot(na.omit(scree) ~ num.nonzero[1:length(na.omit(scree))], type='o',ylab="Proportion of Increased Explained Variance",xlab="The number of non-penalized features",cex.lab=2,lwd=2,cex.axis=1.5)
abline(h=0.1,col='red',lwd=2)
optimal.lambda <-  sequence[sum(scree> 0.1, na.rm=TRUE)]
## optimal.lambda = 8
}
