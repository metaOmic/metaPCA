## wong06
## Caleb
## 2016/08/02
## 1, perform MetaSparseKmeans

rm(list=ls())

install.packages('/home06/MetaOmics/pkg/MetaPCA/metaPCA_1.0.tar.gz',repos=NULL,type="source")
library(metaPCA)

load('/home06/MetaOmics/dat/internal/leukemia/Leukemia_v2.rda')
load('/home06/MetaOmics/dat/internal/leukemia/leukemiaLabel3.Rdata')

sapply(Leukemia,dim)

Leukemia <- lapply(Leukemia,function(x) t(scale(t(x))))
sapply(Leukemia,dim)


## ask sungHwan:
## 1, make PMA a dependent package
## 2, how to tuning parameter
## 3，is.auto.Dim: what is an arbitrary variance quantile 
## 4，is there a way to scilence you output for sparsePCA? Also print out what are the printed numbers?

## argument on GUI:
## clinical dependent: may dependent on clinical variable. (e.g group label)
## 1, method. 2, Meta.DIM. 3, logical: whether sparse. 4, if 3==TRUE, tuning parameter


## faster
library(PMA)
start <- Sys.time()
res6 <- meta.pca(DList=Leukemia, method="SSC", Meta.Dim=2, is.auto.Dim = TRUE)
end <- Sys.time()
time6 <- end - start ## 7.350209 mins
coord <- res6$coord

png('metaPCA_Leukemia_res6.png')
par(mfrow=c(2,2))
for(i in 1:length(coord)){
	acoord <- coord[[i]]
	alabel <- as.factor(label[[i]])
    plot(acoord[,1], acoord[,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n"
         ,ylim=c(min(acoord[,2])-1.5, max(acoord[,2])+1.5)
         ,xlim=c(min(acoord[,1])-1, max(acoord[,1])+1),main=names(Leukemia)[i])
    points(acoord[,1], acoord[,2], col=as.numeric(alabel), cex=1)	
	legend('topright',legend=levels(alabel),col=unique(as.numeric(alabel)),pch=1)
}
dev.off()


## will take more time

library(PMA)
start <- Sys.time()
res4 <- meta.pca(DList=Leukemia, method="SSC", Meta.Dim=2, is.auto.Dim = TRUE, is.sparse=TRUE, Lambda=6)
end <- Sys.time()
time4 <- end - start ## 58.84823 mins

coord <- res4$coord

png('metaPCA_Leukemia_res4.png')
par(mfrow=c(2,2))
for(i in 1:length(coord)){
	acoord <- coord[[i]]
	alabel <- as.factor(label[[i]])
    plot(acoord[,1], acoord[,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n"
         ,ylim=c(min(acoord[,2])-1.5, max(acoord[,2])+1.5)
         ,xlim=c(min(acoord[,1])-1, max(acoord[,1])+1),main=names(Leukemia)[i])
    points(acoord[,1], acoord[,2], col=as.numeric(alabel), cex=1)	
	legend('topright',legend=levels(alabel),col=unique(as.numeric(alabel)),pch=1)
}
dev.off()


