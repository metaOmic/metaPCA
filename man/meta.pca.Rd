\name{meta.pca}
\alias{meta.pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
META-ANALYTIC PRINCIPAL COMPONENT ANALYSIS IN INTEGRATIVE OMICS APPLICATION
}
\usage{
meta.pca(DList, method, Meta.Dim, is.auto.Dim = TRUE, is.equal.Dim = FALSE, e.Dim, is.weight = TRUE, .var.quantile = 0.8, .scaleAdjust = FALSE, is.sparse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DList}{
Input data set matrix (a list of multiple datasets; row=features, column=samples).
}
  \item{method}{
"SSC" = Sum of Squared Cosine, "SV" = Sum of variance
}
  \item{Meta.Dim}{
Dimension size of meta-eigenvector matrix
}
  \item{is.auto.Dim}{
Logical value whether dimension size of each study's eigenvector matrix (SSC) is determined 
by an arbitrary variance quantile
}
  \item{is.equal.Dim}{
Logical value whether dimension size of each study's eigenvector matrix (SSC) is equal across studies
}
  \item{e.Dim}{
Dimension size of each study's eigenvector matrix (SSC) when is.equal.Dim = TRUE
}
  \item{is.weight}{
Logical value whether the reciprocal of the largest eigenvalue is mutiplied to covariance matrix
}
  \item{.var.quantile}{
A threshold indicating the minimum variance of individual study, when is.auto.Dim = TRUE
}
  \item{.scaleAdjust}{
Logical value whether the PC projection is scaled to mean of zero and SD of 1
}
  \item{is.sparse}{
Logical value whether meta-eigenvector matrix is penalized to encourage sparseness
}
\item{Lambda}{
Tuning parameter for penalty constant
}
}
\value{
res = a list of meta analytic PCA \cr
res$v = Meta eigenvector matrix  \cr
res$coord = Meta PC projection from input data
}
\author{
SungHwan Kim \email{swiss747@gmail.com}
}
\examples{
library(yeastCC)
data(yeastCC)
data<-Biobase::exprs(yeastCC)
library(impute)
data.na<-is.na(data)
data.na.length<-apply(data.na, 1, sum)
data.sd<-apply(as.matrix(data), 1, sd, na.rm=TRUE)

new.data<-data[data.na.length<77*0.1 & data.sd>0.45,]    
Spellman <- list(alpha=impute.knn(new.data[,5:22])$data,
                cdc15=impute.knn(new.data[,23:46])$data,
                cdc28=impute.knn(new.data[,47:63])$data,
                elu=impute.knn(new.data[,64:77])$data)

## Standard Principal Component Analysis ##
Spellman.prcomp <- list(alpha=svd(t(Spellman$alpha))$v,
                       cdc15=svd(t(Spellman$cdc15))$v,
                       cdc28=svd(t(Spellman$cdc28))$v,
                       elu=svd(t(Spellman$elu))$v)

# res1 <- meta.pca(DList=Spellman, method="SV", Meta.Dim=2, is.auto.Dim = TRUE)
res2 <- meta.pca(DList=Spellman, method="SSC", Meta.Dim=2, is.auto.Dim = TRUE)

par(mfrow=c(6,4), cex=1 ,oma=c(0,6,2,0), mar=c(0.2,0.2,0.2,0.2))
for(i in 1:4) {
 for(j in 1:4) {
   .coord <- scale(t(Spellman[[j]]), scale=FALSE) %*% Spellman.prcomp[[i]][,1:2]
   #.coord <- t(Spellman[[j]]-pcaPP::l1median(t(Spellman[[j]]),trace=-1)) %*% Spellman.prcomp[[i]]$loadings[,1:2]
   plot(.coord[,1], .coord[,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n"
        ,ylim=c(min(.coord[,2])-1.5, max(.coord[,2])+1.5)
        ,xlim=c(min(.coord[,1])-1, max(.coord[,1])+1))
   text(.coord[,1], .coord[,2], 1:nrow(.coord), cex=1)
   lines(.coord[,1], .coord[,2], col='grey55')
 }
}

## SV
#for(j in 1:4) {
#  .coord <- res1$coord[[j]]
#  plot(.coord[,1], .coord[,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n")
#  text(.coord[,1], .coord[,2], 1:nrow(.coord), cex=1)
#  lines(.coord[,1], .coord[,2])
#  }

## SSC
for(j in 1:4) {
 .coord <- res2$coord[[j]]
 plot(.coord[,1], .coord[,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n"
      ,ylim=c(min(.coord[,2])-1.5, max(.coord[,2])+1.5)
      ,xlim=c(min(.coord[,1])-1, max(.coord[,1])+1))
 text(.coord[,1], .coord[,2], 1:nrow(.coord), cex=1)
 lines(.coord[,1], .coord[,2], col='grey55')
}


## Sparse MetaPCA (SSC)
library(PMA)
library(doMC)
res4 <- meta.pca(DList=Spellman, method="SSC", Meta.Dim=2, is.auto.Dim = TRUE, is.sparse=TRUE, Lambda=6)
res4 <- res4$v
coord <- list()
for(i in 1:4) {
 coord[[i]] <- t(Spellman[[i]]) %*% res4
}
for(j in 1:4) {
 .coord <- coord[[j]]
 plot(.coord[,1], .coord[,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n"
      ,ylim=c(min(.coord[,2])-1.5, max(.coord[,2])+1.5)
      ,xlim=c(min(.coord[,1])-1, max(.coord[,1])+1))
 text(.coord[,1], .coord[,2], 1:nrow(.coord), cex=1)
 lines(.coord[,1], .coord[,2], col='grey55')     
}

####################################################################################################
## Searching the optimal tuning parameter based on the proportion of increased explained variance
####################################################################################################

optimal.lambda <- meta.pca.cv(DList=Spellman, method="SSC", Meta.Dim=2, CV_lambda = seq(1,10,1), is.plot=TRUE)
## optimal.lambda = 8
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
