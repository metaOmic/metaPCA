\name{meta.pca.cv}
\alias{meta.pca.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A algorithm of searching the optimal tuning parameter for meta-analytic principal component analysis
}
\usage{
meta.pca.cv(DList, method, Meta.Dim, is.auto.Dim = TRUE, is.equal.Dim = FALSE, e.Dim, is.weight = TRUE, .var.quantile = 0.8, .scaleAdjust = FALSE, is.sparse = FALSE, CV_lambda = seq(1, 10, 1), is.plot = TRUE)
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
  \item{CV_lambda}{
A set of candidate tuning parameters in which the best tuning parameter is chosen
}
  \item{is.plot}{
Logical value whether visual scree plot is created
}
}
\value{
The value of best tuning parameter selected among considered parameters
}
\author{
SungHwan Kim \email{swiss747@gmail.com}
}
\examples{
library(yeastCC)
data(yeastCC)
data<-Biobase::exprs(yeastCC)
library(impute)
library(doMC)
data.na<-is.na(data)
data.na.length<-apply(data.na, 1, sum)
data.sd<-apply(as.matrix(data), 1, sd, na.rm=TRUE)

new.data<-data[data.na.length<77*0.1 & data.sd>0.45,]    
Spellman <- list(alpha=impute.knn(new.data[,5:22])$data,
                cdc15=impute.knn(new.data[,23:46])$data,
                cdc28=impute.knn(new.data[,47:63])$data,
                elu=impute.knn(new.data[,64:77])$data)
                
####################################################################################################
## Searching the optimal tuning parameter based on the proportion of increased explained variance
####################################################################################################

optimal.lambda <- meta.pca.cv(DList=Spellman, method="SSC", Meta.Dim=2, CV_lambda = seq(1,10,1), is.plot=TRUE)
## optimal.lambda = 8
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.