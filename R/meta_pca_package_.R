################################################################################
################################################################################
################################################################################
#
# META-ANALYTIC PRINCIPAL COMPONENT ANALYSIS IN INTEGRATIVE OMICS APPLICATION
# Version : 1.1.0
# Authors : SungHwan Kim, Dongwan Kang and George C. Tseng
# latest update : 08/09/2016
#
################################################################################
################################################################################
################################################################################

meta.pca <- function(DList, 
                     method, 
                     Meta.Dim, 
                     is.auto.Dim = TRUE, 
                     is.equal.Dim=FALSE, 
                     e.Dim, 
                     is.weight=TRUE,
                     .var.quantile= 0.8,
                     .scaleAdjust=FALSE,
                     is.sparse=FALSE,
                     Lambda)
  { 
  # Purpose 
  #     Dimension reduction and visualizaiton via metaPCA
  # Usage
  #     S = meta.pca(DList,  Method = "Fisher", K.max, is.VO=TRUE, Para.num.gene = 10000)
  # Input
  #     DList = Input data set matrix (a list of multiple datasets; row=features, column=samples).
  #     method = "SSC" = Sum of Squared Cosine, "SV" = Sum of variance.
  #     Meta.Dim = Dimension size of meta-eigenvector matrix.
  #     is.auto.Dim = Logical value whether dimension size of each study's eigenvector matrix (SSC) is determined 
  #                   by the pre-defined level of variance quantile.
  #     is.equal.Dim = Logical value whether dimension size of each study's eigenvector matrix (SSC) is equal across studies.
  #     e.Dim = Dimension size of each study's eigenvector matrix (SSC) when is.equal.Dim = TRUE
  #     is.weight = Logical value whether the reciprocal of the largest eigenvalue is mutiplied to covariance matrix.
  #     .var.quantile = Threshold indicating the minimum variance of individual study, when is.auto.Dim = TRUE
  #     .scaleAdjust = Logical value whether the PC projection is scaled to mean of zero and SD of 1.
  #     is.sparse = Logical value whether meta-eigenvector matrix is penalized to encourage sparseness.
  #     Lambda = Tuning parameter for sparsity
  # Output
  #     res = a list of meta analytic PCA
  #         res$v = Meta eigenvector matrix 
  #         res$coord = Meta PC projection from input data
  
  .X <-list()
  for(i in 1:length(DList)){
    .X[[i]] <- scale(t(DList[[i]]),scale=FALSE)}     
  L.eigen <-list()
  L.value <-list()
  Perturbed.Cov.X <-list()  
  meta <- list()
  
  if(is.auto.Dim & method=="SSC"){
    .L.Dim <- rep(1, length(DList))
    for(i in 1: length(DList)){
      repeat{
        if((sum(((svd(t(DList[[i]]), nu = 0)$d)^2)[1:.L.Dim[i]]) / sum(((svd(t(DList[[i]]), nu = 0)$d)^2))) >= .var.quantile ){break}
        if( .L.Dim[i] == ncol(DList[[i]]) ){break}
        .L.Dim[i] <- .L.Dim[i] + 1
      }
    }
    L.Dim <- .L.Dim
  } 
  if(is.equal.Dim) L.Dim <- rep(e.Dim, length(DList))
  if(method=="SSC"){
    for(J in 1:length(DList)) {    
      L.eigen[[J]] <- svd(.X[[J]], nu = 0)$v[ ,1:L.Dim[J]]
      L.value[[J]] <- (svd(.X[[J]], nu = 0)$d[1])^2
    }  
    Sum <- matrix(0,ncol(.X[[1]]),ncol(.X[[1]]))
    for(i in 1: length(.X) ){
      if(is.weight){  
        Sum <- Sum + (L.eigen[[i]] %*% t(as.matrix(L.eigen[[i]])))*(1/L.value[[i]])}else{
          Sum <- Sum + (L.eigen[[i]] %*% t(as.matrix(L.eigen[[i]])))
        }
    }
    if(is.sparse){
      v <- SPC(Sum, sumabsv=Lambda, K=Meta.Dim, niter=20, orth=TRUE,center=TRUE)$v
    }else{
      v <- eigen(Sum, symmetric = TRUE)$vectors[ ,1:Meta.Dim] 
    }
    } else { if(method=="SV") {
    
    for(J in 1:length(DList)) {    ### Caution the direction for scaling but we don't perform scaling since it's alread done.
      Perturbed.Cov.X[[J]] <- cov(.X[[J]])  
      L.value[[J]] <- (svd(.X[[J]])$d[1])^2
    }
    all.Sum.sv <- matrix(0,ncol(.X[[1]]),ncol(.X[[1]]))
    for(i in 1:length(.X) ){
      all.Sum.sv <- all.Sum.sv + Perturbed.Cov.X[[i]]*(1/L.value[[J]])
    }
    
    if(is.sparse){
      v <- SPC(all.Sum.sv, sumabsv=Lambda, K=Meta.Dim, niter=20, orth=TRUE,center=TRUE)$v
    }else{
      v <- eigen(all.Sum.sv, symmetric = TRUE)$vectors[,1:Meta.Dim]
     }
    }
  }
  
  coord <- list()
  for( i in 1:length(.X)){
    coord[[i]] <- .X[[i]] %*% v }
  
  if(.scaleAdjust){
    .coord <- list()
    for(i in 1:length(coord)) {.coord[[i]] <- scale(coord[[i]],scale=TRUE)}
    coord <- .coord
  }

  names(coord) <- names(DList)
  res<-list()
  res$v <- v 
  res$coord <- coord
  return(res)
}

