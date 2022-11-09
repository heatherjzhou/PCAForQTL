#' Group predict
#'
#' Given dataResponse and dataPredictors, each a matrix, this function predicts each variable in dataResponse with all the variables in dataPredictors using a multiple linear regression.
#'
#' @param dataResponse The matrix containing the response variables (must be observation by feature).
#' @param dataPredictors The matrix containing the predictors (must be observation by feature as well). The rows (i.e., observations) in dataResponse and dataPredictors must match.
#' @param R2Type Either "unadjusted" or "adjusted".
#'
#' @return This function returns a vector of unadjusted or adjusted squares, each corresponding to a response variable.
#'
#' @references
#' Heather J. Zhou, Lei Li, Yumei Li, Wei Li, and Jingyi Jessica Li. PCA outperforms popular hidden variable inference methods for molecular QTL mapping. Genome Biology, 23(1):210, 2022.
#'
#' @export





groupPredict<-function(dataResponse,dataPredictors,R2Type=c("unadjusted","adjusted")){
  # set.seed(1)
  # n<-100
  # dataResponse<-matrix(data=rnorm(n=n*3),nrow=n) #100*3.
  # dataPredictors<-matrix(data=rnorm(n=n*10),nrow=n) #100*10.
  # R2Type<-"unadjusted"

  R2s<-rep(NA,ncol(dataResponse))
  for(j in 1:ncol(dataResponse)){
    # j<-1
    mod<-lm(dataResponse[,j]~as.matrix(dataPredictors))
    # summary(mod)
    if(R2Type=="unadjusted"){
      R2s[j]<-summary(mod)$r.squared
    }else{
      R2s[j]<-summary(mod)$adj.r.squared
    }
  }

  return(R2s)
}

