#' Get residuals
#'
#' Given dataResponse and dataPredictors, each a matrix, this function residualizes each variable in dataResponse against all the variables in dataPredictors using a multiple linear regression.
#'
#' @param dataResponse The matrix containing the response variables (must be observation by feature).
#' @param dataPredictors The matrix containing the predictors (must be observation by feature as well). The rows (i.e., observations) in dataResponse and dataPredictors must match.
#'
#' @return This function returns the residual matrix, which has the same dimensions as dataResponse.
#'
#' @references
#' Heather J. Zhou, Lei Li, Yumei Li, Wei Li, and Jingyi Jessica Li. PCA outperforms popular hidden variable inference methods for molecular QTL mapping. Genome Biology, 23(1):210, 2022.
#'
#' @export





getResiduals<-function(dataResponse,dataPredictors){
  # set.seed(1)
  # n<-100
  # dataResponse<-matrix(data=rnorm(n=n*20000),nrow=n) #100*20,000.
  # dataPredictors<-matrix(data=rnorm(n=n*10),nrow=n) #100*10.

  vecOfOnes<-rep(1,nrow(dataResponse)) #Vector of length 100.
  X<-cbind(vecOfOnes,as.matrix(dataPredictors)) #100*11.
  betaEstimate<-solve(t(X)%*%X)%*%t(X)%*%as.matrix(dataResponse) #11*20,000
  # dim(betaEstimate)

  toReturn<-as.matrix(dataResponse)-X%*%betaEstimate #100*20,000. Residuals.
  # dim(toReturn)

  #Check with lm(). Success.
  # mod<-lm(dataResponse[,1]~X)
  # temp<-mod$residuals

  return(toReturn)
}

