#' The elbow method for choosing K in PCA (based on distance to the diagonal line)
#'
#' This function implements an automatic elbow detection method (based on distance to the diagonal line).
#'
#' A detailed description of this method can be found in our 2022 paper.
#'
#' @param X The data matrix for running PCA (must be observation by feature).
#' @param prcompResult The output from running \code{prcomp()} on \code{X}. Either \code{X} or \code{prcompResult} must be provided. If both are provided, then only \code{prcompResult} will be used. We recommend only providing \code{prcompResult} because it's faster.
#'
#' @return This function returns the number of PCs selected by maximizing the distance to the diagonal line.
#'
#' @references
#' Heather J. Zhou, Lei Li, Yumei Li, Wei Li, and Jingyi Jessica Li. PCA outperforms popular hidden variable inference methods for molecular QTL mapping. Genome Biology, 23(1):210, 2022.
#'
#' @export





#Given X, the data matrix (must be observation by feature),
#and/or prcompResult, the output from running prcomp() on X,
#select the number of PCs by maximizing the distance to the diagonal line (see our 2022 paper or below for details).
#Either X or prcompResult must be provided. If both are provided, then only prcompResult will be used.
#We recommend only providing prcompResult because it's faster.
runElbow<-function(X=NULL,
                   prcompResult=NULL){

  #Obtain prcompResult.
  if(is.null(prcompResult)){
    if(is.null(X)){
      stop("Please input X or prcompResult.")
    }else{
      cat("Running PCA...\n")
      prcompResult<-prcomp(X,center=TRUE,scale.=TRUE) #This may take a moment.
    }
  }else{
    if(class(prcompResult)!="prcomp"){
      stop("prcompResult must be a prcomp object returned by the function prcomp().")
    }
  }

  importanceTable<-summary(prcompResult)$importance
  x<-1:ncol(importanceTable) #PC indices.
  y<-importanceTable[2,] #PVEs.

  #Given x  and y, calculate the distance between each point and the diagonal line (the line connecting the first and last points).
  x1<-x[1] #First point.
  y1<-y[1] #First point.

  x2<-x[length(x)] #Last point.
  y2<-y[length(y)] #Last point.

  x0<-x
  y0<-y

  distancesDenominator<-sqrt((x2-x1)^2+(y2-y1)^2)
  distancesNumerator<-abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))
  distances<-distancesNumerator/distancesDenominator
  # plot(distances)

  numOfPCsChosen<-which.max(distances) #12.
  names(numOfPCsChosen)<-NULL
  return(numOfPCsChosen)
}











