#' The Buja and Eyuboglu (BE) algorithm for choosing K in PCA
#'
#' This function implements the BE algorithm, a permutation-based approach for choosing K in PCA. Intuitively, the BE algorithm retains PCs that explain more variance in the data than by random chance and discards those that do not.
#'
#' The BE algorithm was originally proposed by Buja and Eyuboglu (1992). A detailed description of this method can be found in preprint. For reproducibility, make sure to change the RNG type (necessary unless \code{mc.cores} is 1) and set the seed before using this function (see README file for details).
#'
#' @param X The data matrix for running PCA (must be observation by feature).
#' @param B The number of permutations.
#' @param alpha The significance level.
#' @param mc.cores The number of cores used for parallel computing.
#' @param verbose Logical. If \code{TRUE}, then some detailed messages are printed when the function is run.
#'
#' @return This function returns a list consisting of the following items:
#' \item{pValues}{The p-values corresponding to the PCs.}
#' \item{alpha}{The significance level as inputted.}
#' \item{numOfPCsChosen}{The number of PCs chosen via BE.}
#'
#' @references
#' Andreas Buja and Nermin Eyuboglu. 1992. “Remarks on Parallel Analysis.” Multivariate Behavioral Research 27 (4): 509–40.
#'
#' Preprint.
#'
#' @export





#Given X, the data matrix (must be observation by feature),
#B, the number of permutations (default is 20),
#and alpha, the significance level (default is 0.05),
#run BE and return p-values, alpha, and numOfPCsChosen.
#For reproducibility, make sure to change the RNG type (necessary unless mc.cores is 1) and set the seed before using this function.
runBE<-function(X,B=20,alpha=0.05,
                mc.cores=min(B,parallel::detectCores()-1),
                verbose=FALSE){

  if(alpha<0 || alpha>1){
    stop("alpha must be between 0 and 1.")
  }

  n<-nrow(X) #Number of observations.
  p<-ncol(X) #Number of features.
  d<-min(n,p) #This is the total number of PCs.

  if(verbose) cat("Running PCA on permuted data...\n")

  # testStatsPerm<-matrix(data=NA,nrow=d,ncol=B) #PC by permutation.
  # for(b in 1:B){
  #   # b<-3
  #   if(verbose) cat("b=",b," out of ",B," permutations...\n",sep="")
  #
  #   #Permute each column of X. That is, permute the observations in each feature.
  #   XPermuted<-matrix(data=NA,nrow=n,ncol=p)
  #   for(j in 1:p){
  #     # j<-7
  #     XPermuted[,j]<-sample(x=X[,j],size=n,replace=FALSE)
  #   }
  #
  #   prcompResultPerm<-prcomp(x=XPermuted,center=TRUE,scale.=TRUE) #Key step.
  #   importanceTablePerm<-summary(prcompResultPerm)$importance
  #   testStatsPerm[,b]<-importanceTablePerm[2,] #The second row is PVE.
  # }

  results<-parallel::mclapply(1:B,FUN=function(b){
    # b<-3
    if(verbose) cat("b=",b," out of ",B," permutations...\n",sep="")

    #Permute each column of X. That is, permute the observations in each feature.
    XPermuted<-matrix(data=NA,nrow=n,ncol=p)
    for(j in 1:p){
      # j<-7
      XPermuted[,j]<-sample(x=X[,j],size=n,replace=FALSE)
    }

    prcompResultPerm<-prcomp(x=XPermuted,center=TRUE,scale.=TRUE) #Key step.
    importanceTablePerm<-summary(prcompResultPerm)$importance
    PVEsPerm<-importanceTablePerm[2,]
    return(PVEsPerm)
  },mc.cores=mc.cores) #results is a list of vectors.
  temp<-unlist(results)
  testStatsPerm<-matrix(data=temp,nrow=d,byrow=FALSE) #PC by permutation.

  if(verbose) cat("Running PCA on the unpermuted data...\n")
  prcompResult<-prcomp(x=X,center=TRUE,scale.=TRUE) #Key step.
  importanceTable<-summary(prcompResult)$importance
  PVEs<-importanceTable[2,]
  # Compare PVEs to testStatsPerm.
  # temp<-(testStatsPerm>=PVEs) #temp is calculated as desired.
  pValues<-(rowSums(testStatsPerm>=PVEs)+1)/(B+1) #The p-value for the jth PC is calculated as, roughly speaking, the proportion of permutations where the PVE of the jth PC is greater than or equal to PVE_j.

  for(j in 2:d){ #Enforce monotone increase of the p-values.
    if(pValues[j]<pValues[j-1]){
      pValues[j]<-pValues[j-1]
    }
  }

  numOfPCsChosen<-sum(pValues<=alpha)
  toReturn<-list(pValues=pValues,alpha=alpha,numOfPCsChosen=numOfPCsChosen)
  return(toReturn)
}
















