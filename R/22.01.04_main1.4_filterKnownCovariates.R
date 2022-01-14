#' Filter known covariates
#'
#' This function filters out the known covariates that are captured well by the inferred covariates (unadjusted R squared greater than or equal to 0.9 by default). We use unadjusted R squared instead of adjusted R squared because we do not want to penalize for model complexity here.
#'
#' @param knownCovariates The known covariate matrix (must be observation by feature).
#' @param inferredCovariates The inferred covariate matrix, e.g., the top PCs (must be observation by feature as well). The rows (i.e., observations) in knownCovariates and inferredCovariates must match.
#' @param unadjustedR2_cutoff The cutoff value for unadjusted R squared.
#' @param verbose Logical. If \code{TRUE}, then some detailed messages are printed when the function is run.
#'
#' @return This function returns the filtered known covariates, i.e., the known covariates that should be kept.
#'
#' @references
#' Preprint.
#'
#' @export





filterKnownCovariates<-function(knownCovariates,inferredCovariates,unadjustedR2_cutoff=0.9,
                                verbose=FALSE){

  if(nrow(knownCovariates)!=nrow(inferredCovariates)){ #Check whether the sample sizes are equal in case the rownames are empty.
    stop("knownCovariates and inferredCovariates must have the same number of rows (i.e., observations).")
  }
  if(!identical(rownames(knownCovariates),rownames(inferredCovariates))){ #identical(NULL,NULL) returns TRUE.
    stop("The rownames of knownCovariates and inferredCovariates must match.")
  }
  if(unadjustedR2_cutoff<0 || unadjustedR2_cutoff>1){
    stop("unadjustedR2_cutoff must be between 0 and 1.")
  }

  R2s<-groupPredict(dataResponse=knownCovariates,dataPredictors=inferredCovariates,R2Type="unadjusted") #Use unadjusted because we don't want to penalize for model complexity here.
  indicesOfKnownCovariatesToKeep<-which(R2s<unadjustedR2_cutoff) #This may be integer(0).
  if(verbose){
    cat(ncol(knownCovariates)-length(indicesOfKnownCovariatesToKeep)," out of the ",ncol(knownCovariates)," known covariates has/have been filtered out.\n")
  }

  toReturn<-knownCovariates[,indicesOfKnownCovariatesToKeep] #This syntax works even when all of the known covariates are filtered out.
  return(toReturn)
}






