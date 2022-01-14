




#Given dataResponse and dataPredictors, each a matrix,
#predict each variable in dataResponse with all the variables in dataPredictors.
#Return a vector of R2's.
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

