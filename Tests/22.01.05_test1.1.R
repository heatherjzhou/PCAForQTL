setwd("~/2020.07.13_PEER/PCAForQTL")

#Clean and rebuild.

#Test runElbow().
if(TRUE){
  #Scenario 1 (doesn't work):
  X<-NULL
  prcompResult<-NULL
  resultRunElbow<-runElbow(X=X,prcompResult=prcompResult) #Error.

  #Scenario 2 (only provide prcompResult; recommended approach because it's fast)
  X<-NULL
  prcompResult<-readRDS("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2020.11.10_geneExpressionPCA/Colon_Transverse_prcompResult.rds")
  resultRunElbow<-runElbow(X=X,prcompResult=prcompResult) #12.

  #Scenario 3 (only provide X; slower):
  dataGeneExpressionFP<-readRDS("./Example_Data/Colon_Transverse.v8.normalized_expression.rds") #25,379*372. The first four columns are chr, start, end, and gene_id.
  X<-t(dataGeneExpressionFP[,-(1:4)]) #368*25,379.
  prcompResult<-NULL
  resultRunElbow<-runElbow(X=X,prcompResult=prcompResult) #12.

  #Scenario 4 (both are provided, but only prcompResult will be used):
  dataGeneExpressionFP<-readRDS("./Example_Data/Colon_Transverse.v8.normalized_expression.rds") #25,379*372. The first four columns are chr, start, end, and gene_id.
  X<-t(dataGeneExpressionFP[,-(1:4)]) #368*25,379.
  prcompResult<-readRDS("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2020.11.10_geneExpressionPCA/Colon_Transverse_prcompResult.rds")
  resultRunElbow<-runElbow(X=X,prcompResult=prcompResult) #12.
}


#Test runBE().
if(TRUE){
  dataGeneExpressionFP<-readRDS("./Example_Data/Colon_Transverse.v8.normalized_expression.rds") #25,379*372. The first four columns are chr, start, end, and gene_id.
  expr<-t(dataGeneExpressionFP[,-(1:4)]) #368*25,379.

  #Linux and Mac:
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)
  resultRunBE<-runBE(X=expr,B=20,alpha=0.05,
                     verbose=TRUE) #29.

  #Windows:
  set.seed(1)
  resultRunBE<-runBE(X=expr,B=20,alpha=0.05,
                     mc.cores=1,verbose=TRUE) #29.
}

#Is RNGkind("L'Ecuyer-CMRG") necessary when mc.cores=1? No.
if(FALSE){
  RNGkind("default")

  set.seed(1)
  result1<-parallel::mclapply(1:3,FUN=function(i){
    return(rnorm(n=1))
  },mc.cores=1)

  set.seed(1)
  result2<-parallel::mclapply(1:3,FUN=function(i){
    return(rnorm(n=1))
  },mc.cores=1)

  identical(result1,result2) #TRUE is good.
}


#Test makeScreePlot().
if(TRUE){
  prcompResult<-readRDS("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2020.11.10_geneExpressionPCA/Colon_Transverse_prcompResult.rds")

  makeScreePlot(prcompResult,labels=c("Elbow","BE","GTEx"),values=c(12,29,60),
                titleText="Colon - Transverse")

  makeScreePlot(prcompResult,labels=c("Elbow","BE","GTEx"),values=c(12,29,60),
                titleText="Colon - Transverse",subtitleText="BE - Elbow = 29 - 12 = 17")
}

#Test filterKnownCovariates().
if(TRUE){
  dataCovariates<-readRDS("./Example_Data/Colon_Transverse.v8.covariates.rds") #368*68. The columns are 5 genotype PCs, 60 PEER factors, pcr, platform, and sex.
  knownCovariates<-dataCovariates[,c(1:5,66:68)] #368*8.

  prcompResult<-readRDS("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2020.11.10_geneExpressionPCA/Colon_Transverse_prcompResult.rds")
  PCs<-prcompResult$x #368*368. The columns are PCs.
  PCsTop<-PCs[,1:12]

  knownCovariatesFiltered<-filterKnownCovariates(knownCovariates,PCsTop,unadjustedR2_cutoff=0.9,
                                                 verbose=TRUE)
}




