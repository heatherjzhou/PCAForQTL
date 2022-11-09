#' Make scree plot for PCA
#'
#' This function makes a scree plot showing the numbers of PCs chosen via different methods.
#'
#' @param prcompResult The output from \code{prcomp()}.
#' @param labels A string or vector of strings to be used as the label(s) of the numbers of PCs in \code{values}.
#' @param values An integer or vector of integers representing the numbers of PCs chosen via different methods. \code{labels} and \code{values} must be nonempty, and the entries must match between the two vectors.
#' @param titleText Title text.
#' @param subtitleText Subtitle text.
#' @param maxNumOfPCsToPlot The maximum number of PCs to plot in the scree plot.
#' @param colors A string or vector of strings representing the colors to plot the numbers of PCs in \code{values} in.
#'
#' @return This function returns the produced scree plot.
#'
#' @references
#' Heather J. Zhou, Lei Li, Yumei Li, Wei Li, and Jingyi Jessica Li. PCA outperforms popular hidden variable inference methods for molecular QTL mapping. Genome Biology, 23(1):210, 2022.
#'
#' @import ggplot2
#' @export





makeScreePlot<-function(prcompResult,labels,values,
                        titleText=NULL,subtitleText=NULL,
                        maxNumOfPCsToPlot=max(100,values),colors=NULL){

  if(class(prcompResult)!="prcomp") stop("prcompResult must be a prcomp object returned by the function prcomp().")
  if(!is.character(labels)) stop("labels must be a string or vector of strings.")
  if(!is.numeric(values)) stop("values must be an integer or vector of integers.")
  if(length(labels)!=length(values)) stop("labels and values must have the same length.")

  importanceTable<-summary(prcompResult)$importance
  dataPlot<-as.data.frame(t(importanceTable))
  dataPlot$PCIndex<-1:nrow(dataPlot)
  colnames(dataPlot)<-c("sd","PVE","cumulativePVE","PCIndex")
  # dataPlot<-dataPlot%>%filter(PCIndex<=maxNumOfPCsToPlot)
  dataPlot<-dataPlot[which(dataPlot$PCIndex<=maxNumOfPCsToPlot),] #Avoid using dplyr.

  dataPlotVline<-data.frame(matrix(nrow=length(labels),ncol=2))
  colnames(dataPlotVline)<-c("Category","Value")
  dataPlotVline$Category<-labels
  dataPlotVline$Category<-factor(dataPlotVline$Category,levels=dataPlotVline$Category)
  dataPlotVline$Value<-values
  if(is.null(colors)){
    colors<-RColorBrewer::brewer.pal(n=8,name="Accent")[1:length(labels)]
  }

  p<-ggplot(data=dataPlot,aes(x=PCIndex,y=PVE))+
    geom_point()+
    geom_vline(data=dataPlotVline, #Key.
               aes(xintercept=Value,color=Category),
               size=0.75)+ #Default size is 0.5.
    scale_x_continuous(breaks=seq(0,nrow(dataPlot),by=10))+
    scale_y_continuous(labels=scales::number_format(accuracy=0.01))+
    scale_color_manual(values=colors)+
    labs(x="PC index")+
    theme(
      text=element_text(size=15),

      plot.title=element_text(hjust=0.5), #Plot title.
      plot.subtitle=element_text(hjust=0.5), #Plot subtitle.

      panel.background=element_blank(), #Panel background.
      panel.grid.major.y=element_blank(), #Horizontal grid lines.
      panel.grid.major.x=element_blank(), #Vertical grid lines.

      axis.line=element_line(), #Axis lines (default size is 0.5),

      legend.title=element_blank(), #Legend title.
      legend.key=element_blank(), #Legend key background.
      legend.position="bottom", #Legend position.
    )

  if(!is.null(titleText)) p<-p+labs(title=titleText)
  if(!is.null(subtitleText)) p<-p+labs(subtitle=subtitleText)

  # p
  return(p)
}
