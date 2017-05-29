#' Funnel Plots for catmap
#'
#' \code{catmap.funnel} creates a funnel plot of the individual ORs against the
#' se and saves the file to the current working directory (use getwd() to
#' view).  The plots are not created in the R Graphics device.
#'
#' \code{catmap.funnel} creates a .pdf file of the funnel plot. Plots are not
#' created in the R Graphics device window, but are instead saved to a .pdf
#' file in the current working directory.
#'
#' @param catmapobject The catmap object created by a previous call to catmap
#' @param funnel Logical.  Should a funnel plot be produced?  Funnel plots plot
#' standard error of the log ORs against the ORs along with a solid line at 1.0
#' and a dotted line at the overall pooled OR.  Used to assess publication
#' bias.  Output plot file is saved with the default name of
#' \bold{dataset.funnel.plot.pdf} where dataset is the name of the file given
#' as the first argument to catmap.
#' @author Kristin K. Nicodemus, \email{kristin.nicodemus@@well.ox.ac.uk}
#' @seealso \code{\link{catmap}}, \code{\link{catmap.sense}},
#' \code{\link{catmap.cumulative}}, \code{\link{catmap.forest}}.
#' @keywords methods
#' @examples
#' data(catmapdata)
#' catmapobject <- catmap(catmapdata, 0.95, TRUE)
#' catmap.funnel(catmapobject, TRUE)
#' @export
catmap.funnel<-function(catmapobject, funnel){
  if(funnel==TRUE){
    if(catmapobject$dataset!=catmapdata){
      pdf(file=(paste(catmapobject$dataset, "funnel.plot.pdf", sep=".")))
    }
    if(catmapobject$dataset==catmapdata){
      pdf(file=(paste("catmapdata", "funnel.plot.pdf", sep=".")))
    }
    OddsRatio<-exp(catmapobject$logOR)
    fun.se<-sqrt(catmapobject$comvarlogOR)
    maxse<-max(fun.se)
    minse<-min(fun.se)
    plot(OddsRatio, fun.se, ylim=rev(c(minse,maxse)),log="x", ylab="Standard Error Log OR")
    abline(v=1)
    abline(v=catmapobject$combinedOR, lty=2)
    graphics.off()
  }
}
