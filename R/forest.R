#' Leave-One-Out Sensitivity Analyses and Plots using either Fixed- or
#' Random-Effects Estimates
#'
#' \code{catmap.sense} conducts leave-one-out sensitivity analyses and creates
#' plots of the ORs and confidence intervals using either fixed- or
#' random-effects analyses, which are saved to the current working directory
#' (use getwd() to view) but does not create the plots in the R Graphics
#' device.
#'
#' \code{catmap.sense} conducts leave-one-out sensitivity analyses and creates
#' .pdf files of plots of the ORs and CIs using either the fixed-effect or the
#' random-effect estimates. Plots are not created in the R Graphics device
#' window, but are instead saved to a .pdf file in the current working
#' directory. Likewise the .txt files of results will be saved to the current
#' working directory. To find the current working directory, use getwd()
#'
#' @param catmapobject The catmap object created by a previous call to catmap
#' @param fe.sense Logical.  Should a leave-one-out sensitivity analysis be
#' performed using fixed-effects estimates?  Automatic output result files are
#' saved with the default name of \bold{dataset.fixed.effects.sensitivity.txt},
#' where dataset is the name of the file given as the first argument to catmap.
#' Note that repeated runs of the same input file will be appended to the
#' default output file names.
#' @param re.sense Logical.  Should a leave-one-out sensitivity analysis be
#' performed using random-effects estimates?  Automatic output result files are
#' saved with the default name of
#' \bold{dataset.random.effects.sensitivity.txt}, where dataset is the name of
#' the file given as the first argument to catmap.  Note that repeated runs of
#' the same input file will be appended to the default output file names.
#' @param fe.senseplot Logical.  Should a .pdf plot of the ORs and CIs from the
#' sensitivity analysis using fixed-effects be output?  Can be TRUE only if
#' fe.sense=TRUE.  Output plot file is saved with the default name of
#' \bold{dataset.fixed.effects. sensitivity.plot.pdf} where dataset is the name
#' of the file given as the first argument to catmap.
#' @param re.senseplot Logical.  Should a .pdf plot of the ORs and CIs from the
#' sensitivity analysis using random-effects be output?  Can be TRUE only if
#' re.sense=TRUE.  Output plot file is saved with the default name of
#' \bold{dataset.random.effects. sensitivity.plot. pdf} where dataset is the
#' name of the file given as the first argument to catmap.
#' @author Kristin K. Nicodemus, \email{kristin.nicodemus@@well.ox.ac.uk}
#' @seealso \code{\link{catmap}}, \code{\link{catmap.forest}},
#' \code{\link{catmap.cumulative}}, \code{\link{catmap.funnel}}.
#' @keywords methods
#' @examples
#'
#' \dontrun{
#' data(catmapdata)
#' catmapobject1<-catmap(catmapdata, 0.95, TRUE)
#' catmap.sense(catmapobject1, TRUE, TRUE, TRUE, TRUE)}
#' @export
catmap.sense<-function(catmapobject, fe.sense, re.sense, fe.senseplot, re.senseplot){

  #fixed-effects sensitivity

  if (fe.sense==TRUE){
    sfplot<-matrix(0,nrow(catmapobject$a1),3)
    for(f in 1:nrow(catmapobject$a1)){
      sf.weight<-catmapobject$weight[-f]
      sf.logOR<-catmapobject$logOR[-f]

      sf.combinedLogOR<-((sum(sf.weight*sf.logOR))/sum(sf.weight))
      sf.combinedOR<-exp(sf.combinedLogOR)
      sf.combinedSeLogOR<-(sqrt(1/sum(sf.weight)))
      sf.combinedVarLogOR<-(1/sum(sf.weight))
      sf.combinedChisq<-(((sf.combinedLogOR-0)^2)/sf.combinedVarLogOR)
      sf.combinedValue<-pchisq(sf.combinedChisq, df=1)
      sf.combinedPvalue<-(1-sf.combinedValue)

      #get qnorm values
      alpha<-(1-((1-catmapobject$ci)/2))
      quantilenorm<-qnorm(alpha, 0, 1)

      sf.lbci<-exp(sf.combinedLogOR-(quantilenorm*sf.combinedSeLogOR))
      sf.ubci<-exp(sf.combinedLogOR+(quantilenorm*sf.combinedSeLogOR))
      sf.combinedCI<-c(sf.lbci, sf.ubci)
      sf.SeLogOR<-sqrt(sf.combinedVarLogOR)
      sf.lbci<-exp(sf.logOR-(quantilenorm*sf.SeLogOR))
      sf.ubci<-exp(sf.logOR+(quantilenorm*sf.SeLogOR))

      #calculate heterogeneity
      sf.chisqHet<-(sum(sf.weight*(((sf.logOR-sf.combinedLogOR)^2))))
      sf.df<-(nrow(catmapobject$a1)-2)
      sf.combinedHetValue<-pchisq(sf.chisqHet, df=sf.df)
      sf.heterogeneityPvalue<-(1-sf.combinedHetValue)

      study.removed<-paste("Study Removed =", catmapobject$studyname[f], sep=" ")
      sftable.header<-c("Inverse Variance Fixed-Effects OR", "Inverse Variance Fixed-Effects Lower Bound CI", "Inverse Variance Fixed-Effects Upper Bound CI", "Inverse Variance Fixed-Effects Chi-Square", "Inverse Variance Fixed-Effects p-value", "Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value")
      sftable.fill<-c(sf.combinedOR, sf.combinedCI, sf.combinedChisq, sf.combinedPvalue, sf.chisqHet, sf.heterogeneityPvalue)
      sf.results<-rbind(sftable.header, round(sftable.fill, digits=5))
      sfvalues<-c(sf.combinedOR, sf.combinedCI)
      sfplot[f,]<-sfvalues
      cat("Fixed Effects Sensitivity Analysis\n")
      cat(study.removed, sf.results, sep="\n")
      cat("\n")

      #print results to file dataset.fixed.effects.sensitivity.txt

      if(catmapobject$dataset!=catmapdata){
        sink(paste(catmapobject$dataset, "fixed.effects.sensitivity.txt", sep="."), append=TRUE)
        cat("Fixed Effects Sensitivity Analysis\n")
        cat(study.removed, sf.results, sep="\n")
        cat("\n")
        sink()
      }

      if(catmapobject$dataset==catmapdata){
        sink(paste("catmapdata", "fixed.effects.sensitivity.txt", sep="."), append=TRUE)
        cat("Fixed Effects Sensitivity Analysis\n")
        cat(study.removed, sf.results, sep="\n")
        cat("\n")
        sink()
      }
    }
  }

  #random-effects sensitivity

  if (re.sense==TRUE & catmapobject$tau2 <=0){
    cat("NOTICE: tau2 is less than or equal to 0;\n no random effects estimates will be calculated.\n")
  }

  if (re.sense==TRUE & catmapobject$tau2 > 0){
    srplot<-matrix(0,nrow(catmapobject$a1),3)
    for(r in 1:nrow(catmapobject$a1)){
      sr.weight<-catmapobject$weight[-r]
      sr.logOR<-catmapobject$logOR[-r]
      sr.comvarlogOR<-catmapobject$comvarlogOR[-r]

      sr.combinedLogOR<-((sum(sr.weight*sr.logOR))/sum(sr.weight))
      sr.combinedOR<-exp(sr.combinedLogOR)
      sr.combinedSeLogOR<-(sqrt(1/sum(sr.weight)))
      sr.combinedVarLogOR<-(1/sum(sr.weight))
      sr.combinedChisq<-(((sr.combinedLogOR-0)^2)/sr.combinedVarLogOR)
      sr.combinedValue<-pchisq(sr.combinedChisq, df=1)
      sr.combinedPvalue<-(1-sr.combinedValue)

      #get qnorm values
      alpha<-(1-((1-catmapobject$ci)/2))
      quantilenorm<-qnorm(alpha, 0, 1)

      sr.lbci<-exp(sr.combinedLogOR-(quantilenorm*sr.combinedSeLogOR))
      sr.ubci<-exp(sr.combinedLogOR+(quantilenorm*sr.combinedSeLogOR))
      sr.combinedCI<-c(sr.lbci, sr.ubci)
      sr.SeLogOR<-sqrt(sr.combinedVarLogOR)
      sr.lbci<-exp(sr.logOR-(quantilenorm*sr.SeLogOR))
      sr.ubci<-exp(sr.logOR+(quantilenorm*sr.SeLogOR))

      #calculate heterogeneity
      sr.df<-(nrow(catmapobject$a1)-2)
      sr.chisqHet<-(sum(sr.weight*(((sr.logOR-sr.combinedLogOR)^2))))
      sr.combinedHetValue<-pchisq(sr.chisqHet, df=sr.df)
      sr.heterogeneityPvalue<-(1-sr.combinedHetValue)

      #DerSimonian and Laird random-effects sensitivity analysis

      sr.tau2<-((sr.chisqHet-sr.df)/(sum(sr.weight)-(sum(sr.weight^2)/(sum(sr.weight)))))

      if (sr.tau2 <=0){
        sr.tau2<-0
      }
      srweight.dsl<-(1/(sr.comvarlogOR+sr.tau2))
      srlogOR.dsl<-((sum(srweight.dsl*sr.logOR))/(sum(srweight.dsl)))
      srOR.dsl<-exp(srlogOR.dsl)
      srseLogOR.dsl<-(1/(sqrt(sum(srweight.dsl))))
      srvarLogOR.dsl<-(1/sum(srweight.dsl))
      srlbci.dsl<-exp(srlogOR.dsl-(quantilenorm*srseLogOR.dsl))
      srubci.dsl<-exp(srlogOR.dsl+(quantilenorm*srseLogOR.dsl))
      srci.dsl<-c(srlbci.dsl, srubci.dsl)
      srchisq.dsl<-(((srlogOR.dsl-0)^2)/srvarLogOR.dsl)
      srvalue.dsl<-pchisq(srchisq.dsl, df=1)
      srpvalue.dsl<-(1-srvalue.dsl)

      srstudy.removed<-paste("Study Removed =", catmapobject$studyname[r], sep=" ")
      srtable.header<-c("Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value", "DerSimonian & Laird Random-Effects OR", "DerSimonian & Laird Random-Effects Lower Bound CI", "DerSimonian & Laird Random-Effects Upper Bound CI", "DerSimonian & Laird Random-Effects Chi-Square", "DerSimonian & Laird Random-Effects p-value")
      srtable.fill<-c(sr.chisqHet, sr.heterogeneityPvalue, srOR.dsl, srlbci.dsl, srubci.dsl, srchisq.dsl, srpvalue.dsl)
      sr.results<-rbind(srtable.header, round(srtable.fill, digits=5))
      srvalues<-c(srOR.dsl, srci.dsl)
      srplot[r,]<-srvalues
      cat("Random Effects Sensitivity Analysis\n")
      cat(srstudy.removed, sr.results, sep="\n")
      cat("\n")

      #print results to file dataset.random.effects.sensitivity.txt

      if(catmapobject$dataset!=catmapdata){
        sink(paste(catmapobject$dataset, "random.effects.sensitivity.txt", sep="."), append=TRUE)
        cat("Random Effects Sensitivity Analysis\n")
        cat(srstudy.removed, sr.results, sep="\n")
        cat("\n")
        sink()
      }

      if(catmapobject$dataset==catmapdata){
        sink(paste("catmapdata", "random.effects.sensitivity.txt", sep="."), append=TRUE)
        cat("Random Effects Sensitivity Analysis\n")
        cat(srstudy.removed, sr.results, sep="\n")
        cat("\n")
        sink()
      }
    }
  }

  ci100<- catmapobject$ci*100

  #create fixed-effects sensitivity plots
  if (fe.sense==FALSE & fe.senseplot==TRUE){
    cat("NOTICE: No fixed-efffects sensitivity analyses have been performed. \n Please set fe.sense==TRUE and try again.\n")
  }

  if (fe.senseplot==TRUE){
    if(catmapobject$dataset!=catmapdata){
      pdf(file=(paste(catmapobject$dataset, "fixed.effects.sensitivity.plot.pdf", sep=".")))
    }
    if(catmapobject$dataset==catmapdata){
      pdf(file=(paste("catmapdata", "fixed.effects.sensitivity.plot.pdf", sep=".")))
    }
    sfpstudy<-c(1:nrow(catmapobject$a1))
    sfplx<-max((min(sfplot[,2])-0.25),0)
    sfphx<-max(sfplot[,3])+0.25
    sfply<-min(sfpstudy)-0.5
    sfphy<-max(sfpstudy)+0.5
    mar1<-c(5.1, 7.1, 4.1, 2.1)
    las1<-1
    par("las"=las1)
    par("mar"=mar1)
    sfpdummy<-c(rep(NA, length(sfpstudy)))
    sfpdummy[1]<- sfphx
    sfpdummy[2]<- sfplx
    sfpydummy<-c(rep(NA, length(sfpstudy)))
    sfpydummy[1]<- sfphy
    sfpydummy[2]<- sfply
    xtitle<-paste("OR(", ci100, "% CI)")
    plot(sfpdummy, sfpydummy, type="n", log="x", ylab="", ylim=rev(c(sfply, sfphy)), yaxt="n", xlab=xtitle, main="Sensitivity Analysis: \n Inverse Variance (Fixed-Effects) ORs")
    abline(v=1.0)
    for(z in 1:nrow(catmapobject$a1)){
      points(sfplot[z,1],sfpstudy[z], pch=22, bg="black", col="black")
      segments(sfplot[z,2],sfpstudy[z],sfplot[z,3],sfpstudy[z], bg="black", col="black")
      sfpstudyname<-paste("Study Removed:", catmapobject$studyname[z], sep="\n")
      mtext(paste(sfpstudyname),side=2, at=sfpstudy[z], cex=0.80)
    }
    #dev.off()
  }

  ci100<- catmapobject$ci*100

  #create random-effects sensitivity plots
  if (re.sense==FALSE & re.senseplot==TRUE){
    cat("NOTICE: No random-efffects sensitivity analyses have been performed. \n Please set re.sense==TRUE and try again.\n")
  }
  if(fe.senseplot==TRUE){
    dev.off()
  }
  if (re.senseplot==TRUE & catmapobject$tau2 > 0){
    if(catmapobject$dataset!=catmapdata){
      pdf(file=(paste(catmapobject$dataset, "random.effects.sensitivity.plot.pdf", sep=".")))
    }
    if(catmapobject$dataset==catmapdata){
      pdf(file=(paste("catmapdata", "random.effects.sensitivity.plot.pdf", sep=".")))
    }
    srpstudy<-c(1:nrow(catmapobject$a1))
    srplx<-max((min(srplot[,2])-0.25),0)
    srphx<-max(srplot[,3])+0.25
    srply<-min(srpstudy)-0.5
    srphy<-max(srpstudy)+0.5
    mar1<-c(5.1, 7.1, 4.1, 2.1)
    las1<-1
    par("las"=las1)
    par("mar"=mar1)
    srpdummy<-c(rep(NA, length(srpstudy)))
    srpdummy[1]<- srphx
    srpdummy[2]<- srplx
    srpydummy<-c(rep(NA, length(srpstudy)))
    srpydummy[1]<- srphy
    srpydummy[2]<- srply
    xtitle<-paste("OR(", ci100, "% CI)")
    plot(srpdummy, srpydummy, type="n", log="x", ylab="", ylim=rev(c(srply, srphy)), yaxt="n", xlab=xtitle, main="Sensitivity Analysis: \n DerSimonian & Laird (Random-Effects) ORs")
    abline(v=1.0)
    for(y in 1:nrow(catmapobject$a1)){
      points(srplot[y,1],srpstudy[y], pch=22, bg="black", col="black")
      segments(srplot[y,2],srpstudy[y],srplot[y,3],srpstudy[y], bg="black", col="black")
      srpstudyname<-paste("Study Removed:", catmapobject$studyname[y], sep="\n")
      mtext(paste(srpstudyname),side=2, at=srpstudy[y], cex=0.80)
    }
    graphics.off()
  }
}

#' Cumulative Meta-Analyses and Plots using either Fixed- or Random-Effects
#'
#' \code{catmap.cumulative} conducts cumulative meta-analyses and creates plots
#' of the ORs and confidence intervals using either fixed- or random-effects
#' analyses, and saves text files and plot files to the current working
#' directory (use getwd() to obtain the current working directory). The plots
#' are not created in the R Graphics device.
#'
#' \code{catmap.cumulative} conducts cumulative meta-analyses and creates .pdf
#' files of plots of the ORs and CIs using either the fixed-effect or the
#' random-effect estimates.  \bold{NOTE: The studies should be listed in
#' chronological order in the input file.  \code{catmap.cumulative} does not
#' re-order studies by publication year.}
#'
#' @param catmapobject The catmap object created by a previous call to catmap
#' @param fe.cumulative Logical.  Should a cumulative meta-analysis be
#' performed using fixed-effects estimates?  catmap assumes the order in which
#' the studies are listed is the chronological ordering.  Automatic output
#' result file is saved with the default name of
#' \bold{dataset.fixed.effects.cumulative .txt}, where dataset is the name of
#' the file given as the first argument to catmap.  Note that repeated runs of
#' the same input file will be appended to the default output file.
#' @param re.cumulative Logical.  Should a cumulative meta-analysis be
#' performed using random-effects estimates?  catmap assumes the order in which
#' the studies are listed is the chronological ordering.  Automatic output
#' result file is saved with the default name of
#' \bold{dataset.random.effects.cumulative .txt}, where dataset is the name of
#' the file given as the first argument to catmap.  Note that repeated runs of
#' the same input file will be appended to the default output file.  Also note
#' that random-effects estimates are undefined for a single study.  Therefore,
#' the first study will have 0s and NaNs in the output, although the resulting
#' plot shows the fixed-effect estimate for the first study.  The OR and CI for
#' the first study may be found using the fixed effects estimates.
#' @param fe.cumplot Logical.  Should a .pdf plot of the ORs and CIs from the
#' cumulative meta-analysis using fixed-effects be output?  Can be TRUE only if
#' fe.cumulative =TRUE.  Output plot file is saved with the default name of
#' \bold{dataset.fixed. effects.cumulative.plot.pdf} where dataset is the name
#' of the file given as the first argument to catmap.
#' @param re.cumplot Logical.  Should a .pdf plot of the ORs and CIs from the
#' cumulative meta-analysis using random-effects be output?  Can be TRUE only
#' if re.cumulative= TRUE.  Output plot file is saved with the default name of
#' \bold{dataset.random. effects .cumulative.plot.pdf} where dataset is the
#' name of the file given as the first argument to catmap.
#' @author Kristin K. Nicodemus, \email{kristin.nicodemus@@well.ox.ac.uk}
#' @seealso \code{\link{catmap}}, \code{\link{catmap.forest}},
#' \code{\link{catmap.sense}}, \code{\link{catmap.funnel}}.
#' @keywords methods
#' @examples
#'
#' \dontrun{
#' data(catmapdata)
#' catmapobject1<-catmap(catmapdata, 0.95, TRUE)
#' catmap.cumulative(catmapobject1, TRUE, TRUE, TRUE, TRUE)}
#' @export
catmap.cumulative<-function(catmapobject, fe.cumulative, re.cumulative, fe.cumplot, re.cumplot){

  #fixed-effect cumulative meta analyses
  ci100<- catmapobject$ci*100

  if (fe.cumulative==TRUE){

    cfplot<-matrix(0, nrow(catmapobject$a1),3)
    for(c in 1:nrow(catmapobject$a1)){
      cf.weight<- catmapobject$weight[1:c]
      cf.logOR<- catmapobject$logOR[1:c]

      cf.combinedLogOR<-((sum(cf.weight*cf.logOR))/sum(cf.weight))
      cf.combinedOR<-exp(cf.combinedLogOR)
      cf.combinedSeLogOR<-(sqrt(1/sum(cf.weight)))
      cf.combinedVarLogOR<-(1/sum(cf.weight))
      cf.combinedChisq<-(((cf.combinedLogOR-0)^2)/cf.combinedVarLogOR)
      cf.combinedValue<-pchisq(cf.combinedChisq, df=1)
      cf.combinedPvalue<-(1-cf.combinedValue)

      #get qnorm values
      alpha<-(1-((1-catmapobject$ci)/2))
      quantilenorm<-qnorm(alpha, 0, 1)

      cf.lbci<-exp(cf.combinedLogOR-(quantilenorm*cf.combinedSeLogOR))
      cf.ubci<-exp(cf.combinedLogOR+(quantilenorm*cf.combinedSeLogOR))
      cf.combinedCI<-c(cf.lbci, cf.ubci)
      cf.SeLogOR<-sqrt(cf.combinedVarLogOR)
      cf.lbci<-exp(cf.logOR-(quantilenorm*cf.SeLogOR))
      cf.ubci<-exp(cf.logOR+(quantilenorm*cf.SeLogOR))

      #calculate heterogeneity
      cf.df<-(c-1)
      cf.chisqHet<-(sum(cf.weight*(((cf.logOR-cf.combinedLogOR)^2))))
      cf.combinedHetValue<-pchisq(cf.chisqHet, df=cf.df)
      cf.heterogeneityPvalue<-(1-cf.combinedHetValue)

      study.added<-paste("Study Added =", catmapobject$studyname[c], sep=" ")
      cftable.header<-c("Inverse Variance Fixed-Effects OR", "Inverse Variance Fixed-Effects Lower Bound CI", "Inverse Variance Fixed-Effects Upper Bound CI", "Inverse Variance Fixed-Effects Chi-Square", "Inverse Variance Fixed-Effects p-value", "Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value")
      cftable.fill<-c(cf.combinedOR, cf.combinedCI, cf.combinedChisq, cf.combinedPvalue, cf.chisqHet, cf.heterogeneityPvalue)
      cf.results<-rbind(cftable.header, round(cftable.fill, digits=5))
      cfvalues<-c(cf.combinedOR, cf.combinedCI)
      cfplot[c,]<-cfvalues
      cat("Fixed Effects Cumulative Meta-Analysis\n")
      cat(study.added, cf.results, sep="\n")
      cat("\n")

      #print results to file dataset.fixed.effects.cumulative.txt

      if(catmapobject$dataset!=catmapdata){
        sink(paste(catmapobject$dataset, "fixed.effects.cumulative.txt", sep="."), append=TRUE)
        cat("Fixed Effects Cumulative Meta-Analysis\n")
        cat(study.added, cf.results, sep="\n")
        cat("\n")
        sink()
      }

      if(catmapobject$dataset==catmapdata){
        sink(paste("catmapdata", "fixed.effects.cumulative.txt", sep="."), append=TRUE)
        cat("Fixed Effects Cumulative Meta-Analysis\n")
        cat(study.added, cf.results, sep="\n")
        cat("\n")
        sink()
      }
    }
  }

  #random-effect cumulative meta analyses

  if(re.cumulative==TRUE & catmapobject$tau2 <= 0){
    cat("NOTICE: tau2 is less than or equal to 0;\n no random effects estimates will be calculated.\n")
  }

  if (re.cumulative==TRUE & catmapobject$tau2 > 0){
    crplot<-matrix(0,nrow(catmapobject$a1),3)
    for(u in 1:nrow(catmapobject$a1)){
      #v<-u-1
      cr.weight<- catmapobject$weight[1:u]
      cr.logOR<- catmapobject$logOR[1:u]
      cr.comvarlogOR<- catmapobject$comvarlogOR[1:u]

      cr.combinedLogOR<-((sum(cr.weight*cr.logOR))/sum(cr.weight))
      cr.combinedOR<-exp(cr.combinedLogOR)
      cr.combinedSeLogOR<-(sqrt(1/sum(cr.weight)))
      cr.combinedVarLogOR<-(1/sum(cr.weight))
      cr.combinedChisq<-(((cr.combinedLogOR-0)^2)/cr.combinedVarLogOR)
      cr.combinedValue<-pchisq(cr.combinedChisq, df=1)
      cr.combinedPvalue<-(1-cr.combinedValue)

      #get qnorm values
      alpha<-(1-((1- catmapobject$ci)/2))
      quantilenorm<-qnorm(alpha, 0, 1)

      cr.lbci<-exp(cr.combinedLogOR-(quantilenorm*cr.combinedSeLogOR))
      cr.ubci<-exp(cr.combinedLogOR+(quantilenorm*cr.combinedSeLogOR))
      cr.combinedCI<-c(cr.lbci, cr.ubci)
      cr.SeLogOR<-sqrt(cr.combinedVarLogOR)

      #calculate heterogeneity
      cr.df<-(u-1)
      cr.chisqHet<-(sum(cr.weight*(((cr.logOR-cr.combinedLogOR)^2))))
      cr.combinedHetValue<-pchisq(cr.chisqHet, df=cr.df)
      cr.heterogeneityPvalue<-(1-cr.combinedHetValue)

      #DerSimonian and Laird random-effects cumulative analysis

      cr.tau2c<-(cr.chisqHet-cr.df)/(sum(cr.weight)-(sum(cr.weight^2)/(sum(cr.weight))))

      cr.tau2<-max(0,cr.tau2c)

      #if (cr.tau2 <=0){
      #cr.tau2<-0
      #}

      crweight.dsl<-(1/(cr.comvarlogOR+cr.tau2))
      crlogOR.dsl<-((sum(crweight.dsl*cr.logOR))/(sum(crweight.dsl)))
      crOR.dsl<-exp(crlogOR.dsl)
      crseLogOR.dsl<-(1/(sqrt(sum(crweight.dsl))))
      crvarLogOR.dsl<-(1/sum(crweight.dsl))
      crlbci.dsl<-exp(crlogOR.dsl-(quantilenorm*crseLogOR.dsl))
      crubci.dsl<-exp(crlogOR.dsl+(quantilenorm*crseLogOR.dsl))
      crci.dsl<-c(crlbci.dsl, crubci.dsl)
      crchisq.dsl<-(((crlogOR.dsl-0)^2)/crvarLogOR.dsl)
      crvalue.dsl<-pchisq(crchisq.dsl, df=1)
      crpvalue.dsl<-(1-crvalue.dsl)
      #added } here

      #crOR.dsl[1]<-exp(catmapobject$logOR[1])
      #crlbci.dsl[1]<-catmapobject$lbci.fe[1]
      #crubci.dsl[1]<-catmapobject$ubci.fe[1]

      crstudy.added<-paste("Study Added =", catmapobject$studyname[u], sep=" ")
      crtable.header<-c("Q Statistic (Heterogeneity) Chi-Square", "Q Statistic (Heterogeneity) p-value", "DerSimonian & Laird Random-Effects OR", "DerSimonian & Laird Random-Effects Lower Bound CI", "DerSimonian & Laird Random-Effects Upper Bound CI", "DerSimonian & Laird Random-Effects Chi-Square", "DerSimonian & Laird Random-Effects p-value")
      crtable.fill<-c(cr.chisqHet, cr.heterogeneityPvalue, crOR.dsl, crlbci.dsl, crubci.dsl, crchisq.dsl, crpvalue.dsl)
      cr.results<-rbind(crtable.header, round(crtable.fill, digits=5))
      crvalues<-c(crOR.dsl, crci.dsl)
      crplot[u,]<-crvalues
      crplot[1,1]<-exp(catmapobject$logOR[1])
      crplot[1,2]<-catmapobject$lbci.fe[1]
      crplot[1,3]<-catmapobject$ubci.fe[1]

      cat("Random Effects Cumulative Meta-Analysis\n")
      cat("The first study will NOT have estimates - this is because the random\neffects estimates are not defined for a single study\n")

      cat(crstudy.added, cr.results, sep="\n")
      cat("\n")

      #print results to file dataset.random.effects.cumulative.txt

      if(catmapobject$dataset!=catmapdata){
        sink(paste(catmapobject$dataset, "random.effects.cumulative.txt", sep="."), append=TRUE)
        cat("Random Effects Cumulative Meta-Analysis\n")
        cat("The first study will NOT have estimates - this is because the random\neffects estimates are not defined for a single study\n")
        cat(crstudy.added, cr.results, sep="\n")
        cat("\n")
        sink()
      }

      if(catmapobject$dataset==catmapdata){
        sink(paste("catmapdata", "random.effects.cumulative.txt", sep="."), append=TRUE)
        cat("Random Effects Cumulative Meta-Analysis\n")
        cat("The first study will NOT have estimates - this is because the random\neffects estimates are not defined for a single study\n")
        cat(crstudy.added, cr.results, sep="\n")
        cat("\n")
        sink()
      }
    }
  }

  #create fixed-effects cumulative plots

  if(fe.cumulative==FALSE & fe.cumplot==TRUE){
    cat("NOTICE: No fixed-efffects cumulative analyses have been performed.\n Please set fe.cumulative==TRUE and try again.\n")
  }

  if (fe.cumplot==TRUE){
    if(catmapobject$dataset!=catmapdata){
      pdf(file=(paste(catmapobject$dataset, "fixed.effects.cumulative.plot.pdf", sep=".")))
    }
    if(catmapobject$dataset==catmapdata){
      pdf(file=(paste("catmapdata", "fixed.effects.cumulative.plot.pdf", sep=".")))
    }
    cfpstudy<-c(1:nrow(catmapobject$a1))
    cfplx<-max((min(cfplot[,2])-0.25),0.1)
    cfphx<-max(cfplot[,3])+0.25
    cfply<-1-0.5
    cfphy<-max(cfpstudy)+0.5
    mar1<-c(5.1, 7.1, 4.1, 2.1)
    las1<-1
    par("las"=las1)
    par("mar"=mar1)
    cfpdummy<-c(rep(NA, length(cfpstudy)))
    cfpdummy[1]<- cfphx
    cfpdummy[2]<- cfplx
    cfpydummy<-c(rep(NA, length(cfpstudy)))
    cfpydummy[1]<- cfphy
    cfpydummy[2]<- cfply
    xtitle<-paste("OR(", ci100, "% CI)")
    plot(cfpdummy, cfpydummy, type="n", log="x", ylab="", ylim=rev(c(cfply, cfphy)), yaxt="n", xlab=xtitle, main="Cumulative Meta-Analysis: \n Inverse Variance (Fixed-Effects) ORs")
    abline(v=1.0)
    for(x in 1:nrow(catmapobject$a1)){
      points(cfplot[x,1],cfpstudy[x], pch=22, bg="black", col="black")
      segments(cfplot[x,2],cfpstudy[x],cfplot[x,3],cfpstudy[x], bg="black", col="black")
      cfpstudyname<-paste("Study Added:", catmapobject$studyname[x], sep="\n")
      mtext(paste(cfpstudyname),side=2, at=cfpstudy[x], cex=0.80)
    }
    dev.off()
  }

  #create random-effects cumulative plots
  if(re.cumulative==FALSE & re.cumplot==TRUE){
    cat("NOTICE: No random-efffects cumulative analyses have been performed. \n Please set re.cumulative==TRUE and try again.\n")
  }
  if (re.cumplot==TRUE & catmapobject$tau2 > 0){
    if(catmapobject$dataset!=catmapdata){
      pdf(file=(paste(catmapobject$dataset, "random.effects.cumulative.plot.pdf", sep=".")))
    }
    if(catmapobject$dataset==catmapdata){
      pdf(file=(paste("catmapdata", "random.effects.cumulative.plot.pdf", sep=".")))
    }
    crpstudy<-c(1:nrow(catmapobject$a1))
    crplx<-max((min(crplot[,2])-0.25),0.1)
    crphx<-max(crplot[,3])+0.25
    #crply<-min(crpstudy)-0.5
    crply<-1-0.5
    crphy<-max(crpstudy)+0.5
    mar1<-c(5.1, 7.1, 4.1, 2.1)
    las1<-1
    par("las"=las1)
    par("mar"=mar1)
    crpdummy<-c(rep(NA, length(crpstudy)))
    crpdummy[1]<- crphx
    crpdummy[2]<- crplx
    crpydummy<-c(rep(NA, length(crpstudy)))
    crpydummy[1]<- crphy
    crpydummy[2]<- crply
    xtitle<-paste("OR (", ci100, "% CI)")
    plot(crpdummy, crpydummy, type="n", log="x", ylab="", ylim=rev(c(crply, crphy)), yaxt="n", xlab=xtitle, main="Cumulative Meta-Analysis: \n DerSimonian & Laird (Random-Effects) ORs")
    abline(v=1.0)
    for(w in 1:nrow(catmapobject$a1)){
      points(crplot[w,1],crpstudy[w], pch=22, bg="black", col="black")
      segments(crplot[w,2],crpstudy[w],crplot[w,3],crpstudy[w], bg="black", col="black")
      crpstudyname<-paste("Study Added:", catmapobject$studyname[w], sep="\n")
      mtext(paste(crpstudyname),side=2, at=crpstudy[w], cex=0.80)
    }
    graphics.off()
  }
}

#' Forest Plots using either Fixed- or Random-Effects Pooled ORs and CIs
#'
#' \code{catmap.forest} creates forest plots of the individual study ORs and
#' CIs and the fixed or random effects pooled OR and CI and saves them to the
#' current working directory (use getwd() to view cwd).  The plots are not
#' created in the R Graphics device.
#'
#' \code{catmap.forest} creates forest plots of individual study ORs and CIs
#' plus the pooled estimate of the fixed- or random-effects pooled OR and CI.
#' Plots are not created in the R Graphics device window, but are instead saved
#' to a .pdf file in the current working directory, which can be found using
#' getwd()
#'
#' @param catmapobject The catmap object created by a previous call to catmap
#' @param fe.forest Logical.  Should a forest plot be created using the
#' fixed-effects estimates?  Plots are saved with the default name of
#' \bold{dataset.fixed.effects.forest.pdf}, where dataset is the name of the
#' file given as the first argument to catmap.
#' @param re.forest Logical.  Should a forest plot be created using the
#' random-effects estimates?  Plots are saved with the default name of
#' \bold{dataset.random.effects. forest.pdf} where dataset is the name of the
#' file given as the first argument to catmap.
#' @author Kristin K. Nicodemus, \email{kristin.nicodemus@@well.ox.ac.uk}
#' @seealso \code{\link{catmap}}, \code{\link{catmap.cumulative}},
#' \code{\link{catmap.sense}}, \code{\link{catmap.funnel}}.
#' @keywords methods
#' @examples
#'
#' \dontrun{
#' data(catmapdata)
#' catmapobject1<-catmap(catmapdata, 0.95, TRUE)
#' catmap.cumulative(catmapobject1, TRUE, TRUE, TRUE, TRUE)}
#' @export
catmap.forest<-function(catmapobject, fe.forest, re.forest){
  ci100<- catmapobject$ci*100
  if (fe.forest==TRUE){
    if(catmapobject$dataset!=catmapdata){
      pdf(file=(paste(catmapobject$dataset, "fixed.effects.plot.pdf", sep=".")))
    }
    if(catmapobject$dataset==catmapdata){
      pdf(file=(paste("catmapdata", "fixed.effects.plot.pdf", sep=".")))
    }
    prop.weight<-(catmapobject$weight/(sum(catmapobject$weight))*10)
    OR<-c(catmapobject$combinedOR, exp(catmapobject$logOR))
    study<-1:(nrow(catmapobject$a1)+1)
    lx<-max((min(catmapobject$lbci.fe)-0.25),0.1)
    hx<-max(catmapobject$ubci.fe)+0.25
    ly<-(1-0.5)
    hy<-(length(catmapobject$study)+ 1.5)
    mar1<-c(5.1, 7.1, 4.1, 2.1)
    las1<-1
    par("las"=las1)
    par("mar"=mar1)
    dummy<-c(rep(NA, length(catmapobject$study)))
    dummy[1]<-hx
    dummy[2]<-lx
    ydummy<-c(rep(NA, length(catmapobject$study)))
    ydummy[1]<-hy
    ydummy[2]<-ly
    xtitle<-paste("OR(", ci100, "% CI)")
    plot(dummy, ydummy, type="n", log="x", ylab="", yaxt="n", xlab=xtitle, main="Inverse Variance (Fixed-Effects) ORs")
    points(OR[1], 1, pch=22, cex=1.0, bg="red", col="red")
    segments(catmapobject$lbci, 1, catmapobject$ubci, 1, col="red")
    mtext("Pooled OR",side=2,at=1, cex=0.80)
    abline(v=1.0)
    limit<-nrow(catmapobject$a1)+1
    counts<-c(1:limit)
    for(i in 2:(nrow(catmapobject$a1)+1)){
      j<-i-1
      points(OR[i], counts[i], pch=22, cex=prop.weight[j], bg="black", col="black")
      segments(catmapobject$lbci.fe[j], counts[i], catmapobject$ubci.fe[j], counts[i], bg="black", col="black")
      mtext(paste(catmapobject$studyname[j]),side=2, at=counts[i], cex=0.80)
    }
    dev.off()
  }

  #create random-effects forest plots

  if(re.forest==TRUE & catmapobject$tau2 <= 0){
    cat("NOTICE: tau2 is less than or equal to 0;\n no random effects estimates will be calculated.\n")
  }

  if (re.forest==TRUE & catmapobject$tau2 > 0){
    if(catmapobject$dataset!=catmapdata){
      pdf(file=(paste(catmapobject$dataset, "random.effects.plot.pdf", sep=".")))
    }
    if(catmapobject$dataset==catmapdata){
      pdf(file=(paste("catmapdata", "random.effects.plot.pdf", sep=".")))
    }
    dslprop.weight<-(catmapobject$weight.dsl/(sum(catmapobject$weight.dsl))*10)
    dslOR<-c(catmapobject$OR.dsl, exp(catmapobject$logOR))
    OR<-c(catmapobject$combinedOR, exp(catmapobject$logOR))
    rstudy<-1:(nrow(catmapobject$a1)+1)
    rlx<-max((min(catmapobject$lbci.fe)-0.25),0.1)
    rhx<-max(catmapobject$ubci.fe)+0.25
    rly<-(1-0.5)
    rhy<-(length(catmapobject$study)+1.5)
    mar1<-c(5.1, 7.1, 4.1, 2.1)
    las1<-1
    par("las"=las1)
    par("mar"=mar1)
    rdummy<-c(rep(NA, length(rstudy)))
    rdummy[1]<-rhx
    rdummy[2]<-rlx
    rydummy<-c(rep(NA, length(rstudy)))
    rydummy[1]<-rhy
    rydummy[2]<-rly
    xtitle<-paste("OR(", ci100, "% CI)")
    plot(rdummy, rydummy, type="n", log="x", ylab="", yaxt="n", xlab=xtitle, main="DerSimonian & Laird (Random-Effects) ORs")
    points(dslOR[1], 1, pch=22, cex=1.0, bg="red", col="red")
    segments(catmapobject$lbci.dsl, 1, catmapobject$ubci.dsl, 1, col="red")
    mtext("Pooled OR",side=2,at=1, cex=0.80)
    abline(v=1.0)
    limit<-nrow(catmapobject$a1)+1
    counts<-c(1:limit)
    for(e in 2:(nrow(catmapobject$a1)+1)){
      f<-e-1
      points(OR[e],counts[e], pch=22, cex=dslprop.weight[f], bg="black", col="black")
      segments(catmapobject$lbci.fe[f],counts[e], catmapobject$ubci.fe[f],counts[e], bg="black", col="black")
      mtext(paste(catmapobject$studyname[f]),side=2, at=counts[e], cex=0.80)
    }
    graphics.off()
  }
}
