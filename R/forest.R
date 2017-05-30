#' Make Forest Plot
#'
#' A back-end wrapper function used to make forest plots.
#'
#' @inheritParams catmap.funnel
#' @param summary A character string. The kind of summary statistic
#'  to plot. Select from "fixed" or "random".
#' @param main A character string. The figure title.
makeForest <- function(catmapobject, summary = "", main = "Main Title", mean = exp(catmapobject$logOR),
                       lower = catmapobject$lbci.fe, upper = catmapobject$ubci.fe,
                       studyname = catmapobject$studyname){

  dat <- data.frame(
    lower = c(NA, lower),
    mean  = c(NA, mean),
    upper = c(NA, upper)
  )

  study <- c("Study", sub(",", " ", studyname))
  summaryv <- c(TRUE, rep(FALSE, length(mean)))

  if(summary == "fixed"){

    main <- "Inverse Variance (Fixed-Effects) ORs"
    dat <- rbind(dat, c(NA, NA, NA))
    dat <- rbind(dat, c(catmapobject$lbci, catmapobject$combinedOR, catmapobject$ubci))
    study <- c(study, NA, "Summary")
    summaryv <- c(summaryv, FALSE, TRUE)

  }else if(summary == "random"){

    main <- "DerSimonian & Laird (Random-Effects) ORs"
    dat <- rbind(dat, c(NA, NA, NA))
    dat <- rbind(dat, c(catmapobject$lbci.dsl, catmapobject$OR.dsl, catmapobject$ubci.dsl))
    study <- c(study, NA, "Summary")
    summaryv <- c(summaryv, FALSE, TRUE)
  }

  tex <- cbind(
    study,
    c(as.character(round(dat$mean, 2))),
    c(as.character(round(dat$lower, 2))),
    c(as.character(round(dat$upper, 2)))
  )

  colnames(tex) <- NULL
  tex[1, ] <- c("Study", "OR", paste0("[", catmapobject$ci, "%"), "CI]")

  f <- forestplot::forestplot(
    labeltext = tex, mean = dat$mean, lower = dat$lower, upper = dat$upper,
    graph.pos = 2, zero = 1, is.summary = summaryv, xlog = TRUE,
    col = forestplot::fpColors(box = "black", line = "grey", summary = "red"),
    title = main
  )

  return(f)
}

#' catmap: Forest Plot
#'
#' The \code{catmap.forest} creates forest plots of the individual study
#'  Odds Ratios (OR) and Confidence Intervals (CI). It then summarizes the
#'  data using a fixed-effects or random-effects pooled OR and CI.
#'
#' @inheritParams catmap.funnel
#' @param fe.forest A boolean. Toggles whether the forest plot should get saved
#'  to the current working directory.
#' @param re.forest A boolean. Toggles whether the forest plot should get saved
#'  to the current working directory.
#'
#' @author Algorithm designed and implemented by Kristin K. Nicodemus.
#'  Code modified and updated by Thom Quinn.
#' @seealso \code{\link{catmap}}, \code{\link{catmap.forest}},
#'  \code{\link{catmap.sense}}, \code{\link{catmap.cumulative}},
#'  \code{\link{catmap.funnel}}
#'
#' @examples
#' data(catmapdata)
#' catmapobject <- catmap(catmapdata, 0.95, TRUE)
#' catmap.forest(catmapobject, TRUE, TRUE)
#' @export
catmap.forest <- function(catmapobject, fe.forest = FALSE, re.forest = FALSE){

  if(!fe.forest & !re.forest) par(ask = TRUE)

  # FE
  if(fe.forest) pdf(file = paste0(catmapobject$dataset, ".fixed.effects.plot.pdf"))
  f1 <- makeForest(catmapobject, summary = "fixed")
  if(fe.forest) graphics.off()

  # RE
  if(catmapobject$tau2 <= 0){
    message("NOTICE: tau2 is less than or equal to 0.\n",
            " No random effects estimates calculated.\n")
  }else{
    if(re.forest) pdf(file = paste0(catmapobject$dataset, ".random.effects.plot.pdf"))
    f2 <- makeForest(catmapobject, summary = "random")
    if(re.forest) graphics.off()
    f1 <- list(f1, f2)
  }

  if(!fe.forest & !re.forest) par(ask = FALSE)
  return(f1)
}
