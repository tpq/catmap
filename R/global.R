#' An example data frame for use with catmap
#'
#' A file for use with \code{catmap} containing simulated data.
#'
#' The header must be part of the file and either 0 or NA must be included for
#' data not relevant for the particular study.  For example, using a TDT study
#' the caserisk, controlrisk, casenotrisk and controlnotrisk must have values
#' of either 0 or NA.
#'
#' @name catmapdata
#' @docType data
#' @format A data frame with 5 observations on the following 8 variables.
#' \describe{ \item{name}{a factor with study name and optionally year
#' of publication.  NOTE: if year of publication is included there must be no
#' space between study name and year.  A comma or underscore work nicely.
#' Example: \code{Abrams,2001} \code{Peter,2002} \code{Todd,2003}
#' \code{Wei,2007} \code{Yu,2007}} \item{study}{a numeric vector
#' containing 1 if study is TDT and 2 if case-control} \item{t}{a
#' numeric vector containing counts of alleles transmitted in TDT study}
#' \item{nt}{a numeric vector containing counts of alleles not
#' transmitted in TDT study} \item{caserisk}{a numeric vector
#' containing counts of risk alleles in cases} \item{controlrisk}{a
#' numeric vector containing counts of risk alleles in controls}
#' \item{casenotrisk}{a numeric vector containing counts of non-risk
#' alleles in cases} \item{controlnotrisk}{a numeric vector containing
#' counts of non-risk alleles in controls} }
#' @author Kristin K. Nicodemus, \email{kristin.nicodemus@@well.ox.ac.uk}
#' @seealso \code{\link{catmap.forest}}, \code{\link{catmap.sense}},
#' \code{\link{catmap.cumulative}}, \code{\link{catmap.funnel}}.
#' @keywords datasets
"catmapdata"

#' @importFrom grDevices dev.off graphics.off pdf
#' @importFrom graphics abline mtext par plot points segments
#' @importFrom stats pchisq qnorm
#' @importFrom utils read.table
NULL