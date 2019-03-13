\name{get.biom}
\alias{get.biom}
\alias{print.BMark}
\alias{summary.BMark}
\alias{coef.BMark}
\title{Get biomarkers discriminating between two classes}
\description{Biomarkers can be identified in several ways: the classical way
  is to look at those variables with large model coefficients or large t
  statistics. One other is based on the higher criticism approach (HC),
  and the third possibility assesses the stability of these coefficients
  under subsampling of the data set.
}
\usage{
get.biom(X, Y, fmethod = "all", type = c("stab", "HC", "coef"),
         ncomp = 2, biom.opt = biom.options(), 
         \dots)
\method{coef}{BMark}(object, \dots)
\method{print}{BMark}(x, \dots)
\method{summary}{BMark}(object, \dots)
}
\arguments{
  \item{X}{Data matrix. Usually the number of columns (variables) is
    (much) larger than the number of rows (samples).}
  \item{Y}{Class indication. For classification with two or more factors
    a factor; a numeric vector will be interpreted as a regression
    situation, which can only be tackled by \code{fmethod = "lasso"}.}
  \item{fmethod}{Modelling method(s) employed. The default is to use
    \code{"all"}, which will test all methods in the current
    \code{biom.options$fmethods} list. Note that from version 0.4.0,
    \code{"plsda"} and \code{"pclda"} are no longer in the list of
    methods - they have been replaced by \code{"pls"} and \code{"pcr"},
    respectively. For compatibility reasons, using the old terms will
    not lead to an error but only a warning.} 
  \item{type}{Whether to use coefficient size as a criterion
    (\code{"coef"}), \code{"stab"} or \code{"HC"}.}
  \item{ncomp}{Number of latent variables to use in PCR and PLS (VIP)
    modelling. In function \code{get.biom} this may be a vector; in all
    other functions it should be one number. Default: 2.}
  \item{biom.opt}{Options for the biomarker selection - a list with
    several named elements. See \code{\link{biom.options}}.}
  \item{object, x}{A BMark object.}
  \item{\dots}{Further arguments for modelling functions. Often used to
    catch unused arguments.}
}
\value{Function \code{get.biom} returns an object of class "BMark", a
  list containing an element 
  for every \code{fmethod} that is selected, as well as an element
  \code{info}. The individual elements contain information depending on
  the type chosen: for \code{type == "coef"}, the only element returned
  is a matrix containing coefficient sizes. For \code{type == "HC"}
  and \code{type == "stab"}, a list is returned containing elements
  \code{biom.indices}, and either \code{pvals} (for \code{type == "HC"})
  or \code{fraction.selected} (for \code{type == "stab"}).
  Element \code{biom.indices} contains the indices of
  the selected variables, and can be extracted using function
  \code{selection}. Element \code{pvals} contains the p values
  used to perform HC thresholding; these are presented in the original
  order of the variables, and can be obtained directly from e.g. t
  statistics, or from permutation sampling. Element
  \code{fraction.selected} indicates in what fraction of the
  stability selection iterations a particular variable has been
  selected. The more often it has been selected, the more stable it is
  as a biomarker. Generic function \code{coef.biom} extracts model
  coefficients, p values or stability fractions for types \code{"coef"},
  \code{"HC"} and \code{"stab"}, respectively.}
\author{Ron Wehrens}
\seealso{\code{\link{biom.options}}, \code{\link{get.segments}},
  \code{\link{selection}}}
\examples{
## Real apple data (small set)
data(spikedApples)
X <- sqrt(spikedApples$dataMatrix)
Y <- factor(rep(0:1, each = 10))
 
apple.coef <- get.biom(X = X, Y = Y, ncomp = 2:3, type = "coef")
coef.sizes <- coef(apple.coef) 
sapply(coef.sizes, range)

## stability-based selection
set.seed(17)
apple.stab <- get.biom(X = X, Y = Y, ncomp = 2:3, type = "stab")
selected.variables <- selection(apple.stab)
unlist(sapply(selected.variables, function(x) sapply(x, length)))
## Ranging from approx 70 for pcr, approx 40 for pls and t statistics,
## to 0-29 for the lasso

## check agreement with "true" biomarkers: here the number of true positives
unlist(
  sapply(selected.variables,
         function(x)
         lapply(x, function(xx) sum(xx \%in\% spikedApples$biom))))

## False positives	 
unlist(
  sapply(selected.variables,
         function(x)
         lapply(x, function(xx) sum(!(xx \%in\% spikedApples$biom)))))

## regression example
data(gasoline) ## from the pls package
gasoline.stab <- get.biom(gasoline$NIR, gasoline$octane,
                          fmethod = c("pcr", "pls", "lasso"), type = "stab")


\dontrun{
## Same for HC-based selection, takes some time
apple.HC <- get.biom(X = X, Y = Y, ncomp = 2:3, type = "HC")

## TP
unlist(
  sapply(apple.HC[names(apple.HC) != "info"],
         function(x)
         lapply(x,
	        function(xx) sum(xx$biom.indices \%in\% spikedApples$biom))))
## FP
unlist(
  sapply(apple.HC[names(apple.HC) != "info"],
         function(x) 
         lapply(x,
	        function(xx) sum(!(xx$biom.indices \%in\% spikedApples$biom)))))
}
}
\keyword{multivariate}
\keyword{models}
\keyword{htest}
