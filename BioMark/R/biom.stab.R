lasso.stab <- function(X, Y, sgmnts = NULL, variables = NULL, ...)
{
  lasso.opt <- biom.options()$lasso
  
  ## do one run to obtain lambda sequence
  ## this one also will give a warning in case of an inappropriate family
  lambdas <- as.numeric(colnames(lasso.coef(X, Y, ...)))
  if (is.null(lasso.opt$lambda))
    lasso.opt$lambda <- lambdas
  
  ## get all coefficients - note that they usually are calculated for only
  ## a part of the variables in each iteration
  all.coefs <- lapply(1:ncol(sgmnts),
                      function(i, xx, ss, vv, yy) {
                        huhn <- lasso.coef(xx[-ss[,i], vv[,i]],
                                           yy[-ss[,i]], 
                                           lasso.opt = lasso.opt, ...)
                        (huhn != 0) + 0 ## + 0 to convert to nrs
                      },
                      X, sgmnts, variables, Y)

  x.coef <- matrix(0, ncol(X), length(lambdas))
  dimnames(x.coef) <- list(NULL, lambdas)
  ## count how often non-zero - this can probably be written more elegantly
  for (i in 1:ncol(sgmnts)) 
    x.coef[variables[,i],] <- x.coef[variables[,i],] + all.coefs[[i]]
  
  ## finally, correct for the number of times every variable has had the
  ## chance to be selected. Maximum value in x.coef can be 1.
  noccur <- tabulate(variables, nbins = ncol(X))
  x.coef.sc <- sweep(x.coef, 1, noccur, FUN = "/")

  x.coef.sc
}


########################################################################
### Functions below select by taking the ntop biggest coefficients,
### or, alternatively the ntop fraction of biggest coefficients.
########################################################################

## New selection function for the ntop highest coefs. The last
## dimension will always disappear, and the result is a matrix,
## possibly with only one column. The number of rows is always equal
## to the number of variables

select.aux <- function(object, variables) {
  ntop <- biom.options()$ntop
  nvar <- dim(object)[2]

  ## if ntop is a fraction between 0 and 1, it is taken to mean the
  ## fraction of variables to be selected. Typically 0.1.
  if (ntop > 0 & ntop < 1)
      ntop <- round(ntop * nvar)

  if (is.matrix(object))
    object <- array(object, c(nrow(object), ncol(object), 1))
  
  gooduns <- apply(object, c(1,3),
                   function(x) sort.list(abs(x), decreasing = TRUE)[1:ntop])
  x.coef <- apply(gooduns, 3, function(x) tabulate(x, nbins = nvar))

  ## finally, correct for the number of times every variable has had the
  ## chance to be selected. Maximum value in the result can be 1.
  noccur <- tabulate(variables, nbins = nvar)
  sweep(x.coef, 1, noccur, FUN = "/")
}


pcr.stab <- function(X, Y, ncomp = 2, sgmnts = NULL, variables = NULL, ...)
{
  x.coef <- array(NA, c(ncol(sgmnts), ncol(X), length(ncomp)))
  for (i in 1:ncol(sgmnts))
    x.coef[i,variables[,i],] <-
      pcr.coef(X[-sgmnts[,i], variables[,i]], Y[-sgmnts[,i]],
                 ncomp = ncomp, ...)

  select.aux(x.coef, variables)
}

pls.stab <- function(X, Y, ncomp = 2, sgmnts = NULL, variables = NULL, ...)
{
  x.coef <- array(NA, c(ncol(sgmnts), ncol(X), length(ncomp)))
  for (i in 1:ncol(sgmnts))
    x.coef[i,variables[,i],] <-
      pls.coef(X[-sgmnts[,i], variables[,i]], Y[-sgmnts[,i]],
                 ncomp = ncomp, ...)
  
  select.aux(x.coef, variables)
}

vip.stab <- function(X, Y, ncomp = 2, sgmnts = NULL, variables = NULL, ...)
{
  x.coef <- array(NA, c(ncol(sgmnts), ncol(X), length(ncomp)))
  for (i in 1:ncol(sgmnts)) 
    x.coef[i,variables[,i],] <-
      vip.coef(X[-sgmnts[,i], variables[,i]], Y[-sgmnts[,i]],
               ncomp = ncomp, ...)
  
  select.aux(x.coef, variables)
}

### the dots in the shrinkt.stab and studentt.stab functions are
### necessary to catch extra arguments to other functions.
shrinkt.stab <- function(X, Y, sgmnts = NULL, variables = NULL, ...)
{
  x.coef <- matrix(NA, ncol(sgmnts), ncol(X))
  cat("\n")
  for (i in 1:ncol(sgmnts)) {
    x.coef[i,variables[,i]] <- shrinkt.coef(X[-sgmnts[,i], variables[,i]], 
                                            Y[-sgmnts[,i]], ...)
  }
  
  select.aux(x.coef, variables)
}

studentt.stab <- function(X, Y, sgmnts = NULL, variables = NULL, ...)
{
  x.coef <- matrix(NA, ncol(sgmnts), ncol(X))
  for (i in 1:ncol(sgmnts))
    x.coef[i,variables[,i]] <- studentt.coef(X[-sgmnts[,i], variables[,i]], 
                                             Y[-sgmnts[,i]], ...)
  
  select.aux(x.coef, variables)
}

