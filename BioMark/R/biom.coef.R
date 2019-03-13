pcr.coef <- function(X, Y, ncomp, ...)
{
  if (nlevels(Y) > 2)
    stop("multi-class discrimination not implemented for PCR")

  Y <- as.numeric(Y)
  matrix(svdpc.fit(X, Y, ncomp = max(ncomp),
                   stripped = TRUE)$coefficients[, 1, ncomp],
         ncol(X), length(ncomp))
}


## Changed to widekernelpls.fit because this probably is the most
## relevant situation  
pls.coef <- function(X, Y, ncomp, ...)
{
  if (nlevels(Y) > 2)
    stop("multi-class discrimination not implemented for PLS")

  Y <- as.numeric(Y)
  matrix(widekernelpls.fit(X, Y, ncomp = max(ncomp),
                           stripped = TRUE)$coefficients[, 1, ncomp],
         ncol(X), length(ncomp))
}

vip.coef <- function(X, Y, ncomp, ...)
{
  if (nlevels(Y) > 2)
    stop("multi-class discrimination not implemented for VIP")

  Y <- as.numeric(Y)
  plsmod <- plsr(Y ~ X, ncomp = max(ncomp), method = "widekernelpls")
  ww <- loading.weights(plsmod)

  result <- matrix(NA, ncol(X), length(ncomp))
  for (i in 1:length(ncomp)) {
    var.exp <- diff(c(0, R2(plsmod, estimate = "train",
                            ncomp = 1:ncomp[i], intercept = FALSE)$val))

    result[,i] <- sqrt(ncol(X) * ww[,1:ncomp[i],drop = FALSE]^2 %*%
                       var.exp / sum(var.exp))
  }

  result
}

studentt.coef <- function(X, Y, ...)
{
  if (nlevels(Y) > 2)
    stop("only two-class discrimination implemented for studentt")
  
  TFUN <- studentt.fun(Y)
  
  matrix(TFUN(X), ncol = 1)
}

shrinkt.coef <- function(X, Y, ...)
{
  if (nlevels(Y) > 2)
    stop("only two-class discrimination implemented for shrinkt")
  
  TFUN <- shrinkt.fun(L =  Y, var.equal = FALSE, verbose = FALSE)
  
  matrix(TFUN(X), ncol = 1)
}

## Nov 21, 2011: inclusion of the lasso. For classification, Y should
## be a factor!
lasso.coef <- function(X, Y, lasso.opt = biom.options()$lasso, ...)
{
  ## check whether family and character of Y agree
  fam <- lasso.opt$family
  if (!is.null(fam)) {
    if (!is.factor(Y)) {
      if (fam != "gaussian")
        stop("Attempt of regression with a family different than 'gaussian'")
    } else {
      if (fam != "binomial")
        stop("Attempt of binary classification with a family different than 'binomial'")
    }
  } else {
    if (!is.factor(Y)) {
      lasso.opt$family <- "gaussian"
    } else {
      lasso.opt$family <- "binomial"
    }
  }

  glmargs <- c(list(x = X, y = Y),
                    lasso.opt)
  
  huhn <- do.call(glmnet, glmargs)
  x.coef <- as.matrix(huhn$beta)
  colnames(x.coef) <- huhn$lambda
  
  x.coef
}

