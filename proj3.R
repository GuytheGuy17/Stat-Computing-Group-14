# Project 4

#' Group Members: 
#' Guy McClennan - s2036567
#' Alexandru Girban - s2148980
#' Louis Bennett - s2744241

#' LMMsetup
#' 
#' @description Sets up the matrices involved in the mixed linear model
#' @param form The model formula.
#' @param dat data frame containing all the variables needed in the model.
#' @param ref a list of vectors of variable names specifying the random effects for the Zb part of the model
#' @return The initialized parameters and necessary variables for the mixed linear model, in a list

LMMsetup <- function(form, dat, ref) {
  # loop over ref, creating the blocks in Z
  blocks <- lapply(ref, \(x) {
    # if length is 1, just take the term
    if(length(x) == 1) {
      block_form <- as.formula(paste0(' ~ ', x, ' - 1'))
    } else {
      block_form <- as.formula(paste0(' ~ ', paste(x, collapse = ':'), ' - 1'))
    }
    
    model.matrix(block_form, dat = dat)
  })
  
  # bind into single matrix Z
  Z <- do.call('cbind', blocks)
  
  # lengths of blocks - need this to form the psi matrix later
  lengths <- sapply(blocks, ncol)
  
  # handle the no random effects case
  if(is.null(Z)) {
    Z <- matrix(nrow = 0, ncol = 0)
    lengths <- 0
  }
  
  X <- model.matrix(form, dat = dat)
  
  # simulate random theta initial guess
  theta <- rnorm(1 + length(ref))
  
  # extracting y
  y <- model.response(model.frame(form, dat))
  # y <-  dat[[all.vars(form)[1]]]
  
  return(list(Z = Z, X = X, theta = theta, y = y, lengths = lengths))
}

#' LMMprof
#' 
#' @description evaluates the negative log-likelihood for a given theta
#' @param theta 
#' @param y response part of model
#' @param Z model matrix associated with the random effects
#' @param X predictor part of model
#' @param lengths length of blocks, needed to form the psi matrix
#' @return The initialized parameters and necessary variables for the mixed linear model, in a list

LMMprof <- function(theta, y, Z, X, lengths) {
  # extract basic info from the inputs
  sigma <- exp(2 * theta[1])
  n <- length(y)
  p <- ncol(Z)
  
  # handle the case of no random effects
  if(p == 0) {
    # if p = 0, simplifies to just the linear model case
    det_term <- n * log(sigma)
    
    # use Cholesky factor to compute multiplication by (t(X) %*% X) ^ -1
    L <- chol(t(X) %*% X)
    beta <- backsolve(L, forwardsolve(t(L), t(X) %*% y))
    
    loglik <- as.numeric(t(y - X %*% beta) %*% (y - X %*% beta) + det_term) / 2
  } else {
    # perform QR decomposition - from our tests even running 1000 times took less than 1 second so fine to do inside LMMprof
    qr <- qr(Z)
    R <- qr.R(qr)
    
    # if p = 1 then psi should just be a 1-1 matrix
    if(p > 1) {
      psi <- diag(rep(exp(2 * theta[-1]), lengths))
    } else {
      psi <- exp(2 * theta[-1])
    }
    
    # form Cholesky factor
    L <- chol(R %*% psi %*% t(R) + diag(p) * sigma)
    
    # calculate the determinant part of likelihood
    det_term <- 2 * sum(log(diag(L))) + (n - p) * log(sigma)
    
    ## compute Wy
    # multiplication by t(Q)
    qty <- qr.qty(qr, y)
    
    # block multiplication - 2nd block is simply dividing by sigma
    qty[(p + 1):n] <- qty[(p + 1):n] / sigma
    
    # using the Cholesky decomposition to speed up this multiplication for the 1st block
    qty[1:p] <- backsolve(L, forwardsolve(t(L), qty[1:p]))
    
    # finally, multiplication by Q
    Wy <- qr.qy(qr, qty)
    
    ## compute WX
    # multiplication by t(Q)
    qtX <- qr.qty(qr, X)
    
    # block multiplication - 2nd block
    qtX[(p + 1):n, ] <- qtX[(p + 1):n, ] / sigma
    
    # same use of Cholesky decomposition as before
    qtX[1:p, ] <- backsolve(L, forwardsolve(t(L), qtX[1:p, ]))
    
    # multiplication by Q
    WX <- qr.qy(qr, qtX)
    
    # use another Cholesky factor to speed up the calculation of beta
    L2 <- chol(t(X) %*% WX)
    beta <- backsolve(L2, forwardsolve(t(L2), t(X) %*% Wy))
    
    # final computation of log-likelihood - note we take the positive as we are minimising!
    loglik <- as.numeric(t(y - X %*% beta) %*% (Wy - WX %*% beta) + det_term) / 2
  }
  
  attr(loglik, 'beta') <- as.vector(beta)
  
  return(loglik)
}

lmm <- function(form, dat, ref = list()) {
  ## validating inputs 
  # all elements of ref must be length >= 1
  if(!all(sapply(ref, length)) >= 1) {
    stop('All the random effects must have length >= 1')
  }
  
  # all elements in the formula must be in the data
  vars <- all.vars(form)
  if(!all(vars %in% colnames(dat))) {
    stop(paste0('Please ensure all variables are present in the data: `', vars[which(!vars %in% colnames(dat))[1]], '` wasn\'t found.'))
  }
  
  # make sure ref is a list and dat is a data frame
  ref <- as.list(ref)
  dat <- as.data.frame(dat)
  
  # check all random effect variables are present
  ref_vars <- unique(unlist(ref))
  if(!all(ref_vars %in% colnames(dat))) {
    stop(paste0('Please ensure all variables are present in the data: `', ref_vars[which(!ref_vars %in% colnames(dat))[1]], '` wasn\'t found.'))
  }
  
  # remove NA values from dat & check that there are some non NA observations
  dat <- dat[unique(c(ref_vars, vars))]
  dat <- na.omit(dat)
  
  if(length(dat) == 0) {
    stop('No non-NA cases found in `dat`')
  }
  
  # add any others !
  
  inits <- LMMsetup(form, dat, ref)
  
  
  ## Fit optimise function for single variable case
  if((length(inits$theta)) > 1){
    # optimise over theta
    run <- optim(
      par = inits$theta,
      fn = LMMprof,
      y = inits$y,
      Z = inits$Z,
      X = inits$X,
      lengths = inits$lengths,
      method = 'Nelder-Mead'
    )
    theta <- run$par
    loglik_value <- run$value
  } else {
    # Single-parameter optimization using optimise
    run <- optimise(
      f = LMMprof,
      interval = c(1, 1000),
      y = inits$y,
      Z = inits$Z,
      X = inits$X,
      lengths = inits$lengths
    )
    theta <- run$minimum
    loglik_value <- run$objective
  }
  
  # output the log-likelihood for checking with tests
  print(paste0('The best loglikelihood found was: ', - round(run$value + length(inits$y) / 2 * log(2 * pi), digits = 3)))
  
  # extract the beta vector too
  beta <- attr(LMMprof(theta, y = inits$y, Z = inits$Z, X = inits$X, lengths = inits$lengths), 'beta')
  
  # add labels on both beta and theta
  names(beta) <- colnames(inits$X)
  
  names(theta) <- c('Residual', sapply(ref, \(x) paste(x, collapse = ':')))
  
  return(list(theta = theta, beta = beta, stdevs = exp(theta)))
}
