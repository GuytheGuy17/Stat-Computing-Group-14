# Project 4

#' Group Members: 
#' Guy McClennan - s2036567
#' Alexandru Girban - s2148980
#' Louis Bennett - s2744241

LMMsetup <- function(form, dat, ref) {
  # lopp over ref, creating the blocks in Z
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
  
  X <- model.matrix(form, dat = dat)
  
  # simulate random theta initial guess
  theta <- rnorm(1 + length(ref))
  
  # extracting y
  y <- dat[[all.vars(form)[1]]]
  
  return(list(Z = Z, X = X, theta = theta, y = y, lengths = lengths))
}

LMMprof <- function(theta, y, Z, X, lengths) {
  # extract basic info from the inputs
  sigma <- exp(2 * theta[1])
  p <- ncol(Z)
  n <- length(y)
  
  # perform qr decomposition - from our tests even running 1000 times took less than 1 second so fine to do inside LMMprof
  qr <- qr(Z)
  R <- qr.R(qr)
  
  psi <- diag(rep(exp(2 * theta[-1]), lengths))
  
  # form cholesky factor
  L <- chol(R %*% psi %*% t(R) + diag(p) * sigma)
  
  # calculate the determinant part of likelihood
  det_term <- 2 * sum(log(diag(L))) + (n - p) * log(sigma)
  
  ## compute Wy
  # multiplication by t(Q)
  qty <- qr.qty(qr, y)
  
  # block multiplication - 2nd block is simply dividing by sigma
  qty[(p + 1):n] <- qty[(p + 1):n] / sigma
  
  # using the Cholesky decomposition to speed up this multiplication for the 1st block
  qty[1:p] <- forwardsolve(L, backsolve(L, qty[1:p], transpose = TRUE))
  
  # finally, multiplication by Q
  Wy <- qr.qy(qr, qty)
  
  ## compute WX
  # multiplication by t(Q)
  qtX <- qr.qty(qr, X)
  
  # block multiplication - 2nd block
  qtX[(p + 1):n, ] <- qtX[(p + 1):n, ] / sigma
  
  # same use of Cholesky decomposition as before
  qtX[1:p, ] <- forwardsolve(L, backsolve(L, qtX[1:p, ], transpose = TRUE))
  
  # multiplication by Q
  WX <- qr.qy(qr, qtX)
  
  # use another Cholesky factor to speed up the calculation of beta
  L2 <- chol(t(X) %*% WX)
  beta <- forwardsolve(L2, backsolve(L2, (t(X) %*% Wy), transpose = TRUE))
  
  # final computation of log likelihood - note we take the positive as we are minimising!
  loglik <- as.numeric(t(y - X %*% beta) %*% (Wy - WX %*% beta) + det_term) / 2
  
  attr(loglik, 'beta') <- beta
  
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
  
  # make sure dat is a data frame
  dat <- as.data.frame(dat)
  
  # add any others !
  
  inits <- LMMsetup(form, dat, ref)
  
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
  
  # extract the beta vector too
  beta <- attr(LMMprof(run$par, y = inits$y, Z = inits$Z, X = inits$X, lengths = inits$lengths), 'beta')
  
  return(list(theta = run$par, beta = beta))
}
