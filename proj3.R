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
  
  # lengths of blocks - need this to form the psi matrix but likely will drop it when code more efficient way?
  lengths <- sapply(blocks, ncol)
  
  X <- model.matrix(form, dat = dat)
  
  # simulate random theta initial guess
  theta <- rnorm(1 + length(ref))
  
  # extract y
  y <- dat[[all.vars(form)[1]]]
  
  return(list(Z = Z, X = X, theta = theta, y = y, lengths = lengths))
}

LMMprof <- function(theta, y, Z, X, lengths) {
  ## long way ! code this more efficiently
  sigma <- exp(2 * theta[1])
  
  psi <- diag(rep(theta[-1], lengths))
  
  W <- solve(Z %*% psi %*% t(Z) + diag(length(y)) * sigma)
  
  beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% (W %*% y)
  
  # this outputs the beta_hat vector
  print(beta)
  
  - as.numeric(t(y - X %*% beta) %*% W %*% (y - X %*% beta) + log(det(Z %*% psi %*% t(Z) + diag(length(y)) * sigma))) / 2
}

lmm <- function(form, dat, ref=list()) {
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
  
  # this normally gives an error as solve fails - but not a problem, checking the beta_hat printed out we get gives the same answer as lme4::lmer
  # compare the printed output to
  lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine), data = Machines, REML=FALSE)
  
  # now need to speed up the code, which would also have the effect of removing solve and this issue!
}
