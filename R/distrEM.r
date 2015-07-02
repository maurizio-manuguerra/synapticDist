# EM implementation

# M is number of mixture components - probably best to rename it
# x_expanded is an NxM matrix where column 1 is x_i and column 2 is (x_i - theta_i)
invgaussmixEM <- function(x_expanded, num.components=2, initials=NULL, max.iters=100, epsilon=0.000001) {
  # TODO check that arguments are valid

  # init
  iterations <- 0
  N <- nrow(x_expanded)

  # Set initial conditions
  # FIXME: update this, doesn't work for variable M
  mu <- matrix(rnorm(num.components, 1, 0.1), nrow=1)
  lambda <- matrix(rnorm(num.components, 1, 0.1), nrow=1)
  alpha <- matrix(rep(1/num.components, num.components), nrow=1)  # mixing components

  #if (!is.null(initials)) {
  #    for (m in 1:num.components) {
  #        mu[1, m] <- initials[[m]]$mu
  #        lambda[1, m] <- initials[[m]]$lambda
  #        alpha[1, m] <- initials[[m]]$alpha
  #    }
  #}
  member.prob <- matrix(alpha, nrow=N, ncol=num.components, byrow=TRUE)
  #Avoid starting with all 0.5 as we are going to round for Iz.
  #noise <- runif(N,0,0.01) #needs to be updated to M>2 (but anyway why using the hard complete data likelihood?)
  #member.prob <- member.prob + cbind(noise, -noise)
  diff <- 1
  log.lik <- -10e9

  # x_expanded <- matrix(x, nrow=N, ncol=num.components, byrow=FALSE)

  # algorithm based off http://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm#Gaussian_mixture

  while (diff > epsilon && iterations <= max.iters) {
    iterations <- iterations + 1

    # print(sprintf("iteration %d mu0=%f mu1=%f la0=%f la1=%f al0=%f al1=%f", iterations, mu[1], mu[2], lambda[1], lambda[2], alpha[1], alpha[2]))

    # TODO: this might be faster/better expressed using rep_len
    mu_expanded <- matrix(mu, nrow=N, ncol=num.components, byrow=TRUE)
    lambda_expanded <- matrix(lambda, nrow=N, ncol=num.components, byrow=TRUE)
    alpha_expanded <- matrix(alpha, nrow=N, ncol=num.components, byrow=TRUE)
    # E-step: calculate Q

    # calculate p(l|x_i, Theta^g)
    x.prob <- matrix(dinvgauss(x_expanded, mean=mu_expanded, shape=lambda_expanded), nrow=N, ncol=num.components, byrow=FALSE)  # N x M matrix
    x.prob[is.na(x.prob)] <- 0
    log.lik.old <- log.lik
    #log.lik <- sum(log(rowSums(alpha_expanded * x.prob))) #Ian: incompl. data likelihood
    log.lik <- sum(log(rowSums(member.prob * x.prob))) #Maurizio: compl. data likelihood (soft)
    diff <- (log.lik - log.lik.old)
    if (diff < epsilon) {
      break
    }

    # then come up with new parameter estimates using the equations that you derived

    # membership probabilities
    weighted.prob <- alpha_expanded * x.prob  # per-component weighted sum of probabilities (NxM)
    sum.prob <- rowSums(weighted.prob)  # second line of the T function from wikipedia (Nx1)
    member.prob <- weighted.prob / sum.prob

    # Elements of sum.prob can go to 0 if the density is low enough.
    # This causes divide-by-zero in the member.prob calculation. Replace
    # any NaNs with 0.
    member.prob[is.nan(member.prob)] <- 0

    # we've got components across columns and observations down rows, so we do all of the summations simultaneously on both
    member.prob.sum <- colSums(member.prob)
    alpha.new <- member.prob.sum / N  # should be 1xM matrix
    mu.new <- colSums(x_expanded * member.prob) / member.prob.sum  # should be 1xM matrix
    ##Check next couple lines
    aa=((x_expanded - mu_expanded) ^ 2 * member.prob) / (mu_expanded ^ 2 * x_expanded)
    aa[is.nan(aa)] <- 0
    lambda.new <- member.prob.sum / colSums(aa)

    mu <- mu.new
    lambda <- lambda.new
    alpha <- alpha.new

    # LATER: it might make sense to try this algorithm using RGPU or something - it would be ideal to be able to plug in R code to your CUDA code. tHis might be a halfway option taht is fast enough
  }
  ii = which(member.prob>0)
  Q <- sum(member.prob[ii] * log(x.prob[ii]))
  result <- list(x=x_expanded, alpha=alpha, mu=mu, lambda=lambda, member.prob=member.prob, llik=log.lik, Q=Q, fit_success=TRUE)
  #result <- list(x=x_expanded, alpha=alpha, mu=mu, lambda=lambda, member.prob=round(member.prob, 3), llik=log.lik, Q=Q, fit_success=TRUE)
  class(result) <- "mixEM"
  result
}

invgaussexpneuromixEM <- function(x_expanded, num.components=2, initials=NULL, max.iters=100, epsilon=0.000001) {
  # TODO check that arguments are valid

  # init
  iterations <- 0
  N <- nrow(x_expanded)

  # Set initial conditions
  # FIXME: update this, doesn't work for variable M
  mu <- matrix(rnorm(num.components, 1, 0.1), nrow=1)
  lambda <- matrix(rnorm(num.components, 1, 0.1), nrow=1)
  alpha <- matrix(rep(1/num.components, num.components), nrow=1)  # mixing components

  if (!is.null(initials)) {
    for (m in 1:num.components) {
      mu[1, m] <- initials[[m]]$mu
      lambda[1, m] <- initials[[m]]$lambda
      alpha[1, m] <- initials[[m]]$alpha
    }
  }
  member.prob <- matrix(alpha, nrow=N, ncol=num.components, byrow=TRUE)
  #Avoid starting with all 0.5 as we are going to round for Iz.
  #noise <- runif(N,0,0.01) #needs to be updated to M>2 (but anyway why using the hard complete data likelihood?)
  #member.prob <- member.prob + cbind(noise, -noise)
  diff <- 1
  log.lik <- -10e9

  # x_expanded <- matrix(x, nrow=N, ncol=num.components, byrow=FALSE)

  # algorithm based off http://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm#Gaussian_mixture

  while (diff > epsilon && iterations <= max.iters) {
    iterations <- iterations + 1

    # print(sprintf("iteration %d mu0=%f mu1=%f la0=%f la1=%f al0=%f al1=%f", iterations, mu[1], mu[2], lambda[1], lambda[2], alpha[1], alpha[2]))

    # TODO: should not use same initials for IG and gamma distributions
    mu_expanded <- matrix(mu, nrow=N, ncol=num.components, byrow=TRUE)
    lambda_expanded <- matrix(lambda, nrow=N, ncol=num.components, byrow=TRUE)
    alpha_expanded <- matrix(alpha, nrow=N, ncol=num.components, byrow=TRUE)
    # E-step: calculate Q

    # calculate p(l|x_i, Theta^g)
    x.prob <- matrix(c(dinvgauss(x_expanded[,1], mean=mu_expanded[1], shape=lambda_expanded[1]), dexp(x_expanded[,2], rate=1/mu_expanded[2])), nrow=N, ncol=num.components, byrow=FALSE)  # N x M matrix

    # have we converged to within epsilon?
    # we do this up here and not at the end as we're doing parts of the calculation anyway
    log.lik.old <- log.lik
    #log.lik <- sum(log(rowSums(alpha_expanded * x.prob))) #Ian: incompl. data likelihood
    log.lik <- sum(log(rowSums(member.prob * x.prob))) #Maurizio: compl. data likelihood (soft)
    #Iz <- round(member.prob)
    #log.lik <- sum(log(rowSums(Iz * x.prob))) #Maurizio: compl. data likelihood (hard)
    # TODO this abs isn't ideal - shouldn't be necessary
    # TODO also catch if log-lik increases, which shouldn't happen (supposed to flag non-convergence and then try again with new initial conditions)
    #diff <- abs(log.lik.old - log.lik)
    diff <- (log.lik - log.lik.old)
    if (diff < epsilon) {
      break
    }

    # then come up with new parameter estimates using the equations that you derived

    # membership probabilities
    weighted.prob <- alpha_expanded * x.prob  # per-component weighted sum of probabilities (NxM)
    sum.prob <- rowSums(weighted.prob)  # second line of the T function from wikipedia (Nx1)
    member.prob <- weighted.prob / sum.prob

    # Elements of sum.prob can go to 0 if the density is low enough.
    # This causes divide-by-zero in the member.prob calculation. Replace
    # any NaNs with 0.
    member.prob[is.nan(member.prob)] <- 0

    # we've got components across columns and observations down rows, so we do all of the summations simultaneously on both
    member.prob.sum <- colSums(member.prob)
    alpha.new <- member.prob.sum / N  # should be 1xM matrix
    #IG and exp
    mu.new <- colSums(x_expanded * member.prob) / member.prob.sum  # should be 1xM matrix
    #IG
    aa.1=matrix(((x_expanded[,1] - mu_expanded[,1]) ^ 2 * member.prob[,1]) / (mu_expanded[,1] ^ 2 * x_expanded[,1]),nrow=N, ncol=1)
    aa.1[is.nan(aa.1)] <- 0
    lambda.new.1 <- member.prob.sum[1] / colSums(aa.1)
    lambda.new <- c(lambda.new.1, 0)

    mu <- mu.new
    lambda <- lambda.new
    alpha <- alpha.new

    # LATER: it might make sense to try this algorithm using RGPU or something - it would be ideal to be able to plug in R code to your CUDA code. tHis might be a halfway option taht is fast enough
  }

  result <- list(x=x_expanded, alpha=alpha, mu=mu, lambda=lambda, member.prob=member.prob, llik=log.lik, fit_success=TRUE)
  class(result) <- "neuromixEM"
  result
}


invgaussneuromixEM <- function(x_expanded, num.components=2, initials=NULL, max.iters=100, epsilon=0.000001) {
  # TODO check that arguments are valid

  # init
  iterations <- 0
  N <- nrow(x_expanded)

  # Set initial conditions
  # FIXME: update this, doesn't work for variable M
  mu <- matrix(rnorm(num.components, 1, 0.1), nrow=1)
  lambda <- matrix(rnorm(num.components, 1, 0.1), nrow=1)
  alpha <- matrix(rep(1/num.components, num.components), nrow=1)  # mixing components

  if (!is.null(initials)) {
    for (m in 1:num.components) {
      mu[1, m] <- initials[[m]]$mu
      lambda[1, m] <- initials[[m]]$lambda
      alpha[1, m] <- initials[[m]]$alpha
    }
  }
  member.prob <- matrix(alpha, nrow=N, ncol=num.components, byrow=TRUE)
  #Avoid starting with all 0.5 as we are going to round for Iz.
  #noise <- runif(N,0,0.01) #needs to be updated to M>2 (but anyway why using the hard complete data likelihood?)
  #member.prob <- member.prob + cbind(noise, -noise)
  diff <- 1
  log.lik <- -10e9

  # x_expanded <- matrix(x, nrow=N, ncol=num.components, byrow=FALSE)

  # algorithm based off http://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm#Gaussian_mixture

  while (diff > epsilon && iterations <= max.iters) {
    iterations <- iterations + 1

    # print(sprintf("iteration %d mu0=%f mu1=%f la0=%f la1=%f al0=%f al1=%f", iterations, mu[1], mu[2], lambda[1], lambda[2], alpha[1], alpha[2]))

    # TODO: this might be faster/better expressed using rep_len
    mu_expanded <- matrix(mu, nrow=N, ncol=num.components, byrow=TRUE)
    lambda_expanded <- matrix(lambda, nrow=N, ncol=num.components, byrow=TRUE)
    alpha_expanded <- matrix(alpha, nrow=N, ncol=num.components, byrow=TRUE)
    # E-step: calculate Q

    # calculate p(l|x_i, Theta^g)
    x.prob <- matrix(dinvgauss(x_expanded, mean=mu_expanded, shape=lambda_expanded), nrow=N, ncol=num.components, byrow=FALSE)  # N x M matrix

    # have we converged to within epsilon?
    # we do this up here and not at the end as we're doing parts of the calculation anyway
    log.lik.old <- log.lik
    #log.lik <- sum(log(rowSums(alpha_expanded * x.prob))) #Ian: incompl. data likelihood
    log.lik <- sum(log(rowSums(member.prob * x.prob))) #Maurizio: compl. data likelihood (soft)
    #Iz <- round(member.prob)
    #log.lik <- sum(log(rowSums(Iz * x.prob))) #Maurizio: compl. data likelihood (hard)
    # TODO this abs isn't ideal - shouldn't be necessary
    # TODO also catch if log-lik increases, which shouldn't happen (supposed to flag non-convergence and then try again with new initial conditions)
    #diff <- abs(log.lik.old - log.lik)
    diff <- (log.lik - log.lik.old)
    if (diff < epsilon) {
      break
    }

    # then come up with new parameter estimates using the equations that you derived

    # membership probabilities
    weighted.prob <- alpha_expanded * x.prob  # per-component weighted sum of probabilities (NxM)
    sum.prob <- rowSums(weighted.prob)  # second line of the T function from wikipedia (Nx1)
    member.prob <- weighted.prob / sum.prob

    # Elements of sum.prob can go to 0 if the density is low enough.
    # This causes divide-by-zero in the member.prob calculation. Replace
    # any NaNs with 0.
    member.prob[is.nan(member.prob)] <- 0

    # we've got components across columns and observations down rows, so we do all of the summations simultaneously on both
    member.prob.sum <- colSums(member.prob)
    alpha.new <- member.prob.sum / N  # should be 1xM matrix
    mu.new <- colSums(x_expanded * member.prob) / member.prob.sum  # should be 1xM matrix
    ##Check next couple lines
    aa=((x_expanded - mu_expanded) ^ 2 * member.prob) / (mu_expanded ^ 2 * x_expanded)
    aa[is.nan(aa)] <- 0
    lambda.new <- member.prob.sum / colSums(aa)

    mu <- mu.new
    lambda <- lambda.new
    alpha <- alpha.new

    # LATER: it might make sense to try this algorithm using RGPU or something - it would be ideal to be able to plug in R code to your CUDA code. tHis might be a halfway option taht is fast enough
  }

  result <- list(x=x_expanded, alpha=alpha, mu=mu, lambda=lambda, member.prob=member.prob, llik=log.lik, fit_success=TRUE)
  class(result) <- "neuromixEM"
  result
}

invgaussgammaneuromixEM <- function(x_expanded, num.components=2, initials=NULL, max.iters=100, epsilon=0.000001) {
  require(distr)
  # TODO check that arguments are valid

  # init
  iterations <- 0
  N <- nrow(x_expanded)

  # Set initial conditions
  # FIXME: update this, doesn't work for variable M
  mu <- matrix(rnorm(num.components, 1, 0.1), nrow=1)
  lambda <- matrix(rnorm(num.components, 1, 0.1), nrow=1)
  alpha <- matrix(rep(1/num.components, num.components), nrow=1)  # mixing components

  if (!is.null(initials)) {
    for (m in 1:num.components) {
      mu[1, m] <- initials[[m]]$mu
      lambda[1, m] <- initials[[m]]$lambda
      alpha[1, m] <- initials[[m]]$alpha
    }
  }
  member.prob <- matrix(alpha, nrow=N, ncol=num.components, byrow=TRUE)
  #Avoid starting with all 0.5 as we are going to round for Iz.
  #noise <- runif(N,0,0.01) #needs to be updated to M>2 (but anyway why using the hard complete data likelihood?)
  #member.prob <- member.prob + cbind(noise, -noise)
  diff <- 1
  log.lik <- -10e9

  # x_expanded <- matrix(x, nrow=N, ncol=num.components, byrow=FALSE)

  # algorithm based off http://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm#Gaussian_mixture

  while (diff > epsilon && iterations <= max.iters) {
    iterations <- iterations + 1

    # print(sprintf("iteration %d mu0=%f mu1=%f la0=%f la1=%f al0=%f al1=%f", iterations, mu[1], mu[2], lambda[1], lambda[2], alpha[1], alpha[2]))

    # TODO: should not use same initials for IG and gamma distributions
    mu_expanded <- matrix(mu, nrow=N, ncol=num.components, byrow=TRUE)
    lambda_expanded <- matrix(lambda, nrow=N, ncol=num.components, byrow=TRUE)
    alpha_expanded <- matrix(alpha, nrow=N, ncol=num.components, byrow=TRUE)
    # E-step: calculate Q

    # calculate p(l|x_i, Theta^g)
    x.prob <- matrix(c(dinvgauss(x_expanded[,1], mean=mu_expanded[1], shape=lambda_expanded[1]),
                       dgamma(x_expanded[,2], shape=lambda_expanded[2], scale=mu_expanded[2]/lambda_expanded[2])),
                     nrow=N, ncol=num.components, byrow=FALSE)  # N x M matrix

    # have we converged to within epsilon?
    # we do this up here and not at the end as we're doing parts of the calculation anyway
    log.lik.old <- log.lik
    #log.lik <- sum(log(rowSums(alpha_expanded * x.prob))) #Ian: incompl. data likelihood
    log.lik <- sum(log(rowSums(member.prob * x.prob))) #Maurizio: compl. data likelihood (soft)
    #Iz <- round(member.prob)
    #log.lik <- sum(log(rowSums(Iz * x.prob))) #Maurizio: compl. data likelihood (hard)
    # TODO this abs isn't ideal - shouldn't be necessary
    # TODO also catch if log-lik increases, which shouldn't happen (supposed to flag non-convergence and then try again with new initial conditions)
    #diff <- abs(log.lik.old - log.lik)
    diff <- (log.lik - log.lik.old)
    if (diff < epsilon) {
      break
    }

    # then come up with new parameter estimates using the equations that you derived

    # membership probabilities
    weighted.prob <- alpha_expanded * x.prob  # per-component weighted sum of probabilities (NxM)
    sum.prob <- rowSums(weighted.prob)  # second line of the T function from wikipedia (Nx1)
    member.prob <- weighted.prob / sum.prob

    # Elements of sum.prob can go to 0 if the density is low enough.
    # This causes divide-by-zero in the member.prob calculation. Replace
    # any NaNs with 0.
    member.prob[is.nan(member.prob)] <- 0

    # we've got components across columns and observations down rows, so we do all of the summations simultaneously on both
    member.prob.sum <- colSums(member.prob)
    alpha.new <- member.prob.sum / N  # should be 1xM matrix
    #IG and gamma
    mu.new <- colSums(x_expanded * member.prob) / member.prob.sum  # should be 1xM matrix
    #IG
    aa.1=matrix(((x_expanded[,1] - mu_expanded[,1]) ^ 2 * member.prob[,1]) / (mu_expanded[,1] ^ 2 * x_expanded[,1]),nrow=N, ncol=1)
    aa.1[is.nan(aa.1)] <- 0
    lambda.new.1 <- member.prob.sum[1] / colSums(aa.1)
    #gamma
    aa.2= sum(log(x_expanded[,2]) * member.prob[,2])/member.prob.sum[2] - log(mu[2]/lambda[2])
    aa.2[is.nan(aa.2)] <- 0
    lambda.new.2 <- igamma(aa.2)

    lambda.new <- c(lambda.new.1, lambda.new.2)

    mu <- mu.new
    lambda <- lambda.new
    alpha <- alpha.new

    # LATER: it might make sense to try this algorithm using RGPU or something - it would be ideal to be able to plug in R code to your CUDA code. tHis might be a halfway option taht is fast enough
  }

  result <- list(x=x_expanded, alpha=alpha, mu=mu, lambda=lambda, member.prob=member.prob, llik=log.lik, fit_success=TRUE)
  class(result) <- "neuromixEM"
  result
}

