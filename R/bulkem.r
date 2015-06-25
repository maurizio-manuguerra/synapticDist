#' @useDynLib bulkem

# Fit a number of finite mixture models using the Expectation Maximisation algorithm

# Maximum likelihood estimate of model parameters for data x
invgaussMaximumLikelihood <- function(x) {
  # from http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution#Maximum_likelihood
  mu <- mean(x)
  lambda <- 1 / (1 / length(x) * sum(1 / x - 1 / mu))

  result <- list(mu=mu, lambda=lambda)

  return(result)
}

#' @export
bulkem <- function(datasets, num.components=2, max.iters=100, random.inits=1, use.gpu=TRUE, epsilon=0.000001, verbose=FALSE) {
  # TODO perhaps an interface where you pass NA to use.gpu means 'do your best'; if you pass TRUE or FALSE, require that setting

  if (use.gpu == TRUE) {
    # TODO: check input argument types
    #' @useDynLib bulkem bulkem_gpu_
       fits <- .Call(bulkem_host, datasets, num.components, max.iters, random.inits, epsilon, verbose)

       # TODO: fits$gpu <- TRUE

       for (index in 1:range(length(datasets))) {
         names(fits[[index]]) <- c('lambda', 'mu', 'alpha', 'init_lambda', 'init_mu', 'init_alpha', 'loglik', 'num_iterations', 'fit_success')
       }

       # TODO if GPU failed, fall back to CPU
       # use.gpu <- FALSE
  }

  if (use.gpu == FALSE) {
    if (verbose) {
      print("Using CPU datapath")
    }

    fits <- list()

    # for each dataset...
    # print(sprintf('len %d', length(datasets)))
    for (index in 1:length(datasets)) {
      if (verbose) {
        print(sprintf('index %d', index))
      }
      x <- datasets[[index]]

      # RANDOMISED INITIALISATION
      # For j attempts, sample a few items from the dataset. Calculate the
      # maximum likelihood parameters and use those as initial parameters for a
      # mixture component attempt.

      best.fit <- NULL

      # for j replicates
      for (j in 1:random.inits) {
        if (verbose) {
          print(paste('    attempt', j))
        }
        initials <- list()

        # come up with initial conditions
        for (m in 1:num.components) {
          # take a subset of the data
          # modified from http://stackoverflow.com/a/19861866/591483
          items <- sample(1:length(x), 3)
          randomsubset <- x[items]

          ml <- invgaussMaximumLikelihood(randomsubset)
          ml$alpha <- 1 / num.components  # start with equal mixing proportions

          initials[[m]] <- ml
        }

        # perform the fit
        fit <- invgaussmixEM(x, initials=initials, num.components=num.components, max.iters=max.iters, epsilon=epsilon)
        # print(paste0('fit llik: ', fit$llik, ', alpha: ', fit$alpha, ', lambda=', fit$lambda, ', mu=', fit$mu))

        if (is.null(best.fit)) {
          best.fit <- fit
        } else if (!is.null(fit) && (fit$llik > best.fit$llik)) {
          # print(paste0('found better ', fit$llik))
          best.fit <- fit
        }
      }

      # save results
      fits[[index]] <- best.fit
    }
  }

  return(fits)
}

#' @export
bulkem2 <- function(datasets, num.components=2, max.iters=100, random.inits=1, use.gpu=TRUE, epsilon=0.000001, verbose=FALSE) {
  # TODO perhaps an interface where you pass NA to use.gpu means 'do your best'; if you pass TRUE or FALSE, require that setting

  if (use.gpu == TRUE) {
    # TODO: check input argument types
    #' @useDynLib bulkem bulkem_gpu_
       fits <- .Call(bulkem_host, datasets, num.components, max.iters, random.inits, epsilon, verbose)

       # TODO: fits$gpu <- TRUE

       for (index in 1:range(length(datasets))) {
         names(fits[[index]]) <- c('lambda', 'mu', 'alpha', 'init_lambda', 'init_mu', 'init_alpha', 'loglik', 'num_iterations', 'fit_success')
       }

       # TODO if GPU failed, fall back to CPU
       # use.gpu <- FALSE
  }

  if (use.gpu == FALSE) {
    if (verbose) {
      print("Using CPU datapath")
    }

    fits <- list()

    # for each dataset...
    # print(sprintf('len %d', length(datasets)))
    for (index in 1:length(datasets)) {
      if (verbose) {
        print(sprintf('index %d', index))
      }
      xmat <- datasets[[index]]

      # check that x is an NxM matrix
      # FIXME: this check isn't working
      if (!is.matrix(xmat) || dim(xmat)[2] != num.components) {
        error('each dataset must be a matrix with number of columns equal to num.components')
      }

      # RANDOMISED INITIALISATION
      # For j attempts, sample a few items from the dataset. Calculate the
      # maximum likelihood parameters and use those as initial parameters for a
      # mixture component attempt.

      best.fit <- NULL

      # for j replicates
      for (j in 1:random.inits) {
        if (verbose) {
          print(paste('    attempt', j))
        }
        initials <- list()

        # come up with initial conditions
        # only sample from the corresponding column component
        for (m in 1:num.components) {
          # take a subset of the data
          # modified from http://stackoverflow.com/a/19861866/591483
          x <- xmat[,m]
          items <- sample(1:length(x), 5)
          randomsubset <- x[items]
          randomsubset <- randomsubset[randomsubset!=0]

          ml <- invgaussMaximumLikelihood(randomsubset)
          ml$alpha <- 1 / num.components  # start with equal mixing proportions

          initials[[m]] <- ml
        }
        #initials[[1]]$alpha=0.95
        #initials[[1]]$mu=200
        #initials[[1]]$lambda=800
        #initials[[2]]$alpha=0.05
        #initials[[2]]$mu=10
        #initials[[2]]$lambda=10
        # perform the fit
        fit <- invgaussmixEM(xmat, initials=initials, num.components=num.components, max.iters=max.iters, epsilon=epsilon)
        # print(paste0('fit llik: ', fit$llik, ', alpha: ', fit$alpha, ', lambda=', fit$lambda, ', mu=', fit$mu))

        if (is.null(best.fit)) {
          best.fit <- fit
        } else if (!is.null(fit) && (fit$llik > best.fit$llik)) {
          # print(paste0('found better ', fit$llik))
          best.fit <- fit
        }
      }

      # save results
      fits[[index]] <- best.fit
    }
  }

  return(fits)
}

#' InverseGaussian-Exp mixture
#' @export
bulkem3 <- function(datasets, num.components=2, max.iters=100, random.inits=1, use.gpu=TRUE, epsilon=0.000001, verbose=FALSE) {
  # TODO perhaps an interface where you pass NA to use.gpu means 'do your best'; if you pass TRUE or FALSE, require that setting

  if (use.gpu == TRUE) {
    # TODO: check input argument types
    #' @useDynLib bulkem bulkem_gpu_
       fits <- .Call(bulkem_host, datasets, num.components, max.iters, random.inits, epsilon, verbose)

       # TODO: fits$gpu <- TRUE

       for (index in 1:range(length(datasets))) {
         names(fits[[index]]) <- c('lambda', 'mu', 'alpha', 'init_lambda', 'init_mu', 'init_alpha', 'loglik', 'num_iterations', 'fit_success')
       }

       # TODO if GPU failed, fall back to CPU
       # use.gpu <- FALSE
  }

  if (use.gpu == FALSE) {
    if (verbose) {
      print("Using CPU datapath")
    }

    fits <- list()

    # for each dataset...
    # print(sprintf('len %d', length(datasets)))
    for (index in 1:length(datasets)) {
      if (verbose) {
        print(sprintf('index %d', index))
      }
      xmat <- datasets[[index]]

      # check that x is an NxM matrix
      # FIXME: this check isn't working
      if (!is.matrix(xmat) || dim(xmat)[2] != num.components) {
        error('each dataset must be a matrix with number of columns equal to num.components')
      }

      # RANDOMISED INITIALISATION
      # For j attempts, sample a few items from the dataset. Calculate the
      # maximum likelihood parameters and use those as initial parameters for a
      # mixture component attempt.

      best.fit <- NULL

      # for j replicates
      for (j in 1:random.inits) {
        if (verbose) {
          print(paste('    attempt', j))
        }
        initials <- list()

        # come up with initial conditions
        # only sample from the corresponding column component
        for (m in 1:num.components) {
          # take a subset of the data
          # modified from http://stackoverflow.com/a/19861866/591483
          x <- xmat[,m]
          items <- sample(1:length(x), 30)
          randomsubset <- x[items]
          randomsubset <- randomsubset[randomsubset!=0]

          ml <- invgaussMaximumLikelihood(randomsubset)
          ml$alpha <- 1 / num.components  # start with equal mixing proportions

          initials[[m]] <- ml
        }
        initials[[1]]$alpha=0.95
        initials[[1]]$mu=400
        initials[[1]]$lambda=400
        initials[[2]]$alpha=0.05
        initials[[2]]$mu=0.1
        initials[[2]]$lambda=0
        # perform the fit
        fit <- invgaussexpmixEM(xmat, initials=initials, num.components=num.components, max.iters=max.iters, epsilon=epsilon)
        # print(paste0('fit llik: ', fit$llik, ', alpha: ', fit$alpha, ', lambda=', fit$lambda, ', mu=', fit$mu))

        if (is.null(best.fit)) {
          best.fit <- fit
        } else if (!is.null(fit) && (fit$llik > best.fit$llik)) {
          # print(paste0('found better ', fit$llik))
          best.fit <- fit
        }
      }

      # save results
      fits[[index]] <- best.fit
    }
  }

  return(fits)
}

#' @export
bulkem2 <- function(datasets, num.components=2, max.iters=100, random.inits=1, use.gpu=TRUE, epsilon=0.000001, verbose=FALSE) {
  # TODO perhaps an interface where you pass NA to use.gpu means 'do your best'; if you pass TRUE or FALSE, require that setting

  if (use.gpu == TRUE) {
    # TODO: check input argument types
    #' @useDynLib bulkem bulkem_gpu_
       fits <- .Call(bulkem_host, datasets, num.components, max.iters, random.inits, epsilon, verbose)

       # TODO: fits$gpu <- TRUE

       for (index in 1:range(length(datasets))) {
         names(fits[[index]]) <- c('lambda', 'mu', 'alpha', 'init_lambda', 'init_mu', 'init_alpha', 'loglik', 'num_iterations', 'fit_success')
       }

       # TODO if GPU failed, fall back to CPU
       # use.gpu <- FALSE
  }

  if (use.gpu == FALSE) {
    if (verbose) {
      print("Using CPU datapath")
    }

    fits <- list()

    # for each dataset...
    # print(sprintf('len %d', length(datasets)))
    for (index in 1:length(datasets)) {
      if (verbose) {
        print(sprintf('index %d', index))
      }
      xmat <- datasets[[index]]

      # check that x is an NxM matrix
      # FIXME: this check isn't working
      if (!is.matrix(xmat) || dim(xmat)[2] != num.components) {
        error('each dataset must be a matrix with number of columns equal to num.components')
      }

      # RANDOMISED INITIALISATION
      # For j attempts, sample a few items from the dataset. Calculate the
      # maximum likelihood parameters and use those as initial parameters for a
      # mixture component attempt.

      best.fit <- NULL

      # for j replicates
      for (j in 1:random.inits) {
        if (verbose) {
          print(paste('    attempt', j))
        }
        initials <- list()

        # come up with initial conditions
        # only sample from the corresponding column component
        for (m in 1:num.components) {
          # take a subset of the data
          # modified from http://stackoverflow.com/a/19861866/591483
          x <- xmat[,m]
          items <- sample(1:length(x), 30)
          randomsubset <- x[items]
          randomsubset <- randomsubset[randomsubset!=0]

          ml <- invgaussMaximumLikelihood(randomsubset)
          ml$alpha <- 1 / num.components  # start with equal mixing proportions

          initials[[m]] <- ml
        }
        #initials[[1]]$alpha=0.95
        #initials[[1]]$mu=200
        #initials[[1]]$lambda=800
        #initials[[2]]$alpha=0.05
        #initials[[2]]$mu=10
        #initials[[2]]$lambda=10
        # perform the fit
        fit <- invgaussmixEM(xmat, initials=initials, num.components=num.components, max.iters=max.iters, epsilon=epsilon)
        # print(paste0('fit llik: ', fit$llik, ', alpha: ', fit$alpha, ', lambda=', fit$lambda, ', mu=', fit$mu))

        if (is.null(best.fit)) {
          best.fit <- fit
        } else if (!is.null(fit) && (fit$llik > best.fit$llik)) {
          # print(paste0('found better ', fit$llik))
          best.fit <- fit
        }
      }

      # save results
      fits[[index]] <- best.fit
    }
  }

  return(fits)
}

#' InverseGaussian-Gamma mixture
#' @export
bulkem4 <- function(datasets, num.components=2, max.iters=100, random.inits=1, use.gpu=TRUE, epsilon=0.000001, verbose=FALSE) {
  # TODO perhaps an interface where you pass NA to use.gpu means 'do your best'; if you pass TRUE or FALSE, require that setting

  if (use.gpu == TRUE) {
    # TODO: check input argument types
    #' @useDynLib bulkem bulkem_gpu_
       fits <- .Call(bulkem_host, datasets, num.components, max.iters, random.inits, epsilon, verbose)

       # TODO: fits$gpu <- TRUE

       for (index in 1:range(length(datasets))) {
         names(fits[[index]]) <- c('lambda', 'mu', 'alpha', 'init_lambda', 'init_mu', 'init_alpha', 'loglik', 'num_iterations', 'fit_success')
       }

       # TODO if GPU failed, fall back to CPU
       # use.gpu <- FALSE
  }

  if (use.gpu == FALSE) {
    if (verbose) {
      print("Using CPU datapath")
    }

    fits <- list()

    # for each dataset...
    # print(sprintf('len %d', length(datasets)))
    for (index in 1:length(datasets)) {
      if (verbose) {
        print(sprintf('index %d', index))
      }
      xmat <- datasets[[index]]

      # check that x is an NxM matrix
      # FIXME: this check isn't working
      if (!is.matrix(xmat) || dim(xmat)[2] != num.components) {
        error('each dataset must be a matrix with number of columns equal to num.components')
      }

      # RANDOMISED INITIALISATION
      # For j attempts, sample a few items from the dataset. Calculate the
      # maximum likelihood parameters and use those as initial parameters for a
      # mixture component attempt.

      best.fit <- NULL

      # for j replicates
      for (j in 1:random.inits) {
        if (verbose) {
          print(paste('    attempt', j))
        }
        initials <- list()

        # come up with initial conditions
        # only sample from the corresponding column component
        for (m in 1:num.components) {
          # take a subset of the data
          # modified from http://stackoverflow.com/a/19861866/591483
          x <- xmat[,m]
          items <- sample(1:length(x), 30)
          randomsubset <- x[items]
          randomsubset <- randomsubset[randomsubset!=0]

          ml <- invgaussMaximumLikelihood(randomsubset)
          ml$alpha <- 1 / num.components  # start with equal mixing proportions

          initials[[m]] <- ml
        }
        #initials[[1]]$alpha=0.95
        #initials[[1]]$mu=400
        #initials[[1]]$lambda=400
        #initials[[2]]$alpha=0.05
        #initials[[2]]$mu=5
        #initials[[2]]$lambda=2
        # perform the fit
        fit <- invgaussgammamixEM(xmat, initials=initials, num.components=num.components, max.iters=max.iters, epsilon=epsilon)
        # print(paste0('fit llik: ', fit$llik, ', alpha: ', fit$alpha, ', lambda=', fit$lambda, ', mu=', fit$mu))

        if (is.null(best.fit)) {
          best.fit <- fit
        } else if (!is.null(fit) && (fit$llik > best.fit$llik)) {
          # print(paste0('found better ', fit$llik))
          best.fit <- fit
        }
      }

      # save results
      fits[[index]] <- best.fit
    }
  }

  return(fits)
}



.onUnload <- function (libpath) {
  library.dynam.unload("bulkem", libpath)
}

