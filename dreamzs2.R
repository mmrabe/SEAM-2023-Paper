DREAMzs2 <- function (bayesianSetup, settings = list(iterations = 10000,
                                         nCR = 3, gamma = NULL, eps = 0, e = 0.05, pCRupdate = FALSE,
                                         updateInterval = 10, burnin = 0, thin = 1, adaptation = 0.2,
                                         parallel = NULL, Z = NULL, ZupdateFrequency = 10, pSnooker = 0.1,
                                         DEpairs = 2, consoleUpdates = 10, startValue = NULL, currentChain = 1,
                                         message = FALSE))
{
  message("This is the stochastic sampler!")
  if ("bayesianOutput" %in% class(bayesianSetup)) {
    restart <- TRUE
  }
  else restart <- FALSE
  if (restart) {
    if (is.null(settings))
      settings <- bayesianSetup$settings
    else settings <- applySettingsDefault(settings = settings,
                                          sampler = "DREAMzs")
    settings$adaptation <- 0
  }
  else {
    settings <- applySettingsDefault(settings = settings,
                                     sampler = "DREAMzs")
  }
  if (!restart) {
    setup <- bayesianSetup
  }
  else setup <- bayesianSetup$setup
  setup <- checkBayesianSetup(setup, parallel = settings$parallel)
  if (is.null(settings$parallel))
    settings$parallel = setup$parallel
  if (!restart) {
    if (is.null(settings$startValue)) {
      parLen = length(bayesianSetup$prior$sampler(1))
      X = bayesianSetup$prior$sampler(3)
    }
    if (is.function(settings$startValue)) {
      X = settings$startValue()
    }
    if (class(settings$startValue)[1] == "numeric") {
      X = bayesianSetup$prior$sampler(settings$startValue)
    }
    if (is.matrix(settings$startValue))
      X <- settings$startValue
    if (is.null(settings$Z)) {
      parLen = length(bayesianSetup$prior$sampler(1))
      Z = bayesianSetup$prior$sampler(parLen * 10)
    }
    if (is.function(settings$Z)) {
      Z = settings$Z()
    }
    if (class(settings$Z)[1] == "numeric") {
      Z = bayesianSetup$prior$sampler(settings$Z)
    }
    if (is.matrix(settings$Z))
      Z <- settings$Z
  }
  else {
    X <- bayesianSetup$X
    Z <- bayesianSetup$Z
    if (is.vector(Z))
      Z = as.matrix(Z)
  }
  if (!is.matrix(X))
    stop("wrong starting values")
  if (!is.matrix(Z))
    stop("wrong Z values")
  FUN <- setup$posterior$density
  pCRupdate <- settings$pCRupdate
  nCR <- settings$nCR
  Npar <- ncol(X)
  Nrejected <- 0
  Nrejected_limit <- 3
  Npar12 <- (Npar - 1)/2
  parallel <- settings$parallel
  if (!is.null(parallel)) {
    if (is.numeric(parallel) | parallel == "external")
      parallel <- TRUE
  }
  else parallel <- FALSE
  pCRupdate <- settings$pCRupdate
  nCR <- settings$nCR
  Npar <- ncol(X)
  Npop <- nrow(X)
  if (settings$adaptation < 1)
    settings$adaptation <- settings$adaptation * settings$iterations
  n.iter <- ceiling(settings$iterations/Npop)
  if (n.iter < 2)
    stop("The total number of iterations must be greater than 3")
  settings$burnin <- settings$burnin/Npop
  lChain <- ceiling((n.iter - settings$burnin)/settings$thin) +
    1
  pChain <- array(NA, dim = c(lChain, Npar + 3, Npop))
  M <- nrow(Z[complete.cases(Z), , drop = FALSE])
  Zold <- Z[complete.cases(Z), , drop = FALSE]
  Z <- matrix(NA, nrow = M + floor((n.iter)/settings$ZupdateFrequency) *
                Npop, ncol = Npar)
  Z[1:M, ] <- Zold
  colnames(pChain) <- c(setup$names, "LP", "LL", "LPr")
  logfitness_X <- FUN(X, returnAll = T)
  pChain[1, , ] <- t(cbind(X, logfitness_X))
  counter <- 1
  iseq <- 1:Npop
  gamma <- 2.38/sqrt(settings$DEpairs * Npar)
  delta <- rep(0, settings$nCR)
  funevals <- 0
  if (!restart) {
    pCR = rep(1/nCR, nCR)
    lCR <- rep(0, nCR)
    CR <- matrix(1/nCR, nrow = Npop, ncol = settings$updateInterval)
  }
  else {
    pCR <- bayesianSetup$pCR
    CR <- BayesianTools:::generateCRvalues(pCR, settings, Npop)
  }
  counter_update <- 0
  omega <- numeric()
  eps <- settings$eps
  e <- settings$e
  for (iter in 2:n.iter) {
    xOld <- X
    if (parallel == TRUE) {
      x_prop <- matrix(NA, nrow = Npop, ncol = Npar)
      r_extra <- numeric(Npop)
      for (i in 1:Npop) {
        if (runif(1) > settings$pSnooker) {
          selectedChains1 <- sample((1:M), settings$DEpairs,
                                    replace = FALSE)
          selectedChains2 <- numeric(settings$DEpairs)
          for (k in 1:settings$DEpairs) {
            selectedChains2[k] <- sample((1:M)[-c(selectedChains1[k],
                                                  selectedChains2[1:k])], 1)
          }
          rn <- runif(Npar)
          indX <- which(rn > (1 - CR[i]))
          if (length(indX) == 0)
            indX <- sample(1:Npar, 1)
          x_prop[i, ] <- X[i, ]
          if (runif(1) > 4/5) {
            gamma <- 1
          }
          else {
            gamma <- 2.38/sqrt(settings$DEpairs * length(indX))
          }
          x_prop[i, indX] <- X[i, indX] + (1 + e) *
            gamma * (apply(as.matrix(Z[selectedChains1,
                                       indX]), 2, sum) - apply(as.matrix(Z[selectedChains2,
                                                                           indX]), 2, sum)) + eps * rnorm(length(indX),
                                                                                                          0, 1)
          r_extra[i] <- 0
        }
        else {
          selectSnooker <- sample((1:M), replace = FALSE,
                                  3)
          z <- Z[selectSnooker[1], ]
          x_z <- X[i, ] - z
          D2 <- max(sum(x_z * x_z), 1e-300)
          projdiff <- sum((Z[selectSnooker[1], ] - Z[selectSnooker[2],
          ]) * x_z)/D2
          gamma_snooker <- runif(1, min = 1.2, max = 2.2)
          x_prop[i, ] <- X[i, ] + gamma_snooker * projdiff *
            x_z
          x_z <- x_prop[i, ] - z
          D2prop <- max(sum(x_z * x_z), 1e-300)
          r_extra[i] <- Npar12 * (log(D2prop) - log(D2))
        }
      }
      logfitness_x_prop <- FUN(x_prop, returnAll = T)
      for (i in 1:Npop) {
        if(is.finite(logfitness_x_prop[i, 1])) {
          logfitness_X[i, ] <- FUN(xOld[i, ], returnAll = T)
        }
        if (!is.na(logfitness_x_prop[i, 1] - logfitness_X[i,
                                                          1])) {
          if ((logfitness_x_prop[i, 1] - logfitness_X[i,
                                                      1] + r_extra[i]) > log(runif(1))) {
            X[i, ] <- x_prop[i, ]
            logfitness_X[i, ] <- logfitness_x_prop[i,
            ]
            Nrejected <- 0
          } else {
            Nrejected <- Nrejected + 1
          }
        }
      }
    }
    else {
      for (i in 1:Npop) {
        if (runif(1) > settings$pSnooker) {
          selectedChains1 <- sample((1:M), settings$DEpairs,
                                    replace = FALSE)
          selectedChains2 <- numeric(settings$DEpairs)
          for (k in 1:settings$DEpairs) {
            selectedChains2[k] <- sample((1:M)[-c(selectedChains1[k],
                                                  selectedChains2[1:k])], 1)
          }
          rn <- runif(Npar)
          indX <- which(rn > (1 - CR[i]))
          if (length(indX) == 0)
            indX <- sample(1:Npar, 1)
          x_prop <- X[i, ]
          if (runif(1) > 4/5) {
            gamma <- 1
          }
          else {
            gamma <- 2.38/sqrt(settings$DEpairs * length(indX))
          }
          x_prop[indX] <- X[i, indX] + (1 + e) * gamma *
            (apply(as.matrix(Z[selectedChains1, indX]),
                   2, sum) - apply(as.matrix(Z[selectedChains2,
                                               indX]), 2, sum)) + eps * rnorm(length(indX),
                                                                              0, 1)
          r_extra <- 0
        }
        else {
          selectSnooker <- sample((1:M), replace = FALSE,
                                  3)
          z <- Z[selectSnooker[1], ]
          x_z <- X[i, ] - z
          D2 <- max(sum(x_z * x_z), 1e-300)
          projdiff <- sum((Z[selectSnooker[1], ] - Z[selectSnooker[2],
          ]) * x_z)/D2
          gamma_snooker <- runif(1, min = 1.2, max = 2.2)
          x_prop <- X[i, ] + gamma_snooker * projdiff *
            x_z
          x_z <- x_prop - z
          D2prop <- max(sum(x_z * x_z), 1e-300)
          r_extra <- Npar12 * (log(D2prop) - log(D2))
        }
        logfitness_x_prop <- FUN(x_prop, returnAll = T)
        if(Nrejected > Nrejected_limit) {
          if(is.finite(logfitness_x_prop[1])) {
            logfitness_X[i, ] <- FUN(xOld[i, ], returnAll = T)
          }
          Nrejected <- 0
        }
        if (!is.na(logfitness_x_prop[1] - logfitness_X[i,
                                                       1])) {
          if ((logfitness_x_prop[1] - logfitness_X[i,
                                                   1] + r_extra) > log(runif(1))) {
            X[i, ] <- x_prop
            logfitness_X[i, ] <- logfitness_x_prop
            Nrejected <- 0
          } else {
            Nrejected <- Nrejected + 1
          }
        }
      }
    }
    if ((iter > settings$burnin) && (iter%%settings$thin ==
                                     0)) {
      counter <- counter + 1
      pChain[counter, , ] <- t(cbind(X, logfitness_X))
    }
    if (counter%%settings$ZupdateFrequency == 0) {
      Z[(M + 1):(M + Npop), ] <- X
      M <- M + Npop
    }
    if (iter < settings$adaptation) {
      if (pCRupdate) {
        sdX <- apply(X[, 1:Npar, drop = FALSE], 2, sd)
        delta_Norm <- rowSums(((xOld - X[, 1:Npar, drop = FALSE])/sdX)^2)
        for (k in 1:settings$nCR) {
          ind <- which(abs(CR[, k] - (k/nCR)) < 1e-05)
          delta[k] <- delta[k] + sum(delta_Norm[ind])
        }
      }
      if (iter%%settings$updateInterval == 0) {
        if (pCRupdate) {
          tmp <- AdaptpCR(CR, delta, lCR, settings,
                          Npop)
          pCR <- tmp$pCR
          lCR <- tmp$lCR
        }
        for (out in 1:Npop) {
          omega[out] <- mean(pChain[((counter/2):counter),
                                    Npar + 1, out])
        }
        if (NaN %in% omega) {
          outlierChain <- NULL
        }
        else {
          IQR <- quantile(omega, probs = c(0.25, 0.75))
          outlierChain <- which(omega < IQR[1] - 2 *
                                  (IQR[2] - IQR[1]))
        }
        if (length(outlierChain) > 0) {
          best <- which.max(pChain[counter, Npar + 1,
          ])
          pChain[counter, , outlierChain] <- pChain[counter,
                                                    , best]
        }
      }
    }
    if (iter%%settings$updateInterval == 0) {
      counter_update <- 0
      CR <- BayesianTools:::generateCRvalues(pCR, settings, Npop)
    }
    if (settings$message) {
      if ((iter%%settings$consoleUpdates == 0) | (iter ==
                                                  n.iter))
        cat("\r", "Running DREAM-MCMC, chain ", settings$currentChain,
            "iteration", iter * Npop, "of", n.iter * Npop,
            ". Current logp ", logfitness_X[, 1], "Please wait!",
            "\r")
      flush.console()
    }
  }
  iterationsOld <- 0
  pChain <- pChain[1:counter, , ]
  if (restart) {
    newchains <- array(NA, dim = c((counter + nrow(bayesianSetup$chain[[1]])),
                                   (Npar + 3), Npop))
    for (i in 1:Npop) {
      for (k in 1:(Npar + 3)) {
        newchains[, k, i] <- c(bayesianSetup$chain[[i]][,
                                                        k], pChain[, k, i])
      }
    }
    pChain <- newchains
  }
  pChain <- coda::as.mcmc.list(lapply(1:Npop, function(i) coda::as.mcmc(pChain[,
                                                                               1:(Npar + 3), i])))
  list(chains = pChain, X = as.matrix(X[, 1:Npar]), Z = Z,
       pCR = pCR)
}



