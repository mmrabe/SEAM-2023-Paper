
library(BayesianTools)
library(parallel)
library(dplyr)
library(tidyr)



runs.fitting <- 4

set.seed(1234)


sim.folder <- "./FIT"


jobs <- expand.grid(
  subject = 1:61,
  model = c("swift","seam")
) %>% mutate(seed2 = as.integer(runif(n(), 1, 1000000)))



if(length(commandArgs(TRUE)) >= 1) {
  ind <- as.integer(commandArgs(TRUE)[1])
} else {
  cat("Which job?")
  ind <- scan(what = integer(), n = 1, quiet = TRUE)
}

print(jobs[ind,])


swift.dyn.dir <- file.path(toupper(jobs$model[ind]), "MCMC", "swiftstat7_r.so")




if(jobs$model[ind] == "seam") {
  likprof_ranges <- list(
    delta0 = seq(4, 15, length.out = 50),
    F = seq(0.01, 2.0, length.out = 50),
    log_d = seq(-4, 2, length.out = 50),
    log_misfac = seq(-1, 1, length.out = 50),
    log_refix = seq(-1, 1, length.out = 50),
    log_rfrac = seq(-1.5, 1.5, length.out = 50),
    mu2 = seq(0, 1, length.out = 50),
    msac = seq(1, 6, length.out = 50)
  )
} else if(jobs$model[ind] == "swift") {
  likprof_ranges <- list(
    delta0 = seq(4, 15, length.out = 50),
    log_misfac = seq(-1, 1, length.out = 50),
    log_refix = seq(-1, 1, length.out = 50),
    msac = seq(1, 6, length.out = 50)
  )
}

parnames <- names(likprof_ranges)
parrange <- vapply(likprof_ranges, range, double(2))
prior <- createUniformPrior(lower = parrange[1,], upper = parrange[2,])



load.model <- function (corpus, dir = "../DATA", parfile = file.path(dir, "swpar_default.par"), seed = runif(1, 0, 2**31)) {
  .Call("swiftr_loadmodel", dir, parfile, file.path(dir, sprintf("corpus_%s.dat", corpus)), seed)
}

load.data <- function (seqid, dir = "../DATA") .Call("swiftr_loaddata", file.path(dir, sprintf("fixseqin_%s.dat", seqid)))

transform.parameters <- function(param_value, param_name = names(param_value)) {
  for(i in seq_along(param_value)) {
    if(startsWith(param_name[i], "log_")) {
      param_name[i] <- substr(param_name[i], 5, nchar(param_name[i]))
      param_value[i] <- exp(param_value[i])
    }
  }
  names(param_value) <- param_name
  param_value
}


loglik <- function(model=1L, data=1L, threads=0L) {
  logliks <- .Call("swiftr_eval", model, data, threads)
  system(paste("echo",shQuote(paste0(Sys.getpid()," < ",paste(logliks, collapse = ", "))), ">> pars.log"))
  file.remove(sprintf("pars_%s.log", Sys.getpid()))
  logliks
}

update.parameter <- function(param_name, param_value, model) {
  stopifnot(length(param_name)==length(param_value))
  for(i in seq_along(param_name)) {
    if(startsWith(param_name[i], "log_")) {
      param_name[i] <- substr(param_name[i], 5, 100)
      param_value[i] <- exp(param_value[i])
    }
    if(is.null(.Call("swiftr_update", model, param_name[i], param_value[i]))) {
      stop(sprintf("Setting %s to %f failed.", param_name[i], param_value[i]))
    }
  }
  system(paste("echo",shQuote(paste0(Sys.getpid()," > ",paste(sprintf("%s=%g", param_name, param_value), collapse = ","))), ">> pars.log"))
}

bayesianSetup <- createBayesianSetup(
  likelihood = function(pars) {
    update.parameter(parnames, pars, swift_model)
    loglik(swift_model, swift_data)[1]
  },
  prior = prior,
  names = parnames,
  parallel = FALSE,
  catchDuplicates = FALSE
)


dyn.load(swift.dyn.dir)

source("dreamzs2.R")
assignInNamespace("DREAMzs", DREAMzs2, ns = "BayesianTools")



for(j in ind) {

  rds_file <- sprintf("%s/mcmc_%s_%d.rds", sim.folder, jobs$model[j], jobs$subject[j])

  if(file.exists(rds_file)) {
    mcmc <- readRDS(rds_file)
  } else {


    swift_model <<- load.model(sprintf("RETRO-EN-%s", jobs$model[j]), dir = "./DATA", parfile = sprintf("./DATA/swpar_%s_default.par", jobs$model[j]), seed = jobs$seed2[j])

    update.parameter("nsims", runs.fitting, swift_model)

    swift_data <<- load.data(sprintf("%d-1abcd-train", jobs$subject[j]), "./DATA")

    mcmc <- runMCMC(bayesianSetup, "DREAMzs", settings = list(iterations = 30000*5, nCR = 5, pSnooker = 0.3, adaptation = 0.2, consoleUpdates = 1, startValue = 5, parallel = FALSE))

    saveRDS(mcmc, rds_file)
  }

  summary(mcmc)
}

