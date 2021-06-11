#' Estimate the Frechet Mean of a list of F-matrices using Simulated Annealing
#'
#' This function estimates the Frechet Mean of a given list of F-matrices using Simulated Annealing
#'
#' @param F.list                  list of F-matrices
#' @param weights                 weights for F.list
#' @param temp.schedule           temperature schedule, one of "exp", "linear", or "log"
#' @param track                   boolean, returns tracked stuff if TRUE
#'
#' @export
sa_mean <-
  function(F.list,
           weights = NULL,
           temp.schedule = c("exp", "lin", "log"),
           chain = c("encod", "matching"),
           track = FALSE,
           init.temp = 1000,
           alpha = .9995,
           max.iter = 50e3) {

    temp.schedule <- match.arg(temp.schedule)

    chain <- match.arg(chain)

    qmat <- meanF(F.list)
    # matrix_F <- matrix_list(F.list)




    probab <- function(e1, e2, temp, eps = 0.0001){
      if(e2 < e1) return(1)
      # if(e1 == e2) cat("!")
      return(exp(-(e2-e1+eps)/temp))
    }


    n <- ncol(F.list[[1]]) + 1

    # set.seed(4564)

    if(chain == "matching"){

      energy <- function(state, matr, d="l2"){

        matr <- reduce_to_vector(matr)

        state.F <- reduce_to_vector(phylodyn:::gen_Fmat(tree_from_matching(state)))

        if(d == "l1") distances <- sum(abs(matr - state.F))/n
        if(d == "l2") distances <- sqrt(sum((matr - state.F)**2)/n)

        return(sum(distances))
      }

      proposal <- proposal_matching

      Fmat_from_config <- function(config){
        phylodyn:::gen_Fmat(tree_from_matching(config))
      }

      init <- rcoal(n, br = 1)
      config <- matching(init)

    }
    else if(chain == "encod"){

      energy <- function(state, matr, d="l2"){

        matr <- reduce_to_vector(matr)

        state.F <- reduce_to_vector(Fmat_from_myencod(state))

        if(d == "l1") distances <- sum(abs(matr - state.F))/n
        if(d == "l2") distances <- sqrt(sum((matr - state.F)**2)/n)

        return(sum(distances))
      }

      proposal <- proposal_myencod

      Fmat_from_config <- function(config){
        Fmat_from_myencod(config)
      }

      init <- rcoal(n, br = 1)
      config <- phylodyn:::gen_Fmat(init)
      config <- my_encod(config)


    } else{
      stop("Invalid chain.")
    }

    e1 <- energy(config, matr = qmat) # note sample has been set to default

    # plot(init, direction = "downwards")
    # tiplabels()
    # nodelabels

    min.config <- config
    min.energy <- e1


    temp <- init.temp
    if(is.null(alpha)) alpha <- 0.995

    if(track){
      track.config <- matrix(NA, nrow = max.iter, ncol = length(config))
      track.proposal <- matrix(NA, nrow = max.iter, ncol = length(config))
      track.energy.config <- rep(NA, max.iter)
      track.energy.proposal <- rep(NA, max.iter)
      track.temp <- rep(NA, max.iter)
      track.acceptance <- rep(NA, max.iter)
      track.prob <- rep(NA, max.iter)
    }

    cat(sprintf("Initialising chain,
                schedule = %s,
                alpha = %f,
                initial temperature = %f,
                chain = %s \n"
                , temp.schedule, alpha, init.temp, chain))

    start.time <- Sys.time()
    for(i in 1:max.iter){
      new <- proposal(config)
      e2 <- energy(new, matr = qmat)

      prob <- probab(e1, e2, temp)

      if(track){
        track.config[i,] <- config
        track.proposal[i,] <- new
        track.energy.config[i] <- e1
        track.energy.proposal[i] <- e2
        track.temp[i] <- temp
        track.prob[i] <- prob
      }

      if(runif(1) <= prob){
        config <- new
        e1 <- e2
        if(track){
          track.acceptance[i] <- 1
        }
      } else{
        if(track){
          track.acceptance[i] <- 0
        }
      }

      if(e1 < min.energy) min.config <- config

      if(track & (i > 500)) if(sum(track.acceptance[(i-100):i])==0) break

      if(temp.schedule == "exp") temp <- alpha*temp
      if(temp.schedule == "log")temp <- init.temp/(1+log(1+i))
      if(temp.schedule == "lin") temp <- init.temp/(1+i)

      if(i %% 100 == 0) cat(sprintf("\r Temperature: %f", temp))

    }

    time2 <- Sys.time() - start.time

    cat(sprintf("\n Finished! \n"))
    cat(sprintf("Time elapsed: %f seconds \n", as.double(time2, units = "secs")))


    ret <- Fmat_from_config(min.config)

    if(track){
      plot(track.temp[1:i])
      plot(track.energy.config[1:i])
      plot(track.acceptance[1:i])
    }


    if(track) return(list(ret = ret,
         track_list = list(config = track.config,
                           proposal = track.proposal,
                           energy.config = track.energy.config,
                           energy.proposal = track.energy.proposal,
                           temp = track.temp,
                           acceptance = track.acceptance,
                           prob = track.prob)
         ))
    ret
  }
