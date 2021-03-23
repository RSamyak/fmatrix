

find_sa_mean <- function(input,
                         preprocess = FALSE,
                         track = TRUE,
                         schedule = c("exp", "lin", "log"),
                         init = NULL,
                         max.iter = 1000, init.temp = 1000, alpha = .995){


  if(preprocess){

    cat("\nBeginning preprocessing\n")
    filename <- input$filename
    subsample.size <- ifelse(is.null(input$subsample.size), 1000, input$subsample.size)

    trees_list <- ape::read.nexus(filename)
    cat(sprintf("\nFinished reading %s\n", filename))
    set.seed(1)
    trees_sublist <- trees_list[-1][sample(1:length(trees_list[-1]), 1000)]

    Fmat_sublist <- lapply(trees_sublist, phylodyn:::gen.tr.data)
    cat("\nFinished computing F-matrices\n")
    Finfo_sublist <- lapply(Fmat_sublist, function(u) u[[2]])
    Fmat_sublist <- lapply(Fmat_sublist, function(u) u[[1]])
    Fmat_sublist <- lapply(Fmat_sublist, function(mat) mat[-1, -1])
    hEncod_sublist <- lapply(trees_sublist, gen_hEncod)
    cat("\nFinished computing Encodings\n")

    if(! length(table(sapply(Fmat_sublist, dim))) == 1) {
      ret.list <- list(
        data_list = list(
          trees_list = trees_list,
          trees_sublist = trees_sublist,
          Fmat_sublist = Fmat_sublist,
          Finfo_sublist = Finfo_sublist,
          hEncod_sublist = hEncod_sublist
        )
      )
      warning("error 01")
      ret.list$flag <- "error 01"
      return(ret.list)
    }

    # matrix_F <- matr_list(Fmat_sublist, Finfo_sublist)
  } else{
    trees_sublist <- input$trees_sublist
    # matrix_F <- input$matrix_F
    Fmat_sublist <- input$Fmat_sublist
    Finfo_sublist <- input$Finfo_sublist
    hEncod_sublist <- input$hEncod_sublist
  }

  Mmat <- meanF(Fmat_sublist)


  #########


  schedule <- match.arg(schedule)

  # dist_type <- 1

  if(is.null(init)) init <- trees_sublist[[3]]

  config <- gen_hEncod(init)

  e1 <- energy(config, Mmat = Mmat, d = "l2") # note sample has been set to default


  min.energy <- Inf

  temp <- init.temp
  alpha <- 0.995

  if(track){
    # track.config <- matrix(NA, nrow = max.iter, ncol = length(config))
    # track.proposal <- matrix(NA, nrow = max.iter, ncol = length(config))
    track.energy.config <- rep(NA, max.iter)
    track.energy.proposal <- rep(NA, max.iter)
    track.temp <- rep(NA, max.iter)
    track.acceptance <- rep(NA, max.iter)
    track.prob <- rep(NA, max.iter)
  }

  cat("Beginning chain\n")
  start.time <- Sys.time()
  for(i in 1:max.iter){
    if(i %% 10 == 0) cat(i, ".. ")
    new <- proposal_hEncod(config)
    e2 <- energy(new, Mmat)

    prob <- probab(e1, e2, temp)

    if(e2 < min.energy){
      min.energy <- e2
      min.config <- new
    }

    if(track){
      # track.config[i,] <- config
      # track.proposal[i,] <- new
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

    if(track & (i > 500)) if(sum(track.acceptance[(i-100):i])==0) break

    if(schedule == "exp") temp <- alpha*temp
    if(schedule == "log") temp <- init.temp/(1+log(1+i))
    if(schedule == "lin") temp <- init.temp/(1+i)
    # cat(i,"\n")
  }
  cat("\nFinished! \n")
  time2 <- Sys.time() - start.time
  time2


  attributes_list <- list(time = time2,
                          max.iter = max.iter,
                          final.iter = i,
                          schedule = schedule,
                          alpha = alpha,
                          init.temp = init.temp,
                          final.temp = temp,
                          init = init,
                          final.energy = e1)

  output <- tree_from_hEncod(min.config)

  ret.obj <- list(
    output = output,
    output_config = min.config,
    attributes_list = attributes_list
  )

  if(track){
    track_list <- list(
      track.energy.config = track.energy.config,
      track.energy.proposal = track.energy.proposal,
      track.temp = track.temp,
      track.prob = track.prob,
      track.acceptance = track.acceptance
    )
    ret.obj$track_list <- track_list
  }

  if(preprocess){
    data_list <- list(
      trees_list = trees_list,
      trees_sublist = trees_sublist,
      Fmat_sublist = Fmat_sublist,
      Finfo_sublist = Finfo_sublist,
      hEncod_sublist = hEncod_sublist,
      Mmat = Mmat,
      filename = filename
    )
    ret.obj$data_list <- data_list
  }

  return(ret.obj)
}
