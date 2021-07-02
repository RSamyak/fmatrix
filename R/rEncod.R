#' Generate a random encoding
#'
#' Generate a random encoding according to Blum-Francoid beta splitting model
#'
#' @param n      n+1 is the number of tips of the tree
#' @param b      beta parameter in BF model of SV2016.
#'                 Default b = 0 gives the Yule / Kingman model
#' @param a      alpha parameter in BF model of SV 2016.
#'                 Default is a = b for BF model.
#' @param type   flag to switch between isochronous and hetero
#'                 chronous versions of encodings
#'
#' @export rEncod.single.SV2016
rEncod.single.SV2016 <- function(n, b = 0, a = b,
                   type = c("encod", "hEncod")){

  #TODO: implement b = Inf correctly

  ### This function generates encodings
  ### as per the Blum-Francois beta model
  ### as described in Sainudiin & Veber (2016)

  type <- match.arg(type)

  if(! (b >= -1)) stop("b must be greater than -1")

  B <- rbeta(n,
             shape1 = a + 1,
             shape2 = b + 1,
             ncp = 0)
  U <- runif(n)


  eval_splitpoint <- function(x, i1, i2){
    i1 + x*(i2-i1)
  }


  SplitPoints <- c(0, 1)
  NodeParents <- c(1)

  get_index <- function(b){max(which(SplitPoints < b))}

  encod <- c()

  # Intervals <- list(c(0,1))

  for(i in 1:(n-1)){

    this_index <- get_index(U[i])

    encod <- c(encod, NodeParents[this_index])

    this_splitpoint <- eval_splitpoint(B[i],
                    SplitPoints[this_index],
                    SplitPoints[this_index + 1])

    SplitPoints <- c(SplitPoints[0:this_index],
                     this_splitpoint,
                     SplitPoints[(this_index + 1):length(SplitPoints)])

    suffix <- 0
    if(this_index + 1 <= length(NodeParents)) suffix <- (this_index + 1):length(NodeParents)

    NodeParents <- c(NodeParents[0:(this_index - 1)],
                    i + 1, i + 1,
                    NodeParents[suffix])

  }

  ret <- NA
  if(type == "encod"){
    ret <- encod[-1]
  }
  if(type == "hEncod"){
    ret <- fmatrix:::preprocess_hEncod(encod)
  }

  return(ret)
}


#' Generate random encodings
#'
#' Generate random encodings according to either BF model or uniform
#'
#' @param m      number of trees to generate
#' @param n      n+1 is the number of tips of the tree
#' @param b      beta parameter in BF model of SV2016.
#'                 Default b = 0 gives the Yule / Kingman model
#' @param a      alpha parameter in BF model of SV 2016.
#'                 Default is a = b for BF model.
#' @param type   flag to switch between isochronous and hetero
#'                 chronous versions of encodings
#'
#' @export rEncod
rEncod <- function(m, n,
                   distr = c("BF", "uniform"),
                   b = 0, a = b,
                   type = c("encod", "hEncod"),
                   ...){
  distr <- match.arg(distr)
  type <- match.arg(type)

  if(distr == "BF"){
    ret <- lapply(1:m, function(i){rEncod.single.SV2016(n = n,
                                        a = a,
                                        b = b,
                                        type = type, ...) })
  } else if(distr == "uniform"){
    ret <- lapply(1:m, function(i){rEncod.single.uniform(n = n,
                                                  type = type)})
  }
  return(ret)
}

