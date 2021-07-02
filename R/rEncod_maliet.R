
#' Generate a random encoding
#'
#' Generate a random encoding according to Maliet-Aldous beta splitting model
#'
#' @param n      n+1 is the number of tips of the tree
#' @param b      beta parameter in beta model of Maliet2018.
#'                 Default b = 0 gives the Yule / Kingman model
#' @param alpha  alpha parameter in Maliet2018, age-richness index
#'                 Default is alpha = 1.
#' @param eta    eta parameter in Maliet2018, abundance-richness index (Not implemented!)
#'                 Default would be eta = 1.
#' @param type   flag to switch between isochronous and hetero
#'                 chronous versions of encodings
#'
#' @export rEncod.single.Maliet2018
rEncod.single.Maliet2018 <- function(n, b = 0, alpha = 1, eta = 1,
                                 type = c("encod", "hEncod")){

  #TODO: implement b = Inf correctly
  #TODO: implement eta for nonrandom extinctions



  ### This function generates encodings
  ### as per the Blum-Francois beta model
  ### as described in Sainudiin & Veber (2016)

  type <- match.arg(type)

  #TODO: Implement b in [-2, -1)

  if(! (b >= -1)) stop("b must be greater than -2 for Maliet 2018; only implemented b >= -1")


  U <- runif(n)


  myrbeta <- function(k) {
    rbeta(k,
          shape1 = b + 1,
          shape2 = b + 1,
          ncp = 0)
  }

  eval_splitintervals <- function(x, I){
    list(c(I[1], I[1] + x*(I[2]-I[1])),
         c(I[1] + x*(I[2]-I[1]), I[2]))
  }


  # SplitPoints <- c(0, 1)

  NodeParents <- c(1)

  get_index <- function(b){max(which(SplitPoints < b))}

  encod <- c()

  Intervals <- list(c(0,1))
  Widths <- sapply(Intervals, function(I){I[2] - I[1]})
  Counts <- sapply(Intervals, function(I){sum(U >= I[1] & U <= I[2])})
  Flags <- Counts >= 2

  B <- c()

  counter <- 1
  encod <- c()

  while(TRUE){

    this_index <- fmatrix:::sample.vec((1:length(Intervals))[Flags],
           1,
           prob = Widths[Flags]**alpha)

    thisB <- myrbeta(1)
    B <- c(B, thisB)

    this_splitintervals <- eval_splitintervals(thisB, Intervals[[this_index]])

    this_Counts <- sapply(this_splitintervals, function(I){sum(U >= I[1] & U <= I[2])})



    if(all(this_Counts >= 1)) {
      counter <- counter + 1
      encod <- c(encod, NodeParents[this_index])

      parents <- c(counter, counter)

    } else{
      parents <- NodeParents[this_index]*(this_Counts > 0) + counter*(this_Counts == 0)
    }



    suffix <- 0
    if(this_index + 1 <= length(Intervals)){suffix <- (this_index + 1):length(Intervals)}
    Intervals <- c(Intervals[0:(this_index - 1)],
                   this_splitintervals,
                   Intervals[suffix])
    NodeParents <- c(NodeParents[0:(this_index - 1)],
                     parents,
                     NodeParents[suffix])


    Widths <- sapply(Intervals, function(I){I[2] - I[1]})
    Counts <- sapply(Intervals, function(I){sum(U >= I[1] & U <= I[2])})
    Flags <- Counts >= 2


    if(all(Counts <= 1)) break


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
