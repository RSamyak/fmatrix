#' Get weight matrix from coalescent times
#'
#' Get matrix of weights corresponding to F-matrix for a given
#' vector of coalescent times
#'
#' @param coal_times
#' @param type
#'
#' @export wt_from_times
wt_from_times <- function(coal_times, type = c("2", "1")){
  type <- match.arg(type)

  if(is.unsorted(rev(coal_times))) stop("coal_times not sorted")

  n <- length(coal_times) + 1

  ret <- matrix(0, n - 1, n - 1)

  if(type == "1"){
    appended_coal_times <- append(coal_times, 0)
    ret <- (-1)*diff(appended_coal_times) %*% t(rep(1, n-1))
    ret[upper.tri(ret)] <- 0
  }

  if(type == "2"){
    appended_coal_times <- append(coal_times, 0)

    for(i in 1:(n-1)){
      ret[i:(n-1), i] <- coal_times[i] - appended_coal_times[-(1:i)]
    }
  }


  return(ret)
}


#' Get coalescent times from weight matrix
#'
#' Get coalescent times from matrix of weights
#' corresponding to F-matrix
#'
#' @param wt
#'
#' @export times_from_wt
times_from_wt <- function(wt, type = c("2", "1")){
  type <- match.arg(type)
  if(type == "1"){
    ret <- cumsum(rev(wt[, 1]))
  }

  if(type == "2"){
    ret <- c(wt[nrow(wt), 1], wt[nrow(wt), 1] - wt[-nrow(wt), 1])
  }
  return(ret)
}
