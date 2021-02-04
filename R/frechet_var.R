#' frechet variance of sample
#'
#' frechet variance of sample
#'
#' @param F.list
#' @param dist
#'
#' @export frechet_var_sample
frechet_var_sample <- function(F.list, dist = c("l1", "l2")){
  dist <- match.arg(dist)

  sample_mean_Fmat <- F.list[[brute.mean(F.list, F.list, dist = dist)]]

  distance_list <- sapply(F.list, function(fmat){
    distance_Fmat(fmat, sample_mean_Fmat, dist = dist)
  })

  return(list(var = mean(distance_list**2),
              meandist = mean(distance_list),
              fmean = sample_mean_Fmat))
}

#' frechet variance of population
#'
#' frechet variance of population
#'
#' @param F.list
#' @param probs
#' @param dist
#'
#' @export frechet_var_pop
frechet_var_pop <- function(F.list, probs, dist = c("l1", "l2")){
  dist <- match.arg(dist)

  sample_mean_Fmat <- F.list[[brute.mean.weighted(F.list, F.list, dist = dist, weights = probs)]]

  distance_list <- sapply(F.list, function(fmat){
    distance_Fmat(fmat, sample_mean_Fmat, dist = dist)
  })

  distance_list <- distance_list

  return(list(var = sum(probs * distance_list**2),
              meandist = sum(probs * distance_list),
              fmean = sample_mean_Fmat))
}

#' Entropy of a discrete distribution
#'
#' Entropy of a discrete distribution
#'
#' @param probs
#'
#' @export entropy_pop
entropy_pop <- function(probs){
  stopifnot(all(probs >= 0 & probs <= 1))
  logprobs <- ifelse(probs > 0, log(probs), 0)
  return(sum(- probs * logprobs))
}

