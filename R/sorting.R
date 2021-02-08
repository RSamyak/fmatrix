#' @param list
#' @param `%lessthan%`
#'
#' @export bubbleSort

bubbleSort <- function(list, `%lessthan%` = `%<=%`){
  n <- length(list)

  if(n == 1){return(list)}

  `%morethan%` <- function(x, y) `%lessthan%`(y, x)

  swapped <- TRUE

  while(swapped){
    swapped <- FALSE
    for(i in 2:n){
      # If this pair is out of order
      if(list[[i-1]] %morethan% list[[i]]){
        temp <- list[[i]]
        list[[i]] <- list[[i-1]]
        list[[i-1]] <- temp
        swapped <- TRUE
      }
    }
  }
  return(list)
}
