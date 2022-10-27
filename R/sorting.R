#' @param list
#' @param `%lessthan%`
#'
#' @export bubbleSort

bubbleSort <- function(list, `%lessthan%` = `%<=%`){
  n <- length(list)

  order <- 1:n

  if(n == 1){return(1)}

  `%morethan%` <- function(x, y) `%lessthan%`(y, x)

  swapped <- TRUE

  while(swapped){
    swapped <- FALSE
    for(i in 2:n){
      # If this pair is out of order
      if(list[[ order[i-1] ]] %morethan% list[[ order[i] ]]){
        temp <- order[i]
        order[i] <- order[i-1]
        order[i-1] <- temp
        swapped <- TRUE
      }
    }
  }
  return(order)
}
