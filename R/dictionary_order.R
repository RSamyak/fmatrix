#' @export `%<=%`
`%<=%` <- function(x, y){
  n <- min(length(x), length(y))

  for(i in 1:n){
    if(x[i] > y[i]) return(FALSE)
    if(x[i] < y[i]) return(TRUE)
  }

  if(length(y) < length(x)) return(FALSE)
  return(TRUE)
}

#' @export `%>=%`
`%>=%` <- function(x, y){
  `%<=%`(y, x)
}
