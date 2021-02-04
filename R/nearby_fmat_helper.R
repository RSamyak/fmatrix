#' @param current_cell
#' @param up
#' @param rd
#' @param tol
nearby_Fmat_check_neighbour_col1 <- function(current_cell, up, rd, tol){


  if(all(current_cell %in% c(up, up - 1),
         current_cell >= 0
  )) {
    flag_current_cell <- 1
    Fmat_current_cell <- current_cell

  } else{

    d.same <- abs(current_cell - up)
    d.drop <- abs(current_cell - up + 1)

    if(d.same < tol){
      Fmat_current_cell <- up
      flag_current_cell <- 1
    }
    else if(d.drop < tol){
      Fmat_current_cell <- max(up - 1, 0)
      flag_current_cell <- 1
    }
    else if(d.same < rd){
      Fmat_current_cell <- up
      flag_current_cell <- 0
    }
    else if(d.drop < rd){
      Fmat_current_cell <- max(up - 1, 0)
      flag_current_cell <- 0
    }
    else if(d.same < d.drop) {
      Fmat_current_cell <- up
      flag_current_cell <- -1
    }
    else {
      Fmat_current_cell <- max(up - 1, 0)
      flag_current_cell <- -1
    }

  }
  return(list(Fmat_current_cell = Fmat_current_cell,
              flag_current_cell = flag_current_cell))
}

#' @param current_cell
#' @param up
#' @param left
#' @param upnleft
#' @param leftflag
#' @param rd
#' @param tol
nearby_Fmat_check_neighbour <- function(current_cell, up, left, upnleft, leftflag, rd, tol){
  if(all(current_cell %in% c(up, up - 1),
         current_cell >= left,
         up - current_cell >= upnleft - left)) {

    flag_current_cell <- 1
    Fmat_current_cell <- current_cell
  }
  else {
    d.same <- abs(current_cell - up)
    d.drop <- abs(current_cell - up + 1)

    if (abs(upnleft - left - 1) < tol) {
      Fmat_current_cell <- max(up - 1, left, 0)
      flag_current_cell <- leftflag
    }

    else if (d.same < tol) {
      Fmat_current_cell <- up
      flag_current_cell <- 1
    }
    else if (d.drop < tol) {
      Fmat_current_cell <- max(up - 1, left, 0)
      flag_current_cell <- 1
    }
    else if (d.same < rd) {
      Fmat_current_cell <- up
      flag_current_cell <- 0
    }
    else if (d.drop < rd) {
      Fmat_current_cell <- max(up - 1, left, 0)
      flag_current_cell <- 0
    }
    else if (d.same < d.drop) {
      Fmat_current_cell <- up
      flag_current_cell <- -1
    }
    else {
      Fmat_current_cell <- max(up - 1, left, 0)
      flag_current_cell <- -1
    }

  }
  return(list(Fmat_current_cell = Fmat_current_cell,
              flag_current_cell = flag_current_cell))
}




