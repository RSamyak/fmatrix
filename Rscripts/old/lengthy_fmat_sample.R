n <- 8

times <- rev(sort(sample(1:50, n-1)))

wt <- wt_from_times(times)
times_from_wt(wt)
lengthy_Fmat <- list(f = mat,
                     w = wt)

lengthy_Fmat

View(distance_Fmat)
