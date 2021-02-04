

str1 <- function(n){
  2:(n-1)
}

str2 <- function(n){
  rep(2:n, each = 2)[1:(n-2)]
}

str3 <- function(n){
  rep(2:n, each = 2)[c(2:(n-2), 1)]
}


N <- 300

ret1 <- sapply(3:N, function(u) sum(abs(str1(u) - str2(u))))
ret2 <- sapply(3:N, function(u) sum(abs(str1(u) - str3(u))))
ret3 <- sapply(3:N, function(u) sum(abs(str2(u) - str3(u))))

plot(ret1)
points(ret2, col = "red")
points(ret3, col = "blue")

any(ret1 < ret2)
any(ret2 < ret1)
