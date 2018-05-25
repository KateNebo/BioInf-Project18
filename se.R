se.bm <- function(x, g = function(x) x) {
  n <- length(x)
  b <- floor(sqrt(n))
  a <- floor(n / b)
  
  y <- sapply(1:a, function(k) mean(g(x[((k - 1) * b + 1):(k * b)])))
  gx <- g(x)
  mu.hat <- mean(gx)
  var.hat <- b * sum((y - mu.hat)^2) / (a - 1)
  se <- sqrt(var.hat / n)
  
  list(mu = mu.hat, se.mean = se, var = var.hat)
}

se.obm <- function(x, g = function(x) x) {
  n <- length(x)
  b <- floor(sqrt(n))
  a <- n - b + 1
  
  gx <- g(x)
  cgx <- c(0, cumsum(gx))
  
  y <- diff(cgx, lag = b) / b
  mu.hat <- mean(gx)
  var.hat <- n * b * sum((y - mu.hat)^2) / (a - 1) / a
  se <- sqrt(var.hat / n)
  lambda <- var(gx)
  
  list(mu = mu.hat, se.mean = se, var = var.hat, lambda = lambda)
}


se.tukey <- function(x, g = function(x) x) {
  n <- length(x)
  b <- floor(sqrt(n))
  a <- floor(n / b)
  
  alpha <- 1:b
  alpha <- (1 + cos(pi * alpha / b)) / 2 * (1 - alpha / n)
  gx <- g(x)
  mu.hat <- mean(gx)
  # FIXME: FFT
  R <- sapply(0:b, function(j) mean((gx[1:(n - j)] - mu.hat) * (gx[(j + 1):n] - mu.hat)))
  var.hat <- R[1] + 2 * sum(alpha * R[-1])
  se <- sqrt(var.hat / n)
  
  list(mu = mu.hat, se.mean = se, var = var.hat)
}


se.bartlett <- function(x, g = function(x) x) {
  n <- length(x)
  b <- floor(sqrt(n))
  a <- floor(n / b)
  
  alpha <- 1:b
  alpha = (1 - abs(alpha) / b) * (1 - alpha / n)
  gx <- g(x)
  mu.hat <- mean(gx)
  # FIXME: FFT
  R <- sapply(0:b, function(j) mean((gx[1:(n - j)] - mu.hat) * (gx[(j + 1):n] - mu.hat)))
  var.hat <- R[1] + 2 * sum(alpha * R[-1])
  se <- sqrt(var.hat / n)
  lambda <- var(gx)

  list(mu = mu.hat, se.mean = se, var = var.hat, lambda = lambda)
}
