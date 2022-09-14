rm(list = ls())
gc()
# setwd(r"(C:\Users\user\Dropbox\MyFolder\Topics\Sparse Kernel clustering\R code)")
setwd(r"(C:\Users\Beom\GitHub\SparseKK-means)")

generateMultiorange = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1)
{
  set.seed(seed)
  X = matrix(nrow = n, ncol = p)
  y = numeric(n)
  k = 1
  while (k <= n) {
    x = rnorm(p, sd = 2)
    sx = sum(x^2)
    if (sx <= 0.5) {
      y[k] = 1
      X[k, ] = x
      k = k + 1
    }
    else if (2.0 < sx & sx <= 3) {
      y[k] = 2
      X[k, ] = x
      k = k + 1
    }
    else if (5 < sx & sx <= 7) {
      y[k] = 3
      X[k, ] = x
      k = k + 1
    }
  }
  if (with_noise) {
    noise_dat = matrix(rnorm(n * noise_p, sd = 2), n, noise_p)
    X = cbind(X, noise_dat)
  }
  return(list(x = X, y = y))
}


generateTwoorange = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1)
{
  set.seed(seed)
  X = matrix(nrow = n, ncol = p)
  y = numeric(n)
  k = 1
  while (k <= n) {
    x = rnorm(p, sd = 2)
    sx = sum(x^2)
    if (sx <= 0.5) {
      y[k] = 1
      X[k, ] = x
      k = k + 1
    }else if (5 < sx & sx <= 7) {
      y[k] = 2
      X[k, ] = x
      k = k + 1
    }
  }
  if (with_noise) {
    noise_dat = matrix(rnorm(n * noise_p, sd = 2), n, noise_p)
    X = cbind(X, noise_dat)
  }
  return(list(x = X, y = y))
}


generateMultiMoon = function(each_n = 100, sigma = 1, noise_p = 4, noise_sd = 3, seed = NULL)
{
  set.seed(seed)
  x = runif(each_n, 0, pi)
  c1 = cbind(5 * cos(x) - 3.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  x = runif(each_n, pi, 2 * pi)
  c2 = cbind(5 * cos(x) + 3.5 + rnorm(each_n) * sigma, 10 * sin(x) +
               0.5 + rnorm(each_n) * sigma)
  x = runif(each_n, 0, pi)
  c3 = cbind(5 * cos(x) + 10.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  X = rbind(c1, c2, c3)
  noise_X = matrix(rnorm(3 * each_n * noise_p, 0, noise_sd), nrow = 3 * each_n, ncol = noise_p)
  X = cbind(X, noise_X)
  y = rep(c(1, 2, 3), each = each_n)
  return(list(x = X, y = y))
}

generateTwoMoon = function(each_n = 100, sigma = 1, noise_p = 4, noise_sd = 3, seed = NULL)
{
  set.seed(seed)
  x = runif(each_n, 0, pi)
  c1 = cbind(5 * cos(x) - 3.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  x = runif(each_n, pi, 2 * pi)
  c2 = cbind(5 * cos(x) + 1.5 + rnorm(each_n) * sigma, 10 * sin(x) +
               0.5 + rnorm(each_n) * sigma)
  
  X = rbind(c1, c2)
  noise_X = matrix(rnorm(2 * each_n * noise_p, 0, noise_sd), nrow = 2 * each_n, ncol = noise_p)
  X = cbind(X, noise_X)
  y = rep(c(1, 2), each = each_n)
  return(list(x = X, y = y))
}


require(sparcl)
require(kernlab)
require(caret)
source(r"(.\R\subfuncs.R)")
source(r"(.\R\main.R)")


n = 200
p = 2
# dat = generateMultiorange(n = n, p = p, seed = 2, with_noise = TRUE, noise_p = 5)
dat = generateTwoorange(n = n, p = p, seed = 2, with_noise = TRUE, noise_p = 0)
# dat = generateMultiMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 5, noise_sd = 3)
# dat = generateTwoMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 5, noise_sd = 3)

# sigma = kernlab::sigest(scale(dat$x), scale = FALSE)[3]
sigma = 0.8

# Sparse kernel k-means algorithm
tuned_res = tune.skkm(x = dat$x, nCluster = 2, s = NULL, ns = 10, nPerms = 20,
                      nStart = 1, kernel = "gaussian-2way", kparam = sigma, opt = TRUE)
tuned_res$optModel$opt_theta

tuned_res = skkm(x = dat$x, nCluster = 2, s = 3, 
                nStart = 1, kernel = "gaussian-2way", kparam = sigma, 
                 opt = TRUE, eps = 1e-9)
tuned_res$opt_theta
# tuned_res$max_bcd
# tuned_res$res[[1]]$td
tuned_res$res[[1]]$wcd

aa = make_anovaKernel(dat$x, dat$x, kernel = "gaussian-2way", kparam = sigma)
# theta = rep(1 / sqrt(3), 3)
theta = tuned_res$opt_theta
K = combine_kernel(aa, theta)
res = kkmeans2(K, centers = 2)
# res@.Data
# GetWCD(aa, theta, res@.Data, weights = rep(1, nrow(dat$x)))
sum(theta * GetWCD(aa, theta = theta, clusters = res@.Data, weights = rep(1, nrow(dat$x))))

# a = sapply(lapply(tuned_res$res, "[[", "bcd"), max)
# which(a == max(a))
# tuned_res$res[[29]]$theta
plot(dat$x[, 1:2], col = tuned_res$opt_clusters, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Proposed method")
plot(dat$x[, 1:2], col = res@.Data, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Proposed method")
plot(dat$x[, 1:2], col = tuned_res$optModel$opt_clusters, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Proposed method")


# tuned_res$opt_s
# tuned_res$optModel$opt_clusters
# tuned_res$optModel$opt_theta

# Sparse k-means algorithm
tuned_scl = KMeansSparseCluster.permute(x = dat$x, K = 2)
# opt_scl = KMeansSparseCluster(x = dat$x, K = 2, nstart = 1, wbounds = tuned_scl$bestw)
opt_scl = KMeansSparseCluster(x = dat$x, K = 2, nstart = 1, wbounds = 2.5)
opt_scl[[1]]$crit
# Kernel k-means algorithm
kkm_res = kkmeans(dat$x, centers = 2, kernel = "rbfdot", kpar = list(sigma = sigma))

par(mfrow = c(1, 4))
plot(dat$x[, 1:2], col = dat$y, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "True clusters")
plot(dat$x[, 1:2], col = tuned_res$optModel$opt_clusters, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Proposed method")
plot(dat$x[, 1:2], col = opt_scl[[1]]$Cs, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Sparse k-means clustering")
plot(dat$x[, 1:2], col = kkm_res@.Data, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Kernel k-means clustering")
