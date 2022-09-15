rm(list = ls())
gc()

# setwd(r"(C:\Users\user\Dropbox\MyFolder\Topics\Sparse Kernel clustering\R code)")
setwd(r"(C:\Users\Beom\GitHub\SparseKK-means)")

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
# dat = generateMultiMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 0, noise_sd = 3)
# dat = generateTwoMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 0, noise_sd = 3)
plot(dat$x[, 1:2], col = dat$y, xlim = c(-15, 20))


sigma = kernlab::sigest(dat$x, scale = FALSE)[3]
# sigma = kernlab::sigest(dat$x, scale = FALSE)[2]
# sigma = 0.01
# sigma = 1

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


tuned_res = tuned_res$optModel

tuned_res = skkm(x = dat$x, nCluster = 2, s = 4.2, nStart = 20,
                 kernel = "gaussian-2way", kparam = sigma, opt = TRUE)
tuned_res2 = skkm2(x = dat$x, nCluster = 2, s = 4.2, nStart = 20,
                 kernel = "gaussian-2way", kparam = sigma, opt = TRUE)
tuned_res$opt_theta
tuned_res2$opt_theta

tuned_res2$res[[1]]

tuned_res$max_bcd
tuned_res2$max_bcd


# Sparse k-means algorithm
tuned_scl = KMeansSparseCluster.permute(x = dat$x, K = 2)
# opt_scl = KMeansSparseCluster(x = dat$x, K = 2, nstart = 1, wbounds = tuned_scl$bestw)
opt_scl = KMeansSparseCluster(x = dat$x, K = 2, nstart = 1, wbounds = 2.5)
opt_scl[[1]]$crit

tuned_scl = KMeansSparseCluster(x = dat$x, K = 2, wbounds = 4.5, nstart = 1)
tuned_scl[[1]]$crit
# tuned_scl = KMeansSparseCluster.permute(x = dat$x, K = 3)
# opt_scl = KMeansSparseCluster(x = dat$x, K = 3, nstart = 1, wbounds = tuned_scl$bestw)

# Kernel k-means algorithm
kkm_res = kkmeans(dat$x, centers = 2, kernel = "rbfdot", kpar = list(sigma = sigma))
kkm_res@withinss

par(mfrow = c(1, 4))
plot(dat$x[, 1:2], col = dat$y, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "True clusters")
# plot(dat$x[, 1:2], col = tuned_res$optModel$opt_clusters, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Proposed method")
plot(dat$x[, 1:2], col = tuned_res$opt_clusters, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Proposed method")
plot(dat$x[, 1:2], col = tuned_res2$opt_clusters, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Proposed method")
# plot(dat$x[, 1:2], col = opt_scl[[1]]$Cs, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Sparse k-means clustering")
plot(dat$x[, 1:2], col = tuned_scl[[1]]$Cs, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Sparse k-means clustering")
plot(dat$x[, 1:2], col = kkm_res@.Data, pch = 16, cex = 1.5, xlab = "x1", ylab = "y1", main = "Kernel k-means clustering")
