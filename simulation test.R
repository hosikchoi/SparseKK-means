rm(list = ls())
gc()

# setwd(r"(C:\Users\Beom\GitHub\SparseKK-means)")
setwd(r"(C:\Users\user\GitHub\SparseKK-means)")

require(sparcl)
require(kernlab)
require(caret)
require(fossil)
source(r"(.\R\subfuncs.R)")
source(r"(.\R\main.R)")
source(r"(.\R\simfuncs.R)")

n = 100
p = 2
# dat = generateMultiorange(n = n, p = p, seed = 2, with_noise = TRUE, noise_p = 5)
dat = generateTwoorange(n = n, p = p, seed = 2, with_noise = TRUE, noise_p = 5)
# dat = generateMultiMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 5, noise_sd = 3)
# dat = generateTwoMoon(each_n = n, sigma = 0.5, seed = 1, noise_p = 5, noise_sd = 3)

# sigma = kernlab::sigest(scale(dat$x), scale = FALSE)[3]
sigma = 1

# Sparse kernel k-means algorithm
tuned_skkm = tune.skkm(x = dat$x, nCluster = 2, s = NULL, ns = 10, nPerms = 25,
                      nStart = 20, kernel = "gaussian-2way", kparam = sigma, opt = TRUE)
skkm_clusters = tuned_skkm$optModel$opt_clusters
adj_skkm = adj.rand.index(dat$y, skkm_clusters)
# plot(dat$x[, 1:2], col = skkm_clusters,
#      pch = 16, cex = 1.5, 
#      xlab = "x1", ylab = "y1", main = "Proposed method")



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
