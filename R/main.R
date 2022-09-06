skkm = function(x, nCluster, nStart = 10, s = 1.5,
               kernel = c("linear", "gaussian", "spline-t",
                          "gaussian-2way", "spline-t-2way"),
               kparam = 1, opt = TRUE, ...) 
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel)
  
  x = as.matrix(x)
  n = nrow(x)
  d = ncol(x)
  
  res = vector("list", length = nStart)
  
  seeds = seq(1, nStart, by = 1)
  for (j in 1:length(seeds)) {
    # initialization
    set.seed(seeds[j])
    clusters0 = sample(1:nCluster, size = n, replace = TRUE)
    res[[j]] = skkm_core(x = x, clusters0 = clusters0, theta0 = NULL, s = s,
                         kernel = kernel, kparam = kparam, ...)
  }
  if (opt) {
    bcd_list = sapply(res, function(x) {
      bcd = x$bcd[x$iteration]
    })
    
    opt_ind = which(bcd_list == max(bcd_list))
    if (length(opt_ind) > 1) {
      warning("")
      opt_ind = opt_ind[1]
    }
    
    out$opt_clusters = res[[opt_ind]]$clusters
    out$opt_theta = res[[opt_ind]]$theta
    out$max_bcd = bcd_list[opt_ind]
  }
  out$res = res
  return(out)
}

skkm_core = function(x, clusters0 = NULL, theta0 = NULL, s = 1.5,
               kernel = "linear", kparam = 1, maxiter = 100, eps = 1e-5) 
{
  call = match.call()
  n = nrow(x)
  d = ncol(x)
  
  # initialization
  init_clusters = clusters0
  anovaKernel = make_anovaKernel(x, x, kernel, kparam)
  if (is.null(theta0)) {
    theta0 = rep(1 / sqrt(anovaKernel$numK), anovaKernel$numK)
  }
  bcd_vec = c()
  
  for (i in 1:maxiter) {
    
    # Update clusters
    clusters = updateCs(anovaKernel = anovaKernel, theta = theta0, clusters = clusters0)$clusters
    # RKHS_d2 = RKHS_dist2(x, theta = theta0, g = clusters0, kernel = kernel, kparam = kparam)
    # clusters = apply(RKHS_d, MARGIN = 1, which.min)
    
    # Update theta
    wcd = GetWCD(anovaKernel, theta = theta0, clusters)
    td = GetWCD(anovaKernel, theta = theta0, rep(1, length(clusters)))
    bcd = td - wcd
    
    delta = BinarySearch(coefs = bcd, s = s)
    
    theta_tmp = soft_threshold(bcd, delta = delta)
    theta = normalization(theta_tmp)
    
    # print((sum(abs(theta - theta0)) / sum(theta0)))
    bcd_vec[i] = sum(theta * (td - GetWCD(anovaKernel, theta = theta, clusters = clusters)))
    
    if ((sum(abs(theta - theta0)) / sum(theta0)) < eps) {
      break
    } else {
      theta0 = theta
      clusters0 = clusters
    }
  }
  return(list(clusters = clusters, theta = theta, iteration = i, bcd = bcd_vec, init_clusters = init_clusters))
}
