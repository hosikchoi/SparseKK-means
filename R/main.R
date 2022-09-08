tune.skkm = function(x, nCluster, nPerms = 20, s = NULL, ns = 100, nStart = 10, weights = NULL, 
                     kernel = "linear", kparam = 1, opt = TRUE, ...)
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel, c("linear", "gaussian", "spline-t",
                               "gaussian-2way", "spline-t-2way"))
  
  if (!is.matrix(x)) {
    x = as.matrix(x)
  }
  
  p = ncol(x)
  
  if (is.null(s)) {
    if (kernel %in% c("gaussian-2way", "spline-t-2way")) {
      nv = p + p * (p - 1) / 2
    } else {
      nv = p
    }
    s = exp(seq(log(1), log(sqrt(nv)), length.out = ns))
  }
  
  perm_list = vector("list", nPerms)
  for (i in 1:nPerms) {
    perm_list[[i]] = sapply(1:p, function(j) sample(x[, j]))
  }
    
  org_bcd = numeric(ns)
  for (j in 1:ns) {
    org_fit = skkm(x, nCluster = nCluster, nStart = nStart, s = s[j], weights = weights,
                   kernel = kernel, kparam = kparam, opt = TRUE, ...)
    org_bcd[j] = org_fit$max_bcd
  }
    
  perm_bcd_list = matrix(0, nrow = nPerms, ncol = ns)
  for (b in 1:nPerms) {
    perm_bcd = numeric(ns)
    for (j in 1:ns) {
      perm_fit = skkm(x = perm_list[[b]], nCluster = nCluster, nStart = nStart, s = s[j], weights = weights,
                      kernel = kernel, kparam = kparam, opt = TRUE, ...)
      perm_bcd[j] = perm_fit$max_bcd
    }
    perm_bcd_list[b, ] = perm_bcd
  }
    
  out$org_bcd = org_bcd
  out$perm_bcd = perm_bcd_list
  out$gaps = log(org_bcd) - colMeans(log(perm_bcd_list))
  out$opt_ind = which.max(out$gaps)
  out$opt_s = s[out$opt_ind]
    
  if (opt) {
    opt_fit = skkm(x = x, nCluster = nCluster, nStart = nStart, s = out$opt_s, weights = weights,
                  kernel = kernel, kparam = kparam, opt = TRUE, ...)  
    out$optModel = opt_fit
  }
  
  out$call = call
  return(out)
}

skkm = function(x, nCluster, nStart = 10, s = 1.5, weights = NULL,
               kernel = "linear", kparam = 1, opt = TRUE, ...) 
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel, c("linear", "gaussian", "spline-t",
                               "gaussian-2way", "spline-t-2way"))
  
  x = as.matrix(x)
  n = nrow(x)
  d = ncol(x)
  
  if (is.null(weights)) {
    weights = rep(1, n)
    attr(weights, "type") = "auto"
  }
  
  res = vector("list", length = nStart)
  
  seeds = seq(1, nStart, by = 1)
  for (j in 1:length(seeds)) {
    # initialization
    # set.seed(seeds[j])
    clusters0 = sample(1:nCluster, size = n, replace = TRUE)
    res[[j]] = skkm_core(x = x, clusters0 = clusters0, theta0 = NULL, s = s, weights = weights,
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

skkm_core = function(x, clusters0 = NULL, theta0 = NULL, s = 1.5, weights = NULL,
               kernel = "linear", kparam = 1, maxiter = 100, eps = 1e-5) 
{
  call = match.call()
  n = nrow(x)
  d = ncol(x)
  
  if (is.null(weights)) {
    weights = rep(1, n)
  }
  
  # initialization
  init_clusters = clusters0
  anovaKernel = make_anovaKernel(x = x, y = x, kernel = kernel, kparam = kparam)
  
  if (is.null(theta0)) {
    theta0 = rep(1 / sqrt(anovaKernel$numK), anovaKernel$numK)
  }
  
  td_vec = wcd_vec = bcd_vec = c()
  
  for (i in 1:maxiter) {
    
    # Update clusters
    # clusters = kkk(combine_kernel(anovaKernel, theta = theta0), 2)
    # clusters0 = sample(1:nCluster, size = n, replace = TRUE)
    
    if (attr(weights, "type") == "auto") {
      weights0 = weights
      NbyC = table(clusters0)
      weights = 1 / as.numeric(NbyC[clusters0])
      attr(weights, "type") = attr(weights0, "type")
    }
    
    clusters = updateCs(anovaKernel = anovaKernel, theta = theta0, 
                        clusters = clusters0, weights = weights)$clusters
    
    # plot(dat$x[, 1:2], col = clusters)
    # RKHS_d2 = RKHS_dist2(x, theta = theta0, g = clusters0, kernel = kernel, kparam = kparam)
    # clusters = apply(RKHS_d, MARGIN = 1, which.min)
    
    # Update theta
    wcd = GetWCD(anovaKernel, theta = theta0, clusters = clusters, weights = weights)
    td = GetWCD(anovaKernel, theta = theta0, rep(1, length(clusters)), weights = weights)
    bcd = td - wcd
    
    delta = BinarySearch(coefs = bcd, s = s)
    
    theta_tmp = soft_threshold(bcd, delta = delta)
    theta = normalization(theta_tmp)
    
    td_new = GetWCD(anovaKernel, theta = theta, rep(1, length(clusters)), weights = weights)
    wcd_new = GetWCD(anovaKernel, theta = theta, clusters = clusters, weights = weights)
    
    td_vec[i] = sum(theta * td_new)
    wcd_vec[i] = sum(theta * wcd_new)
    bcd_vec[i] = sum(theta * (td_new - wcd_new))
    
    # print((sum(abs(theta - theta0)) / sum(theta0)))
    if ((sum(abs(theta - theta0)) / sum(theta0)) < eps) {
      break
    } else {
      theta0 = theta
      clusters0 = clusters
    }
  }
  return(list(clusters = clusters, theta = theta, weights = weights, 
              iteration = i, td = td_vec, wcd = wcd_vec, bcd = bcd_vec, init_clusters = init_clusters))
}
