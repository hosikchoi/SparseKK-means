skc = function(x, n_cluster, s = 1.5,
               kernel = "gaussian", kparam = 1,
               maxiter = 100, init_delta = 15, eps = 1e-5) 
{
  n = nrow(x)
  d = ncol(x)
  
  # initialization
  init_clusters = sample(1:n_cluster, size = n, replace = TRUE)
  clusters0 = init_clusters
  anovaKernel = make_anovaKernel(x, x, kernel, kparam)
  theta0 = rep(1 / sqrt(anovaKernel$numK), anovaKernel$numK)
  bcd_vec = c()
  
  for (i in 1:maxiter) {
    
    # Update clusters
    # cluster update를 한번만 하는 것이 아니라 여러번 해야함. update_Cs 함수 수정
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
    if ((sum(abs(theta - theta0)) / sum(theta0)) < eps) {
      break
    } else {
      theta0 = theta
      clusters0 = clusters
    }
    bcd_vec[i] = sum(theta * (td - GetWCD(anovaKernel, theta = theta, clusters = clusters)))
  }
  return(list(clusters = clusters, theta = theta, iteration = i, bcd = bcd_vec, init_clusters = init_clusters))
}
