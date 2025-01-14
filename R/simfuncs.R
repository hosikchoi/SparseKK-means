generateMultiorange = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1)
{
  set.seed(seed)
  X = matrix(nrow = n, ncol = p)
  y = numeric(n)
  nclass = numeric(3)
  each_n = floor(n / 3)
  k = 1
  while (k <= n) {
    x = rnorm(p, sd = 2)
    sx = sum(x^2)
    if (sx <= 0.5 & nclass[1] < each_n) {
      nclass[1] = nclass[1] + 1
      y[k] = 1
      X[k, ] = x
      k = k + 1
    }
    else if (2.0 < sx & sx <= 3 & nclass[2] < each_n) {
      nclass[2] = nclass[2] + 1
      y[k] = 2
      X[k, ] = x
      k = k + 1
    }
    else if (5 < sx & sx <= 7 & nclass[3] < each_n) {
      nclass[3] = nclass[3] + 1
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
  nclass = numeric(2)
  each_n = floor(n / 2)
  k = 1
  while (k <= n) {
    x = rnorm(p, sd = 2)
    sx = sum(x^2)
    if (sx <= 0.5 & nclass[1] < each_n) {
      nclass[1] = nclass[1] + 1
      y[k] = 1
      X[k, ] = x
      k = k + 1
    }else if (5 < sx & sx <= 7 & nclass[2] < each_n) {
      nclass[2] = nclass[2] + 1
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
  # c1 = cbind(5 * cos(x) - 3.5 + rnorm(each_n) * sigma, 10 * sin(x) -
  #              2.5 + rnorm(each_n) * sigma)
  # x = runif(each_n, pi, 2 * pi)
  # c2 = cbind(5 * cos(x) + 3.5 + rnorm(each_n) * sigma, 10 * sin(x) +
  #              0.5 + rnorm(each_n) * sigma)
  # x = runif(each_n, 0, pi)
  # c3 = cbind(5 * cos(x) + 10.5 + rnorm(each_n) * sigma, 10 * sin(x) -
  #              2.5 + rnorm(each_n) * sigma)
  c1 = cbind(7.5 * cos(x) - 5.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  x = runif(each_n, pi, 2 * pi)
  c2 = cbind(7.5 * cos(x) + 3.5 + rnorm(each_n) * sigma, 10 * sin(x) +
               3.5 + rnorm(each_n) * sigma)
  x = runif(each_n, 0, pi)
  c3 = cbind(7.5 * cos(x) + 12.5 + rnorm(each_n) * sigma, 10 * sin(x) -
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
  c1 = cbind(7.5 * cos(x) - 3.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  x = runif(each_n, pi, 2 * pi)
  c2 = cbind(7.5 * cos(x) + 3.5 + rnorm(each_n) * sigma, 10 * sin(x) +
               2.5 + rnorm(each_n) * sigma)
  
  X = rbind(c1, c2)
  noise_X = matrix(rnorm(2 * each_n * noise_p, 0, noise_sd), nrow = 2 * each_n, ncol = noise_p)
  X = cbind(X, noise_X)
  y = rep(c(1, 2), each = each_n)
  return(list(x = X, y = y))
}