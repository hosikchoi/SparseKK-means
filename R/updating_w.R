S <- 1.1
a <- c(1.9628, 8.4325, 
       9.2427, 0.8907, 3.0255, 2.2164, 0.8837, 2.6337, 1.3604, 2.1935)

S <- 1.5
a <- c(2.0715, 8.3090, 
       9.2584, 1.5142, 2.9718, 2.1594, 0.3180, 1.6977, 1.4135, 2.2254)

sa <- sort(a,decreasing = TRUE)

p <- length(a)  
lam1 <- rep(Inf,p)
l1norm <- rep(Inf,p)
l2norm <- rep(Inf,p)
cnt = 0 
tot <- c()
for(l in p:1){
  cnt <- cnt + 1
  t1 <- sum(sa[1:l])
  t2 <- sum(sa[1:l]^2)

  bb <- t1^2/l-(t1^2-t2*S^2)/(l-S^2)
  b <- t1*(l-S^2)
  aval <- l*(l-S^2)
  cval <- t1^2-t2*S^2
  dett <- b^2- aval*cval
  #print(c(bb, dett))
  if(bb<0){
    next
  }
  lam1[cnt] <- t1/l - 1/sqrt(l)*sqrt(bb)
  w_temp = pmax(a-lam1[cnt],0)
  w <- w_temp/sqrt(sum(w_temp^2))
  l1norm[cnt] <- sum(w)
  l2norm[cnt] <- sum(w^2)
  tot <- rbind(tot, c(lam1[cnt], l2norm[cnt], l1norm[cnt], round(w,4)))
}
tot <- data.frame(tot)
colnames(tot)[1:3] <- c("lambda", "l2", "l1")
colnames(tot)[4:(3+p)] <- paste0("x", seq(1,p))
print(tot)

