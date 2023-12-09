#############################################################
## There are functions used in homework                    ##
#############################################################


#' @import stats
#' @import boot
#' @import bootstrap
#' @import Rcpp
#' @import DAAG
#' @import microbenchmark
#' @useDynLib SA23204179
NULL





# Homework1---------------------------------------------------------------------

#' @title My sample function 
#' @description Use the inverse transform method to realize sample function when replace equals True.
#' @param x the data.
#' @param size the number of samples.
#' @param prob the sample probability given to the data \code{x}.
#' @return a random sample of size \code{size}
#' @examples
#' \dontrun{
#'     n<-1000
#'     X1<-my.sample(1:5,n)
#'     X2<-my.sample(1:5,n,prob = c(1,2,3,4,5))
#'     X3<-my.sample(c(3,2,1,4,5),n,prob = c(5,4,3,2,1))
#' }
#' @export
my.sample<-function(x,size,prob=NULL){
  stopifnot(exprs = {
    is.numeric(x)
    is.finite(x)
    is.numeric(size)
  })
  if(is.null(prob)){
    prob<-rep(1,length(x))
  }else{
    stopifnot(is.numeric(prob))
  }
  if(length(x)!=length(prob)) stop("x should have the same length with prob")
  
  #main algorithm
  
  prob<-prob/sum(prob)
  prob[order(x)]<-prob
  x<-sort(x)
  cp <- cumsum(prob)
  U <- runif(size)
  r <- x[findInterval(U,cp,left.open = T)+1]
  return(r)
}


#' @title My Laplace function
#' @description Use the inverse transform method to generate a random sample from Laplace distribution
#' @param n the number of samples.
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#'     X<-my.laplace(1000)
#' }
#' @export
my.laplace<-function(n){
  #Use the inverse transform method to generate a random sample of size n from Laplace distribution
  U<-runif(n,0,1)
  mylaplace<-ifelse(U>0.5,-log(2*(1-U)),log(2*U))
  return(mylaplace)
}


#' @title My random beta function
#' @description Use the acceptance-rejection method to generate a random sample from beta distribution
#' @param n the number of samples.
#' @param a,b parameters of beta distribution
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#'     X<-my.rbeta(1000,3,2)
#' }
#' @export
my.rbeta<-function(n,a,b){
  #generate a random sample of size n from the Beta(a, b) distribution by the acceptance-rejection method
  stopifnot(a>0&b>0)
  rn<-numeric(n)
  c<-factorial(a+b-1)*(a-1)^(a-1)*(b-1)^(b-1)/factorial(a-1)/factorial(b-1)/(a+b-2)^(a+b-2)
  k<-1
  while(k<=n){
    x<-runif(1)
    u<-runif(1)
    if(u<x^(a-1)*(1-x)^(b-1)*(a+b-2)^(a+b-2)/(a-1)^(a-1)/(b-1)^(b-1)){
      rn[k]<-x
      k<-k+1
    }
  }
  return(rn)
}

#' @title My random rescaled Epanechnikov kernel function
#' @description Use the algorithm in E3.9 to generate a random sample from the rescaled Epanechnikov kernel.
#' @param n the number of samples.
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#'     X<-rEpa(1000)
#' }
#' @export
rEpa<-function(n){
  #the rescaled Epanechnikov kernel
  stopifnot(n>0)
  rn<-numeric(n)
  u1<-runif(n,-1,1)
  u2<-runif(n,-1,1)
  u3<-runif(n,-1,1)
  rn<-ifelse(abs(u3)>=abs(u2)&abs(u3)>=abs(u1),u2,u3)
  return(rn)
}


# Homework2---------------------------------------------------------------------


#' @title Estimate the value of pi
#' @description Use needle dropping methods to estimate the value of pi.
#' @param rho the ratio of l to d, shoud be a numeric between 0 and 1.
#' @return an estimated value of pi.
#' @examples
#' \dontrun{
#'     my_pi(0.5)
#'     my_pi(0.8)
#'     my_pi(1)
#' }
#' @export
my_pi<-function(rho){
  d <- 1
  l <- rho
  m <- 1e6
  X <- runif(m,0,d/2)
  Y <- runif(m,0,pi/2)
  return(2*l/d/mean(l/2*sin(Y)>X))
}


#' @title Use MC to estimate integration by simple MC and antithetic variate approach
#' @description Use the antithetic variate approach and the simple Monte Carlo method to estimate the value of theta.
#' @param K the number of repeated trials.
#' @return an list includes the results by the antithetic bariate approach and by the simple Monte Carlo method under \code{K} times.
#' @examples
#' \dontrun{
#'     result<-mc(1000)
#'     mean(result$simple)
#'     mean(result$anti)
#'     (var(result$simple)-var(result$anti))/var(result$simple)
#' }
#' @export
mc<-function(K){
simple<-numeric(K)
anti<-numeric(K)
for(i in 1:K){
  U1<-runif(1e5)
  U2<-runif(5e4)
  simple[i]<-mean(exp(U1))
  anti[i]<-mean((exp(U2)+exp(1-U2))/2)
}
return(list(simple=simple,anti=anti))
}



# Homework3---------------------------------------------------------------------

#' @title Use MC to estimate integration by importance sampling
#' @description Use the importance sampling method to estimate the value of theta.
#' @param n the number of sampling.
#' @return an list includes mean and est by two different importance funcion.
#' @examples
#' \dontrun{
#'     mc1(1e6)
#' }
#' @export
mc1<-function(n){
  g <- function(x) x^2*exp(-0.5*x^2)/sqrt(2*pi)
  f1 <- function(x) exp(-0.5*x^2)/sqrt(2*pi)/(1-pnorm(1))
  f2 <- function(x) x^2*exp(-0.5*x)/16/(1-pgamma(1,shape = 3,scale = 2))
  m <- n
  est <- sd <- numeric(2)
  u <- runif(m) #f1, inverse transform method
  x <- qnorm(u*(1-pnorm(1))+pnorm(1))
  fg <- g(x) / f1(x)
  est[1] <- mean(fg)
  sd[1] <- sd(fg)/sqrt(m)
  u <- runif(m) #f4, inverse transform method
  x <- qgamma(u*(1-pgamma(1,shape = 3,scale = 2))+pgamma(1,shape = 3,scale = 2),shape = 3,scale = 2)
  fg <- g(x) / f2(x)
  est[2] <- mean(fg)
  sd[2] <- sd(fg)/(sqrt(m))
  return(list(est=est,sd=sd))
}

#' @title Use MC to estimate integration by the stratified importance sampling
#' @description Use the stratified importance sampling method to estimate the value of theta.
#' @param n the number of sampling.
#' @return an list includes mean and est by the stratified importance sampling.
#' @examples
#' \dontrun{
#'     mc2(2000)
#' }
#' @export
mc2<-function(n){
  f<-function(j,x) (exp(-(j-1)/5)-exp(-j/5))/(1+x^2)
  F_inv<-function(j,x) -log(exp(-(j-1)/5)-x*(exp(-(j-1)/5)-exp(-j/5)))
  m<-n
  theta<-numeric(5)
  v<-numeric(5)
  set.seed(1234)
  for(j in 1:5){
    u<-runif(m)
    x<-F_inv(j,u)
    int<-f(j,x)
    theta[j]<-mean(int)
    v[j]<-var(int)
  }
  return(list(est=sum(theta),sd=sqrt(mean(v))))
}


#' @title Test the empirical Type 1 error under MC estimation.
#' @description according to E6.A.
#' @param n the number of random numbers.
#' @param k the number of trials.
#' @return an list includes empirical Type 1 error under three different distribution.
#' @examples
#' \dontrun{
#'     mc3(20,1e4)
#' }
#' @export
mc3<-function(n,k){
  n<-n
  k<-k
  result<-matrix(nrow=k,ncol=3)
  set.seed(123)
  for(i in 1:k){
    r1<-rchisq(n,df=1)
    result[i,1]<-((mean(r1)-qt(0.975,n-1)*sd(r1)/sqrt(n)>1)|(mean(r1)+qt(0.975,n-1)*sd(r1)/sqrt(n)<1))
    r2<-runif(n,0,2)
    result[i,2]<-((mean(r2)-qt(0.975,n-1)*sd(r2)/sqrt(n)>1)|(mean(r2)+qt(0.975,n-1)*sd(r2)/sqrt(n)<1))
    r3<-rexp(n)
    result[i,3]<-((mean(r3)-qt(0.975,n-1)*sd(r3)/sqrt(n)>1)|(mean(r3)+qt(0.975,n-1)*sd(r3)/sqrt(n)<1))
  }
  return(apply(result,2,mean))
}

# Homework4---------------------------------------------------------------------

#' @title Compare the effect of Bonferroni and B-H method.
#' @description the evaluating metrics are FWER, FDR and TPR.
#' @param M the number of simulations.
#' @return the results of the two methods in terms of FWER, FDR and TPR.
#' @examples
#' \dontrun{
#'     B_B.H(1000)
#' }
#' @export
B_B.H<-function(M){
  M<-M
  bonf_fwer<-numeric(M)
  bonf_fdr<-numeric(M)
  bonf_tpr<-numeric(M)
  B_H_fwer<-numeric(M)
  B_H_fdr<-numeric(M)
  B_H_tpr<-numeric(M)
  for (i in 1:M){
    H_0<-runif(950)
    H_1<-rbeta(50,0.1,1)
    p<-c(H_0,H_1)
    p.adj1<-p.adjust(p,method = "bonferroni")
    p.adj2<-p.adjust(p,method="fdr")
    bonf_fwer[i]<-sum(p.adj1[1:950]<0.1)>0
    bonf_fdr[i]<-sum(p.adj1[1:950]<0.1)/sum(p.adj1<0.1)
    bonf_tpr[i]<-mean(p.adj1[951:1000]<0.1)
    B_H_fwer[i]<-sum(p.adj2[1:950]<0.1)>0
    B_H_fdr[i]<-sum(p.adj2[1:950]<0.1)/sum(p.adj2<0.1)
    B_H_tpr[i]<-mean(p.adj2[951:1000]<0.1)
  }
  bonf<-round(c(mean(bonf_fwer),mean(bonf_fdr),mean(bonf_tpr)),4)
  B_H<-round(c(mean(B_H_fwer),mean(B_H_fdr),mean(B_H_tpr)),4)
  return(list(bonf=bonf,B_H=B_H))
}

lambda<-function(x,i) 1/mean(x[i])

#' @title My bootstrap function to estimate the parameter of exponential distribution.
#' @description My bootstrap function to estimate the parameter of exponential distribution.
#' @param n the size of sample.
#' @return compare the mean bootstrap bias and bootstrap standard error woth the thepretical ones.
#' @examples
#' \dontrun{
#'     set.seed(1234)
#'     n5<-boots(5)
#'     n10<-boots(10)
#'     n20<-boots(20)
#'     round(rbind(n5,n10,n20),4)
#' }
#' @export
my.boots<-function(n){
  m<-1000
  bias<-numeric(m)
  se<-numeric(m)
  for(i in 1:m){
    x<-rexp(n,2)
    obj<-boot(x,lambda,1000)
    bias[i]<-mean(obj$t)-obj$t0
    se[i]<-sd(obj$t)
  }
  return(c(Bias=mean(bias),Bias.Theo=2/(n-1),Se=mean(se),Se.Theo=2*n/(n-1)/sqrt(n-2)))
}



# Homework5---------------------------------------------------------------------

#' @title Jackknife estimate
#' @description Refer to Exercise 7.7. Obtain the jackknife estimates of bias and standard error of theta
#' @param data a dataframe.
#' @return an vector of length 2 including the jackknife estimates of bias and standard error.
#' @examples
#' \dontrun{
#'     attach(scor)
#'     jackn9(scor)
#'     detach(scor)
#' }
#' @export
jackn<-function(data){
n<-nrow(data)
#首先计算theta
decom<-eigen(cov(data))
theta<-decom$values[1]/sum(decom$values)
jack_theta<-numeric(n)
for(i in 1:n){
  decom<-eigen(cov(data[-i,]))
  jack_theta[i]<-decom$values[1]/sum(decom$values)
}
return(c(jack.bias=(n-1)*(mean(jack_theta)-theta),jack.se=sqrt((n-1)^2/n*var(jack_theta))))
}




# Homework6---------------------------------------------------------------------



#' @title Compute W-2 statistics 
#' @description Implement the two-sample Cram`er-von Mises test for equal distributions as a permutation test.
#' @param z data
#' @param i index, a vector.
#' @param sizes sizes.
#' @return a numeric.
#' @examples
#' \dontrun{
#'     attach(chickwts)
#'     x <- as.vector(weight[feed == "soybean"]);y <- as.vector(weight[feed == "linseed"])
#'     detach(chickwts)
#'     z <- c(x, y)
#'     N <- c(length(x), length(y))
#'     boot.obj <- boot(data = z, statistic = W2, R = 9999,sim = "permutation", sizes = N)
#'     ts <- c(boot.obj$t0,boot.obj$t)
#'     p.value <- mean(ts>=ts[1])
#' }
#' @export
W2<-function(z,i,sizes){
  jieguo<-0
  z<-z[i]
  n<-sizes[1]
  m<-sizes[2]
  x0<-z[1:n]
  y0<-z[-(1:n)]
  Fn<-function(x) sum(x0<=x)/n
  Gm<-function(x) sum(y0<=x)/m
  for(j in 1:(m+n)) jieguo<-jieguo+(Fn(z[j])-Gm(z[j]))^2
  return(m*n/(m+n)^2*jieguo)
}

count5<-function(z,i,sizes){
  z <- z[i]
  x <- z[1:sizes[1]]
  y <- z[-(1:sizes[1])]
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx,outy)))
}

#' @title  The Count 5 test for equal variances
#' @description The Count 5 test for equal variances in Section 6.4 is based on the maximum number of extreme points. Example 6.15 shows that the Count 5 criterion is not applicable for unequal sample sizes. Implement a permutation test for equal variance based on the maximum number of extreme points that applies when sample sizes are not necessarily equal.
#' @param x data.
#' @param y data.
#' @return logical number.
#' @examples
#' \dontrun{
#'     n1 <- 20
#'     n2 <- 30
#'     mu1 <- mu2 <- 0
#'     sigma1<-sigma2 <- 1
#'     m <- 1000
#'     set.seed(1234)
#'     alphahat1 <- mean(replicate(m, expr={
#'       x <- rnorm(n1, mu1, sigma1)
#'       y <- rnorm(n2, mu2, sigma2)
#'       perm_count5test(x, y)
#'     }))
#'     alphahat1
#' }
#' @export
perm_count5test <- function(x, y) {
  z <- c(x,y)
  N <- c(length(x), length(y))
  boot.obj <- boot(data = z, statistic = count5, R = 500,
                   sim = "permutation", sizes = N)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  # return 1 (reject) or 0 (do not reject H0)
  return(p.value<0.05)
}


# Homework7---------------------------------------------------------------------

#' @title Estimate a in Logistic model 
#' @description a function that takes as input values N,b1,b2,b3 and f0, and produces the output a.
#' @param N the number of random numbers.
#' @param b1,b2,b3 the parameters in logistic model.
#' @param f0 a parameter.
#' @return a numeric.
#' @examples
#' \dontrun{
#'     solve_alpha(1e6,0,1,-1,0.1)
#' }
#' @export
solve_alpha<-function(N,b1,b2,b3,f0){
  x1 <- rpois(N,1)
  x2 <- rexp(N,1)
  x3 <- rbinom(N,1,0.5)
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
  }
  return(uniroot(g,c(-100,100))$root)
}

#' @title Random walk Metropolis sampler for generating the standard Laplace distribution 
#' @description a function that generates the standard Laplace distribution random numbers by random walk Metropolis sampling.
#' @param sigma standard variance of proposal distribution N(x_t,sigma).
#' @param x0 initial value.
#' @param N size of random numbers required.
#' @return a list including a vector of random numbers and repetition times.
#' @examples
#' \dontrun{
#'     rw.Metropolis(.5, 100, 10000)
#' }
#' @export
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(-abs(y)) /exp(-abs(x[i-1])) ){
      x[i] <- y  
      k <- k + 1
    }
    else {
      x[i] <- x[i-1]
    }
  }
  return(list(x=x, k=k))
}


#' @title  Gibbs sampler to generate a bivariate normal chain
#' @description Implement a Gibbs sampler to generate a bivariate normal chain with zero means, unit standard deviations, and correlation 0.9.
#' @param N the length of chain.
#' @param burn burn-in length.
#' @return a random vector with the burn-in period removed.
#' @examples
#' \dontrun{
#'     my.gibbs(15000,5000)
#' }
#' @export
my.gibbs<-function(N,burn){
  N <- N #length of chain
  burn <- burn #burn-in length
  X <- matrix(0, N, 2) #the chain, a bivariate sample
  rho <- 0.9 #correlation
  mu1 <- 0
  mu2 <- 0
  sigma1 <- 1
  sigma2 <- 1
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  ###### generate the chain #####
  set.seed(1234)
  X[1, ] <- c(mu1, mu2) #initialize
  for (i in 2:N) {
    x2 <- X[i-1, 2]
    m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X[i, 1] <- rnorm(1, m1, s1)
    x1 <- X[i, 1]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] <- rnorm(1, m2, s2)
  }
  b <- burn + 1
  x <- X[b:N, ]
  return(x)
}

#' @title   Gelman-Rubin method to monitor convergence of the chain
#' @description Refer to Example 9.1. Use the Gelman-Rubin method to monitor convergence of the chain, and run the chain until the chain has converged approximately to the target distribution.
#' @param psi psi[i,j] is the statistic psi(X[i,1:j]) for chain in i-th row of X.
#' @return hat_r
#' @export
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

# Homework8---------------------------------------------------------------------

mle<-function(x,u,v){
  sum((-u*exp(-x*u)+v*exp(-x*v))/(exp(-x*u)-exp(-x*v)))
}



#' @title  An EM algorithm to estimate the lambda of exponential distribution.
#' @description An EM algorithm to estimate the lambda of exponential distribution.
#' @param x0 initial value for the iterative procedure.
#' @param u,v the known interval.
#' @param max.it the maximum iterative times, the default is 500.
#' @param eps the threshold condition for the end of the iteration, the default is 1e-7.
#' @return an estimator for lambda by EM algorithm.
#' @examples
#' \dontrun{
#'     u<-c(11,8,27,13,16,0,23,10,24,2)
#'     v<-c(12,9,28,14,17,1,24,11,25,3)
#'     em(1,u,v)
#' }
#' @export
em<-function(x0,u,v,max.it=500,eps=1e-7){
  n<-length(u)
  x<-c(x0,-100)
  for(i in 1:max.it){#EM算法迭代过程
    x[2]<-n/(n/x[1]-mle(x[1],u,v))
    if(abs(x[1]-x[2])<eps) break
    x[1]<-x[2]
  }
  return(list(lambda=x[1],error=abs(x[1]-x[2])))
}


#' @title  The Morra game's solution
#' @description In the Morra game, the set of optimal strategies are not changed if a constant is subtracted from every entry of the payoff matrix, or a positive constant is multiplied times every entry of the payoff matrix. However, the simplex algorithm may terminate at a different basic feasible point (also optimal).
#' @param A a matrix.
#' @return a list to store the solution of the Morra game.
#' @examples
#' \dontrun{
#'     A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
#'     2,0,0,0,-3,-3,4,0,0,
#'     2,0,0,3,0,0,0,-4,-4,
#'     -3,0,-3,0,4,0,0,5,0,
#'     0,3,0,-4,0,-4,0,5,0,
#'     0,3,0,0,4,0,-5,0,-5,
#'     -4,-4,0,0,0,5,0,0,6,
#'     0,0,4,-5,-5,0,0,0,6,
#'     0,0,4,0,0,5,-6,-6,0), 9, 9)
#'     result <- solve_game(A)
#' }
#' @export
solve_game <- function(A) {
  #solve the two player zero-sum game by simplex method
  #optimize for player 1, then player 2
  
  min.A <- min(A)
  A <- A - min.A #so that v >= 0
  max.A <- max(A)
  A <- A / max(A)
  m <- nrow(A)
  n <- ncol(A)
  it <- n^3
  a <- c(rep(0, m), 1) #objective function
  A1 <- -cbind(t(A), rep(-1, n)) 
  b1 <- rep(0, n)
  A3 <- t(as.matrix(c(rep(1, m), 0))) 
  b3 <- 1
  sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
                maxi=TRUE, n.iter=it)
  
  #minimize v subject to ...
  #let y strategies 1:n, with v as extra variable
  a <- c(rep(0, n), 1) #objective function
  A1 <- cbind(A, rep(-1, m)) #constraints <=
  b1 <- rep(0, m)
  A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
  b3 <- 1
  sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
                maxi=FALSE, n.iter=it)
  soln <- list("A" = A * max.A + min.A,
               "x" = sx$soln[1:m],
               "y" = sy$soln[1:n],
               "v" = sx$soln[m+1] * max.A + min.A)
  return(soln)
}



# Homework9---------------------------------------------------------------------

#' @title  The Gibbs sampler for the binomial-beta distribution.
#' @description The target bivariate density can be shown that for fixed a,b,n, the conditional distributions are B(n,y) and Beta(x+a,n−x+b).
#' @param num number of observations.
#' @param a,b,n the parameters of the density, keep a close watch on the description for details.
#' @return a vector of target random numbers.
#' @examples
#' \dontrun{
#'     r_gibbs(1000,2,3,10)
#' }
#' @export
r_gibbs<-function(num,a,b,n){
  num<-num+500
  random<-matrix(nrow = num,ncol=2)
  random[1,]<-c(0,0)
  for(i in 2:num){
    random[i,1]<-rbinom(1,n,random[i-1,2])
    random[i,2]<-rbeta(1,random[i,1]+a,n-random[i,1]+b)
  }
  return(random[-c(1:500),])
}


