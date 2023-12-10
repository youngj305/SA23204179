## -----------------------------------------------------------------------------
library(SA23204179)
n<-1000
X1<-my.sample(1:5,n)
X2<-my.sample(1:5,n,prob = c(1,2,3,4,5))
X3<-my.sample(c(3,2,1,4,5),n,prob = c(5,4,3,2,1))
X4<-my.laplace(1000)
X5<-my.rbeta(1000,3,2)
X6<-rEpa(1000)

## -----------------------------------------------------------------------------
my_pi(0.5)
my_pi(0.8)
my_pi(1)
result<-mc(1000)
mean(result$simple)
mean(result$anti)
(var(result$simple)-var(result$anti))/var(result$simple)

## -----------------------------------------------------------------------------
a<-mc1(1e4)
b<-mc2(2000)
c<-mc3(20,1e4)

## -----------------------------------------------------------------------------
a1<-B_B.H(1000)
set.seed(1234)
n5<-my.boots(5)
n10<-my.boots(10)
n20<-my.boots(20)
round(rbind(n5,n10,n20),4)

## -----------------------------------------------------------------------------
library(bootstrap)
jackn(scor)

## -----------------------------------------------------------------------------
library(boot)
set.seed(1234)
attach(chickwts)
x <- as.vector(weight[feed == "soybean"]);y <- as.vector(weight[feed == "linseed"])
detach(chickwts)
z <- c(x, y)
N <- c(length(x), length(y))
boot.obj <- boot(data = z, statistic = W2, R = 9999,sim = "permutation", sizes = N)
ts <- c(boot.obj$t0,boot.obj$t)
mean(ts>=ts[1])

## -----------------------------------------------------------------------------
solve_alpha(1e6,0,1,-1,0.1)
a7<-rw.Metropolis(.5, 100, 10000)
b7<-my.gibbs(15000,5000)

## -----------------------------------------------------------------------------
u<-c(11,8,27,13,16,0,23,10,24,2)
v<-c(12,9,28,14,17,1,24,11,25,3)
em(1,u,v)
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
result <- solve_game(A)

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
a<-2
b<-3
n<-10
num<-1000
ts <- microbenchmark(gibbR=r_gibbs(num,a,b,n),gibbC=rcpp_gibbs(num,a,b,n))
ts

