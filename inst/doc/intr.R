## -----------------------------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
summary(lm.D9)$coef
knitr::kable(head(iris))
par(mar=c(1,1,1,1))
plot(lm.D9)

## -----------------------------------------------------------------------------
n<-1e4
for (deta in 1:4) {
   x1<-rnorm(n,0,1)
   x2<-rnorm(n,0,1)
   y<-deta*sqrt(x1^2+x2^2)
   hist(y,freq = F,main=paste("Histogram of y with deta=",deta))
   z<-seq(0,15,.1)
   lines(z,(z/deta^2)*exp((-z^2)/(2*deta^2)))
}


## -----------------------------------------------------------------------------
n <- 1e3
X1 <- rnorm(n,0,1)
X2 <- rnorm(n,3,1)
for (p in c(0.25,0.5,0.75)) {
  r <- sample(c(1,0),n,replace=TRUE,prob = c(p,1-p))
  Z <- r*X1+(1-r)*X2
  hist(Z,freq = F,main=paste("Histogram of z with p1=",p))
}


## -----------------------------------------------------------------------------
n<-10
d<-8
A<-matrix(nrow=d,ncol=d)
for (i in 1:d){
  A[i,i]<-rchisq(1,n-i+1)
  j<-i+1
 while (j<=d){
   A[i,j]=0
   j=j+1
 }                      
  k=i-1
 while (k<i&&k>0) {
   A[i,k]<-rnorm(1,0,1)
   k=k-1
  }                     
}
V<-matrix(nrow = d,ncol = d)
for (i in 1:d) {
  for (j in 1:d) {
    if(i==j){V[i,j]=1}
    else V[i,j]=0
  }
}
L<-chol(V)             
X<-(L%*%A)%*%t(L%*%A)

## -----------------------------------------------------------------------------
set.seed(12345)
m<-1e4
x<-runif(m,min = 0,max = pi/3)
theta.hat<-mean(sin(x))*3/pi
print(c(theta.hat,1-cos(pi/3)))

## -----------------------------------------------------------------------------
set.seed(12345)
m<-1e4
k<-10 
r<-m/k 
n<-50 
T<-numeric(k)
est<-matrix(0, n, 2)
g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)
for (i in 1:n) {
est[i, 1] <- mean(g(runif(m)))
for(j in 1:k)T[j]<-mean(g(runif(r,(j-1)/k,j/k)))
est[i, 2] <- mean(T)
}
SD<-as.array(apply(est,2,sd))
knitr::kable(rbind(apply(est,2,mean),apply(est,2,sd)),format = 'html')
print((SD[1]-SD[2])/SD[1])

## -----------------------------------------------------------------------------
set.seed(12345)
m<-1e4
k<-5 
r<-m/k 
n<-100
T<-numeric(k)
est<-matrix(0, n, 2)
g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)
for (i in 1:n) {
  u<-runif(m)
  x<- -log(1-u*(1-exp(-1)))
  fg<-g(x)/(exp(-x)/(1-exp(-1)))
  est[i, 1] <- mean(fg)
  for(j in 1:k){
    u<-runif(r,(j-1)/k,j/k)
    z<- -log(1-(u*(1-exp(-1)))/5)
    fg1<-g(z)/((5*(exp(-z))/(1-exp(-1))))
    T[j]<-mean(fg1)
    }
  est[i, 2]<-sum(T)
}
knitr::kable(rbind(apply(est,2,mean),apply(est,2,sd)),format='html')

## -----------------------------------------------------------------------------
n<-20
m<-1e4
set.seed(1234)
mean.hat<-t.val<-numeric(m)
for (i in 1:m) {
  x<-rnorm(1)
  y<-rchisq(20,2)
  mean.hat[i]<-abs(x/sqrt(mean(y)))
  t.val[i]<-qt(0.975,40)
}
print(mean(mean.hat<=t.val))


## -----------------------------------------------------------------------------
n<-1e3
m<-1e4
set.seed(1234)
ske<-numeric(m)
for (i in 1:m) {
  x<-rnorm(n)
  ske[i]<-(mean((x-mean(x))^3))/(mean((x-mean(x))^2))^(3/2)
}
est<-quantile(ske,prob=c(0.025,0.05,0.95,0.975))
j=1
cv<-Var<-numeric(4)
for (p in c(0.025,0.05,0.95,0.975)) {
  cv[j] <- qnorm(p, 0, sqrt(6/n))
  Var[j]<-(p*(1-p))/n*(dnorm(est[j],0,sqrt(6/n))^2)
  j=j+1
}
res.compare<-data.frame(est,cv,Var)
knitr::kable(res.compare)


## -----------------------------------------------------------------------------
beta<-0.1
n <- 20
m <- 1e4
set.seed(1234)
epsilon <- c(seq(1,20,1))
N <- length(epsilon)
pwr <- numeric(N)
cv <- qnorm(1-beta/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
sk <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
for (j in 1:N) { 
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { 
    x<-rbeta(n,e,e)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type='b',xlab = 'alpha', ylim = c(0,.1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m)
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)


## -----------------------------------------------------------------------------
n <- 100
alpha <- 0.1
mu0 <- 1 
m <- 1e4
set.seed(1234)
p <- numeric(m) 
for (i in 1:m) {
  x <- rchisq(n, mu0)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0)
  p[i] <- ttest$p.value
}

p.hat <- mean(p<alpha)
se.hat <- sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat, se.hat))

## -----------------------------------------------------------------------------
n <- 100
alpha <- 0.1
mu0 <- 1 
m <- 1e4
set.seed(1234)
p <- numeric(m) 
for (i in 1:m) {
  x <- runif(n,0,2)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0)
  p[i] <- ttest$p.value
}
p.hat <- mean(p<alpha)
se.hat <- sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat, se.hat))

## -----------------------------------------------------------------------------
n <- 100
alpha <- 0.1
mu0 <- 1 
m <-1e4
set.seed(1234)
p <- numeric(m) 
for (i in 1:m) {
  x <- rexp(n)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0)
  p[i] <- ttest$p.value
}

p.hat <- mean(p<alpha)
se.hat <- sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat, se.hat))

## -----------------------------------------------------------------------------
set.seed(12345)
library(bootstrap)
library(boot)
scor = data.matrix(scor)
pairs(scor, main='the scatter plots for each pair of test scores', pch = 18)
b.cor = function(x,i){  
  r = cor(x[i,])
  c(r[1,2], r[3,4], r[3,5], r[4,5])
}
sc.boot = boot(data = scor, statistic = b.cor, R=2000)
apply(sc.boot$t, 2, FUN = sd)

## -----------------------------------------------------------------------------
set.seed(12345)
library(boot)
# normal
skew = 0
sk = function(x,i) {
  xbar = mean(x[i])
  m3 = mean((x[i] - xbar)^3)
  m2 = mean((x[i] - xbar)^2)
  return( m3 / m2^1.5 )
}
n=10
m=1000
ci.norm = ci.basic = ci.perc = matrix(0,m,2)

for(i in 1:m){
 x = rnorm(n)
 b = boot(data = x, statistic = sk, R = 2000)
 ci = boot.ci(b, type=c("norm","basic","perc"))
 ci.norm[i,] = ci$norm[2:3]
 ci.basic[i,] = ci$basic[4:5]
 ci.perc[i,] = ci$percent[4:5]
}
# output
output1 = output2 = output3 = matrix(0,2,3)
rownames(output1)=rownames(output2)=rownames(output3)=c("N(0,1)","X^2(5)")
colnames(output1)=colnames(output2)=colnames(output3)=c("norm","basic","perc")

output1[1,] = c(mean(ci.norm[,1]<=skew & ci.norm[,2]>=skew), mean(ci.basic[,1]<=skew & ci.basic[,2]>=skew), mean(ci.perc[,1]<=skew & ci.perc[,2]>=skew))
output2[1,] = c(mean(ci.norm[,1]>skew), mean(ci.basic[,1]>skew), mean(ci.perc[,1]>skew)) 
output3[1,] = c(mean(ci.norm[,2]<skew), mean(ci.basic[,2]<skew), mean(ci.perc[,2]<skew)) 
# chi-square
skew = sqrt(8/5)
ci.norm = ci.basic = ci.perc = matrix(0,m,2)
for(i in 1:m){
 x = rchisq(n,5)
 b = boot(data = x, statistic = sk, R = 2000)
 ci = boot.ci(b, type=c("norm","basic","perc"))
 ci.norm[i,] = ci$norm[2:3]
 ci.basic[i,] = ci$basic[4:5]
 ci.perc[i,] = ci$percent[4:5]
}
# output
output1[2,] = c(mean(ci.norm[,1]<=skew & ci.norm[,2]>=skew), mean(ci.basic[,1]<=skew & ci.basic[,2]>=skew), mean(ci.perc[,1]<=skew & ci.perc[,2]>=skew))
output2[2,] = c(mean(ci.norm[,1]>skew), mean(ci.basic[,1]>skew), mean(ci.perc[,1]>skew)) 
output3[2,] = c(mean(ci.norm[,2]<skew), mean(ci.basic[,2]<skew), mean(ci.perc[,2]<skew)) 
output1
output2
output3


## -----------------------------------------------------------------------------
data(scor, package = "bootstrap")
n<-nrow(scor)
f<-function(x){
  mat<-cov(x)
  evals<-eigen(mat)$values
  return(evals[1]/sum(evals))
}
theta.hat<-f(scor)
theta.jack<-numeric(n)
for (i in 1:n) {
  theta.jack[i]=f(scor[-i,])
}
bias<-(n-1)*(mean(theta.jack)-theta.hat)
print(bias)
se <- sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2))
print(se)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic) 
e1 <- e2 <- e3 <- e4 <- numeric(n)
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2)+I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k]+J4$coef[3]*chemical[k]^2+J4$coef[4]*chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
summary(J1)
summary(J2)
summary(J3)
summary(J4)

## -----------------------------------------------------------------------------
count5<- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}
set.seed(1234)
n1 <- 20
n2 <- 30
K=1:(n1+n2)
mu1 <- mu2 <- 0
sigma1 <-2
sigma2 <-1
R<-1e4
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
z<-c(x,y)
n<-length(x)
reps<-numeric(R)
t0<-count5(x,y)
for (i in 1:R) {
  k<-sample(K,size = n,replace = F)
  x1<-z[k]
  y1<-z[-k]
  reps[i]=count5(x1,y1)
}
p<-mean(c(t0,reps)>t0)
round(p,3)

## -----------------------------------------------------------------------------
distcor.test<-function(x,y){
  dCov <- function(x, y) {
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d); M <- mean(d)
    a <- sweep(d, 1, m); b <- sweep(a, 2, m)
    b + M
  }
A<- Akl(x); B <- Akl(y)
sqrt(mean(A * B))
}
ndCov2 <- function(z, ix, dims) {
  p <- dims[1]
  q <- dims[2]
  d <- p + q
  x <- z[ , 1:p] 
  y <- z[ix, -(1:p)] 
  return(nrow(z) * dCov(x, y)^2)
}
set.seed(1234)
z <-cbind(x,y)
boot.obj <- boot(data = z, statistic = ndCov2, R = 999,
sim = "permutation", dims = c(2, 2))
tb <- c(boot.obj$t0, boot.obj$t)
p.cor <- mean(tb>=tb[1])
}

## -----------------------------------------------------------------------------
library('MASS')
library('Ball')
library('boot')
n<-30
mu<-matrix(c(0,0),2,1)
sigma<-matrix(c(1,0,0,1),2,2)
x<-mvrnorm(n,mu,sigma)
e<-mvrnorm(n,mu,sigma)
y<-(x/4+e)
p.dist<-distcor.test(x,y)
p.ball <- bcov.test(x,y,R=999)$p.value
round(c(p.dist,p.ball),4)

## -----------------------------------------------------------------------------
library('MASS')
library('Ball')
library('boot')
n<-30
mu<-matrix(c(0,0),2,1)
sigma<-matrix(c(1,0,0,1),2,2)
x<-mvrnorm(n,mu,sigma)
e<-mvrnorm(n,mu,sigma)
y<-(x/4)*e
p.dist<-distcor.test(x,y)
p.ball <- bcov.test(x,y,R=999)$p.value
round(c(p.dist,p.ball),4)

## -----------------------------------------------------------------------------
library(GeneralizedHyperbolic)
Laplace.Metropolis<-function(sigma,N){
  x<-numeric(N)
  x[1]<-rnorm(1,0,sigma)
  u<-runif(N)
  k<-0
  for (i in 2:N) {
    y<-rnorm(1,x[i-1],sigma)
    if(u[i]<=dskewlap(y)/dskewlap(x[i-1]))
      x[i]=y
    else{
      x[i]=x[i-1]
      k=k+1
    }
  }
  return(list(x=x,k=k))
}
N<-1e4
sigma<-c(1,4,9,16)
set.seed(1234)
laplace1<-Laplace.Metropolis(sigma[1],N)
laplace2<-Laplace.Metropolis(sigma[2],N)
laplace3<-Laplace.Metropolis(sigma[3],N)
laplace4<-Laplace.Metropolis(sigma[4],N)
print(c(1-laplace1$k/N,1-laplace2$k/N,1-laplace3$k/N,1-laplace4$k/N))

## -----------------------------------------------------------------------------
exp(log(20)) == log(exp(20))
all.equal(exp(log(20)),log(exp(20)))

## -----------------------------------------------------------------------------
gammaratio <- function(k){
  if (k %% 2==0) 
    return(exp(sum(log((k/2-1):1))-sum(log(seq(1,k-3,2)/2))-log(gamma(1/2))))
  else
    return(exp(sum(log(seq(1,k-2,2)/2))+log(gamma(1/2))-sum(log(((k-1)/2-1):1))))
}   ### Gamma ratio computing for large k
g<- function(a,k){
  f <- function(x, N) (1 + x^2/N)^(-(N + 1)/2)
  I1 <- integrate(f,lower = 0,upper = sqrt(a^2*(k-1)/(k-a^2)),rel.tol=.Machine$double.eps^0.25,N = k-1 )$value
  I2 <- integrate(f,lower = 0,upper = sqrt(a^2*k/(k+1-a^2)),rel.tol=.Machine$double.eps^0.25,N = k )$value
  return(2*gammaratio(k)*I1/sqrt(k-1)-2*gammaratio(k+1)*I2/sqrt(k))
}   # the equation in 11.5
f1 <- function(x,k) pt(sqrt(x^2*(k-1)/(k-x^2)),k-1)-pt(sqrt(x^2*k/(k+1-x^2)),k)
##  computing method in 11.4
a <- seq(0.001,6,0.0001)
y1 <-numeric(length(a))
for ( i in 1:length(a)){
  y1[i] <- g(a[i],100) 
}
plot(cbind(a,y1),type = "l",xlab="a",ylab = "value")

## -----------------------------------------------------------------------------
nk <- c(4:25,100,500,1000)
ak <-numeric(25)
qk <-numeric(25)
for (i in 1:25) {
  ak[i] <- uniroot(g,c(0.00001,ifelse(nk[i] <=7,sqrt(nk[i])-0.01,2.5)),k=nk[i])$root # the solution of the equation
  qk[i] <- uniroot(f1,c(0.00001,ifelse(nk[i] <=7,sqrt(nk[i])-0.01,2.5)),k=nk[i])$root # the solution of 11.4
}
cbind(ak,qk)
plot(cbind(nk,ak),xlab = "k",ylab = "A(k)")
lines(cbind(nk,ak),lwd =2,col = "red")

## -----------------------------------------------------------------------------
loglike <- function(pq,ng){
  p <- pq[1]
  q <- pq[2]
  r <- 1-p -q
  a <- c(p^2,q^2,r^2,2*p*r,2*q*r,2*p*q)
  return(-sum(log(a)*ng))
}
na <- 28
nb <- 24
noo <- 41
nab <- 70
n <-na + nb +noo+nab
eloglike <- function(pq){
  p <- pq[1]
  q <- pq[2]
  r <- 1-p -q
  a <- c(p^2+2*p*r,q^2+2*q*r,r^2,2*p*q)
  ng <- c(na,nb,noo,nab)
  return(-sum(log(a)*ng))
}
r <- sqrt(noo/n)
f2 <- function(x) x^2 +2*x*r - na/n
p <- uniroot(f2,c(0,1))$root
q <- 1-p-r
m <- 200
pqr <- matrix(rep(0,m*3),nrow = m)
ll <- numeric(m)
pqr[1,]<-c(p,q,r)
for (i in 2:m) {
  p <- pqr[i-1,1]
  q <- pqr[i-1,2]
  r <- pqr[i-1,3]
  naa <- p/(p+2*r)*na
  nao <- 2*r/(p+2*r)*na
  nbb <- q/(q+2*r)*nb
  nbo <- 2*r/(q+2*r)*nb
  ng <- c(naa,nbb,noo,nao,nbo,nab)  #estimation step
  res <- optim(c(p,q),loglike,method="Nelder-Mead",ng=ng) #maximum step
  pqr[i,]<-c(res$par,1-res$par[1]-res$par[2])
  ll[i] <- -eloglike(res$par)
}

## -----------------------------------------------------------------------------
pqr[200,]

## -----------------------------------------------------------------------------
res <- optim(c(p,q),eloglike,method="Nelder-Mead")
res$par
plot(ll[-1],type = "l",ylab = "log-likelihood")
abline(a=ll[200],b=0,col="red")

## -----------------------------------------------------------------------------
formulas <- list(mpg ~ disp,mpg ~ I(1 / disp),mpg ~ disp + wt,mpg ~ I(1 / disp) + wt)
out1<-vector("list",length(formulas))
for (i in 1:length(formulas)) {
  out1[[i]]<-lm(formulas[[i]],data=mtcars)
}
out2<-lapply(formulas,lm,data=mtcars)
out1

## -----------------------------------------------------------------------------
bootstraps<-lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
  }
)
out3<-vector("list",10)
for (i in 1:10) {
  out3[[i]]<-lm(mpg~disp,data=bootstraps[[i]])
}
out4<-lapply(bootstraps,lm,formula=mpg~disp)
out3

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
lapply(out1,rsq)
lapply(out3,rsq)

## -----------------------------------------------------------------------------
trials<-replicate(100, t.test(rpois(10, 10), rpois(7, 10)), simplify = FALSE)
p_value<-sapply(trials, function(x) x$p.value)
p_value2<-sapply(trials, '[[',3)
p_value
p_value2

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
library(GeneralizedHyperbolic)
cppFunction('List Cpplaplace(double Sigma,int N){
  NumericVector x(N);
  int k=0;
  x[0]=as<double>(rnorm(1,25,Sigma));
  for(int i=1;i<N;i++){
    double y=as<double>(rnorm(1,x[i-1],Sigma));
    double u=as<double>(runif(1));
    if (u<=exp(-abs(y))/exp(-abs(x[i-1]))) {
    x[i]=y;
    k=k+1;  
    }
    else x[i]=x[i-1];
  }
  return(List::create(Named("x")=x,Named("k")=k));
}')
Laplace.Metropolis<-function(sigma,N){
  x<-numeric(N)
  x[1]<-rnorm(1,0,sigma)
  u<-runif(N)
  k<-0
  for (i in 2:N) {
    y<-rnorm(1,x[i-1],sigma)
    if(u[i]<=dskewlap(y)/dskewlap(x[i-1])){
      x[i]=y
      k=k+1
    }
    else
      x[i]=x[i-1]
  }
  return(list(x=x,k=k))
}
N<-10000
sigma<-c(1,4,9,16)
set.seed(1234)
laplace1<-Laplace.Metropolis(sigma[1],N)
laplace2<-Laplace.Metropolis(sigma[2],N)
laplace3<-Laplace.Metropolis(sigma[3],N)
laplace4<-Laplace.Metropolis(sigma[4],N)
Claplace1<-Cpplaplace(sigma[1],N)
Claplace2<-Cpplaplace(sigma[2],N)
Claplace3<-Cpplaplace(sigma[3],N)
Claplace4<-Cpplaplace(sigma[4],N)
par(mar=c(1,1,1,1))
qqplot(laplace1$x,Claplace1$x)
qqplot(laplace2$x,Claplace2$x)
qqplot(laplace3$x,Claplace3$x)
qqplot(laplace4$x,Claplace4$x)
compare1<-microbenchmark(Laplace.Metropolis(sigma[1],N),Cpplaplace(sigma[1],N)) 
compare2<-microbenchmark(Laplace.Metropolis(sigma[2],N),Cpplaplace(sigma[2],N))
compare3<-microbenchmark(Laplace.Metropolis(sigma[3],N),Cpplaplace(sigma[3],N))
compare4<-microbenchmark(Laplace.Metropolis(sigma[4],N),Cpplaplace(sigma[4],N))
summary(compare1)[c(1,3,5,6)]
summary(compare2)[c(1,3,5,6)]
summary(compare3)[c(1,3,5,6)]
summary(compare4)[c(1,3,5,6)]

