## -----------------------------------------------------------------------------
cor.sprman<-function(z,ix){
  #compute the Spearman rank correlation test statistic
  #z:two samples x and y
  #ix:a permutation of row indices of y
  x<-z[,1]
  y<-z[ix,2]#permutate rows of y
  return(cor(x,y,method = "spearman"))
}

## -----------------------------------------------------------------------------
set.seed(123)
library(boot)
z<-as.matrix(iris[1:50,c(1,3)])#Sepal.length and Petal.length
asl<-numeric(4)
asl<-replicate(4,expr={#compute achieved significance level of the permutation test
  boot.obj<-boot(data=z,statistic=cor.sprman,R=999,sim="permutation")
  tb<-c(boot.obj$t0,boot.obj$t)
  hist(tb,nclass="scott",xlab="cor_spearman",main="",freq=FALSE,)
  points(boot.obj$t0,0,cex=1,pch=16)
  mean(abs(tb)>=boot.obj$t0)
})

## -----------------------------------------------------------------------------
print(asl)

## -----------------------------------------------------------------------------
x<-z[,1]#Sepal.length
y<-z[,2]#Petal.length
c<-cor.test(x,y,alternative="two.sided",method="spearman",exact = FALSE)
print(c$p.value)

## -----------------------------------------------------------------------------
library(MASS)
library(RANN)#for nn2 function 
library(energy)#for eqdist.etest()
library(Ball)
knn<-function(z,ix,sizes,k){
  #compute k-Nearest neighbor statistic
  n1<-sizes[1]
  n2<-sizes[2]
  n<-n1+n2
  if(is.vector(z))z<-data.frame(z)
  z<-z[ix,]
  NN<-nn2(z,k=k+1)
  block1<-NN$nn.idx[1:n1,-1]
  block2<-NN$nn.idx[(n1+1):n,-1]
  i1<-sum(block1<n1+.5)
  i2<-sum(block2>n1+.5)
  return((i1+i2)/(k*n))
}
knn.test<-function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=knn,R=999, sim="permutation", sizes = sizes,k=k)
  tb<- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(tb>=boot.obj$t0)
  list(statistic=boot.obj$t0,p.value=p.value)
}

## -----------------------------------------------------------------------------
mu2<-mu1<-c(0,0,0)
s1<-diag(3)
s2<-matrix(c(3,0,0,0,2,0,0,0,3),nrow=3,ncol=3)
n1=n2=15
n <- n1+n2 
N = c(n1,n2)
k=3
m=100
powr<- numeric(3)#store power for three kinds tests
set.seed(11)
p.values<-replicate(m,expr = {
  ##compute p-values for nn,energy,ball methods
  x<- mvrnorm(n1,mu1,s1)
  y<- mvrnorm(n2,mu2,s2)
  z<- rbind(x,y)
  p1 <- knn.test(z,sizes=N,k)$p.value#NN method
  p2<- eqdist.etest(z,sizes=N,R=999)$p.value#energy method
  p3 <- bd.test(x,y,num.permutations=999)$p.value#ball method
  list(p1,p2,p3)
})
for(i in 1:3){
  powr[i]=mean(as.numeric(p.values[i,])<0.05)##set alpha as 0.05
} 
print(c(power_nn=powr[1],power_energy=powr[2],power_ball=powr[3]))


## -----------------------------------------------------------------------------
mu1<-c(0,0,0);mu2<-c(1,1,1)
s1<-diag(3)
s2<-matrix(c(2,0,0,0,2,0,0,0,3),nrow=3,ncol=3)
n1=n2=15
n <- n1+n2 
N = c(n1,n2)
k=3
m=100
powr<- numeric(3)#store power for three kinds tests
set.seed(12)
p.values<-replicate(m,expr = {
  ##compute p-values for nn,energy,ball methods
  x<- mvrnorm(n1,mu1,s1)
  y<- mvrnorm(n2,mu2,s2)
  z<- rbind(x,y)
  p1 <- knn.test(z,sizes=N,k)$p.value#NN method
  p2<- eqdist.etest(z,sizes=N,R=999)$p.value#energy method
  p3 <- bd.test(x,y,num.permutations=999)$p.value#ball method
  list(p1,p2,p3)
})
for(i in 1:3){
  powr[i]=mean(as.numeric(p.values[i,])<0.05)##set alpha as 0.05
} 
print(c(power_nn=powr[1],power_energy=powr[2],power_ball=powr[3]))


## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
m=100
powr<- numeric(3)#store power for three kinds tests
set.seed(111)
p.values<-replicate(m,expr = {
  ##compute p-values for nn,energy,ball methods
  x<-as.matrix(rt(n1,1,2));y<-as.matrix(rt(n2,2,5))
  z<- rbind(x,y)
  p1 <- knn.test(z,sizes=N,k)$p.value#NN method
  p2<- eqdist.etest(z,sizes=N,R=999)$p.value#energy method
  p3 <- bd.test(x,y,num.permutations=999)$p.value#ball method
  list(p1,p2,p3)
})
for(i in 1:3){
  powr[i]=mean(as.numeric(p.values[i,])<0.05)##set alpha as 0.05
} 
print(c(power_nn=powr[1],power_energy=powr[2],power_ball=powr[3]))


## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
m=100
powr<- numeric(3)#store power for three kinds tests
set.seed(111)
p.values<-replicate(m,expr = {
  ##compute p-values for nn,energy,ball methods
  e1<-sample(0:1,size=20,replace=TRUE,prob = c(0.7,0.3))
  x<-e1*rnorm(n1,0,1)+(1-e1)*rnorm(n1,1,2)#x is data from bimodel distribution
  x<-as.matrix(x)
  e2<-sample(0:1,size=20,replace=TRUE,prob = c(0.5,0.5))
  y<-e2*rnorm(n2,0,2)+(1-e2)*rnorm(n2,1,1)#y is data from bimodel distribution
  y<-as.matrix(y)
  z<- rbind(x,y)
  p1 <- knn.test(z,sizes=N,k)$p.value#NN method
  p2<- eqdist.etest(z,sizes=N,R=999)$p.value#energy method
  p3 <- bd.test(x,y,num.permutations=999)$p.value#ball method
  list(p1,p2,p3)
})
for(i in 1:3){
  powr[i]=mean(as.numeric(p.values[i,])<0.05)##set alpha as 0.05
} 
print(c(power_nn=powr[1],power_energy=powr[2],power_ball=powr[3]))


## -----------------------------------------------------------------------------
mu2<-mu1<-c(0,0,0)
s1<-diag(3)
s2<-matrix(c(3,0,0,0,2,0,0,0,3),nrow=3,ncol=3)
n1=10;n2=100
n <- n1+n2 
N = c(n1,n2)
k=3
m=100
powr<- numeric(3)#store power for three kinds tests
set.seed(11)
p.values<-replicate(m,expr = {
  ##compute p-values for nn,energy,ball methods
  x<- mvrnorm(n1,mu1,s1)
  y<- mvrnorm(n2,mu2,s2)
  z<- rbind(x,y)
  p1 <- knn.test(z,sizes=N,k)$p.value#NN method
  p2<- eqdist.etest(z,sizes=N,R=999)$p.value#energy method
  p3 <- bd.test(x,y,num.permutations=999)$p.value#ball method
  list(p1,p2,p3)
})
for(i in 1:3){
  powr[i]=mean(as.numeric(p.values[i,])<0.05)##set alpha as 0.05
} 
print(c(power_nn=powr[1],power_energy=powr[2],power_ball=powr[3]))


## -----------------------------------------------------------------------------
    set.seed(10)
    m <- 10000
    x <- numeric(m)#初始化随机向量
    x[1] <- rnorm(1)
    #建议分布选择方差为1的正态分布,从标准正态分布中生成X0并储存在x[1]中
    k <- 0
    u <- runif(m)

    for (i in 2:m) {
        xt <- x[i-1]
        y <- rnorm(1,mean=xt,sd=1)#从N(Xt,1)=N(x[i-1],1)中生成y
        num <- dcauchy(y) * dnorm(xt, mean=y,sd=1)
        den <- dcauchy(xt) * dnorm(y, mean=xt,sd=1)
        if (u[i] <= num/den){
          x[i] <- y
        } else {
          x[i] <- xt
          k <- k+1     #y is rejected
        }
    }
    index<-5000:5500
    y1 <- x[index]
    plot(index, y1, type="l", main="", ylab="x")


## -----------------------------------------------------------------------------

    b <- 1001      #discard the burn-in sample
    y <- x[b:m]
    a <- seq(0.1,0.9,0.1)
    QC<- tan(pi*(a-0.5)) #quantiles of Cauchy
    Q <- quantile(x, a) #the deciles of the generated chain
    par(mfrow=c(1,2))
    hist(y, breaks="scott", main="Standard Cauchy", xlab="", freq=FALSE)
    lines(QC, dcauchy(QC),col='blue')
    qqplot(QC, Q, main="",xlab="Cauchy deciles", ylab="Sample deciles")
    abline(0,1,col='red',lwd=2)
    

## -----------------------------------------------------------------------------
    Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
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

    cauchy.chain <- function( sigma,N, X1) {
        #generates a Metropolis chain for stantard cauchy distribution
        #with Normal(X[t], 1) proposal distribution
        #and starting value X1
        x <- rep(0, N)
        x[1] <- X1
        u <- runif(N)

        for (i in 2:N) {
            xt <- x[i-1]
            y <- rnorm(1, mean=xt, sd=sigma)     #candidate point
            r1 <- dcauchy(y) * dnorm(xt, y, sd=sigma)
            r2 <- dcauchy(xt) * dnorm(y, xt, sd=sigma)
            r <- r1 / r2
            if (u[i] <= r) x[i] <- y else
                 x[i] <- xt
            }
        return(x)
        }

  
    k <- 4          #number of chains to generate
    n <- 15000      #length of chains
    b <- 1000       #burn-in length

    #choose overdispersed initial values
    x0 <- c(-10, -5, 5, 10)

    #generate the chains
    set.seed(11)
    sigma<-seq(1.5,2.5,0.1)
    GR<-numeric(length(sigma))
    for(j in 1:length(sigma)){
      X <- matrix(0, nrow=k, ncol=n)
      for (i in 1:k)
        X[i, ] <- cauchy.chain(sigma[j],n, x0[i])

      #compute diagnostic statistics
      psi <- t(apply(X, 1, cumsum))
      for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))
      GR[j]=Gelman.Rubin(psi)
    }
    print(GR)

## -----------------------------------------------------------------------------
    #plot psi for the four chains
    sigma=1.9
    X <- matrix(0, nrow=k, ncol=n)
      for (i in 1:k)
        X[i, ] <- cauchy.chain(sigma,n, x0[i])

      #compute diagnostic statistics
      psi <- t(apply(X, 1, cumsum))
      for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))
    for (i in 1:k)
      if(i==1){
        plot((b+1):n,psi[i, (b+1):n],ylim=c(-2,2), type="l",
            xlab='Index', ylab=bquote(phi))
      }else{
        lines(psi[i, (b+1):n], col=i)
      }
    print(Gelman.Rubin(psi))
    par(mfrow=c(1,1)) #restore default
    
    #plot the sequence of R-hat statistics
    rhat <- rep(0, n)
    for (j in (b+1):n)
        rhat[j] <- Gelman.Rubin(psi[,1:j])
    plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
    abline(h=1.2, lty=2,col="blue")

## -----------------------------------------------------------------------------
#initialize constants and parameters
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
a=2;b=3;n=10#初始化参数
x1=5;x2=0.5
###### generate the chain #####
set.seed(15)
X[1, ] <- c(x1, x2) #initialize
for (i in 2:N) {
x2 <- X[i-1, 2]
X[i, 1] <- rbinom(1, size=n, prob=x2)
x1 <- X[i, 1]
X[i, 2] <- rbeta(1, x1+a, n-x1+b)
}
bu <- burn + 1
x <- X[bu:N, ]

cat('Means: ',round(colMeans(x),2))
cat('Standard errors: ',round(apply(x,2,sd),2))
cat('Correlation coefficients: ', round(cor(x[,1],x[,2]),2))

## -----------------------------------------------------------------------------
plot(x[,1],type='l',col=1,lwd=2,xlab='Index',ylab='Random numbers')
lines(x[,2],col=2,lwd=2)
legend('bottomright',c(expression(X[1]),expression(X[2])),col=1:2,lwd=2)

## -----------------------------------------------------------------------------
    Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
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

        x.chain <- function( N, X1,yvalue) {
        #generates a Metropolis chain for Binormial distribution
        #with poisson(Xt) proposal distribution
        #and starting value X1
        x <- rep(0, N)
        x[1] <- X1
        u <- runif(N)

        for (i in 2:N) {
            xt <- x[i-1]
            y <- rpois(1,lambda=xt)     #candidate point
            r1 <- dbinom(y,size=10,prob=yvalue) * dpois(xt, y)#n=10
            r2 <- dbinom(xt,size=10,prob=yvalue) * dpois(y, xt)
            r <- r1 / r2
            if (u[i] <= r) x[i] <- y else
                 x[i] <- xt
            }
        return(x)
        }
        y.chain <- function( N, X1,xvalue) {
        #generates a Metropolis chain for stantard cauchy distribution
        #with chisq(X[t], 1) proposal distribution
        #and starting value X1
        x <- rep(0, N)
        x[1] <- X1
        u <- runif(N)

        for (i in 2:N) {
            xt <- x[i-1]
            y <- rchisq(1, df=xt)     #candidate point
            r1 <- dbeta(y,xvalue+2,13-xvalue) * dchisq(xt, y)#a=2,b=3
            r2 <- dbeta(xt,xvalue+2,13-xvalue) * dchisq(y, xt)#a=2,b=3
            r <- r1 / r2
            if (u[i] <= r) x[i] <- y else
                 x[i] <- xt
            }
        return(x)
        }

  
    k <- 4          #number of chains to generate
    n <- 15000      #length of chains
    b <- 1000       #burn-in length

    #choose overdispersed initial values
    x0 <- c(1,2,3,4)
    y0<-c(0.2,0.3,0.4,0.5)
    #generate the chains
    set.seed(16)
      X <-Y<- matrix(0, nrow=k, ncol=n)
      for (i in 1:k){
        X[i, ] <- x.chain(n, x0[i],yvalue=0.8)
        Y[i, ] <- y.chain(n, y0[i],xvalue=3)}
      #compute diagnostic statistics
      psi1 <- t(apply(X, 1, cumsum))
      psi2 <- t(apply(Y, 1, cumsum))
      for (i in 1:nrow(psi1)){
        psi1[i,] <- psi1[i,] / (1:ncol(psi1))
        psi2[i,] <- psi2[i,] / (1:ncol(psi2))}
    print(c(Gelman.Rubin(psi1),Gelman.Rubin(psi2)))

## -----------------------------------------------------------------------------
    #plot the sequence of R-hat statistics
    rhat1<- rhat2<- rep(0, n)
    for (j in (b+1):n){
        rhat1[j] <- Gelman.Rubin(psi1[,1:j])
        rhat2[j] <- Gelman.Rubin(psi2[,1:j])
        }
    plot(rhat1[(b+1):n], type="l", xlab="", ylab="R",ylim = c(1,1.3))
    abline(h=1.2, lty=2,col="blue")
    plot(rhat2[(b+1):n], type="l", xlab="", ylab="R")
    abline(h=1.2, lty=2,col="blue")

## -----------------------------------------------------------------------------
options(warn=-1)
##compute terms
fa<-function(k,a){
  ##compute kth term for the series
  d<-length(a)
  i<-rep(c(1,-1),length=k)#1 or -1 for k
  ik<-i[k]
  k<-k-1
  c1=exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
  c2<-factorial(k)*2^k*(2*k+1)*(2*k+2)
  c<-ik*c1*norm(a)^(2*k+2)/c2
  return(c)
}


##compute sum
fb<-function(k,a){
  ##return sum for k terms
  d<-length(a)
  i<-rep(c(1,-1),length=k)#1 or -1 for k
  K<-0:(k-1)
  c1=exp(lgamma((d+1)/2)+lgamma(K+3/2)-lgamma(K+d/2+1))
  c2<-factorial(K)*2^K*(2*K+1)*(2*K+2)
  c<-i*c1*norm(a)^(2*K+2)/c2
  return(sum(c))
}

a<-matrix(c(1,2))
print(c(fa(1,a),fa(2,a),fa(3,a)))#terms for k=0,1,2

f<-k<-0:25
for(i in 1:length(k))
  f[i]<-fb(k[i],a)
print(f[c(5,10,15,20,25)])#sums until k=4,9,14,19,24
knitr::kable(rbind(k[15:25],f[15:25]))

## -----------------------------------------------------------------------------
##function for S_k(a)
sk<-function(a,k){
  x<-sqrt(a^2*k/(k+1-a^2))
  pt(x,df=k,lower.tail=FALSE)

}
f<-function(a,k) sk(a,k)-sk(a,k-1)#xlim is (0,sqrt(k))

K<-c(4 : 25,100, 500, 1000)
root1<-numeric(length(K))
##first to find the approximation of the root for f
a<-seq(0,5,0.01)
plot(a,f(a,100),lty=1,type="l",xlim=c(0,5))
lines(a,f(a,500),lty=2,col="blue",xlim=c(0,5))
lines(a,f(a,1000),lty=3,col="red",xlim=c(0,5))
abline(h=0)
legend("topright", legend=c("k=100", "k=500", "k=1000"), lty=1:3)


## -----------------------------------------------------------------------------
for(i in 1:length(K)){
  k<-K[i]
  out<-uniroot(function(a){
    sk(a,K[i])-sk(a,K[i]-1)},
    c(1,2))
  root1[i]<-out$root
}
knitr::kable(rbind(K[1:12],root1[1:12]))
knitr::kable(rbind(K[-(1:12)],root1[-(1:12)]))

## -----------------------------------------------------------------------------
fck<-function(a,k){
  ##compute integration
  ck<-sqrt(a^2*k/(k+1-a^2))
  v<-integrate(function(u){
    (1+u^2/k)^(-(k+1)/2)},
    lower = 0,upper = ck,
    rel.tol = .Machine$double.eps^0.25)$value
  return(exp(lgamma((k+1)/2)-lgamma(k/2))/sqrt(k)*v)
}

##write function to compute the right minus the left
fun<-function(a,k){
  ##a is a numeric vector
  m<-length(a)
  x<-1:m
  for(i in 1:m)
    x[i]<-fck(a[i],k)-fck(a[i],k-1)
  return(x)
}

K<-s.root<-c(4 : 25,100, 500, 1000)
a<-seq(0,5,0.01)
plot(a,fun(a,100),lty=1,type="l",xlim=c(0,5))
lines(a,fun(a,500),lty=2,col="blue",xlim=c(0,5))
lines(a,fun(a,1000),lty=3,col="red",xlim=c(0,5))
abline(h=0)
legend("topright", legend=c("k=100", "k=500", "k=1000"), lty=1:3)

## -----------------------------------------------------------------------------
root2<-numeric(length(K))
for(i in 1:length(K)){
  k<-K[i]
  out<-uniroot(function(a){
    fck(a,k)-fck(a,k-1)},
    c(1,2))
  root2[i]<-out$root
}
knitr::kable(rbind(K[1:12],root2[1:12],root1[1:12]))
knitr::kable(rbind(K[-(1:12)],root2[-(1:12)],root1[-(1:12)]))
#root1 are the roots of ex 11.4;root2 are the roots of ex 11.5

## ----eval=FALSE---------------------------------------------------------------
#  library(nloptr)
#  y<-c(0.54, 0.48, 0.33, 0.43, 0.91, 0.21, 0.85,1.00, 1.00,1.00)
#  
#  eval.lam<-function(lam,lam0,data=y){
#    c<-10*log(lam)+(sum(data[1:7])+3*lam0+3)/lam
#    return(c)
#  }
#  
#  N<-10000
#  l0<-as.vector(0.8)##initial est. for lambdas
#  tol<-.Machine$double.eps^0.5
#  opts = list("algorithm"="NLOPT_LN_COBYLA",
#               "xtol_rel"=tol)
#  d=1
#  j=1
#  while(abs(d)>tol){
#    out<-nloptr(x0=0.8,
#                eval_f=eval.lam,
#                lb=0,ub=10,
#                opts = opts,
#                lam0=l0[j],data=y)
#    d<-l0[j]-out$solution
#    l0<-c(l0,out$solution)
#    j<-j+1
#  }
#  print(l0)#lambdas
#  k<-length(l0)
#  lam.EM<-(sum(y[1:7])+3)/7
#  knitr::kable(cbind(k,l0[k],lam.EM),col.names = c("iteration times","test result","theoretical values"))

## -----------------------------------------------------------------------------
set.seed(17)
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
unlist(lapply(trims, function(trim) mean(x, trim = trim)))
unlist(lapply(trims, mean, x = x))

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
attach(mtcars)

formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

#method 1 for excercise 3
l3=lapply(formulas, function(x) lm(formula = x, data = mtcars))
rsq3.1=lapply(l3,rsq)
unlist(rsq3.1)
#method 2 for excercise 3
rsq3.2<-lapply(formulas, function(x){
  mod<-lm(formula = x, data = mtcars)
  rsq(mod)}
  )
unlist(rsq3.2)


bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
#method 1 for excercise 4
l4<-lapply(bootstraps, function(x) lm(formula =mpg ~ disp, data = x))
rsq4.1=lapply(l4,rsq)#method 1 for excercise 4
unlist(rsq4.1)

#method 2 for excercise 4
rsq4.2<-lapply(bootstraps, function(x){
  l<-lm(formula =mpg ~ disp, data = x)
  rsq(l)
})
unlist(rsq4.2)


## -----------------------------------------------------------------------------
#(a)Compute sd of every column in a numeric data frame
unlist(lapply(mtcars,class))
#mtcars is a numeric data frame,so we take it as example.
vapply(mtcars,sd,numeric(1))

#(b)Compute sd of every numeric column in a mixed data frame
attach(CO2)#Carbon Dioxide Uptake in Grass Plants;CO2 is a mixed data frame
idx<-vapply(CO2,is.numeric,logical(1))
vapply(CO2[idx],sd,numeric(1))

## ----eval=FALSE---------------------------------------------------------------
#  library(parallel)
#  
#  sapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
#      FUN <- match.fun(FUN)
#      answer <- lapply(X = X, FUN = FUN, ...)
#      if (USE.NAMES && is.character(X) && is.null(names(answer)))
#          names(answer) <- X
#      if (!isFALSE(simplify) && length(answer))
#          simplify2array(answer, higher = (simplify == "array"))
#      else answer
#  }
#  
#  
#  mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
#    FUN <- match.fun(FUN)
#    answer <- mclapply(X = X, FUN = FUN, ...)
#    if (USE.NAMES && is.character(X) && is.null(names(answer)))
#      names(answer) <- X
#    if (!isFALSE(simplify) && length(answer))
#      simplify2array(answer, higher = (simplify == "array"))
#    else answer
#  }
#  
#  
#  boot_df <- function(x) x[sample(nrow(x), rep = T), ]
#  rsquared <- function(mod) summary(mod)$r.square
#  boot_lm <- function(i) {
#    dat <- boot_df(mtcars)
#    rsquared(lm(mpg ~ wt + disp, data = dat))
#  }
#  n <- 1e4
#  system.time(sapply(1:n, boot_lm))
#  system.time(mcsapply(1:n, boot_lm,mc.cores=1))
#  

## -----------------------------------------------------------------------------
gibbsR<- function(N,n=9,a=2,b=3) {
  mat <- matrix(nrow = N, ncol = 2)
  x <-5; y <- 0.5
  for (i in 1:N) {
      x <- rbinom(1, size = n, prob=y)
      y <- rbeta(1, x+a, n-x+b)
      mat[i, ] <- c(x, y)
  }
  mat
}


## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('
NumericMatrix gibbs_C(int N){
  NumericMatrix mat(N, 2);
  int x0=5;
  double y0=0.5;
  int n=9;
  int a=2;
  int b=3;
  for (int i=0; i<N;i++) {
    int x = rbinom(1, n, y0)[0];
    mat(i,0)=x;
    x0=x;
    double y = rbeta(1, x0+a, n-x0+b)[0];
    mat(i,1)=y;
    y0=y;
  }
  return(mat);
}')

## -----------------------------------------------------------------------------
options(warn = -1)
set.seed(19)
x1<-gibbsR(5000)##random numbers with R
x2<-gibbs_C(5000)##random numbers with Rcpp

## -----------------------------------------------------------------------------
xR<-x1[,1]
xC<-x2[,1]
yR<-x1[,2]
yC<-x2[,2]
par(mfrow=c(1,2))
qqplot(xR[-(1:1000)],xC[-(1:1000)])
abline(0,1,col='red',lwd=2)
qqplot(yR[-(1:1000)],yC[-(1:1000)])
abline(0,1,col='red',lwd=2)


## -----------------------------------------------------------------------------
library(microbenchmark)
ts<-microbenchmark(gibR=gibbsR(5000),gibC=gibbs_C(5000))
summary(ts)[,c(1,3,5,6)]

