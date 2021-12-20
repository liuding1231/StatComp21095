## -----------------------------------------------------------------------------
df<-cars
print(df[1:10,])
knitr::kable(summary(df))

## -----------------------------------------------------------------------------
x=df[,"speed"]
y=df[,"dist"]
plot(x,y,xlab="speed",ylab="distance",xlim=c(0,max(df$speed)+1),ylim=c(0,max(df$dist)+1),'p')
model<-lm(y~x,data=df)
abline(model)

## -----------------------------------------------------------------------------
fe<-function(){
  u<-runif(3,-1,1)#生成服从均匀分布U(-1,1)的三个随机数
  uabs=abs(u)#绝对值化
  indx<-which(uabs==max(uabs))#找绝对值化后最大值的下标
  u2=u[-indx]#去掉绝对值最大的数
  u2[2]#输出需要的随机数
}

## -----------------------------------------------------------------------------
library(MASS)
x<-1:1000
for(i in 1:1000){
 x[i]<-fe() 
}
hist(x,prob=TRUE,breaks=20,main=expression(f(x)==3/4*(1-x^2)))
x1<-seq(-1,1,.01)
lines(x1,col="red",3/4*(1-x1^2))

## -----------------------------------------------------------------------------
ppnorm<-function(p1,p2){
  ber<-sample(0:1,size=1000,replace=TRUE,prob=c(p2,p1))
  x11<-1:1000
  for(i in 1:1000){
    x11[i]<-ber[i]*rnorm(1,0,1)+(1-ber[i])*rnorm(1,3,1)
  }
  hist(x11,prob=TRUE,main=expression(f(x)==p1*N(0,1)+p2*N(3,1)))
  y<-seq(floor(min(x11)),ceiling(max(x11)),0.01)
  lines(y,p1*dnorm(y,0,1)+p2*dnorm(y,3,1),col="red")
}
ppnorm(.75,.25)

## -----------------------------------------------------------------------------
par(mfrow=c(1,3))
for(i in 1:9){
  p<-i/10
  ppnorm(p,1-p)
}

## -----------------------------------------------------------------------------
pogam<-function(lmd,alpha=1,rate=1){#模拟混合分布算法
  t0<-10
  tn<-rexp(100,lmd)
  sn<-cumsum(tn)
  n<-min(which(sn>t0))
  xt<-rgamma(n,1,1)
  sum(xt)
}
lmd<-2#给参数lambda赋值
pg<-1:1000
for(i in 1:1000){
  pg[i]<-pogam(lmd,alpha,rate)
}#生成1000个混合分布的随机样本

## -----------------------------------------------------------------------------
tab1<-matrix(0,3,5)
tab2<-matrix(0,3,5)
row1<-c('λ','样本均值','理论均值')
row2<-c('λ','样本方差','理论方差')
dimnames(tab1)<-list(row1,c('','','','',''))
dimnames(tab2)<-list(row2,c('','','','',''))
for(i in 1:5){
  lmd<-i
  p<-pogam(lmd)
  m_p<-mean(p)#样本均值
  s_p<-var(p)#样本方差
  m<-10*lmd
  s<-20*lmd
  tab1[,i]<-c(i,m_p,m)
  tab2[,i]<-c(i,s_p,s)
}

knitr::kable(tab1)
knitr::kable(tab2)

## -----------------------------------------------------------------------------
beta_cdf<-function(x,m=1000){
  u<-runif(m)
  cdf<-numeric(length(x))
  for(i in 1:length(x)){
    g<-u^2*x[i]^3*(1-u*x[i])^2
    cdf[i]<-mean(g)/beta(3,3)
  }
  cdf
}

## -----------------------------------------------------------------------------
set.seed(1)
x<-seq(0.1,0.9,length=9)
CDF<-beta_cdf(x)
beta_real<-pbeta(x,3,3)
knitr::kable(round(rbind(x,CDF,beta_real),3))


## -----------------------------------------------------------------------------
MC.ray<-function(x,s,R=10000,antithetic=TRUE){
  u<-runif(R/2)
  if(!antithetic) v<-runif(R/2)
  else v<-1-u
  u<-c(u,v)
  cdf<-numeric(length(x))
  for(i in 1:length(x)){
    g<-u*x[i]^2*exp(-(u*x[i])^2/(2*s^2))
    cdf[i]<-mean(g)/s^2
  }
  cdf  
  
}

## -----------------------------------------------------------------------------
x<-seq(.1,2.5,length=10)
set.seed(2)
MC1<-MC.ray(x,s=1,antithetic = FALSE)
MC2<-MC.ray(x,s=1)
knitr::kable(round(rbind(x,MC1,MC2),4))

## -----------------------------------------------------------------------------
m<-1000
MC1<-MC2<-numeric(m)
x<-1.5
for(i in 1:m){
  MC1[i]<-MC.ray(x,s=1,R=1000,antithetic = FALSE)
  MC2[i]<-MC.ray(x,s=1,R=1000)
}
print(sd(MC1))
print(sd(MC2))
print((var(MC1)-var(MC2))/var(MC1))

## -----------------------------------------------------------------------------
m<-10000
theta.hat<-se<-numeric(2)
g<-function(x){
  x^2/sqrt(2*pi)*exp(-x^2/2)
}
f1<-function(x){2*dnorm(x,1,1)*(x>=1)}
f2<-function(x){
  exp(1-x)
}
x1<-abs(rnorm(m,0,1))+1
x2<-rexp(m,1)+1
fg1<-g(x1)/f1(x1)
fg2<-g(x2)/f2(x2)
theta.hat[1]<-mean(fg1)
theta.hat[2]<-mean(fg2)
se[1]<-sd(fg1)
se[2]<-sd(fg2)
rbind(theta.hat,se)#得到估计值和标准差

## -----------------------------------------------------------------------------
x<-seq(1,4,0.01)
par(mfrow=c(1,2))
#figue(a)
plot(x,g(x),type="l",main="figure(a)",ylim=c(0,1),lwd=2)
lines(x,f1(x),lty=2,lwd=2,col="red")
lines(x,f2(x),lty=3,lwd=2,col="blue")
legend("topright",legend=c("g",1:2),lty=1:3,lwd=2,inset=0.02)
#figue(b)
plot(x,g(x),type="l",main="figure(b)",ylim=c(0,2),lwd=2)
lines(x,g(x)/f1(x),lty=2,lwd=2,col="red")
lines(x,g(x)/f2(x),lty=3,lwd=2,col="blue")
legend("topright",legend=c("g",1:2),lty=2:3,lwd=2,inset=0.02)


## -----------------------------------------------------------------------------
print(theta.hat[1])

## -----------------------------------------------------------------------------
set.seed(1)
n<-20
I<-replicate(1000,expr={
  x<-rchisq(n,df=2)
  se<-sqrt(var(x))
  UCL1<-mean(x)-se/sqrt(n)*qt(.975,df=n-1)#CI lower
  UCL2<-mean(x)+se/sqrt(n)*qt(.975,df=n-1)#CI upper
  (UCL2>2)*(UCL1<2)# mean 2 of chisq is in CI or not
})
sum(I)
mean(I)


## -----------------------------------------------------------------------------
set.seed(2)
n<-c(10,20,30,50,100,500,1000)#a vector of sample sizes
p1.hat<-p2.hat<-p3.hat<-numeric(length(n))#to store TYPE I ERROR for different n
m<-10000#num. repl. each sim.
for(i in 1:length(n)){
  p3<-p2<-p1<-numeric(m)
  se3<-se2<-se1<-numeric(m)
  for(j in 1:m){
    x1<-rchisq(n[i],df=1)
    x2<-runif(n[i],0,2)
    x3<-rexp(n[i],rate=1)
    ttest1<-t.test(x1,alternative="two.sided",mu=1)
    ttest2<-t.test(x2,alternative="two.sided",mu=1)
    ttest3<-t.test(x3,alternative="two.sided",mu=1)
    p1[j]<-ttest1$p.value
    p2[j]<-ttest2$p.value
    p3[j]<-ttest3$p.value
  }
  p1.hat[i]<-mean(p1<0.05)
  p2.hat[i]<-mean(p2<0.05)
  p3.hat[i]<-mean(p3<0.05)
}
se1<-sqrt(p1.hat*(1-p1.hat)/m)
se2<-sqrt(p2.hat*(1-p2.hat)/m)
se3<-sqrt(p3.hat*(1-p3.hat)/m)
options(scipen=10)
print(round(rbind(n,p1.hat,se1,p2.hat,se2,p3.hat,se3),7))

## -----------------------------------------------------------------------------
n<-c(20,30,50,100,500)#n是试验中样本量组成的向量
d=3
cv<-qchisq(.95,df=d*(d+1)*(d+2)/6)

## -----------------------------------------------------------------------------
library(MASS)
sk<-function(x){
  Sigma.hat<-(nrow(x)-1)*cov(x)/nrow(x)#协方差矩阵的极大似然估计
  sigma.ni<-ginv(Sigma.hat)
  x.new=x
  for(i in 1:d){
    x.new[,i]<-x[,i]-mean(x[,i])
  }
  y<-(x.new%*%sigma.ni%*%t(x.new))^3
  nrow(x)*mean(y)/6
}

## -----------------------------------------------------------------------------
m<-1e4
set.seed(111)
mean<-rep(0,d);sigma<-diag(d)#标准正态分布的均值和协方差
p.reject<-numeric(length(n))
for(i in 1:length(n)){
     p.reject[i]<-mean(replicate(m,expr={
       x<-mvrnorm(n[i],mean,sigma)#抽n[i]个向量样本,服从多元标准正态分布
       as.integer(sk(x)>=cv)}))
} 
print(rbind(n,p.reject))

## -----------------------------------------------------------------------------
set.seed(1)
n<-30
m<-2500
epsilon<-c(seq(0,0.15,0.01),seq(0.15,1,0.05))
N<-length(epsilon)
pwr<-numeric((N))#储存不同epsilon对应的功效值
mu2<-mu1<-c(0,0,0)
s1<-diag(3)
s2<-diag(3)*100
for(i in 1:N){
  e<-epsilon[i]
  pwr[i]<-mean(replicate(m,expr={
    r<-sample(c(0,1),size=n,replace=TRUE,prob=c(1-e,e))
    x<-(1-r)*mvrnorm(1,mu1,s1)+ r*mvrnorm(n,mu2,s2)
    as.integer(sk(x)>=cv)
  }))
}

#plot power(pwr) vs epsilon
plot(epsilon,pwr,type="b",xlab=bquote(epsilon),ylim=c(0,1))
abline(h=0.05,lty=3)
se<-sqrt(pwr*(1-pwr)/m)#add se
lines(epsilon,pwr+se,lty=3)
lines(epsilon,pwr-se,lty=3)
epsilon[which(pwr==max(pwr))]#power最大值对应的epsilon

## -----------------------------------------------------------------------------
set.seed(1)
library(bootstrap)#for the score data
library(MASS)
library(boot)
#write a function to compute the theta
theta.boot<-function(score,i){
  x<-score[i,]
  n<-nrow(score)#sample size
  #cov(x)-unbised estimator,need to *(n-1)/n to get MLE
  lmd<-eigen((n-1)*cov(x)/n)$values#vector of values
  lmd[1]/sum(lmd)#return estimator of theta according to score
}

#bootstrap
B<-1000#number of replicates
n<-nrow(scor)#sample size
lmd.hat<-eigen((n-1)*cov(scor)/n)$values
theta.hat<-lmd.hat[1]/sum(lmd.hat)
R<-numeric(B)
result<-boot(data=scor,statistic=theta.boot,R=1000)
theta_b<-result$t#extract the replicates in $t
bias_b<-mean(theta_b)-theta.hat# Bootstrap Estimation of Bias
se_b<-sd(theta_b)#Bootstrap Estimation of Standard Error
print(cbind(theta.hat,theta_boots=mean(theta_b),bias_boots=bias_b,se_boots=se_b))

## -----------------------------------------------------------------------------
#JACKknife
set.seed(2)
theta.jack<-numeric(n)#store JACKknife estimation of theta
for(i in 1:n){
  x<-scor[-i,]#delet the ith data to get a set of JACKKNIFE samples
  lmd.J<-eigen((n-1)*cov(x)/n)$values#values of MLE of Sigma 
  theta.jack[i]<-lmd.J[1]/sum(lmd.J)
}
theta_j<-mean(theta.jack)
bias_j<-(n-1)*(mean(theta.jack)-theta.hat)
se_j<-sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))
print(cbind(theta.hat,theta_j,bias_j,se_j))

## -----------------------------------------------------------------------------
#result is by bootstrap in 7.7
print(result)
print(boot.ci( result,conf = 0.95,type=c("perc","bca")))

## -----------------------------------------------------------------------------
m=500#Mont carlo study times
n=100#sample size
sk.chi<-sqrt(8/5)#skweness of Chi-square distribution with df=5
sk<-function(x){
  #compute the sample skewness coeff
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  return(m3/m2^1.5)
}
skew.boot<-function(data,i){
  #function to compute the statistic
  x<-data[i,]
  sk(x)
}
x.cp<-y.cp<-numeric(m)#store decision values (0or1) for theta is in CI or not. 
norm_prob<-numeric(9)#store coverage prob and proportion of times of miss on the two sides about N(0,1)
chis_prob<-numeric(9)#store coverage prob and proportion of times of miss on the two sides about N(0,1)
set.seed(3)
x.cp<-replicate(500,expr={
  x<-rnorm(n,0,1)
  x<-matrix(x)
  boot.x<-boot(x,statistic = skew.boot,R=1000)
  ci=boot.ci(boot.x,conf = 0.95,type=c("norm","basic","perc"))
  x.p<-(ci$perc[4]<0)*(0<ci$perc[5])#theta=0 is in perc CI or not
  x.b<-(ci$basic[4]<0)*(0<ci$basic[5])#0 is basic CI or not
  x.n<-(ci$norm[2]<0)*(0<ci$norm[3])#0 is in normal CI or not
  x.pleft<-(ci$perc[4]>=0)*1;x.pright<-(ci$perc[5]<=0)*1
  x.bleft<-(ci$basic[4]>=0)*1;x.bright<-(ci$basic[5]<=0)*1
  x.nleft<-(ci$normal[2]>=0)*1;x.nright<-(ci$normal[3]<=0)*1
  list(x.n,x.b,x.p,x.nleft,x.bleft,x.pleft,x.nright,x.bright,x.pright)
}) 
for(i in 1:9){
  norm_prob[i]<-mean(as.numeric(x.cp[i,]))
}
y.cp<-replicate(500,expr={
  y<-rchisq(n,df=5)
  y<-matrix(y)
  boot.y<-boot(y,statistic = skew.boot,R=2000)
  ci=boot.ci(boot.y,conf = 0.95,type=c("norm","basic","perc"))
  y.b<-(ci$basic[4]<sk.chi)*(sk.chi<ci$basic[5])#theta= is basic CI or not
  y.n<-(ci$norm[2]<sk.chi)*(sk.chi<ci$norm[3])#0 is in normal CI or not
  y.p<-(ci$perc[4]<sk.chi)*(sk.chi<ci$perc[5])#0 is in perc CI or not
  y.pleft<-(ci$perc[4]>=sk.chi)*1;y.pright<-(ci$perc[5]<=sk.chi)*1
  y.bleft<-(ci$basic[4]>=sk.chi)*1;y.bright<-(ci$basic[5]<=sk.chi)*1
  y.nleft<-(ci$normal[2]>=sk.chi)*1;y.nright<-(ci$normal[3]<=sk.chi)*1
  list(y.n,y.b,y.p,y.nleft,y.bleft,y.pleft,y.nright,y.bright,y.pright)
}) 
for(i in 1:9){
  chis_prob[i]<-mean(as.numeric(y.cp[i,]))
}
coverprob=matrix(c(norm_prob[1:3],chis_prob[1:3]),byrow=TRUE,nrow=2)
norm_miss=matrix(norm_prob[4:9],byrow=TRUE,nrow=2)
chis_miss=matrix(chis_prob[4:9],byrow=TRUE,nrow=2)
row.names(coverprob)<-c("normal_distrib","chis_distrib")
row.names(norm_miss)<-row.names(chis_miss)<-c("left_miss_proportion","right_miss_proportion")
colnames(coverprob)<-colnames(norm_miss)<-colnames(chis_miss)<-c("norm CI","basic CI","perc CI")

## -----------------------------------------------------------------------------
print(coverprob)#coverage probability of three kinds CI for two different distributions

## -----------------------------------------------------------------------------
print(norm_miss)#proportion of times of miss on the two sides about N(0,1)
print(chis_miss)#proportion of times of miss on the two sides about chis-distribution with df=5

