#' @title Bootunif Sampling
#' @description Generate new samples based on the existing samples. Different from the Bootsrap method, we do not take the original samples, but generates new random numbers from the uniform distribution between two samples .
#' @param x the original samples with size \code{n} as a numeric vector
#' @examples 
#' \dontrun{
#' x<-rnorm(99,1,5)
#' y<-bootunif.sample(x)
#' mu=mean(y)
#' bias=mu-mean(x)
#' }
#' @importFrom stats rnorm runif
#' @export
bootunif.sample<-function(x){
  if(!is.numeric(x))
    stop("data must be a numeric vector!\n")
  if(!(is.atomic(x)||is.list(x)))
    stop("data must be a numeric vector!\n")
  n<-length(x)
  x0<-sort(x)
  data<-sapply(1:(n-1), function(i){
    runif(1,x0[i],x0[i+1])
  })
  data1<-c(x0[1],data,x0[n])
  return(data1)
}

#' @title Compute statistic with bootunif method
#' @description Generate R replicates of a statistic applied to data. 
#' @param data The data as a numeric vector.
#' @param statistic A function which when applied to data returns a vector containing the statistic(s) of interest.The first argument to statistic must be the data. No need to pass the indcices to the statistic.
#' @param R The number of sampling replicates. Usually this will be a single positive integer. 
#' @return The returned value is a list, containing the following components: \code{t0} The observed value of statistic applied to data;\code{t} A matrix with sum(R) rows each of which is a bootstrap replicate of the result of calling statistic;\code{bootunif.Statistics} A matrix like the Bootstrap Statistics. 
#' @examples 
#' \dontrun{
#' x<-rnorm(19,2,3)
#' y<-bootunif(x,mean,1000)
#' meanx<-function(x,i){
#' x<-x[i]
#' mean(x)
#' }
#' z<-boot(x,meanx,1000)
#' y$bootunif.Statistics
#' z
#' }
#' @importFrom boot boot
#' @importFrom stats rnorm runif sd
#' @export
bootunif<-function(data,statistic,R){
  if(!is.numeric(data))
    stop("data must be a numeric vector!\n")
  if(!(is.atomic(data)||is.list(data)))
    stop("data must be a numeric vector!\n")
  data0<-replicate(R,expr = {
    list(bootunif.sample(data))
  })
  res0<-lapply(list(data), as.function(statistic))
  res<-lapply(data0, as.function(statistic))
  t0<-unlist(res0)
  t<-unlist(res)
  mu=mean(t)
  b=mu-t0
  s=sd(t)
  out<-matrix(c(t0,b,s),nrow = 1,dimnames = list("values",c("original","bias","std.error")))
  t=matrix(t)
  return(list(t0=t0,t=t,bootunif.Statistics=out))
}

#' @title Gauss kernel
#' @description Kernel estimate using Gauss density function and choose appropriate bandwidth.
#' @param x the estimated value
#' @param h bandwith
#' @param data samples
#' @return a vector of Gauss estimators for \code{x}
#' @examples 
#' \dontrun{
#' y<-sample(0:1,size=1000,replace=TRUE,prob = c(0.7,0.3))
#' data1<-y*rnorm(1000,0,1)
#' data2<-(1-y)*rnorm(1000,1,0.3)
#' data<-data1+data2
#' x0<-c(1,2.5,-0.5)
#' fh.hat(x0,0.1,data)
#' }
#' @importFrom stats rnorm dnorm
#' @export
fh.hat<-function(x,h,data){
  if(!is.numeric(data))
    stop("data must be a numeric vector!\n")
  if(!(is.atomic(data)||is.list(data)))
    stop("data must be a numeric vector!\n")
  f<-sapply(1:length(x), function(i)(1/h)*mean(dnorm((x[i]-data)/h)))
  return(f)
}

#' @title Gauss leave-one-out CV
#' @description choose appropriate bandwidth for Gauss kernel estimator.
#' @param h a vector as an initial interval for bandwiths 
#' @param x a vector of samples
#' @return a number of optimal bandwith in \code{h}.
#' @examples 
#' \dontrun{
#' y<-sample(0:1,size=1000,replace=TRUE,prob = c(0.7,0.3))
#' data1<-y*rnorm(1000,0,1)
#' data2<-(1-y)*rnorm(1000,1,0.3)
#' data<-data1+data2
#' h<-c(0.005,0.25)
#' CVh(h,data)
#' }
#' @importFrom stats rnorm
#' @export
CVh<-function(h,x){
  if(!(length(h)==2))
    stop("h must be a vector with length 2 as initial intervals!")
  if(!is.numeric(x))
    stop("data must be a numeric vector!\n")
  if(!(is.atomic(x)||is.list(x)))
    stop("data must be a numeric vector!\n")
  hs<-seq(h[1],h[2],by=0.1)
  wid=0.1
  cv<-sapply(hs,function(h){
    e<-sapply(1:length(x), function(i){
      (x[i]-fh.hat(x[i],h,x[-i]))^2
    })
    mean(e)
  })
  n<-length(hs)
  pos<-which.min(cv)
  while((pos==1)||(pos==n)||(wid>0.005)){
    hs<-seq(hs[pos]-wid,hs[pos]+wid,by=wid/10)
    if(!((pos==1)||(pos==n)))wid=wid/10
    n<-length(hs)
    cv<-sapply(hs,function(h){
      e<-sapply(1:length(x), function(i){
        (x[i]-fh.hat(x[i],h,x[-i]))^2
      })
      mean(e)
    })
    pos<-which.min(cv)
  }
  v<-hs[pos]
  return(v)
}






