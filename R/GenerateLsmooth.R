#rm(list = ls())
##########################################################
##This function generate piecewise linear functions,######
##with random slope and intercept, given number of knots##
##########################################################
piecewise<-function(p,knots=3,len,sd=0.07,a=1,b=2,seed=25)
{
  len=len
  k<-numeric(len)
  knots = knots
  sd<-sd
  a<-a
  b<-b
  l<-ceiling(length(k)/knots)
  c<-0
  for ( i in 1:knots)
  {
    set.seed(seed+i)
    if(i!= knots)
    {
      k[(1+c*l):(l+c*l)]=runif(1,-a,a)+runif(1,-b,b)*(((c*l+1):((c+1)*l))/p)+ rnorm(1:(l), sd=sd)
    }else{
      k[(1+(i-1)*l):(length(k))]=runif(1,-a,a)+runif(1,-b,b)*(((1+(i-1)*l):(length(k)))/p)+ rnorm(((1+(i-1)*l):(length(k))), sd=sd)
    }
    c= c+1
  }
  return(k)
}

########################################################################
### This function generate random Markov process########################
### discribed in the paper##############################################
########################################################################

genmarkov<-function(size,prob=0.9,b=0.5,sigma=0.5,seed=25)
{
  y=v=x=c()
  z= rnorm(size,0,sigma)
  v[1]=runif(1,-b,b)
  x[1]=0
  y[1]=x[1]+z[1]
  for( t in 2:size)
  {
    set.seed(seed+t)
    indic<-rbinom(1,1,prob)
    v[t]=indic*(v[t-1])+(1-indic)*runif(1,-b,b)
    x[t]=x[t-1] + v[t-1]
    y[t]= x[t]+z[t]
  }
  return<-list(y=y,x=x)
  return(return)
}
#############################################################################3
##############################################################################
#############################################################################
fused1d<-function(len,coef1,coef2,coef3)
{
#  len=p-3;seed=25;knots=3
  knots= 3
  k<-numeric(len)
  knots = knots
  half<-ceiling(length(k)/2)
  threequart<-ceiling(3*length(k)/4)
  rest = length(k)- threequart+1
  l<-ceiling(length(k)/knots)
  c=0
  if(l>2)
  {
  for ( i in 1:knots)
  {
    if(i== 1)
    {
      k[(1):(half)]= coef1#rnorm(1:(l), sd=0.02) ##runif(1,-1,1)++runif(1,-2,2)*(((c*l+1):((c+1)*l))/p)
    }else if (i==2){
      k[(1+half):(threequart)]=coef2# runif(1,0.2,0.4)+rnorm(1:(l), sd=0.02) 
    }else{
      k[(threequart+1):(length(k))]=coef3#runif(1,0.5,0.7)+ rnorm((1:(length(k)-(i-1)*l)), sd=0.02) #+runif(1,-2,2)*(((1+(i-1)*l):(length(k)))/p)
    }
    c= c+1
  }
  }
  return(k)
}

###############################################
######This funciton generates piecivise AR model
##############################################3
piecewiseAR<-function(p,a=-0.9,b=0.9,seed=25)
{
  len=p
  k<-numeric(p)
  knots = 3
  a<-a
  b<-b
  half<-ceiling(length(k)/2)
  threequart<-ceiling(3*length(k)/4)
  rest = length(k)- threequart+1
  k[(1):(half)]=arima.sim(list(order=c(1,0,0), ar=a), n=half)
  k[(half+1):(threequart)]=arima.sim(list(order=c(1,0,0), ar=b), n=(threequart-half))
  k[(threequart+1):(length(k))]=arima.sim(list(order=c(1,0,0), ar=a), n=(length(k)-threequart))
  return(k)
}
#################################################################################
#################################################################################
#################################################################################

#' Title
#'
#' @param p    Number of observations
#' @param band Number of band
#' @param scaled If TRUE the scaled \code{L} is returned
#' @param case   case as in Dallakyan and Pourahmadi 2019 paper.
#' @param seed   
#' @param ... 
#'
#' @return Returns Standard Cholesky factor \code{L} and modified Cholesky Factor \code{T}
#' @export
#'
#' @examples
#' generateL(p = 10, band = 4, case = "c")
generateL<-function(p,band=NULL,scaled=FALSE,case=c("b","c","d"),seed=25,...)
{
  p=p
 # knots =knots
  if(is.null(band))
  {
    band =p
  }else{
  band =band
  }
  case = match.arg(case)
  ii <- toeplitz(1:p)
  ii[upper.tri(ii)]<-0
  K <- band + 1
  L <- (ii <= K) &(ii !=0)
  for(i in 2:(band+1))
  {
    if (i==p-1)
    {
      L[p,1]=0#diag(T[-(1),-(p)])[1]
    }else{
      set.seed(seed)
 #   if(case=="a")
#    {
#      L[ii == i]=piecewise(p=p,len=p-i,seed=seed,...)/scale#diag(T[-(1),-(p)])[1:(p-i)]
#    }
#    if(case=="b")
#    {
    #   size= p-i+1
    #   L[ii == i]= (1 - abs((1:size)-size/2) + rnorm(size, sd=sd))/scale#diag(T[-(1),-(p)])[1:(p-i)]
    # }
     if(case=="d")
    {
      size= p-i+1
      L[ii == i]=genmarkov(size,seed=seed+i,...)$y/scale#diag(T[-(1),-(p)])[1:(p-i)]
    }
     if(case=="c")
     {
       L[ii == i]=(2*((1:(p-i+1))/((p)))^2- 0.5)/5 
     }
    if(case=="b")
    {
      if (i==2)
      {
      L[ii == i]=fused1d(p-i+1,0.9,1.69,1.32)
      }else if (i==3)
      {
        L[ii == i]=fused1d(p-i+1,0,-0.81,-0.81)
      }else
      {
        L[ii == i]=0
      }
    }
    # if(case=="e")
    # {
    #   L[ii == i]=fused1d(knots=3,len=p-i+1,seed=seed)
    # }
    if(case=="a")
    {
      size= p-i+1
      sign <- ((runif(1, 0, 1) > 0.5) - 1/2) * 2
      l<-sign*runif(1,0.1,0.3)
      L[ii == i]=l
    } 
    }
  }
  diag(L) <- rep(0, p)
  L <- L - upper.tri(L) * L
  T <- L
  L <- diag(rep(1, p)) - L
  if(isTRUE(scaled))
  {
  L <- diag(rep(1, p)) %*% L
  }else{
     L <- diag(1/log(8*(1:p)/50+2))%*% L ##diag(1/runif(p, 0.5, 2)) %*%L
  }
  output=list(L=L,T=T)
  return(output)
}

