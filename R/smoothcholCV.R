#### Example 1 
DD <- 0:3
PP <- c(0.292, 0.525, 0.175, 0.008)

X11()
plot(DD,PP,type="h",col=2,main="Pmf from user list",xlab="x",ylab="p(x)")
points(DD,PP,col=2);abline(h=0,col=3)


### Binomial Distribution
##Example 4
n = 10
p= 0.5
dbinom(4, n, p)
####
pbinom(5, n, p)
##or
sum(dbinom(0:5, n, p))

## Example 5
pbinom(9,12,1/6) - pbinom(5,12,1/6)
## or you can use
sum(dbinom(6:9, 12 ,1/6))


## Example 7
pbinom(8,10,0.7) - pbinom(4,10,0.7)
## or you can use
sum(dbinom(5:8, 10 ,0.7))

### Hypergeometric
### Example 9
dhyper(3, 4, 8, 6)

### Example 11
1-pnbinom(1, 1, 0.75)
1 - pgeom(3,0.5)
####
#### Illustration of continuous random variables.
####

## Example 15
a <- 0.14
b <-0.71
density_to_integrate <- function(x) { 3*x^2 }
integrate(density_to_integrate, lower = a, upper = b)

## Example 16
a <- 3.4
b <-7.1
density_to_integrate <- function(x) { 3/(x^4) }
integrate(density_to_integrate, lower = a, upper = b)

## Example 18
##a)
pnorm(1.75,1.5,0.5) - pnorm(0.1,1.5,0.5)

## or using standard normals
pnorm(0.5) - pnorm(-2.8)
## or
integrate(dnorm, lower = -2.8, upper = 0.5)
##b)
1- pnorm(2.5,1.5,0.5)
## or
1- pnorm(2)
## or 
1-integrate(dnorm, lower = -Inf, upper = 2)

## Example 19
qnorm(0.95)

## Example 20
qnorm(0.95,8,0.46)
## or
8 + 0.46* qnorm(0.95)


##
## One example of a continuous distribution is the chi-square distribution. This 
## distribution has one parameter, which is called its number of "degrees of freedom". 
## Call this parameter m. Then the mean and variance of a chi-square random variable 
## equal m and 2m, respectively. 
##
## Try varying 'm'.
m <- 5

x_sim <- rchisq(1e5, df = m)
x_sim[1:20]
hist(x_sim)


x<-1:10
var(x)
## The probability of observing a value between a and b can be calculated either by 
## integrating the pdf or by taking differences using the cdf.
a <- 5
b <- 10
density_to_integrate <- function(x) { dchisq(x, m) }
integrate(density_to_integrate, lower = a, upper = b)
pchisq(b, m) - pchisq(a, m)

## Probabilities are *integrals*, areas under the curve, of the pdf.
cc <- curve(dchisq(x, m), from = min(x_sim), to = max(x_sim), xlab = "x", ylab = "f(x)")
polygon(c(a, cc$x[cc$x > a & cc$x <= b], b), c(0, cc$y[cc$x > a & cc$x <= b], 0), 
  col = "red")

## Equal to a corresponding *difference* in cdf values.
curve(pchisq(x, m), from = min(x_sim), to = max(x_sim), xlab = "x", ylab = "F(x)")
points(a, pchisq(a, m), pch = 20, cex = 2, col = "red")
points(b, pchisq(b, m), pch = 20, cex = 2, col = "red")


#### Example 11
density<-function (x) 3*x^2

integrate(density,lower = 0.14,upper = 0.71)
#### Example 12
density<-function (x) 3/x^4

integrate(density,lower = 3.4,upper = 7.1)

##
## Normal distributions.
##

## Compare N(0, 1), the "standard normal" distribution, with N(0, 2) and N(0.5, 0.75).
curve(dnorm(x, mean = 0, sd = 1), from = -3.5, to = 3.5, ylim = c(0, 0.5), lwd = 2, 
      col = "red", ylab = "f(x)")
curve(dnorm(x, mean = 0, sd = sqrt(2)), add = TRUE, lwd = 2, col = "blue")
curve(dnorm(x, mean = 0.5, sd = sqrt(0.75)), add = TRUE, lwd = 2, col = "green")

### Example 14
##1
pnorm(1.75,1.5,0.5)-pnorm(0.1,1.5,0.5)
## or 
pnorm(0.5)-pnorm(-2.8)

1- pnorm(2.5,1.5,0.5)

##Exercise 15
qnorm(0.95,8,46)


##
## Exponential Distribution
##
lambda <- 2

x_sim <- rexp(1e5, rate = lambda)
x_sim[1:20]

hist(x_sim)
plot(hist(x_sim)$density,type="l",col="red")

lambda_pos<-seq(0.1,1.8,0.3)
sim <-function(lambda)
{
  rexp(1e5, rate = lambda)
}
variable<-lapply(lambda_pos,sim)
dens<-lapply(variable,density)
names(dens)<-lambda_pos
plot(density(x_sim),type="l",col="red")
mapply(lines, dens, col=1:length(dens))
legend("topright", legend=names(dens), fill=1:length(dens))

# The probability of observing a value between a and b can be calculated either by 
## integrating the pdf or by taking differences using the cdf.
a <- 5
b <- 10
density_to_integrate <- function(x) { dexp(x, lambda) }
integrate(density_to_integrate, lower = a, upper = b)
pexp(b, lambda) - pchisq(a, lambda)


#### Example 17

density<-function (x) 0.001*exp(-0.001*x)
###a)
-exp(-0.001*1200)+exp(-0.001*900)
##or 
integrate(density,lower = 900,upper = 1200)
#### b) 
1- integrate(density,lower=0,upper=1200)$value


###Example 18
pgamma(1,3,2)






##
## Beta Distribution
##
  

curve(dbeta(x, shape1 = 5, shape2 =4), from = 0, to = 1, ylim = c(0, 6), lwd = 2, 
      col = "red", ylab = "f(x)")
curve(dbeta(x, shape1 = 7, shape2  =4), add = TRUE, lwd = 2, col = "blue")
curve(dbeta(x, shape1 =  7, shape2=5), add = TRUE, lwd = 2, col = "green")

##
## Gamma Distribution
##


curve(dgamma(x, shape = 1, rate =1), from = 0, to = 10, ylim = c(0, 2), lwd = 2, 
      col = "red", ylab = "f(x)")
curve(dgamma(x, shape = 2, rate  =1), add = TRUE, lwd = 2, col = "blue")
curve(dgamma(x, shape = 2, rate= 5), add = TRUE, lwd = 2, col = "green")
