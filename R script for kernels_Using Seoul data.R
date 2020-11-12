library( xlsx )

rzichtjuni2019_v2.xlsx", 
                       sheet=1) %>%
  select( OD = ODc2, serumtype =starts_with("serum"), PCR=starts_with( "PCR") ) %>%
  ##colnames( data.all ) <- c(OD= ODc2, Serum.Type= starts_with("serum"), PCR= starts_with("PCR") ) %>%
  mutate( Serum.Type = as.factor(serumtype), PCR=as.factor(PCR), logOD=log10(OD)) %>%
  filter( PCR != "ND" , !is.na(logOD)) %>%
  droplevels()

##take the logOD values from data set
data.log <- data.all %>% mutate( logOD = log10(OD), 
                                 !is.na(logOD))

## aesthetic mapping (histogram)
p <- ggplot(data.log, aes(x=logOD))
p <- p + geom_histogram(aes(y = ..density..), colour="black", fill="red" )
p



## Supply initial values, mu_pos, mu_neg, sd_pos, sd_neg, prev 
## How best to do this??
## [mu1=-1.5, sd1=0.5, mu2=0.25, sd2=0.1]
## plot(hatf, fx) once initialised later to see this

## Plot density function (distributions)
plot( density(data.log$logOD), 
      main= "Distribution function",
      col="brown",
      ylim=c(0,2.0) ) ## set y-axis limit
##par(bty="n"), ## remove box outline

## Begin plotting the kernel fits to the logOD data
## Find best fit

plot(kdensity(data.log$logOD, start= "gaussian"),
     main= "Parametric starts", ylim= range(0, 2.0), col= "brown")
lines(kdensity(data.log$logOD, start= "gumbel"), col= "green")
lines(kdensity(data.log$logOD, start= "laplace"), col= "blue")
legend(-2.5, 2.0, legend= c("Gaussian", "Gumbel", "Laplace"), fill= c("brown", "green", "blue"))      
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")


## Begin plotting the kernel fits to the data
## Find best fit

## gaussian start
plot(kdensity(data.log$logOD, start="gaussian"),
     main= "Gaussian start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, kernel= "cosine"), col= "gold") 
lines(kdensity(data.log$logOD, kernel= "laplace"), col= "black")
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine", "Laplace"), fill= c("brown", "blue", "green", "gold", "black"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

## gumbel start
plot(kdensity(data.log$logOD, start="gumbel"),
     main= "Gumbel start", ylim= range(0, 2.0), col= "brown") ## gumbel distribution
lines(kdensity(data.log$logOD, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, kernel= "cosine"), col= "gold") 
lines(kdensity(data.log$logOD, kernel= "laplace"), col= "black")
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine", "Laplace"), fill= c("brown", "blue", "green", "gold", "black"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

## laplace start
plot(kdensity(data.log$logOD, start="laplace"),
     main= "Laplace start", ylim= range(0, 2.0), col= "brown") ## gumbel distribution
lines(kdensity(data.log$logOD, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, kernel= "cosine"), col= "gold") 
lines(kdensity(data.log$logOD, kernel= "laplace"), col= "black")
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine", "Laplace"), fill= c("brown", "blue", "green", "gold", "black"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")



## Begin building derivative orders
## For the 'esimate' function
## package ks using dkde as pdf estimate
hatf <- dkde(data.log$logOD, deriv.order = 0, na.rm=TRUE)
hatf1 <- dkde(data.log$logOD, deriv.order = 1, na.rm=TRUE)
hatf2 <- dkde(data.log$logOD, deriv.order = 2, na.rm=TRUE)
hatf3 <- dkde(data.log$logOD, deriv.order = 3, na.rm=TRUE)


## Produce derivative estimators
## Initialise x-value

x <- data.log$logOD ## initiliase x

## Need to correctly initialise mu and sigma
## It seems that two means exist: mu_neg= -1.5 & mu_pos= 0.25
## What about sigma_neg & sigma_pos ??
## ???
fx <- function(x) { 0.5 * dnorm(x,-1.5,0.5) + 0.5 * dnorm(x,0.25,0.1) }
fx1 <- function(x) { 0.5 * (-4*x-6)* dnorm(x,-1.5,0.5) + 0.5 * (-4*x+6) * dnorm(x,0.25,0.1) }
fx2 <- function(x) { 0.5 * ((-4*x-6)^2 - 4) * dnorm(x,-1.5, 0.5) + 0.5 * ((-4*x+6)^2 - 4) * dnorm(x,0.25,0.1) }
fx3 <- function(x) { 0.5 * (-4*x-6) * ((-4*x-6)^2 - 12) * dnorm(x,-1.5,0.5) + 0.5 * (-4*x+6) * ((-4*x+6)^2 - 12) * dnorm(x,0.25,0.1) }


## Plot derivative estimators
## Against the order derivative of function
plot(hatf, fx = fx)
plot(hatf1, fx = fx1)
plot(hatf2, fx = fx2)
plot(hatf3, fx = fx3)


## Find optimal 'h' bandwidth selection
## Guassian kernel
dkde.gaus <- dkde(data.log$logOD, deriv.order = 0, kernel= "gaussian")
dkde.gaus1 <- dkde(data.log$logOD, deriv.order = 1, kernel= "gaussian")
dkde.gaus2 <- dkde(data.log$logOD, deriv.order = 2, kernel= "gaussian")
dkde.gaus3 <- dkde(data.log$logOD, deriv.order = 3, kernel= "gaussian")

## Epanechnikov kernel
dkde.epan <- dkde(data.log$logOD, deriv.order = 0, kernel= "epanechnikov")
dkde.epan1 <- dkde(data.log$logOD, deriv.order = 1, kernel= "epanechnikov")
## Error for epanechnikov kernel derivative >= 2
## dkde.epan2 <- dkde(data.log$logOD, deriv.order = 2, kernel= "epanechnikov")
## dkde.epan3 <- dkde(data.log$logOD, deriv.order = 3, kernel= "epanechnikov")

## Cosine kernel
dkde.cos <- dkde(data.log$logOD, deriv.order = 0, kernel= "cosine")
dkde.cos1 <- dkde(data.log$logOD, deriv.order = 1, kernel= "cosine")
dkde.cos2 <- dkde(data.log$logOD, deriv.order = 2, kernel= "cosine")
dkde.cos3 <- dkde(data.log$logOD, deriv.order = 3, kernel= "cosine")

## Plot derivative functions
## derivative 0
plot(dkde.gaus, col= "blue")
lines(dkde.epan, lty=1, col= "pink", lwd=2)
lines(dkde.cos, lty=1, col= "blue", lwd=1)
lines(hatf, lty=1, col= "green", lwd=2)
## not much is seen of these

## derivative 1
plot(dkde.gaus1)
lines(dkde.epan1, lty=1, col= "pink")
lines(dkde.cos1, lty=1, col= "blue", lwd=1)
lines(hatf1, lty=1, col= "gold", lwd=2)
## epanechnikov with derivative>1 does not plot


## Selecting the optimal bandwidth for kdensity function
## AMISE - optimal bandwidth

## Gaussian kernel
##
h.amise.gaus0 <- h.amise(data.log$logOD, deriv.order = 0, kernel= "gaussian")
h.amise.gaus1 <- h.amise(data.log$logOD, deriv.order = 1, kernel= "gaussian")
h.amise.gaus2 <- h.amise(data.log$logOD, deriv.order = 2, kernel= "gaussian")
h.amise.gaus3 <- h.amise(data.log$logOD, deriv.order = 3, kernel= "gaussian")

## Epanechnikov kernel
##
h.amise.epan0 <- h.amise(data.log$logOD, deriv.order = 0, kernel= "epanechnikov")
h.amise.epan1 <- h.amise(data.log$logOD, deriv.order = 1, kernel= "epanechnikov")
h.amise.epan2 <- h.amise(data.log$logOD, deriv.order = 2, kernel= "epanechnikov")
h.amise.epan3 <- h.amise(data.log$logOD, deriv.order = 3, kernel= "epanechnikov")

## Cosine kernel
##
h.amise.cos0 <- h.amise(data.log$logOD, deriv.order = 0, kernel= "cosine")
h.amise.cos1 <- h.amise(data.log$logOD, deriv.order = 1, kernel= "cosine")
h.amise.cos2 <- h.amise(data.log$logOD, deriv.order = 2, kernel= "cosine")
h.amise.cos3 <- h.amise(data.log$logOD, deriv.order = 3, kernel= "cosine")

## to call on bw= 'h' examples
## h.amise.gaus0$h
## h.amise.epan0$h
## h.amise.cos0$h


## Unbiased cross validation

## Gaussian kernel
##
h.ucv.gaus0 <- h.ucv(data.log$logOD, deriv.order= 0, kernel= "gaussian")
h.ucv.gaus1 <- h.ucv(data.log$logOD, deriv.order= 1, kernel= "gaussian")
h.ucv.gaus2 <- h.ucv(data.log$logOD, deriv.order= 2, kernel= "gaussian")
h.ucv.gaus3 <- h.ucv(data.log$logOD, deriv.order= 3, kernel= "gaussian")

## Epanechnikov kernel
##
h.ucv.epan0 <- h.ucv(data.log$logOD, deriv.order= 0, kernel= "epanechnikov")
h.ucv.epan1 <- h.ucv(data.log$logOD, deriv.order= 1, kernel= "epanechnikov")
h.ucv.epan2 <- h.ucv(data.log$logOD, deriv.order= 2, kernel= "epanechnikov")
h.ucv.epan3 <- h.ucv(data.log$logOD, deriv.order= 3, kernel= "epanechnikov")

## Cosine kernel
##
h.ucv.cos0 <- h.ucv(data.log$logOD, deriv.order= 0, kernel= "cosine")
h.ucv.cos1 <- h.ucv(data.log$logOD, deriv.order= 1, kernel= "cosine")
h.ucv.cos2 <- h.ucv(data.log$logOD, deriv.order= 2, kernel= "cosine")
h.ucv.cos3 <- h.ucv(data.log$logOD, deriv.order= 3, kernel= "cosine")

## Plot the UCV functions
## gaussian
##
for (i in 0:3) 
  plot(h.ucv(data.log$logOD, deriv.order = i, kernel= "gaussian"))

## 'abline' to visually show the minimum value
## That bw 'h' will take
plot(h.ucv.gaus0)
abline(v= h.ucv.gaus0$h, col= "blue")

## epanechnikov
##
for (i in 0:1) 
  plot(h.ucv(data.log$logOD, deriv.order = i, kernel= "epanechnikov"))

## cosine
##
for (i in 0:3) 
  plot(h.ucv(data.log$logOD, deriv.order = i, kernel= "cosine"))


## Begin plotting the distribution function
## Against kernel density
## With NEW adjusted bw='h'

##
## AMISE for derivative order 0

plot(kdensity(data.log$logOD, start="gaussian"),
     main= "Gaussian start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, bw= h.amise.gaus0$h, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, bw= h.amise.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, bw= h.amise.cos0$h, kernel= "cosine"), col= "gold") 
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

plot(kdensity(data.log$logOD, start="gumbel"),
     main= "Gumbel start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, bw= h.amise.gaus0$h, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, bw= h.amise.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, bw= h.amise.cos0$h, kernel= "cosine"), col= "gold") 
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine", "Laplace"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

plot(kdensity(data.log$logOD, start="laplace"),
     main= "Laplace start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, bw= h.amise.gaus0$h, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, bw= h.amise.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, bw= h.amise.cos0$h, kernel= "cosine"), col= "gold") 
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

## The AMISE does not seem to plot new pdf well
## Perhaps to data being belonging to bimodal instead of asymptotic

##
## UCV for derivative order 0

## No start indicated with new 'bw' estimator
##
plot(kdensity(data.log$logOD, start= "gaussian"),
     main= "Gaussian start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, bw= h.ucv.gaus0$h, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, bw= h.ucv.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, bw= h.ucv.cos0$h, kernel= "cosine"), col= "gold") 
## error - will not integrate: lines(kdensity(data.log$logOD, bw= h.ucv.gaus0$h, kernel= "laplace"), col= "black")
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

plot(kdensity(data.log$logOD, start="gumbel"),
     main= "Gumbel start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, bw= h.ucv.gaus0$h, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, bw= h.ucv.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, bw= h.ucv.cos0$h, kernel= "cosine"), col= "gold") 
## error - will not integrate: lines(kdensity(data.log$logOD, bw= h.ucv.gaus0$h, kernel= "laplace"), col= "black")
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

plot(kdensity(data.log$logOD, start="laplace", kernel= "gaussian"),
     main= "Laplace start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, bw= h.ucv.gaus0$h, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, bw= h.ucv.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, bw= h.ucv.cos0$h, kernel= "cosine"), col= "gold") 
## error - will not integrate: lines(kdensity(data.log$logOD, bw= h.ucv.gaus0$h, kernel= "laplace"), col= "black")
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

## Gaussian start chosen
##
plot(kdensity(data.log$logOD, start="gaussian"),
     main= "Gaussian start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, start= "gaussian", bw= h.ucv.gaus0$h, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, start= "gaussian", bw= h.ucv.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, start= "gaussian", bw= h.ucv.cos0$h, kernel= "cosine"), col= "gold") 
## error - will not integrate:  lines(kdensity(data.log$logOD, start= "gaussian", bw= h.ucv.gaus0$h, kernel= "laplace"), col= "black")
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

## Gaussian start chosen
## with bw=0.05
##
plot(kdensity(data.log$logOD, start="gaussian"),
     main= "Gaussian start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, start= "gaussian", bw= 0.05, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, start= "gaussian", bw= h.ucv.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, start= "gaussian", bw= h.ucv.cos0$h, kernel= "cosine"), col= "gold") 
## error - will not integrate:  lines(kdensity(data.log$logOD, start= "gaussian", bw= h.ucv.gaus0$h, kernel= "laplace"), col= "black")
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

## Gumbel start chosen
##
plot(kdensity(data.log$logOD, start="gumbel"),
     main= "Gumbel start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, start="gumbel", bw= h.ucv.gaus0$h, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, start="gumbel", bw= h.ucv.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, start="gumbel", bw= h.ucv.cos0$h, kernel= "cosine"), col= "gold") 
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

## Laplace start chosen
##
plot(kdensity(data.log$logOD, start="laplace"),
     main= "Laplace start", ylim= range(0, 2.0), col= "brown") ## normal distribution
lines(kdensity(data.log$logOD, start="laplace", bw= h.ucv.gaus0$h, kernel= "gaussian"), col= "blue")
lines(kdensity(data.log$logOD, start="laplace", bw= h.ucv.epan0$h, kernel= "epanechnikov"), col= "green")
lines(kdensity(data.log$logOD, start="laplace", bw= h.ucv.cos0$h, kernel= "cosine"), col= "gold") 
## error - will not integrate: lines(kdensity(data.log$logOD, start="laplace", bw= h.ucv.gaus0$h, kernel= "laplace"), col= "black")
legend(-2.5, 2.0, legend= c("Start", "Gaussian", "Epanechnikov", "Cosine"), fill= c("brown", "blue", "green", "gold"))
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")




##further testing
## is this the best adjusted with 'h' kdenstiy distribution ??
## ???

plot(kdensity(data.log$logOD, h.ucv.gaus0$h, start= "gaussian", kernel= "gaussian"),
     main= "bw=ucv.gaussian, deriv=0", ylim= range(0, 2.0), col="red")
lines(kdensity(data.log$logOD, bw=0.05, start= "gaussian", kernel= "gaussian"), lty=2, col= "gold")
lines(kdensity(data.log$logOD, start= "gaussian"), lty=3, col= "blue")
rug(data.log$logOD, ticksize=-0.02, lwd=1, col="black")

## E Step V
param <- c(1,0.5, 1)
param <- as.data.frame(param)
plot(kdensity(data.log$logOD, kernel= "gaussian"))
lines(kdensity(data.log$logOD, h.ucv.gaus0$h, kernel= "gaussian"))
print(kdensity(data.log$logOD, kernel= "gaussian"))
newdata <- data.log$logOD
estep <- estepV(data= newdata, param)


## Dirichlet Process Gaussian
dp <- DirichletProcessGaussian(data.log$logOD)
plot(dp)
dpxtra <- Fit(dp, 100)
plot(dpxtra)
plot(dpxtra, data_method= "hist", xlim= range(-3,1))
print(dpxtra)

## np EM (non parametric, adapts kernel for EM algorithm)
plot(kdensity(data.log$logOD, kernel= "gaussian"), col= "black", xlim= range(-3, 1), ylim= range(0,0.8))
lines(kdensity(data.log$logOD, h.ucv.gaus0$h, kernel= "gaussian"), col="darkblue")

## how many clusters
## run before and after ??
numclusters <- sapply(dpxtra$weightsChain, FUN= length)
qplot(x=numclusters, geom="histogram", binwidth=1)

## From the above 6 plots
## Can see that the UCV bw='h'
## More efficient than the AMISE bw
