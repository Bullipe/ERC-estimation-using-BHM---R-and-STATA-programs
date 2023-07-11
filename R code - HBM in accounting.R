##Hierarchical Bayesian models in accounting: A tutorial 
##Bullipe R. Chintha and Sanjay Kallapur
##This R-code builds hierarchical Bayesian models described in the paper as HBM1, HBM2 and HBM3
##Open this file in the working directory where the input files are created by the accompanied Stata do file are stored

#Load packages 
library(haven)
library(tidyverse)
library(readxl)
library(broom)
library(fuzzyjoin)
library(dplyr)

fulldat <- read.csv("Input_for_HBM2_HBM3.csv")
coeff_est <- read.csv("Input_for_HBM1.csv")

###HBM1########################

N_firms <- nrow(coeff_est) #a scalar data value of the number of firms, 
#equaling 300 in this case. coeff_est is the data table containing OLS estimates, with one row per firm
b_ols <- coeff_est$ERC_ols #a vector of 300 values
s_b_ols <- coeff_est$ERC_ols_sd #vector, 300
df_b_ols <- coeff_est$df_resid #vector, 300
data <- list ("N_firms", "b_ols", "s_b_ols", "df_b_ols")

parameters <- c("b", "mu_b", "sigma_b")

##Model 1
model_hbm1_code <- function() {
  for (j in 1:N_firms) {
    b_ols[j] ~ dt(b[j], tau_b_ols[j], df_b_ols[j]) #t distr for each b_ols 
    #OpenBUGS requires the scale parameter tau = 1/s^2 instead of s or s^2. 
    #Alternatively if normal distribution were assumed for b_ols,
    # we would instead write b_ols[j] ~ dnorm (b[j], tau_b_ols[j])
    b[j] ~ dnorm(mu_b, tau_b) #distribution for unknown b[j]
    tau_b_ols[j] <- pow(s_b_ols[j], -2) #definition of tau_b_ols
  }
  mu_b ~ dnorm (0.0, 1.0E-6) #uninformative prior for mu_b, with 
  #high variance = 1000000, i.e., sd=1000
  #high is defined in the context of expected mu_b values, which in 
  #the case of ERC rarely exceed 10 or 15. The prior gives 68%
  #probability to values between 0-1000 and 0+1000, a wide range, and 
  #therefore uninformative
  tau_b <- pow(sigma_b, -2) #definition of tau_b
  sigma_b ~ dunif(0, 1000) #uninformative uniform prior for sigma_b 
  #in the range 0 to 1000
}

inits <- function(){
  list(b=rnorm(N_firms, 0, 100), mu_b=rnorm(1, 0, 100),
       sigma_b=runif(1, 0, 100))
} #initial values of the parameters need to be specified. Here we are giving   #random values because we want them to be different in each chain (see below)

library(R2OpenBUGS) #load the package
hbm1_bugs <- bugs(data, inits, parameters, model.file=model_hbm1_code,
                  n.chains=3, n.iter=5000, debug = FALSE) 

print(hbm1_bugs)

OLS <- coeff_est$ERC_ols
HBM1 <- hbm1_bugs$mean$b
MSE1 <- sum((coeff_est$ERC_ols_post - HBM1)^2)/length(HBM1)
print(MSE1)
print(sum(HBM1<0))

#################################################
###HBM2########################

y<-as.vector(fulldat$return) #returns data
n<-length(y) #=4500 (300 firms running for 15 years)
x<-as.vector(fulldat$delta_E) #change in earnings data
gvkey<-as.vector(fulldat$gvkey) 
Firm.name<-unique(gvkey) #creating firm identifier
J<-length(Firm.name) #=300 (Number of firms)
Firm<-rep(NA,J)
for (i in 1:J){
  Firm[gvkey==Firm.name[i]]<-i
}
set.seed(149)

##Model 2
model_hbm2_code<- function(){
  for(i in 1:n){
    y[i]~dnorm(y.hat[i],tau.y)
    y.hat[i]<-a[Firm[i]]+b[Firm[i]]*x[i] #We define the varying intercept, varying slope model here
  }
  tau.y<-pow(sigma.y,-2) #definition of tau.y
  sigma.y~dunif (0,100) #uninformative uniform prior
  
  for (j in 1:J){
    a[j]<-xi.a*B.raw[j,1] #Clubbing a[j] and b[j] coefficients into matrix B
    b[j]<-xi.b*B.raw[j,2] 
    B.raw[j,1:2]~dmnorm(B.raw.hat[j,], Tau.B.raw[,]) #Hyper-parameters
    B.raw.hat[j,1]<-mu.a.raw
    B.raw.hat[j,2]<-mu.b.raw
  }
  mu.a<-xi.a*mu.a.raw
  mu.b<-xi.b*mu.b.raw
  mu.a.raw~dnorm(0,0.0001) #uninformative normal 
  mu.b.raw~dnorm(0,0.0001) #uninformative normal 
  
  xi.a~dunif(0,100) #uninformative uniform prior
  xi.b~dunif(0,100) #uninformative uniform prior
  
  Tau.B.raw[1:2,1:2]~dwish(W[,],df) #Define the inverse-Wishart prior distribution for the variance-covariance matrix
  df<-3
  Sigma.B.raw[1:2,1:2]<-inverse(Tau.B.raw[,])
  sigma.a<-xi.a*sqrt(Sigma.B.raw[1,1])
  sigma.b<-xi.b*sqrt(Sigma.B.raw[2,2])
  rho<-Sigma.B.raw[1,2]/sqrt(Sigma.B.raw[1,1]*Sigma.B.raw[2,2])
}

library("MCMCpack")

W<-diag(2)
data <-list ("n","J","y","Firm","x","W")
inits<-function(){
  list(B.raw=array(rnorm(2*J),c(J,2)),mu.a.raw=rnorm(1),
       mu.b.raw=rnorm(1),sigma.y=runif(1), 
       Tau.B.raw=rwish(3,diag(2)),
       xi.a=runif(1),xi.b=runif(1))
}

parameters<-c("a","b","mu.a","mu.b","sigma.y","sigma.a","sigma.b","rho")

library(R2OpenBUGS)
hbm2_bugs<-bugs(data, inits,parameters,model.file=model_hbm2_code,n.chains=3,n.iter=5000, debug=FALSE)

HBM2 <- hbm2_bugs$mean$b
MSE2 <- sum((coeff_est$ERC_ols_post - HBM2)^2)/length(HBM2)
print(MSE2)
print(sum(HBM2<0))

####################################
###HBM3########################

y<-as.vector(fulldat$return) #returns data
n<-length(y) #=4500 (300 firms running for 15 years)
x<-as.vector(fulldat$delta_E) #change in earnings data
x1<- as.vector (fulldat$rf) #load data of risk free rates
gvkey<-as.vector(fulldat$gvkey)
Firm.name<-unique(gvkey) #creating firm identifier
J<-length(Firm.name)
Firm<-rep(NA,J)
for (i in 1:J){
  Firm[gvkey==Firm.name[i]]<-i
}
set.seed(149)

##Model 3
model_hbm3_code<- function(){
  for(i in 1:n){
    y[i]~dnorm(y.hat[i],tau.y)
    y.hat[i]<-a[Firm[i]]+b[Firm[i]]*x[i]+b.0*x1[i]*x[i] #We define the varying intercept, varying slope model here. Note that we have a new coefficient b.0
  }     #b.0 captures the association between ERCs and the risk-free rates. The division by 100 is to convert the percentage rates
  tau.y<-pow(sigma.y,-2) #definition of tau.y
  sigma.y~dunif (0,100) #uninformative uniform prior
  
  for (j in 1:J){
    a[j]<-xi.a*B.raw[j,1] #Clubbing a[j] and b[j] coefficients into matrix B
    b[j]<-xi.b*B.raw[j,2]
    B.raw[j,1:2]~dmnorm(B.raw.hat[j,], Tau.B.raw[,]) #Hyper-parameters
    B.raw.hat[j,1]<-mu.a.raw
    B.raw.hat[j,2]<-mu.b.raw
  }
  mu.a<-xi.a*mu.a.raw
  mu.b<-xi.b*mu.b.raw
  mu.a.raw~dnorm(0,0.0001) #uninformative normal 
  mu.b.raw~dnorm(0,0.0001) #uninformative normal 
  
  b.0~dnorm(0,0.0001) #uninformative normal
  
  xi.a~dunif(0,100)  #uninformative uniform prior
  xi.b~dunif(0,100)  #uninformative uniform prior
  
  Tau.B.raw[1:2,1:2]~dwish(W[,],df) #Define the inverse-Wishart prior distribution for the variance-covariance matrix
  df<-3
  Sigma.B.raw[1:2,1:2]<-inverse(Tau.B.raw[,])
  sigma.a<-xi.a*sqrt(Sigma.B.raw[1,1])
  sigma.b<-xi.b*sqrt(Sigma.B.raw[2,2])
  rho<-Sigma.B.raw[1,2]/sqrt(Sigma.B.raw[1,1]*Sigma.B.raw[2,2])
}

##For calling bugs from R
library("MCMCpack")

W<-diag(2)
data <-list ("n","J","y","Firm","x","x1","W")
inits<-function(){
  list(B.raw=array(rnorm(2*J),c(J,2)),mu.a.raw=rnorm(1),
       mu.b.raw=rnorm(1),sigma.y=runif(1),
       Tau.B.raw=rwish(3,diag(2)), b.0=rnorm(1),
       xi.a=runif(1),xi.b=runif(1))
}

parameters<-c("a","b","mu.a","mu.b","b.0","sigma.y","sigma.a","sigma.b","rho")

library(R2OpenBUGS)
hbm3_bugs<-bugs(data, inits,parameters,model.file=model_hbm3_code,n.chains=3,n.iter=5000, debug=FALSE)

HBM3 <- hbm3_bugs$mean$b
MSE3 <- sum((coeff_est$ERC_ols_post - HBM3)^2)/length(HBM3)
print(MSE3)
print(sum(HBM3<0))
ratecoeff <- hbm3_bugs$mean$b.0
print(ratecoeff)

##Ouput data
df<-data.frame(OLS,HBM1,HBM2,HBM3)
write.csv(df, "BayesoutputR.csv", row.names = FALSE)


