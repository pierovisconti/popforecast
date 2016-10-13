rm(list=ls())


library(plyr)
library(plotrix)
library(MASS)
library(filzbach)

# ##Figure1 ---------------------------------------------------------------
forest=seq(1,0,-0.01)
plot(1:101,100*(1-forest)^1,type='l',lty=1, xlab='% Forest cover',ylab='% Carrying Capacity')
lines(100*(1-forest)^2,type='l',lty=6)
lines(100*(1-forest)^5,type='l',lty=5)
lines(100*(1-forest)^0.8,type='l',lty=4)
lines(100*(1-forest)^0.5,type='l',lty=3)
lines(100*(1-forest)^0.2,type='l',lty=2)


legend(0, 98,c(0.2,0.5,0.8,1,2,5),bty="n",lty=c(2,3,4,1,5,6))
text(6,99,paste(expression(phi), " values"))
dev.copy(pdf,'Phi values carrying capacity.pdf')
dev.off()



# Figures2-4 --------------------------------------------------------------
# this code with this parameters was used to produce fig 2,3,4 in github


# CREATE "REAL" PARAMETERS FOR THE LOGISTIC GROWTH MODEL
No<-10 # initial pop
# k<-runif(15,95,105) # variable carrying capacity
# k
# k<-100 # fixed karrying capacity
real_a=2 # lambda
t<-1:20 # time-step vector
forest=pmin(1,seq(0.76,0.38,-0.02)+ runif(max(t),-0.02,0.02)) #constantly declining forest cover with random fluctuation
theta=1.5 #scaling of carrying capacity respect to forest cover, slightly sublinear
sigma=0.05 ## std of demographic stochasticity
mu=0 # when !=0 the growth rate is deterministically pushed up or down (e.g. other threats exept deforestation, or conservation management)

# define logistic growth model
f_var=function(k,No,a,t){
  N=matrix(nrow=max(t))
  for (i in 1:max(t)){
    if (i==1){
      N[i]=No+(No*a*(1-No/k[i]))
    } else {
      N[i]=N[i-1]+N[i-1]*a*(1-N[i-1]/k[i])
    }
  }
  return(N)
}


# define logistic growth model with k function of forestcover
f_var_for=function(forest,t,No,a,theta){
  N=matrix(nrow=max(t))
  for (i in 1:max(t)){
    k=100*forest[i]^theta
    if (i==1){
      N[i]=No+(No*a*(1-No/k))
    } else {
      N[i]=N[i-1]+N[i-1]*a*(1-N[i-1]/k)
    }
    if (N[i]<=0){
      N[i:max(t)]=0; break
    }
  }
  return(N)
}


# define logistic growth  model with k function of forestcover + noise and bias in growth rate
# this accounts for other factors that may trigger decline that aren't related with variable K
f_var_for_rand=function(forest,t,No,a,theta,mu,sigma){
  N=matrix(nrow=max(t))
  for (i in 1:max(t)){
    k=100*forest[i]^theta
    if (i==1){
      N[i]=No+(No*(a*(1-No/k)+rnorm(1,mean=mu,sd=sigma)))
    } else {
      N[i]=N[i-1]+(N[i-1]*(a*(1-N[i-1]/k)+rnorm(1,mean=mu,sd=sigma)))
    }
    if (N[i]<=0){
      N[i:max(t)]=0; break
    }
  }
  return(N)
}



# define logistic growth  model with fixed k
f_fix=function(k,t,No,a){
  N=matrix(nrow=max(t))
  for (i in 1:max(t)){
    if (i==1){
      N[i]=No+(No*a*(1-No/k))
    } else {
      N[i]=N[i-1]+N[i-1]*a*(1-N[i-1]/k)
    }
  }
  return(N)
}


# create populations data from set parameters
pop1<-f_var_for(forest,t,No,real_a,theta)#create "real" population without stochasticity

pop<-f_var_for_rand(forest,t,No,real_a,theta,mu,sigma) # create "real" population with stochsticity
plot(pop,ylim=c(0,70))
points(pop1,col='red')



# LOG-LIK FUNCTIONS -------------------------------------------------------

## DEFINE LOG-LIK FUNCTIONS AND FUNCTIONS FOR RESIDUALS
# here there is no stochasticity so the residuals rather than the likelihood are minimized 
# function to calculate residuals from observed and modelled pop data with unknown growth rate, carrying cap and initial population
residuals=function(a,k,No){
  res=-sqrt(colSums((pop-f_fix(k,t,No,a))^2))
  return(res)
} 

# function to calculate residuals when both growth rate and theta (relation between forest cover and K) are unknown
residuals2=function(a,theta){
  res=-sqrt(colSums((pop-f_var_for(forest,t,No,a,theta))^2))
  return(res)
} 

## function to obtain loglikelihood when growth rate, theta, and standard deviation of white noise around nominal growth are to be fond
loglik3=function(a,theta,sigma){
  res=matrix(nrow=nobs-1)
  prob=matrix(nrow=nobs-1)
  for (te in 2:nobs){
    k=100*forest[te]^theta
    res[te-1]=(pop[te]/pop[te-1])-1-(a*(1-pop[te-1]/k))
    prob[te-1]=dnorm(res[te-1],mean=0,sd=sigma)
    loglike=sum(log(prob))
  }
  return(loglike)
} 

## function to obtain loglikelihood when growth rate, theta, and both parameter of deviation from nominal growth are to be found
loglik4=function(a,theta,mu,sigma){
  res=matrix(nrow=nobs-1)
  prob=matrix(nrow=nobs-1)
  for (te in 2:nobs){
    k=100*forest[te]^theta
    res[te-1]=(pop[te]/pop[te-1])-1-(a*(1-pop[te-1]/k))
    prob[te-1]=dnorm(res[te-1],mean=mu,sd=sigma)
    loglike=sum(log(prob))
  }
  return(loglike)
} 

## function to obtain loglikelihood when growth rate,, and both parameter of deviation from nominal growth are to be found
loglik5=function(a,mu,sigma,nobs){ # with continuous growth rate r as opposed to lambda and a ricker moel
  res=matrix(nrow=nobs)
  prob=matrix(nrow=nobs-1)
  for (te in 2:nobs){
    res[te-1]=log(pop[te]/pop[te-1])-(log(a)*(1-pop[te-1]))
    prob[te-1]=dnorm(res[te-1],mean=mu,sd=sigma)
    loglike=sum(log(prob))
  }
  return(loglike)
} 



# filzbach parameters list
filz.par=list(a=c(0,4,1,0,-1,1),
              k=c(0,200,150,0,-1,1), # change last value to length(k) if variable carrying cap
              No=c(1,30,3,0,-1,1))

filz.par2=list(a=c(0,4,1,0,-1,1),
              theta=c(0,2,1.5,0,-1,1)) #

filz.par3=list(a=c(0,4,1,0,-1,1),
               theta=c(0,2,1,0,-1,1),
               sigma=c(0,0.6,0.2,0,-1,1)) #

filz.par4=list(a=c(0,4,1,0,-1,1),
               theta=c(0,2,1,0,-1,1),
               mu=c(-0.1,0.1,0,0,-1,1),
               sigma=c(0,0.6,0.2,0,-1,1))#, nobs=c(20,20,20,1,1,1)) #

filz.par5=list(a=c(0,20,1,0,-1,1),
               mu=c(-0.1,0.1,0,0,-1,1),
               sigma=c(0,0.6,0.2,0,-1,1)) #
    

# ARGUMENTS:
#10k burn-in iter, 10k sampling iter, residual function used in place of log-likelihood, 20 observations, list of parameters
# parameters=filzbach(10000,10000,residuals2,20,filz.par2)
# 
# parameters=filzbach(20000,20000,residuals3,19,filz.par3) #use log-likelihood and noise

parameters=matrix(nrow=0,ncol=4)
for (i in 1:5){
parameters=rbind(parameters,filzbach(20000,20000,loglik4,19,filz.par4,thinning=10)) #use log-likelihood and noise and bias in growth rate
}
# you could also add error and bias in count with error in Nt drawn from normal with mean=0 (white noise) or !=0 with bias


summary(parameters) # the median value is extremely close to the "real" value




boxplot(parameters) # boxplot of posterior distribution of parameters
points(1:4,c(real_a,theta,mu,sigma),col='red',pch=2,cex=2) # the actual paraemeters


dev.copy(pdf,'Real and estimated parameters values 20data-points.pdf')
dev.off()




# Figures 5-6 -----------------------------------------------------------------
library(filzbach)
parameters2=matrix(nrow=0,ncol=4)


nobs<-15

for (i in 1:5){
  parameters2=rbind(parameters2,filzbach(20000,20000,loglik4,15,filz.par4,thinning=10)) #use log-likelihood and noise and bias in growth rate
}

summary(parameters2)

boxplot(parameters2) # boxplot of posterior distribution of parameters
points(1:4,c(real_a,theta,mu,sigma),col='red',pch=2,cex=2) # the actual paraemeters
dev.copy(pdf,'Real and estimated parameters values 15data-points.pdf')
dev.off()



## run simulations with the posterior distro using 15 points to fit the model ## do the same with parameters and save 20data-points too
populations=apply(parameters2[1:500,],MARGIN=1,FUN=function(x) f_var_for_rand(forest,t,No,a=x[1],theta=x[2],mu=x[3],sigma=x[4]))
slope=matrix(nrow=500,ncol=2)
plot.new()
plot.window(xlim=c(0,21), ylim=c(0,70))
for (i in 1:500){
  model=rlm(populations[,i]~c(1:20))
  if(min(populations[,i])==0){next}
  else {slope[i,]=summary(model)$coefficients[2,1:2]}
  lines(1:20,populations[,i])
}
axis(side=1);axis(side=2)
points(1:20,pop,col='yellow',pch=19)
title(xlab='years',ylab='population size')
dev.off()

