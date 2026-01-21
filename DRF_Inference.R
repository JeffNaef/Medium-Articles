library(kernlab)
library(drf)
library(Matrix)
library(Hmisc)
library(MASS)
library(ggplot2)


set.seed(2)

n<-2000
beta1<-1
beta2<--1.8

# Model Simulation
X<-mvrnorm(n = n, mu=c(0,0), Sigma=matrix(c(1,0.7,0.7,1), nrow=2,ncol=2))
u<-rnorm(n=n, sd = sqrt(exp(X[,1])))
Y<- matrix(beta1*X[,1] + beta2*X[,2] + u, ncol=1)


# Choose an x that is not too far out
x<-matrix(c(1,1),ncol=2)

# Choose alpha for CIs
alpha<-0.05



## Fit the new DRF framework
drf_fit <- drf(X=X, Y=Y, min.node.size = 5, num.trees=2000, num.features=10, ci.group.size=2000/50)

## predict weights
DRF = predict(drf_fit, newdata=x, estimate.uncertainty = TRUE)
weights <- DRF$weights[1,]
weightsb<-DRF$weights.uncertainty[[1]]


################################
#Example 1: Conditional Expectation
################################

# Estimate the conditional expectation at x:
condexpest<- sum(weights*Y)

# Use the distribution of weights, see below
distofcondexpest<-apply(weightsb,1, function(wb)  sum(wb*Y)  )

# Can either use the above directly to build confidence interval, or can use the normal approximation.
# We will use the latter
varest<-var(distofcondexpest-condexpest)

# build 95%-CI
lower<-condexpest - qnorm(1-alpha/2)*sqrt(varest)
upper<-condexpest + qnorm(1-alpha/2)*sqrt(varest)
c(round(lower,2), round(condexpest,2), round(upper,2))

#(-1.16, -0.59, -0.03)


################################
#Example 2: Conditional Variance
################################

# Estimate the conditional expectation at x:
condvarest<- sum(weights*Y^2) - condexpest^2

distofcondvarest<-apply(weightsb,1,  function(wb)  {
  sum(wb*Y^2) - sum(wb*Y)^2
}  )
  

# Can either use the above directly to build confidence interval, or can use the normal approximation.
# We will use the latter
varest<-var(distofcondvarest-condvarest)

# build 95%-CI
lower<-condvarest - qnorm(1-alpha/2)*sqrt(varest)
upper<-condvarest + qnorm(1-alpha/2)*sqrt(varest)


c(round(lower,2), round(condvarest,2), round(upper,2))

#(1.29, 2.24, 3.18)


################################
# Causal Analysis Example: Witness Function
################################


Witdrf<- function(DRF, x, groupingvar, alpha=0.05, ...){
  
  ### Function to calculate the conditional witness function with
  ### confidence bands from DRF
  ### DRF: DRF object
  ### x: Testpoint
  
  if (is.null(dim(x)) ){
    
    stop("x needs to have dim(x) > 0")
  }
  
  ntest <- nrow(x)
  n <- nrow(DRF$Y)
  coln<-colnames(DRF$Y.orig)
  
  drfpred<- predict(DRF, newdata=x, estimate.uncertainty = TRUE)
  weightsall<-drfpred$weights[1,]
  weightsb<-drfpred$weights.uncertainty[[1]]
  
  #weightsall0<-weightsall[, DRF$Y[, groupingvar]==0, drop=F]
  #weightsall1<-weightsall[,DRF$Y[, groupingvar]==1, drop=F]
  
  
  # Get the weights of the respective classes (need to standardize by propensity!)
  weightsall0<-weightsall*(DRF$Y.orig[, groupingvar]==0)/sum(weightsall*(DRF$Y.orig[, groupingvar]==0))
  weightsall1<-weightsall*(DRF$Y.orig[, groupingvar]==1)/sum(weightsall*(DRF$Y.orig[, groupingvar]==1))
  
  
  bandwidth_Y <- drf:::medianHeuristic(DRF$Y.orig)
  k_Y <- rbfdot(sigma = bandwidth_Y)
  
  K<-kernelMatrix(k_Y, DRF$Y.orig[,coln[coln!=groupingvar]], y = DRF$Y.orig[,coln[coln!=groupingvar]])
  
  
  nulldist <- apply(weightsb,1, function(wb){
    # iterate over class 1
    
    wb0<-wb*(DRF$Y[, groupingvar]==0)/sum(wb*(DRF$Y[, groupingvar]==0))
    wb1<-wb*(DRF$Y[, groupingvar]==1)/sum(wb*(DRF$Y[, groupingvar]==1))
    
    
    diag( ( wb0-weightsall0 - (wb1-weightsall1) )%*%K%*%( wb0-weightsall0 - (wb1-weightsall1) )  )
    
    
  })
  
  # Choose the right quantile
  c<-quantile(nulldist, 1-alpha)
  
  
  return(list(c=c, k_Y=k_Y, Y=DRF$Y[,coln[coln!=groupingvar]], nulldist=nulldist, weightsall0=weightsall0, weightsall1=weightsall1))
  
  
  
}




set.seed(2)

n<-5000
p<-2

X<-matrix(runif(n*p), ncol=p)
W<-rbinom(n,size=1, prob= exp(-X[,2])/(1+exp(-X[,2])))

Y<-(W-0.2)*X[,1] + rnorm(n)
Y<-matrix(Y,ncol=1)



x<-matrix(runif(1*p), ncol=2)
Yall<-cbind(Y,W)
## For the current version of the Witdrf function, we need to give
## colnames to Yall
colnames(Yall) <- c("Y", "W")

## Fit the new DRF framework
drf_fit <- drf(X=X, Y=Yall, min.node.size = 5, num.trees=4000, ci.group.size=4000/40)

Witobj<-Witdrf(drf_fit, x=x, groupingvar="W", alpha=0.05)

hatmun<-function(y,Witobj){
  
  c<-Witobj$c
  k_Y<-Witobj$k_Y
  Y<-Witobj$Y
  weightsall1<-Witobj$weightsall1
  weightsall0<-Witobj$weightsall0
  Ky=t(kernelMatrix(k_Y, Y , y = y))
  
  out<-list()
  out$val <- (Ky%*%weightsall1 - Ky%*%weightsall0)
  out$upper<-  out$val+sqrt(c)
  out$lower<-  out$val-sqrt(c)
  
  return( out )
  
  
  
}

all<-hatmun(sort(Witobj$Y),Witobj)

plot(sort(Witobj$Y),all$val , type="l", col="darkblue", lwd=2, ylim=c(min(all$lower), max(all$upper)),
     xlab="y", ylab="witness function")
lines(sort(Witobj$Y),all$upper , type="l", col="darkgreen", lwd=2 )
lines(sort(Witobj$Y),all$lower , type="l", col="darkgreen", lwd=2 )
abline(h=0)



# Simulate truth for a large number of samples ntest
ntest<-10000
Xtest<-matrix(runif(ntest*p), ncol=2)

Y1<-(1-0.2)*Xtest[,1] + rnorm(ntest)
Y0<-(0-0.2)*Xtest[,1] + rnorm(ntest)

## Plot the test data without adjustment
plotdf = data.frame(Y=c(Y1,Y0), W=c(rep(1,ntest),rep(0,ntest) ))
plotdf$weight=1
plotdf$plotweight[plotdf$W==0] = plotdf$weight[plotdf$W==0]/sum(plotdf$weight[plotdf$W==0])
plotdf$plotweight[plotdf$W==1] = plotdf$weight[plotdf$W==1]/sum(plotdf$weight[plotdf$W==1])

plotdf$W <- factor(plotdf$W)

#plot pooled data
ggplot(plotdf, aes(Y)) +
  geom_density(adjust=2.5, alpha = 0.3, show.legend=TRUE,  aes(fill=W, weight=plotweight)) +
  theme_light()+
  scale_fill_discrete(name = "Group", labels = c('0', "1"))+
  theme(legend.position = c(0.83, 0.66),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.background = element_rect(fill=alpha('white', 0.5)),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=19),
        axis.title.y = element_text(size=19))+
  labs(x='y')

