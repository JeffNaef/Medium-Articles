## Load packages and functions needed
library(drf)
library(mice)

## Set parameters
set.seed(10)
n<-1000

##Simulate Data that experiences both a mean as well as sd shift
# Simulate from X
x1 <- runif(n,-1,1)
x2 <- runif(n,-1,1)
x3 <- x1+ runif(n,-1,1)
X0 <- matrix(runif(7*n,-1,1), nrow=n, ncol=7)
Xfull <- cbind(x1,x2, x3, X0)
colnames(Xfull)<-paste0("X", 1:10)

# Simulate dependent variable Y
Y <- as.matrix(rnorm(n,mean = 0.8*(x1 > 0), sd = 1 + 1*(x2 > 0)))

colnames(Y)<-"Y"

##Also add MAR missing values using ampute from the mice package
X<-ampute(Xfull)$amp

head(cbind(Y,X))


# Obtain Test point
x<-matrix(c(0.2, 0.4, runif(8,-1,1)), nrow=1, ncol=10)

print(x)



# Fit DRF with 50 CI groups, each 20 trees large. This results in 50 uncertainty weights 
DRF <- drf(X=X, Y=Y,num.trees=2000, min.node.size = 5, ci.group.size=2000/50)
DRFpred<-predict(DRF, newdata=x, estimate.uncertainty=TRUE)


## Sample from P_{Y| X=x}
Yxs<-Y[sample(1:n, size=n, replace = T, DRFpred$weights[1,])]
hist(Yxs, prob=T)
z<-seq(-6,7,by=0.01)
d<-dnorm(z, mean=0.8 * (x[1] > 0), sd=(1+(x[2] > 0)))
lines(z,d, col="darkred"  )



# Calculate quantile prediction as weighted quantiles from Y
qx <- quantile(Yxs, probs = c(0.05,0.95))

# Calculate conditional mean prediction
mux <- mean(Yxs)

# True quantiles
q1<-qnorm(0.05, mean=0.8 * (x[1] > 0), sd=(1+(x[2] > 0)))
q2<-qnorm(0.95, mean=0.8 * (x[1] > 0), sd=(1+(x[2] > 0)))
mu<-0.8 * (x[1] > 0)

hist(Yxs, prob=T)
z<-seq(-6,7,by=0.01)
d<-dnorm(z, mean=0.8 * (x[1] > 0), sd=(1+(x[2] > 0)))
lines(z,d, col="darkred"  )
abline(v=q1,col="darkred" )
abline(v=q2, col="darkred" )
abline(v=qx[1], col="darkblue")
abline(v=qx[2], col="darkblue")
abline(v=mu, col="darkred")
abline(v=mux, col="darkblue")




# Calculate uncertainty
alpha<-0.05
B<-nrow(DRFpred$weights.uncertainty[[1]])
qxb<-matrix(NaN, nrow=B, ncol=2)
muxb<-matrix(NaN, nrow=B, ncol=1)
for (b in 1:B){
  Yxsb<-Y[sample(1:n, size=n, replace = T, DRFpred$weights.uncertainty[[1]][b,])]
  qxb[b,] <- quantile(Yxsb, probs = c(0.05,0.95))
  muxb[b] <- mean(Yxsb)
}

CI.lower.q1 <- qx[1] - qnorm(1-alpha/2)*sqrt(var(qxb[,1]))
CI.upper.q1 <- qx[1] + qnorm(1-alpha/2)*sqrt(var(qxb[,1]))

CI.lower.q2 <- qx[2] - qnorm(1-alpha/2)*sqrt(var(qxb[,2]))
CI.upper.q2 <- qx[2] + qnorm(1-alpha/2)*sqrt(var(qxb[,2]))

CI.lower.mu <- mux - qnorm(1-alpha/2)*sqrt(var(muxb))
CI.upper.mu <- mux + qnorm(1-alpha/2)*sqrt(var(muxb))

hist(Yxs, prob=T)
z<-seq(-6,7,by=0.01)
d<-dnorm(z, mean=0.8 * (x[1] > 0), sd=(1+(x[2] > 0)))
lines(z,d, col="darkred"  )
abline(v=q1,col="darkred" )
abline(v=q2, col="darkred" )
abline(v=qx[1], col="darkblue")
abline(v=qx[2], col="darkblue")
abline(v=mu, col="darkred")
abline(v=mux, col="darkblue")
abline(v=CI.lower.q1, col="darkblue", lty=2)
abline(v=CI.upper.q1, col="darkblue", lty=2)
abline(v=CI.lower.q2, col="darkblue", lty=2)
abline(v=CI.upper.q2, col="darkblue", lty=2)
abline(v=CI.lower.mu, col="darkblue", lty=2)
abline(v=CI.upper.mu, col="darkblue", lty=2)






## Variable importance for conditional Quantile Estimation


#' Variable importance for Distributional Random Forests
#'
#' @param X Matrix with input training data.
#' @param Y Matrix with output training data.
#' @param X_test Matrix with input testing data. If NULL, out-of-bag estimates are used.
#' @param num.trees Number of trees to fit DRF. Default value is 500 trees.
#' @param silent If FALSE, print variable iteration number, otherwise nothing is print. Default is FALSE.
#'
#' @return The list of importance values for all input variables.
#' @export
#'
#' @examples
compute_drf_vimp <- function(X, Y, X_test = NULL, num.trees = 500, silent = FALSE){
  
  # fit initial DRF
  bandwidth_Y <- drf:::medianHeuristic(Y)
  k_Y <- rbfdot(sigma = bandwidth_Y)
  K <- kernelMatrix(k_Y, Y, Y)
  DRF <- drf(X, Y, num.trees = num.trees)
  wall <- predict(DRF, X_test)$weights
  
  # compute normalization constant
  wbar <- colMeans(wall)
  wall_wbar <- sweep(wall, 2, wbar, "-")
  I0 <- as.numeric(sum(diag(wall_wbar %*% K %*% t(wall_wbar))))
  
  # compute drf importance dropping variables one by one
  I <- sapply(1:ncol(X), function(j) {
    if (!silent){print(paste0('Running importance for variable X', j, '...'))}
    DRFj <- drf(X = X[, -j, drop=F], Y = Y, num.trees = num.trees) 
    DRFpredj <- predict(DRFj, X_test[, -j])
    wj <- DRFpredj$weights
    Ij <- sum(diag((wj - wall) %*% K %*% t(wj - wall)))/I0
    return(Ij)
  })
  
  # compute retraining bias
  DRF0 <- drf(X = X, Y = Y, num.trees = num.trees)
  DRFpred0 = predict(DRF0, X_test)
  w0 <- DRFpred0$weights
  vimp0 <- sum(diag((w0 - wall) %*% K %*% t(w0 - wall)))/I0
  
  # compute final importance (remove bias & truncate negative values)
  vimp <- sapply(I - vimp0, function(x){max(0,x)})
  
  names(vimp)<-colnames(X)
  
  return(vimp)
  
}

## For the conditional quantiles we use a measure that considers the whole distribution,
## i.e. the MMD based measure of DRF.
MMDVimp <- compute_drf_vimp(X=X,Y=Y)
sort(MMDVimp, decreasing = T)

X2          X1          X8          X6          X3         X10 
0.852954299 0.124110913 0.012194176 0.009578300 0.008191663 0.007517931 
X9          X7          X5          X4 
0.006861688 0.006632175 0.005257195 0.002401974



