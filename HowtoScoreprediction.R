library(dplyr)

#Create some variables:
# Simulate data for 100 individuals
n <- 5000

# Generate age between 20 and 60
age <- round(runif(n, min = 20, max = 60))

# Define education levels
education_levels <- c("High School", "Bachelor's", "Master's")

# Simulate education level probabilities
education_probs <- c(0.4, 0.4, 0.2)

# Sample education level based on probabilities
education <- sample(education_levels, n, replace = TRUE, prob = education_probs)

# Simulate experience correlated with age with some random error
experience <- age - 20 + round(rnorm(n, mean = 0, sd = 3)) 

# Define a non-linear function for wage
wage <- exp((age * 0.1) + (case_when(education == "High School" ~ 1,
                                     education == "Bachelor's" ~ 1.5,
                                     TRUE ~ 2)) + (experience * 0.05) + rnorm(n, mean = 0, sd = 0.5))

hist(wage)



ageDave<-30
educationDave<-"Bachelor's"
experienceDave <- 10

wageDave <- exp((ageDave * 0.1) + (case_when(educationDave == "High School" ~ 1,
                                             educationDave == "Bachelor's" ~ 1.5,
                                             TRUE ~ 2)) + (experienceDave * 0.05) + rnorm(n, mean = 0, sd = 0.5))

hist(wageDave, main="Wage Distribution for Dave", xlab="Wage")



## Generate test set
ntest<-1000

# Generate age between 20 and 60
agetest <- round(runif(ntest, min = 20, max = 60))

# Sample education level based on probabilities
educationtest <- sample(education_levels, ntest, replace = TRUE, prob = education_probs)

# Simulate experience correlated with age with some random error
experiencetest <- agetest - 20 + round(rnorm(ntest, mean = 0, sd = 3))

## Generate ytest that we try to predict:

wagetest <- exp((agetest * 0.1) + (case_when(educationtest == "High School" ~ 1,
                                             educationtest == "Bachelor's" ~ 1.5,
                                             TRUE ~ 2)) + (experiencetest * 0.05) + rnorm(ntest, mean = 0, sd = 0.5))





conditionalmeanest <-
  function(age, education, experience, N = 1000) {
    mean(exp((age * 0.1) + (
      case_when(
        education == "High School" ~ 1,
        education == "Bachelor's" ~ 1.5,
        TRUE ~ 2
      )
    ) + (experience * 0.05) + rnorm(N, mean = 0, sd = 0.5)
    ))
  }

conditionalmedianest <-
  function(age, education, experience, N = 1000) {
    median(exp((age * 0.1) + (
      case_when(
        education == "High School" ~ 1,
        education == "Bachelor's" ~ 1.5,
        TRUE ~ 2
      )
    ) + (experience * 0.05) + rnorm(N, mean = 0, sd = 0.5)
    ))
  }


hist(wageDave, main="Wage Distribution for Dave", xlab="Wage")
abline(v=conditionalmeanest(ageDave, educationDave, experienceDave), col="darkred", cex=1.2)
abline(v=conditionalmedianest(ageDave, educationDave, experienceDave), col="darkblue", cex=1.2)



Xtest<-data.frame(age=agetest, education=educationtest, experience=experiencetest)

meanest<-sapply(1:nrow(Xtest), function(j)  conditionalmeanest(Xtest$age[j], Xtest$education[j], Xtest$experience[j])  )
median<-sapply(1:nrow(Xtest), function(j)  conditionalmedianest(Xtest$age[j], Xtest$education[j], Xtest$experience[j])  )



(MSE1<-mean((meanest-wagetest)^2))
(MSE2<-mean((median-wagetest)^2))

MSE1 < MSE2
### Method 1 (the true mean estimator) is better than method 2!

# but the MAE is actually worse of method 1!
(MAE1<-mean(abs(meanest-wagetest)) )
(MAE2<-mean( abs(median-wagetest)))

MAE1 < MAE2
### Method 2 (the true median estimator) is better than method 1!








library(scoringutils)

## Define conditional quantile estimation
conditionalquantileest <-
  function(probs, age, education, experience, N = 1000) {
    quantile(exp((age * 0.1) + (
      case_when(
        education == "High School" ~ 1,
        education == "Bachelor's" ~ 1.5,
        TRUE ~ 2
      )
    ) + (experience * 0.05) + rnorm(N, mean = 0, sd = 0.5)
    )
    , probs =
      probs)
  }

## Define a very naive estimator that will still have the required coverage
lowernaive <- 0
uppernaive <- max(wage)

# Define the quantile of interest
alpha <- 0.05

lower <-
  sapply(1:nrow(Xtest), function(j)
    conditionalquantileest(alpha / 2, Xtest$age[j], Xtest$education[j], Xtest$experience[j]))
upper <-
  sapply(1:nrow(Xtest), function(j)
    conditionalquantileest(1 - alpha / 2, Xtest$age[j], Xtest$education[j], Xtest$experience[j]))

## Calculate the scores for both estimators

# 1. Score the alpha/2 quantile estimate
qs_lower <- mean(quantile_score(wagetest,
                                predictions = lower,
                                quantiles = alpha / 2))
# 2. Score the alpha/2 quantile estimate
qs_upper <- mean(quantile_score(wagetest,
                                predictions = upper,
                                quantiles = 1 - alpha / 2))

# 1. Score the alpha/2 quantile estimate
qs_lowernaive <- mean(quantile_score(wagetest,
                                     predictions = rep(lowernaive, ntest),
                                     quantiles = alpha / 2))
# 2. Score the alpha/2 quantile estimate
qs_uppernaive <- mean(quantile_score(wagetest,
                                     predictions = rep(uppernaive, ntest),
                                     quantiles = 1 - alpha / 2))

# Construct the interval score by taking the average
(interval_score <- (qs_lower + qs_upper) / 2)
# Score of the ideal estimator: 187.8337

# Construct the interval score by taking the average
(interval_scorenaive <- (qs_lowernaive + qs_uppernaive) / 2)
# Score of the naive estimator: 1451.464


library(scoringRules)

## Ideal "estimate": Simply sample from the true conditional distribution 
## P(Y | X=x) for each sample point x
distributionestimate <-
  function(age, education, experience, N = 100) {
    exp((age * 0.1) + (
      case_when(
        education == "High School" ~ 1,
        education == "Bachelor's" ~ 1.5,
        TRUE ~ 2
      )
    ) + (experience * 0.05) + rnorm(N, mean = 0, sd = 0.5))
  }

## Naive Estimate: Only sample from the error distribution, without including the 
## information of each person.
distributionestimatenaive <-
  function(age, education, experience, N = 100) {
    exp(rnorm(N, mean = 0, sd = 0.5))
  }

scoretrue <- mean(sapply(1:nrow(Xtest), function(j)  {
  wageest <-
    distributionestimate(Xtest$age[j], Xtest$education[j], Xtest$experience[j])
  return(scoringRules::es_sample(y = wagetest[j], dat = matrix(wageest, nrow=1)))
}))

scorenaive <- mean(sapply(1:nrow(Xtest), function(j)  {
  wageest <-
    distributionestimatenaive(Xtest$age[j], Xtest$education[j], Xtest$experience[j])
  return(scoringRules::es_sample(y = wagetest[j], dat = matrix(wageest, nrow=1)))
}))

## scoretrue: 761.026
## scorenaive: 2624.713
