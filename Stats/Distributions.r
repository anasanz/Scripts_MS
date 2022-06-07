
## -------------------------------------------------
##                 Distributions
## ------------------------------------------------- 

## ---- Binomial ----

# Study with N = 10 nest of birds,
# Each nest has a p.detection of 0.5
n <- rbinom(3, size=10, prob=0.5) # Generate 3 binomial outcomes
# Each time we go (J = 3) we observe a different number

# Probability of getting n = 5 when the number of trials is 10 and 
# the probability of success is 0.5 according to the binomial pmf
dbinom(5, 10, 0.5)

# You can plot the whole pmf
plot(0:10, dbinom(0:10, 10, 0.5), type="h", ylab="Probability", 
     xlab="Number of shad caught (X) after 10 casts")

#♣ De 10, hay una probabilidad distinta de cada numero. Esta probabilidad
# depende de la probabilidad de sucess to get the number of trials
dbinom(8, 10, 0.8)
# The sum must be 1
sum(dbinom(0:10, size=10, p=0.5))

# Gives us important summaries: mean and variance
# Mean = Expected value of a binomial distribution of 10 trials and p = 0.5
mean(rbinom(10000, 10, 0.5)) 
# We can get expected values exactly
10*0.5
# Expected value: the average of all possible values of the random variable, weighted by their probabilities.
# When we fit models, we are often modeling changes in the expected value of some random variable

#Variance given by k (trials) and p and the formula
10* 0.5*(1-0.5)

# Example with variable x
x <- rbinom(100000, 10, 0.5)
hist(x)
mean(x)
mean((x-mean (x))^2)


## ---- Bernouilli ----

# Binomial distribution when N = 1
# USeful when individuals have different p
n1 <- rbinom(1, size=1, prob=0.5)
n2 <- rbinom(1, size=1, prob=0.3)
# N is the sum of all bernouilli outcomes n1+n2

# Implementing the pmf: Bernoulli probability of observing n = 1 given p = 0.3
dbinom(1, 1, 0.3)

## ---- Poisson ----

rpois(5, 100)
dpois(107,100) #Probability of getting a realized value of 107 when lambda is 100?

## ---- Uniform ----
# Priors for p and priors for coordinates (independent)

D <- 100 # points per unit area 
A <- 1 # Area of unit square 
N <- rpois(1, D*A)
plot(s <- cbind(runif(N, min = 0, max = 1), runif(N, min = 0, max = 1)))
# s (the coordinates) are distributed following a random distribution

## ---- Normal ----
# For prior of beta coefficcients of covariates

## ---- Beta ----
# Prior for probabilities: express either a lack of knowledge or very precise knowledge about a parameter.

rbeta(1,1,1) # Is the same as runif:
runif(1,0,1)

# It can be use to create informative priors:
curve(dbeta(x, 1, 1), col="black", ylim=c(0,5))
curve(dbeta(x, 10, 10), col="blue", add=TRUE) 
curve(dbeta(x, 10, 20), col="darkgreen", add=TRUE)

curve(dbeta(x, 1, 2), col="black", ylim=c(0,5))
