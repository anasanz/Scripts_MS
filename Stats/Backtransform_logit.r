

## -------------------------------------------------
##   Backtransformation logistic regressiom
##    Understand probability of detection 
## ------------------------------------------------- 

# Examples: Effect of the applicant having honours on 
# job callback rates (probability of being called)

x <- seq(0,1,0.005)
plot(qlogis(x) ~ x, type = "l", xlab = "p", ylab = "logit(p)")


alpha = -2.5
beta.honors = 0.87
honors = rbinom(50,1,0.5)

logitp <- alpha + beta.honors*honors

# lp is not between 0 and 1, it is in the logit scale (qlogis(x))
# To know the probability of an applicant being called with
# and without honours, need to backtransform

# The inverse of the logit is the logistic function which is plogis(x)

p <- plogis(logitp) # Probability between 0 and 1
                      # of getting a callback with and without honors
qlogis(p) # With this we get again the values -inf,+inf


# invlogit (plogis function in r) = exp(x)/(1+exp(x))

## ---- Priors for p ----
beta.norm <- rnorm(1000,0,5)
beta.unif <- runif(1000,0,1)

hist(beta.norm)
hist(beta.unif)

# If pnorm is the values that the beta coefficient can take
# it is important to choose an appropiate prior (even if its not informative)

hist(plogis(beta.norm), breaks = 100)
hist(plogis(beta.unif), breaks = 100)

# So for PROBABILITIES, always choose UNIFORM for informative prior!!!
