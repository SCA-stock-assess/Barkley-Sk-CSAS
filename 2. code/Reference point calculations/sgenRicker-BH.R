## Sgen Ricker and Beverton -Holt 


sGenSolverRicker <- function (loga, b) {
  # Function to estimate Sgen from loga and b Ricker parameters
  # Assuming Ricker form: R=S * exp( loga - b * S), and explicit soln for SMSY
  sMSY <- (1 - gsl::lambert_W0(exp(1 - loga))) / b
  fn <- function(S){ -sum( dnorm ( log(sMSY) - log(S) - loga + b*S, 0, 1, log = T))}
  fit <- optimize(f = fn, interval = c(0, sMSY))
  return(fit$minimum)
}


# or more simply, an explicit solution for Sgen- Ricker:
sGenDirect <- function(loga, b){
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  a <- exp(loga)
  -1/b*gsl::lambert_W0(-b*sMSY/a)
}


sGenSolverBH <- function (a, b) {
  # Function to estimate Sgen from a and b BH parameters
  # Assuming BH form: R = aS/(1 + (a/b) *S)
  # Assuming approximation for SMSY: b*sqrt(1/a) - b/a
  sMSY <- b*sqrt(1/a) - b/a
  fn <- function(S){ -sum( dnorm ( log(sMSY) - log( a*S/ (1+ (a/b) * S) ), 
                                   0, 1, log = T)) }
  fit <- optimize(f = fn, interval = c(0, sMSY))
  return(fit$minimum)
}

# # Check to make sure BH equation works:
# a <- 5
# b <-1000
# Sgen <- sGenSolverBH (a,b)
# S <- 1:1000
# R <- a*S/(1 + (a/b) *S)
# plot(S,R)
# points(1:1000, 1:1000, type="line")
# sMSY <- b*sqrt(1/a) - b/a
# Sgen <- sGenSolverBH (a,b)
# abline(v=sMSY)
# abline(v=Sgen)

