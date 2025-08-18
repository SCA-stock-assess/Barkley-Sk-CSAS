# Libraries used to build functions ---------------------------------------

library(gsl)

# Reference point calculations --------------------------------------------

# Both taken from Dylan Glaser's Fraser Pink analysis:
# https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/analysis/R/functions.R
get_Smsy <- function(a, b){
  Smsy <- (1 - gsl::lambert_W0(exp(1 - a))) / b
  if(Smsy <0){Smsy <- 0.001} #hack for low draws so Smsy doesnt go negative
  return(Smsy)
}

get_Sgen <- function(a, b, int_lower, int_upper, Smsy){
  fun_Sgen <- function(Sgen, a, b, Smsy) {Sgen * a * exp(-b * Sgen) - Smsy}
  Sgen <- uniroot(fun_Sgen, interval=c(int_lower, int_upper), a=a, b=b, Smsy=Smsy)$root
  return(Sgen)
}

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



# Management table rules --------------------------------------------------


# Somass, based on management table brackets
somass_mgt_rule <- function(somass_run) {
  
  if (somass_run < 200000)  {somass_hr <- 0}
  else if (somass_run < 250000)  {somass_hr <- 0.15}
  else if (somass_run < 300000)  {somass_hr <- 0.2}
  else if (somass_run < 350000)  {somass_hr <- 0.229166667}
  else if (somass_run < 400000)  {somass_hr <- 0.25}
  else if (somass_run < 450000)  {somass_hr <- 0.291666667}
  else if (somass_run < 500000)  {somass_hr <- 0.324074074}
  else if (somass_run < 550000)  {somass_hr <- 0.35}
  else if (somass_run < 600000)  {somass_hr <- 0.397727273}
  else if (somass_run < 650000)  {somass_hr <- 0.4375}
  else if (somass_run < 700000)  {somass_hr <- 0.471153846}
  else if (somass_run < 750000)  {somass_hr <- 0.5}
  else if (somass_run < 800000)  {somass_hr <- 0.522222222}
  else if (somass_run < 850000)  {somass_hr <- 0.541666667}
  else if (somass_run < 900000)  {somass_hr <- 0.558823529}
  else if (somass_run < 950000)  {somass_hr <- 0.574074074}
  else if (somass_run < 1000000)  {somass_hr <- 0.587719298}
  else if (somass_run < 1050000)  {somass_hr <- 0.6}
  else if (somass_run < 1100000)  {somass_hr <- 0.618253968}
  else if (somass_run < 1150000)  {somass_hr <- 0.634848485}
  else if (somass_run < 1200000)  {somass_hr <- 0.65}
  else if (somass_run < 1250000)  {somass_hr <- 0.658928571}
  else if (somass_run < 1300000)  {somass_hr <- 0.667142857}
  else if (somass_run < 1350000)  {somass_hr <- 0.674725275}
  else if (somass_run < 1400000)  {somass_hr <- 0.681746032}
  else if (somass_run < 1450000)  {somass_hr <- 0.688265306}
  else if (somass_run < 1500000)  {somass_hr <- 0.694334975}
  else {somass_hr <- 0.7}
 
  return(somass_hr) 
 
}


# Hucuktlis, based on management table brackets
hucuktlis_mgt_rule <- function(hucuktlis_run) {
  
  if (hucuktlis_run < 20000)  {hucuktlis_hr <- 0.15}
  else if (hucuktlis_run < 25000)  {hucuktlis_hr <- 0.18125}
  else if (hucuktlis_run < 33000)  {hucuktlis_hr <- 0.20}
  else if (hucuktlis_run < 35000)  {hucuktlis_hr <- 0.25454545}
  else if (hucuktlis_run < 45000)  {hucuktlis_hr <- 0.26428571}
  else if (hucuktlis_run < 50000)  {hucuktlis_hr <- 0.3}
  else if (hucuktlis_run < 52500)  {hucuktlis_hr <- 0.34}
  else if (hucuktlis_run < 60000)  {hucuktlis_hr <- 0.35714286}
  else if (hucuktlis_run < 67500)  {hucuktlis_hr <- 0.4}
  else if (hucuktlis_run < 75000)  {hucuktlis_hr <- 0.45555556}
  else {hucuktlis_hr <- 0.5}
  
  return(hucuktlis_hr) 
  
}

