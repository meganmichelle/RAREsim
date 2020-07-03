#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param temp A data frame with the first column sample size and the second variants per Kb
#'
#' @return Vector of parameters - phi and omega
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#'
#' @export
#'

#### need to go in and rename the different pieces of the function

Fit_var_per_kb <- function(temp){
  leastsquares <- function(tune){
    a <- tune[1]*(temp[,1]^(tune[2]))
    b <- a - temp[,2]
    c <- b^2
    d <- sum(c)
    return(d)
  }
  
  hin.tune <- function(x) { ### constraints
    h <- numeric(3)
    h[1] <- x[1] ### phi greater than 0
    h[2] <- x[2] ### omega greater than 0
    h[3] <- 1 - x[2] ### omega less than 1
    return(h)
  } 
  ### optimize synonymous
  phi <- temp[which.max(temp[,1]),2]/(temp[which.max(temp[,1]),1])^0.45
  tune <- c(phi, 0.45)
  re_LS <- suppressMessages(slsqp(tune, fn = leastsquares, hin = hin.tune,
                                  control = list(xtol_rel = 1e-12)))
  if(re_LS$value > 1000){
    re_tab1<-c()
    for(omega in c(seq(0.15,0.65, by=0.1))){
      ### optimize synonymous
      phi <- temp[which.max(temp[,1]),2]/(temp[which.max(temp[,1]),1])^omega
      tune <- c(phi, omega)
      re_LS1 <- suppressMessages(slsqp(tune, fn = leastsquares, hin = hin.tune,
                                      control = list(xtol_rel = 1e-12)))
      to_bind1 <- c(re_LS1$par, re_LS1$value )
      re_tab1 <- rbind(re_tab1, to_bind1)
    }
    re_fin <-re_tab1[which.min(re_tab1[,3]),]
    re_fin  <- as.vector(re_fin)
    return(list( phi=re_fin[1], omega=re_fin[2]))
  }else{
    return(list( phi=re_LS$par[1], omega=re_LS$par[2]))

  }
}
