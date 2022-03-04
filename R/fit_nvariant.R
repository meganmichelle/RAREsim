#' Given target data, fit the Number of Variants function
#'
#' This function takes Number of Variants target data and estimates parameters for the Number of Variants function.
#' A dataframe specifying the number of variants per Kb at various sample sizes is required to fit the data
#'
#' @param Observed_variants_per_kb A data frame with the first column sample size and the second variants per Kb, both numeric
#'
#' @return Vector of parameters - phi and omega
#'
#' @export
#'
#' @examples
#' data("nvariant_afr")
#' fit_nvariant(nvariant_afr)
#'
#' @importFrom nloptr slsqp
#'

fit_nvariant <- function(Observed_variants_per_kb){
    
    # Check that each column is numeric
    if(is.numeric(Observed_variants_per_kb[,1]) == FALSE | is.numeric(Observed_variants_per_kb[,2]) == FALSE){
        stop('Error: Columns needs to  be numeric')
    }
    
    # Check that there are not any NA values
    if(anyNA(Observed_variants_per_kb[,2]) == TRUE){
        stop('Number of variants per Kb need to be numeric with no NA values')
    }
    
    # Check that the sample sizes go from smallest to largest
    if(is.unsorted(Observed_variants_per_kb[,1])==TRUE){
        stop('The sample sizes need to be ordered from smallest to largest')
    }
    
    # define the least squares loss function
    leastsquares <- function(tune){ # function of the parameters (phi and omega)
        E <- tune[1]*(Observed_variants_per_kb[,1]^(tune[2])) # calculated the expected number of variants (from the function)
        sq_err <- (E - Observed_variants_per_kb[,2])^2 ## calculate the squared errror of expected - observed
        d <- sum(sq_err) # sum over the sample sizes
        return(d) # return the squared error
    }
    
    hin.tune <- function(x) { ### constraints
        h <- numeric(3)
        h[1] <- x[1] ### phi greater than 0
        h[2] <- x[2] ### omega greater than 0
        h[3] <- 1 - x[2] ### omega less than 1
        return(h)
    }
    
    ### define the starting value for phi so the end of the function matches with omega = 0.45
    phi <- Observed_variants_per_kb[which.max(Observed_variants_per_kb[,1]),2]/(Observed_variants_per_kb[which.max(Observed_variants_per_kb[,1]),1])^0.45
    tune <- c(phi, 0.45) # specify the starting values
    
    ### Use SLSQP to find phi and omega
    re_LS <- suppressMessages(slsqp(tune, fn = leastsquares, hin = hin.tune,
                                    control = list(xtol_rel = 1e-12))) # suppressMessages because SLSQP always give a warning
    
    # If the original starting value resulted in a large loss (>1000), iterate over starting values
    if(re_LS$value > 1000){
        re_tab1<-c() # create to hold the new parameters
        for(omega in c(seq(0.15,0.65, by=0.1))){ ### optimize with different values of omega
            
            # specify phi to  fit the end of the function with the current value of omega
            phi <- Observed_variants_per_kb[which.max(Observed_variants_per_kb[,1]),2]/(Observed_variants_per_kb[which.max(Observed_variants_per_kb[,1]),1])^omega
            tune <- c(phi, omega) # updated  starting values
            
            re_LS1 <- suppressMessages(slsqp(tune, fn = leastsquares, hin = hin.tune,
                                             control = list(xtol_rel = 1e-12))) # estimate parameters with SLSQP
            to_bind1 <- c(re_LS1$par, re_LS1$value ) # record parameters and loss value
            re_tab1 <- rbind(re_tab1, to_bind1) # bind information from each iteration together
        }
        
        re_fin <-re_tab1[which.min(re_tab1[,3]),] # select the minimum least squared error
        re_fin  <- as.vector(re_fin)
        return(list( phi=re_fin[1], omega=re_fin[2])) # return the phi and omega that resulted in the smallest loss
    }
    else{ # if the loss was <1000, bring the parameters forward
        return(list( phi=re_LS$par[1], omega=re_LS$par[2])) # return phi and omega
    }
}
