#' Given target data, fit the AFS function
#'
#' This function takes AFS target data and estimates parameters for the AFS function
#' A dataframe specifying the rare MAC bins and the observed proportion of variants is used to fit the data
#' The proportion of rare variants (p_rv) is by default the sum of the rare allele count bins. The proportion be manually specified if desired
#' 
#' @param  Observed_bin_props data frame with 3 columns, Lower, Upper (of MAC bins)  and proportion of variants in that MAC bin
#'
#' @param p_rv proportion of rare variants - default is the sum of the rare MAC bin proportions
#'
#' @return list of parameters - alpha, beta, and b as well as fitted proportions
#'
#' @examples 
#' data("afs_afr")
#' colnames(afs_afr)[3] <- 'Prop'
#' fit_afs(Observed_bin_props = afs_afr)
#'
#' @export
#' 
#' @importFrom nloptr slsqp
#'

fit_afs <- function(Observed_bin_props, p_rv = NULL){
  
  # Check the column names of Observed_bin_props
  if((colnames(Observed_bin_props)[1] == 'Lower') == FALSE | (colnames(Observed_bin_props)[2]  == 'Upper') == FALSE |
     (colnames(Observed_bin_props)[3] == 'Prop') == FALSE ){
    stop('Observed_bin_props needs to have column names Lower, Upper, and Prop')
  }
  
  # Make sure the Observed_bin_props are numeric
  Observed_bin_props$Prop <- (as.numeric(as.character(Observed_bin_props$Prop)))
  
  # Make sure there are not any NA's in the proportions
  if(anyNA(Observed_bin_props$Prop) == TRUE){
    stop('Proportions in Observed_bin_props need to be numeric with no NA values')
  }
  
  if((is.numeric(Observed_bin_props$Lower) ==  FALSE) | (is.numeric(Observed_bin_props$Upper) == FALSE)){
    stop('Observed_bin_props MAC bins need to be numberic')
  }
  
  # Check the order of the MAC bins
  if(is.unsorted(Observed_bin_props$Upper)==TRUE){
    stop('The MAC bins need to be ordered from smallest to largest')
  }
  
  # Set the default value for p_rv to the sum of the rare variant bins
  if(is.null(p_rv)){
    p_rv <- sum(Observed_bin_props$Prop)
  }
  
  # specify the function to define the constraints
  hin.tune <- function(x) {
    h <- numeric(1)
    h[1] <- x[1] # the first parameter (alpha) must be > 0
    return(h)
  }

  c1 <- 1:Observed_bin_props$Upper[nrow(Observed_bin_props)] #### individual MACs to use in the function
  
  calc_prob_LS<-function(tune){ # define the least squares loss function
    
    ### calculate b
    indivual_prop_no_b <- 1/((c1+tune[2])^(tune[1])) # calculate the function completely without b for each individual MAC
    b <- p_rv/sum(indivual_prop_no_b) #### solve for b
    
    #### calculate the function with b for each individual MAC
    indivual_prop <- b*indivual_prop_no_b
    
    all <- 0
    for(i in 1:nrow(Observed_bin_props)){ # loop over the bins
      E <- sum(indivual_prop[Observed_bin_props$Lower[i]:Observed_bin_props$Upper[i]]) # Calculate expected (from the function)
      O <- Observed_bin_props$Prop[i] # record the observed proportion in the target data
      c <- (E-O)^2 # calculate the squared error
      all <- all+c # sum the squared error over all MAC bins
    }
    
    return(all) # The output here is the sum of the squared error over MAC bins
  }  
  
  tune <- c(1,0) #### start with the function 1/x (alpha = 1, beta = 0)
  
  # Minimize with the SLSQP function using the starting values (tune), the least squares loss function (calc_prob_LS), and constraints (hin.tune)
  S <- suppressMessages(slsqp(tune, fn = calc_prob_LS, hin = hin.tune )) # suppressMessages b/c each produces a warning
  b <- p_rv/sum(1/((c1+S$par[2])^(S$par[1]))) #back calculate b after the parameters have been solved for
  
  
  # and calculate the MAC bin proportions given the parameters:
  
  re <- afs(alpha = S$par[1], beta = S$par[2], b = b, mac_bins = Observed_bin_props[,c(1,2)])
  
  # Return the parameters alpha, beta, and b, as well as the proportions as calculated by the function
  return(list(alpha = S$par[1], beta = S$par[2], b = b, Fitted_results = re))
  
}
