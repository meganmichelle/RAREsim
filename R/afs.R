#' afs
#'
#' The afs function calculates the proportion of variants in minor allele count (MAC) bins given parameters, as described in RAREsim.
#' The default parameters will be used if an ancestrial population is specified: (AFR, EAS, NFE, or SAS) is specified.
#' Otherwise, the parameters alpha, beta, and b need to be provided.
#' Alpha, Beta, and b can be estimated from target data using the *Fit_afs* function.
#' 
#' @details The default parameters will be used if an ancestrial population is specified: (AFR, EAS, NFE, or SAS) is specified.
#' Otherwise, the parameters alpha, beta, and b need to be provided.
#' Alpha, Beta, and b can be estimated from target data using the *Fit_afs* function.
#'
#' @param alpha AFS function parameter alpha
#'
#' @param beta AFS function parameter beta
#'
#' @param b AFS function parameter b
#' 
#' @param mac_bins The MAC bins to use, with lower and upper boundaries defined
#' 
#' @param pop The population - specified when using default parameters
#'
#' @return data frame with the MAC bins provided and proportion of variants in each bin
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#'
#' @export
#' 

afs <- function(alpha=NULL, beta=NULL, b=NULL, mac_bins, pop=NULL){
  
  # Check that the default population is correctly specified
  if(is.null(pop)==FALSE){
    if((pop == 'AFR' | pop == 'EAS' | pop == 'NFE' | pop == 'SAS') ==FALSE){
    stop('Default ancestries must be specified as AFR, EAS, NFE, or SAS.')
  }}
  
  # check that the parameters are numeric
  if((is.null(pop)==TRUE)  & ((is.numeric(alpha) == FALSE) | (is.numeric(beta) == FALSE) |
     (is.numeric(b) == FALSE))){
    stop('Error: at least one parameter is not numeric')
  }

  # Check to the correct column names for the MAC bins
  if((colnames(mac_bins)[1] == 'Lower') == FALSE | (colnames(mac_bins)[2]  == 'Upper') == FALSE){
    stop('mac_bins files needs to have column names Lower and Upper')
  }
  
  # Check that the MAC bins are defined with numeric values
  if((is.numeric(mac_bins$Lower) ==  FALSE) | (is.numeric(mac_bins$Upper) == FALSE)){
    stop('mac_bins columns need to be numberic')
  }
  
  # Check the order of the MAC bins
  if(is.unsorted(mac_bins$Upper)==TRUE){
    stop('The MAC bins need to be ordered from smallest to largest')
  }
  
  ### check that if parameters are null, a population specified
  if((is.null(alpha) | is.null(beta) | is.null(b)) & is.null(pop)){
    stop('Error: either a default population should be specified or
         all three parameters provided')
  }

  # specify the default parameters
  if(is.null(alpha)){
    if(pop == 'AFR'){
      alpha = 1.5882
      beta = -0.3083
      b = 0.2872
    }
    if(pop == 'EAS'){
      alpha = 1.6656
      beta = -0.2951
      b =  0.3137
    }
    if(pop == 'NFE'){
      alpha  = 1.9470
      beta  = 0.1180
      b = 0.6676
    }
    if(pop == 'SAS'){
      alpha = 1.6977
      beta = -0.2273
      b =  0.3564
    }
  }
  
  # create a placeholder matrix with a column for each MAC
  fit <- as.data.frame(matrix(nrow =  1,  ncol = mac_bins$Upper[nrow(mac_bins)]))
  
  # calculate the proportion of variants at each individual MAC
  for(i in 1:mac_bins$Upper[nrow(mac_bins)]){
    fit[,i] <- b/((beta+i)^alpha)
  }
  
  # create the dataframe that will be the results
  re  <- mac_bins # define the bins for the results
  re$Prop <- '.' # create a column for the proportion of variants per bin
  
  # For each bin, sum over the individual MAC to get the proportion for the bin
  for(i in 1:nrow(re)){
    re$Prop[i] <- sum(fit[,c(re$Lower[i]:re$Upper[i])])
  }
  
  # make sure the proportions are numeric
  re$Prop <- as.numeric(as.character(re$Prop))

  return(re)
}
