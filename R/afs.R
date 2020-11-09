#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param alpha AFS function parameter alpha
#'
#' @param beta AFS function parameter beta
#'
#' @param b AFS function parameter b
#' 
#' @param mac The MAC bins to use, with lower and upper boundaries defined
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

afs <- function(alpha=NULL, beta=NULL, b=NULL, mac, pop=NULL){
  
  # if(((pop==NULL)==FALSE)  & ((is.numeric(alpha) == FALSE) | (is.numeric(beta) == FALSE) |
  #    (is.numeric(b) == FALSE))){
  #   stop('All parameters are required to be numeric')
  # } ###  This does not work!s
  
  if((colnames(mac)[1] == 'Lower') == FALSE | (colnames(mac)[2]  == 'Upper') == FALSE){
    stop('mac files needs to have column names Lower and Upper')
  }
  
  if((is.numeric(mac$Lower) ==  FALSE) | (is.numeric(mac$Upper) == FALSE)){
    stop('mac columns need to be numberic')
  }
  
  ### check that ALL parameters are null and a population specified
  if(is.null(alpha) & is.null(pop)){
    stop('a population must be specified if using default parameters')
  }

  ####  make sure the option for the population are just afr, eas,  nfe, and  sas?
  
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
  
  fit <- as.data.frame(matrix(nrow =  1,  ncol = mac$Upper[nrow(mac)]))
  
  for(i in 1:mac$Upper[nrow(mac)]){
    fit[,i] <- b/((beta+i)^alpha)
  }
  
  
  re  <- mac
  re$Prop <- '.'
  
  for(i in 1:nrow(re)){
    re$Prop[i] <- sum(fit[,c(re$Lower[i]:re$Upper[i])])
  }
  
  re$Prop <- as.numeric(as.character(re$Prop))
  
  return(re)
}
