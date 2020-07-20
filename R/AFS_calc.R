#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param alpha Observed probabilities for each MAC bin
#'
#' @param beta Number of individuals in the target data
#'
#' @param b percent of rare variants - just the sum of the proportions
#' 
#' @param mac The MAC bins to use, with lower and upper boundaries defined
#' 
#' @param pop The population -  specified when using default  parameters
#'
#' @return data frame with the MAC bins provided and proportion of variants in  each bin
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#'
#' @export
#' 
#'

AFS_calc <- function(alpha=NULL, beta=NULL, b=NULL, mac, pop=NULL){
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
  re$Fitted_prop <- '.'
  
  for(i in 1:nrow(re)){
    re$Fitted_prop[i] <- sum(fit[,c(re$Lower[i]:re$Upper[i])])
  }
  
  return(re)
}
