#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param phi parameter phi
#'
#' @param omega parameter omega
#'
#' @param n sample size in number of individuals
#' 
#' @param pop population - only needs to be specified if  using default  parameters
#'
#' @return the number of variants per kb
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#'
#' @export
#'
#'

Variants_per_Kb<-function(phi=NULL, omega=NULL, n,  pop=NULL){
  if(n>125,000){
    warning('We currently do not recommend simulating sample sizes over 125,000')
  }
  if(is.null(alpha) & is.null(pop)){
    stop('a population must be specified if using default parameters')
  }
  
  if(is.null(phi)){
    if(pop == 'AFR'){
      phi = 0.1576
      omega = 0.6247
    }
    if(pop == 'EAS'){
      phi = 0.1191
      omega = 0.6369
    }
    if(pop == 'NFE'){
      phi = 0.1073
      omega = 0.6539
    }
    if(pop == 'SAS'){
      phi =0.1249
      omega = 0.6495
    }
  }
  if(n > 125000){
    warning('We do not currently recommend simulating sample sizes over 125,000')
  }
  fn <- phi*(n^omega)
  return(fn)
}



