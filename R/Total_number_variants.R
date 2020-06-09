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

Total_number_variants<-function(phi, omega, n){
  fn <- phi*(n^omega)
  return(fn)
}



