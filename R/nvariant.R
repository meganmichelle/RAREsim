#' Number of Variants function
#'
#' The Number of Variants (nvariant) function calculates the number of variants per kilobase, as described in RAREsim.
#' N is the number of individuals. The Number of Variants function changes with N.
#' The default parameters will be used if an ancestrial population is specified: (AFR, EAS, NFE, or SAS) is specified.
#' Otherwise, the parameters phi and omega need to be provided.
#' Phi and omega can be estimated from target data using the Fit_nvariant function.
#'
#' @param phi parameter phi
#'
#' @param omega parameter omega
#'
#' @param N sample size in number of individuals
#'
#' @param pop population - only needs to be specified if  using default  parameters
#'
#' @return the number of variants per kb
#'
#' @examples
#' nvariant(N=8128, pop = 'AFR')
#' nvariant(phi = 0.1638108, omega = 0.6248848, N = 8128)
#'
#' @keywords RAREsim
#'
#' @export
#'

nvariant<-function(phi=NULL, omega=NULL, N,  pop=NULL){
    
    # Check if the simulation sample size is larger than 125k
    if(N>125000){
        warning('We currently do not recommend simulating sample sizes over 125,000')
    }
    
    if(N<2000){
        warning('To simulate <2000 individuals, use RAREsim to simulate 2000 individuals
            and randomly downsample to the desired size')
    }
    
    # If either both parameters should be specified or the ancestry for default parameters
    if((is.null(phi) | is.null(omega)) & is.null(pop)){
        stop('A population must be specified when using default parameters.')
    }
    
    # Ensure the default parameter ancestry is correctly specified
    if(is.null(pop)==FALSE && (pop == 'AFR' | pop == 'EAS' | pop == 'NFE' | pop == 'SAS') == FALSE){
        stop('Default ancestries must be specified as AFR, EAS, NFE, or SAS.')
    }
    
    # If parameters are specified, they must be numeric
    if((is.null(pop)==TRUE)  & ((is.numeric(phi) == FALSE) | (is.numeric(omega) == FALSE))){
        stop('Error: at least one parameter is not numeric')
    }
    
    # Specify the default parameters
    if(is.null(phi) & is.null(omega)){
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
    
    # calculate the number of variants per Kb at sample size N
    fn <- phi*(N^omega)
    
    # Return the calculated function
    return(fn)
}



