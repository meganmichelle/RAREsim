#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param alpha AFS funciton parameter alpha
#' 
#' @param beta AFS funciton parameter beta
#'  
#' @param b AFS funciton parameter b
#'   
#' @param phi Variants per Kb funciton parameter phi
#' 
#' @param omega Variants per Kb funciton parameter omega
#' 
#' @param mac The MAC bins to use, with lower and upper boundaries defined. Only define for rare variants
#' 
#' @param N sample size to simulate
#' 
#' @param Size Simulation region sample size in Kb
#' 
#' @param pop Ancestry population - needs to be specified ('AFR', 'NFE', 'EAS', or 'SAS') if using default
#'
#' @return data frame with the MAC bins and expected variants
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#' @export


Expected_variants <- function(alpha   = NULL, beta = NULL, b = NULL,
                              phi = NULL, omega = NULL, N, Size, pop = NULL,
                              mac){
  
  # if((is.numeric(alpha) == FALSE) | (is.numeric(beta) == FALSE) |
  #   (is.numeric(b) == FALSE) | (is.numeric(phi) == FALSE) | (is.numeric(omega) == FALSE)){
  #   stop('All parameters are required to be numeric')
  # }
  # 
  ### columns need to be names lower and upper
   if((colnames(mac)[1] == 'Lower') == FALSE | (colnames(mac)[2]  == 'Upper') == FALSE){
     stop('mac files needs to have column names Lower and Upper')
    }
  
  ## mac needs  to have numeric values
  if((is.numeric(mac$Lower) ==  FALSE) | (is.numeric(mac$Upper) == FALSE)){
    stop('mac columns need to be numberic')
  }
  
  ### check that ALL parameters are null and a population specified
  if(is.null(alpha) & is.null(pop)){
    stop('a population must be specified if using default parameters')
  }

  if(N > 125000){
    warning('We do not currently recommend simulating sample sizes over 125,000')
  }
  if(is.null(alpha)){
    if(pop == 'AFR'){
      alpha <- 1.5883
      beta <- -0.3083
      b <- 0.2872
    }
    if(pop == 'EAS'){
      alpha <- 1.6656
      beta <- -0.2951
      b <- 0.3137
    }
    if(pop == 'NFE'){
      alpha <- 1.9470
      beta <- 0.1180
      b <- 0.6676
    }
    if(pop == 'SAS'){
      alpha <- 1.6977
      beta <- -0.2273
      b <-0.3564
    }
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
  
  m <- floor(N*2*0.01) ### RV MAC in the target data
  fit <- as.data.frame(matrix(nrow=1, ncol = m))
  
  for(i in 1:m){
    fit[,i] <- b/((beta+i)^alpha)
  }
  

  mac$Prop <- '.'

  
  for(k in 1:nrow(mac)){
    mac$Prop[k] <- sum(fit[1,c(mac$Lower[k]:mac$Upper[k])])
  }
  mac$Prop <- as.numeric(as.character(mac$Prop))
  
  
  mac$Expected_var <- Size*(phi*N^omega)*mac$Prop
  
  mac <- mac[,c(1,2,4)]
  
  return(mac)
}

