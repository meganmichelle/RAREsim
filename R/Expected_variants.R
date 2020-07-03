#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param alpha parameter from AFS funciton
#' 
#' @param beta another parameter from AFS funciton
#'  
#' @param b final parameter from AFS funciton
#'   
#' @param phi parameter from total variants funciton
#' 
#' @param omega other parameter from total variants funciton
#' 
#' @param Ntar sample size of the target data
#' 
#' @param Size Simulation region sample size in Kb
#' 
#' @param pop Ancestry population - needs to be specified ('AFR', 'NFE', 'EAS', or 'SAS') if using default
#'
#' @return data frame with  the bins and expected variants
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#'
#' @export
#'
#'

Expected_variants <- function(alpha   = NULL, beta = NULL, b = NULL,
                              phi = NULL, omega = NULL, Ntar, Size, pop = NULL){

  if(Ntar > 125000){
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
  
  m <- floor(Ntar*2*0.01) ### RV MAC in the target data
  fit <- as.data.frame(matrix(nrow=1, ncol = m))
  
  for(i in 1:m){
    fit[,i] <- b/((beta+i)^alpha)
  }
  
  #### add the loop matrix here?? that might be easiest
  loop <- as.data.frame(matrix(nrow=7,  ncol = 3))
  colnames(loop) <- c('Lower', 'Upper', 'Prop')
  loop$Lower <- c(1,2,3,6,11,21,m/2)
  loop$Upper <- c(1,2,5,10,20,((m/2)+1),m)
  
  for(k in 1:nrow(loop)){
    loop$Prop[k] <- sum(fit[1,c(loop$Lower[k]:loop$Upper[k])])
  }
  
  loop$Expected_var <- (Size)*(phi*Ntar^omega)*loop$Prop
  
  loop <- loop[,c(1,2,4)]
  
  return(loop)
}

